//
// Created by munk on 22-09-15.
//

#include "simX/parser/DetectionSystemParser.h"
#include "simX/Detection/S3.h"
#include "simX/Detection/W1.h"
#include "simX/Detection/SiPad.h"
#include "simX/Detection/YY1.h"
#include "simX/Detection/Scaler.h"

#include <ausa/util/FileUtil.h>
#include <ausa/json/IO.h>
#include <ausa/json/WrongTypeException.h>
#include <ausa/setup/RoundDSSD.h>
#include <ausa/setup/SquareDSSD.h>
#include <ausa/setup/PadDetector.h>
#include <ausa/setup/YY1.h>
#include <ausa/util/memory>
#include <ausa/json/JSONUtil.h>

#include <rapidjson/document.h>
#include <simX/daq/SimpleDAQ.h>
#include <simX/parser/DAQTriggerParser.h>
#include <simX/parser/AdcThresholdParser.h>
#include <simX/parser/TdcThresholdParser.h>
#include <ausa/util/UnitParser.h>
#include <TF1.h>
#include <TSystem.h>
#include <ausa/util/StringUtil.h>

using namespace std;
using namespace simX::detection;
using namespace simX::parser;
using namespace AUSA;
using namespace AUSA::JSON;
using namespace rapidjson;

namespace {
    auto CONVERSION_FACTOR = 2*sqrt(2*log(2));

    simX::Detector::FluctuationFunction parseDSSDFluctuation(Value& opt, int frontCount) {
        if (!opt.HasMember("resolution")) return [](int, double) { return 0;};

        auto& obj = readObject(opt, "resolution");
        auto type = readString(obj, "type");
        if (type != "GAUSS") throw invalid_argument("Currently we only support GAUSS resolution!");

        UnitParser energyParser('k', "eV");

        auto& res = readObject(obj, "fwhm");
        auto fSigma = energyParser.parse(readString(res, "front")) / CONVERSION_FACTOR;
        auto bSigma = energyParser.parse(readString(res, "back")) / CONVERSION_FACTOR;

        TRandom3 fR;
        TRandom3 bR;

        return [=](int chan, double E) mutable {
            if (chan <= frontCount) return fR.Gaus(0, fSigma);
            else return bR.Gaus(0, bSigma);
        };
    }

    simX::Detector::FluctuationFunction parseSSDFluctuation(Value& opt, int frontCount) {
        if (!opt.HasMember("resolution")) return [](int, double) { return 0;};

        auto& obj = readObject(opt, "resolution");
        UnitParser energyParser('k', "eV");

        auto type = readString(obj, "type");
        if (type != "GAUSS") throw invalid_argument("Currently we only support GAUSS resolution!");

        auto sigma = energyParser.parse(readString(obj, "fwhm")) / CONVERSION_FACTOR;

        TRandom3 random;

        return [=](int chan, double E) mutable {
            return random.Gaus(0, sigma);
        };
    }
}

struct DetectionSystemParser::Impl {
    using Builders = std::unordered_map<std::type_index, DetectorBuilder>;
    using ModuleTrigger = simX::daq::SimpleDAQ::ModuleTrigger;
    using AdcTrigger = ModuleTrigger;
    using TdcTrigger = ModuleTrigger;

    Impl(Builders& b, bool inverseCalibrate) : builders(b), energyParser('k', "eV"), inverseCalibrate(inverseCalibrate) {
        empty.SetObject();
    }

    DetectionSystemParser::Result parse(const string& s, std::string prefix) {
        Document document;
        tryParse(s, document);

        string setupFile;
        shared_ptr<Setup> setup;
        setupFile = readString(document, "ausalib_setup");

        try{
            setup = AUSA::JSON::readSetupFromJSON(setupFile);
        } catch(ios_base::failure e) {
            setupFile = prefix + setupFile;
            setup = AUSA::JSON::readSetupFromJSON(setupFile);
        }

        auto& json = document.HasMember("detectors") ? document["detectors"] : empty;

        auto detSys = buildSetup(*setup, json);
        detSys.setTriggerFunction(buildTriggerFunction(document, detSys));

        return Result{setupFile, detSys};
    }

    DetectionSystem buildSetup(AUSA::Setup& s, Value& json) {
        vector<DetectionSystem::DetectorPtr> detectors;
        size_t dCount = s.dssdCount();
        vector<AdcTrigger> adcTriggers;
        vector<TdcTrigger> tdcTriggers;

        Document d;
        auto& alloc = d.GetAllocator();
        auto& global = json.HasMember("__GLOBAL__") ? json["__GLOBAL__"] : empty;

        for (size_t i = 0; i < dCount + s.singleCount(); i++) {
            AUSA::Detector* detector;
            if (i < dCount) detector = s.getDSSD(i).get();
            else detector = s.getSingleSided(i-dCount).get();

            auto& ref = *detector;
            auto index = std::type_index(typeid(ref));

            // Make deep copy of global json
            Value opt(global, alloc);
            if (json.HasMember(detector->getName().c_str())) merge(opt, json[detector->getName().c_str()], alloc); // Merge local options with global

            if (builders.count(index)) {
                auto simDetector = builders[index](ref, opt);
                if (inverseCalibrate) {
                    if (!ref.isCalibrated())
                        throw std::runtime_error("Inverse calibration of " + ref.getName() + " is requested. But it has no calibration.");

                    simDetector->setCalibration(*ref.getCalibration());
                }
                detectors.emplace_back(simDetector);

                auto adcString = readStringOrDefault(opt, "adc_threshold", "0keV");
                adcTriggers.push_back(buildModuleTrigger(adcString, *detectors.back(), parseAdcThresholdsFromFile, "ADC"));

                auto tdcString = readStringOrDefault(opt, "trigger_threshold", "0keV");
                if (AUSA::endswith(tdcString, ".root"))
                    tdcTriggers.push_back(buildModuleTrigger(tdcString, *detectors.back(), parseTdcThresholdsFromFile, "TDC"));
                else
                    tdcTriggers.push_back(buildModuleTrigger(tdcString, *detectors.back(), parseAdcThresholdsFromFile, "TDC"));
            }
            else
                cerr << "I do not have a simX model for '" << ref.getName() << "' with type " << index.name() << endl;
        }

        std::vector<DetectionSystem::ScalerPtr> scalers;
        for (auto& scalerName : s.getScaler()) {
            scalers.emplace_back(std::make_shared<Scaler>(scalerName, 29912922));
        }

        DetectionSystem detSys{move(detectors), move(scalers)};

        detSys.setADCTrigger([=](int det, int chan, const daq::DAQ::BufferItem& b) {
            return adcTriggers[det](det, chan, b);
        });

        detSys.setTDCTrigger([=](int det, int chan, const daq::DAQ::BufferItem& b) {
            return tdcTriggers[det](det, chan, b);
        });

        return detSys;
    }

    template <class FileConstructor>
    daq::SimpleDAQ::ModuleTrigger buildModuleTrigger(const std::string& input, Detector& d, FileConstructor f, std::string type) {
        try {
            if (gSystem->GetPathInfo(input.c_str(), nullptr, (Long64_t*) nullptr, nullptr, nullptr) == 0) {
                auto thresholds = f(input, d.getNumberOfChannels());
                return [=](int det, int chan, const daq::DAQ::BufferItem& b) {
                    return thresholds.isTriggering(chan, b.energy);
                };
            }
            else {
                auto threshold = UnitParser{'k', "eV"}.parse(input);
                return [=](int det, int chan, const daq::DAQ::BufferItem& b) {
                    return b.energy > threshold;
                };
            }
        } catch (exception& e) {
            cerr << e.what() << endl;
            cerr << "Failed to parse " << type << " threshold entry for " << d.getName() << endl;
            throw;
        }
    }

    daq::SimpleDAQ::TriggerFunction buildTriggerFunction(Value& v, const DetectionSystem& detectionSystem) {
        if (!v.HasMember("daq_trigger")) {
            return [](const daq::SimpleDAQ::TriggerStatus& s) {
                return s.any();
            };
        }

        auto str = readString(v, "daq_trigger");
        return parser::parseTriggerFunction(str, detectionSystem);
    }

    void merge(Value& opt, Value& det, Document::AllocatorType& alloc) {
        for (auto i = det.MemberBegin(); i != det.MemberEnd(); ++i) {
            if (opt.HasMember(i->name)) opt.RemoveMember(i->name);
            opt.AddMember(i->name, i->value, alloc);
        }
    }

    Builders& builders;
    Value empty;
    UnitParser energyParser;
    bool inverseCalibrate;
};

simX::parser::DetectionSystemParser::DetectionSystemParser(bool inverseCalibrate) {
    registerBuilder<RoundDSSD>([](AUSA::Detector& in, Value& v) {
        auto& d = static_cast<RoundDSSD&>(in);

        auto ringGap = d.getRingPitch()-d.getRingWidth();
        auto spokeGap = d.getSpokeGap();

        auto share = readBoolOrDefault(v, "sharing", false);

        auto ionizing = readBoolOrDefault(v, "ionizing", false);

        auto det = std::make_shared<S3>(in.getName(), d.getCenter(), d.getNormal(), d.getOrientation(),
                                    d.frontStripCount(), d.backStripCount(), d.getInnerRadius(), d.getRingPitch(), ringGap, spokeGap,
                                    share, d.getThickness(), d.getFrontDeadLayer(), d.getBackDeadLayer(),
                                    d.isFrontStripsReversed(), d.isBackStripsReversed());

        det->setFluctuationFunction(parseDSSDFluctuation(v, d.frontStripCount()));
        det->detectOnlyIonizingEnergy(ionizing);

        return det;
    });

    registerBuilder<SquareDSSD>([](AUSA::Detector& in, Value& v) {
        auto& d = static_cast<SquareDSSD&>(in);

        auto fP = d.getFrontPitch();
        auto bP = d.getBackPitch();
        auto fG = fP-d.getFrontWidth();
        auto bG = bP-d.getBackWidth();

        auto share = readBoolOrDefault(v, "sharing", false);

        auto ionizing = readBoolOrDefault(v, "ionizing", false);

        auto det = std::make_shared<W1>(in.getName(), d.frontStripCount(), d.backStripCount(),
                                    d.getCenter(), d.getNormal(), d.getOrientation(),
                                    share, d.getThickness(), d.getFrontDeadLayer(), d.getBackDeadLayer(),
                                    d.getGridThickness(), d.getBackContactThickness(),
                                    fP, bP, fG, bG,
                                    d.getGridWidth(), d.isFrontStripsReversed(), d.isBackStripsReversed());

        det->setFluctuationFunction(parseDSSDFluctuation(v, d.frontStripCount()));
        det->detectOnlyIonizingEnergy(ionizing);

        return det;
    });
    registerBuilder<PadDetector>([](AUSA::Detector& in, Value& v) {
        auto& d = static_cast<PadDetector&>(in);

        auto ionizing = readBoolOrDefault(v, "ionizing", false);

        auto det = std::make_shared<SiPad>(in.getName(),
                                    d.getCenter(), d.getNormal(), d.getOrientation(),
                                    d.getThickness(), d.getDeadLayer(), d.getBackDeadLayer(),
                                    d.getTransverseSize(), d.getUpSize());

        det->setFluctuationFunction(parseSSDFluctuation(v, 1));
        det->detectOnlyIonizingEnergy(ionizing);

        return det;
    });

    registerBuilder<YY1>([](AUSA::Detector& in, Value& v) {
        auto& d = static_cast<YY1&>(in);

        auto ionizing = readBoolOrDefault(v, "ionizing", false);

        auto det = std::make_shared<Detection::YY1>(in.getName(),
                                           d.getCenter(), d.getNormal(), d.getOrientation(),
                                           d.getThickness(), d.getDeadLayer(), d.getBackDeadLayer());

        det->setFluctuationFunction(parseSSDFluctuation(v, 1));
        det->detectOnlyIonizingEnergy(ionizing);

        det->reverseStripOrdering(0, d.isReverse());

        return det;
    });

    pimpl = std::make_unique<Impl>(builders, inverseCalibrate);
}

DetectionSystemParser::Result simX::parser::DetectionSystemParser::parseString(const std::string& s, std::string prefix) {
    return pimpl -> parse(s, prefix);
}

DetectionSystemParser::Result DetectionSystemParser::parseFile(const std::string& input) {
    try {
        return parseString(AUSA::asciiFileToString(input), AUSA::extractDirectory(input) + "/");
    } catch (...) {
        cerr << "Failed to parse detection system file: " << input << endl;
        throw;
    }
}

DetectionSystem DetectionSystemParser::buildSimXSetup(AUSA::Setup& s) {
    return pimpl->buildSetup(s, pimpl->empty);
}

DetectionSystemParser::~DetectionSystemParser() {
    // Does nothing...
}
