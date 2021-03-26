//
// Created by munk on 07-08-15.
//

#include "simX/parser/ConfigParser.h"
#include "simX/Layer.h"
#include "simX/Beam.h"
#include "simX/propagator/IonizingPropagator.h"

#include <ausa/json/JSONUtil.h>
#include <ausa/util/Resource.h>
#include <ausa/eloss/Default.h>
#include <ausa/eloss/SRIMTabulation.h>
#include <ausa/eloss/Material.h>
#include <ausa/eloss/RangeInterpolator.h>
#include <ausa/eloss/EnergyLossRangeInverter.h>
#include <ausa/util/memory>
#include <ausa/util/FileUtil.h>
#include <simX/propagator/NonIonizingPropagator.h>
#include <simX/propagator/NoTrackPropagator.h>
#include <cmath>
#include <simX/parser/ElossFactory.h>


using namespace std;
using namespace simX;
using namespace simX::parser;
using namespace AUSA::EnergyLoss;
using namespace AUSA::JSON;
using namespace rapidjson;

using PropagatorFactories = std::map<std::string, ConfigParser::PropagatorFactory>;

namespace {

    Config::PropagatorFactory buildPropagatorFactory(ConfigParser::Options& opt, const char* propagator, PropagatorFactories& m) {

        ConfigParser::Options* options = &opt;

        string name = "ELOSS";
        if(opt.HasMember(propagator)) {
            auto& object = opt[propagator];

            if (object.IsString())
                name = object.GetString();
            else {
                name = readString(object, "type");
                options = &object;
            }
        }

        if (!m.count(name)) throw invalid_argument("Unknown propagator: " + name);
        return m[name](*options);
    }

    boost::optional<string> readOptString(Value& v, const char* member) {
        return v.HasMember(member) ? boost::optional<string>(readString(v, member)) : boost::optional<string>();
    }
}

ConfigParser::ConfigParser() {
    registerFactory("ELOSS", elossFactory);
    registerFactory("NOLOSS", noLossFactory);
    registerFactory("GAUSSSTRAG", gaussianStragglingFactory);
    registerFactory("MCSTRAG", mcStragglingFactory);
}

Config ConfigParser::parse(const std::string& s) {
    Document doc;
    tryParse(s, doc);

    auto targetFactory = buildPropagatorFactory(doc, "target_propagator", propFactories);
    auto detectionFactory = buildPropagatorFactory(doc, "detection_propagator", propFactories);
    auto beamFactory = buildPropagatorFactory(doc, "beam_propagator", propFactories);

    auto target = readOptString(doc, "target");
    auto beam = readOptString(doc, "beam");
    auto reaction = readOptString(doc, "reaction");
    auto detSys = readOptString(doc, "detection_system");

    return Config{targetFactory, detectionFactory, beamFactory, beam, target, reaction, detSys};
}

Config ConfigParser::parseFile(const std::string& file) {
    try {
        return parse(AUSA::asciiFileToString(file));
    } catch (...) {
        cerr << "Failed to parse simX config file " << file << endl;
        throw;
    }
}

void ConfigParser::registerFactory(std::string key, PropagatorFactory f) {
    propFactories[key] = f;
}


