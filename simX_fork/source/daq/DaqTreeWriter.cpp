//
// Created by munk on 23-09-15.
//

#include "simX/daq/DaqTreeWriter.h"
#include "simX/Detection/Scaler.h"
#include "simX/parser/JsonOutputParser.h"
#include "simX/git.h"

#include <TTree.h>
#include <TObjString.h>
#include <ausa/json/IO.h>
#include <ausa/util/FileUtil.h>
#include <ausa/AUSA.h>
#include <simX/Detection/SegmentedDetector.h>

simX::daq::DaqTreeWriter::DaqTreeWriter(std::shared_ptr<AUSA::TFileWrapper> output,
                                        const std::vector<std::shared_ptr<Detector>>& system,
                                        const std::vector<std::shared_ptr<detection::Scaler>>& scalers,
                                        const std::string &ausalibFile)
    : wrapper(output), scalers(scalers)
{
    TDirectory::TContext cxt{nullptr};
    file = &wrapper -> get();
    file -> cd();
    tree = new TTree("h101", "h101");

    simX::parser::JsonOutputParser parser(AUSA::asciiFileToString(ausalibFile));

    // Prereserve space so pointers stay valid.
    fMul.reserve(system.size()); bMul.reserve(system.size());
    fSeg.reserve(system.size()); bSeg.reserve(system.size());
    fEne.reserve(system.size()); bEne.reserve(system.size());
    fTdc.reserve(system.size()); bTdc.reserve(system.size());

    for (size_t i = 0; i < system.size(); ++i) {
        auto& d = system[i];
        fMul.push_back(0);
        fSeg.emplace_back(d->getNumberOfChannels());
        fEne.emplace_back(d->getNumberOfChannels());
        fTdc.emplace_back(d->getNumberOfChannels());

        bMul.push_back(0);
        bSeg.emplace_back(d->getNumberOfChannels());
        bEne.emplace_back(d->getNumberOfChannels());
        bTdc.emplace_back(d->getNumberOfChannels());

        auto abs = std::dynamic_pointer_cast<SegmentedDetector>(d);
        if (!abs) throw std::runtime_error("There is current only support for subclasses of SegmentedDetector");

        if (abs->nSegmentations() == 1) {
            auto& m = parser.getMapping(d->getName());
            mapMapping(m, fMul.data()+i, fSeg[i].data(), fEne[i].data(), fTdc[i].data());
        }
        else if (abs->nSegmentations() == 2){
            auto& m = parser.getDoubleMapping(d->getName());
            mapMapping(m.front, fMul.data()+i, fSeg[i].data(), fEne[i].data(), fTdc[i].data());
            mapMapping(m.back, bMul.data()+i, bSeg[i].data(), bEne[i].data(), bTdc[i].data());
        }
        else {
            throw std::runtime_error("There is currently only support for single or double segmented detectors");
        }
        frontStrips.push_back(abs->nSegments(0));
    }

    // Map scalers to tree.
    for (auto& scaler : scalers) {
        map(scaler->getName(), scaler->getPointer());
    }

    TObjString ausaSha(AUSA::GIT_HASH);
    TObjString ausaBranch(AUSA::GIT_BRANCH);
    TObjString simXSha(simX::GIT_HASH);
    TObjString simXBranch(simX::GIT_BRANCH);

    file -> WriteTObject(&ausaBranch, "AUSALIB_HASH");
    file -> WriteTObject(&ausaSha, "AUSALIB_BRANCH");
    file -> WriteTObject(&simXSha, "SIMX_HASH");
    file -> WriteTObject(&simXBranch, "SIMX_BRANCH");
}

simX::daq::DaqTreeWriter::DaqTreeWriter(std::shared_ptr<AUSA::TFileWrapper> output,
                                        const simX::detection::DetectionSystem& system,
                                        const std::string& ausalibFile)
    : DaqTreeWriter(output, system.getDetectors(), system.getScalers(), ausalibFile)
{

}

simX::daq::DaqTreeWriter::~DaqTreeWriter() {
    if (file != nullptr) {
        file -> WriteTObject(tree);
        file -> Flush();
    }
}

void simX::daq::DaqTreeWriter::map(UInt_t* p, const std::string& b, const std::string& mul) {
    if (b != AUSA::JSON::IGNORE) tree->Branch(b.c_str(), p, (b + "[" + mul + "]/i").c_str());
}

void simX::daq::DaqTreeWriter::map(const std::string& b, UInt_t *p) {
    if (b != AUSA::JSON::IGNORE) tree->Branch(b.c_str(), p);
}

void simX::daq::DaqTreeWriter::mapMapping(simX::parser::JsonOutputParser::Mapping& m, UInt_t *mul, UInt_t *seg,
                                          UInt_t *adc, UInt_t *tdc) {
    if (m.mul != AUSA::JSON::IGNORE) {
        map(m.mul, mul);
        map(seg, m.seg, m.mul);
        map(adc, m.adc, m.mul);
        map(tdc, m.tdc, m.mul);
    }
    else {
        map(m.adc, adc);
        map(m.tdc, tdc);

        manualClear.push_back(adc);
        manualClear.push_back(tdc);
    }
}

void simX::daq::DaqTreeWriter::feed(DAQ::Output& output) {
    clear();
    for (auto& b : output) {

        int i = b.detector;
        auto max = frontStrips[i];
        if (b.channel < max) {
            auto& m = fMul[i];
            fSeg[i][m] = b.channel + 1;
            fEne[i][m] = b.energy;
            fTdc[i][m] = b.time;
            ++m;
        }
        else {
            auto& m = bMul[i];
            bSeg[i][m] = b.channel - max + 1;
            bEne[i][m] = b.energy;
            bTdc[i][m] = b.time;
            ++m;
        }
    }
    tree->Fill();
}

void simX::daq::DaqTreeWriter::clear() {
    for (auto& a : fMul) a = 0;
    for (auto& a : bMul) a = 0;

    for (auto ptr : manualClear) *ptr = 0;
}

