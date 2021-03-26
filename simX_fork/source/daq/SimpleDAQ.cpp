#include "simX/daq/SimpleDAQ.h"
#include <math.h>
#include <iostream>

using namespace std;
using namespace simX;
using namespace simX::daq;


namespace {
    // Define your trigger logic here
    bool globalOR (const SimpleDAQ::TriggerStatus& buf ) {
        return buf.any();
    }

    vector<size_t> computeChannels(std::vector<std::shared_ptr<Detector>>& d) {
        vector<size_t> res;
        for (auto& p : d) {
            res.push_back(p->getNumberOfChannels());
        }
        return res;
    }
}

SimpleDAQ::SimpleDAQ(std::vector<std::shared_ptr<Detector>> detectors, ModuleTrigger tdcTrigger,
                     ModuleTrigger adcTrigger)
    : trigger(globalOR), tdcTrigger(tdcTrigger), adcTrigger(adcTrigger), detectors(detectors)
{
    auto numberOfChannels = computeChannels(detectors);
    // Resize buffer and threshold-vectors (to minimize memory usage)
    auto N = numberOfChannels.size();
    buffer.resize(N);
    status.resize(N);
    for (int i = 0; i < N; i++) {
        auto L = numberOfChannels[i];
        buffer[i].resize(L);
    }

    // Clear buffer
    clearBuffer();
}


SimpleDAQ::~SimpleDAQ() 
{    
    // do nothing    
}

void SimpleDAQ::feed(const simX::DetectionSimulator::Detections& detections) {
    for (auto& hit : detections) {
        feed(hit.detectorID, hit.output);
    }
}

void SimpleDAQ::feed( size_t n, const simX::Detector::DetectorOutput& dO ) {
    // loop over data items
    for (auto& d : dO) {

        // access data
        int cl = d.channel; // channel
        double el = d.energy; // energy

        el = detectors[n]->energyToChannel(cl, el);

        // store data in buffer
        auto& a = buffer[n][cl];
        a.energy += el;

        hits.emplace_back(n, cl);
    }
}

bool SimpleDAQ::getData( DAQ::Output& out ) {

    // Make sure output is empty
    out.clear();
    status.reset(); // Set all to 0

    std::sort(hits.begin(), hits.end());
    hits.erase(std::unique(hits.begin(), hits.end()),hits.end());

    for (auto& h : hits) {
        auto det = h.first;
        auto chan = h.second;

        auto& bij = buffer[det][chan];

        bij.energy += addNoise();   // add electronics noise
        if (tdcTrigger(det, chan, bij)) status[det] |= bij.trigger = true;  // is energy above tdc threshold?
        if (!adcTrigger(det, chan, bij)) bij.energy = 0.;    // is energy below adc threshold?
        // fill output
        if (bij.energy>0) {
             out.emplace_back(OutputItem{det, chan, bij.energy, 0});
        }
    }
    
    // Determine if data triggers DAQ based on user-specified trigger logic
    bool tr = trigger(status);

    // If so, clear buffer
    if (!tr) out.clear();

    return tr;
}


void SimpleDAQ::clearBuffer() {
    for (auto& h : hits) {
        auto& bij = buffer[h.first][h.second];
        bij.energy = 0;
        bij.trigger = false;
    }

    hits.clear();
}

/*
bool SimpleDAQ::defaultTrigger() {
    // Loop over entries in buffer
    for (int i=0; i<buffer.size(); i++) { 
        auto& bi = buffer[i];
        for (int j=0; j<bi.size(); j++) {
            auto& bij = bi[j];            
            if (bij.trigger) return true;
        }
    }
    return false;
}
*/

void SimpleDAQ::clear() {
    clearBuffer();
}


double SimpleDAQ::addNoise() {

    // *** not yet implemented !!!
    return 0;
}


void SimpleDAQ::setCommonTDCThreshold( double thres ) {
    tdcTrigger = [=](int i, int j, const DAQ::BufferItem& b) {
        return b.energy > thres;
    };
}


void SimpleDAQ::setCommonADCThreshold( double thres ) {
    adcTrigger = [=](int i, int j, const DAQ::BufferItem& b) {
        return b.energy > thres;
    };
}

void SimpleDAQ::setTrigger(const TriggerFunction& trigger) {
    SimpleDAQ::trigger = trigger;
}


