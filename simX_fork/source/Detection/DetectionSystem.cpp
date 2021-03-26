//
// Created by munk on 15-12-15.
//

#include "simX/Detection/DetectionSystem.h"

using namespace simX;
using namespace simX::detection;

namespace {
    bool globalOr(const daq::SimpleDAQ::TriggerStatus& s) {
        return s.any();
    }

    bool moduleOR(int det, int chan, const daq::DAQ::BufferItem& b) {
        return b.energy > 0;
    }
}

DetectionSystem::DetectionSystem()
{
    init();
}

DetectionSystem::DetectionSystem(const std::vector<DetectionSystem::DetectorPtr>& detectors,
                                 const std::vector<ScalerPtr>& scalers)
        : detectors(detectors), scalers(scalers)
{
    init();
}

void DetectionSystem::init() {
    triggerFunction = globalOr;
    adc = tdc = moduleOR;
}

const std::vector<DetectionSystem::DetectorPtr>& DetectionSystem::getDetectors() const {
    return detectors;
}

const std::vector<DetectionSystem::ScalerPtr>& DetectionSystem::getScalers() const {
    return scalers;
}

size_t DetectionSystem::size() const {
    return detectors.size();
}

const daq::SimpleDAQ::TriggerFunction& DetectionSystem::getTriggerFunction() const {
    return triggerFunction;
}

void DetectionSystem::setTriggerFunction(const daq::SimpleDAQ::TriggerFunction& triggerFunction) {
    DetectionSystem::triggerFunction = triggerFunction;
}


const daq::SimpleDAQ::ModuleTrigger& DetectionSystem::getADCTrigger() const {
    return adc;
}

void DetectionSystem::setADCTrigger(const daq::SimpleDAQ::ModuleTrigger& trigger) {
    adc = trigger;
}

const daq::SimpleDAQ::ModuleTrigger& DetectionSystem::getTDCTrigger() const {
    return tdc;
}

void DetectionSystem::setTDCTrigger(const daq::SimpleDAQ::ModuleTrigger& trigger) {
    tdc = trigger;
}



