//
// Created by munk on 22-09-15.
//

#ifndef SIMX_DETECTION_DETECTIONSYSTEM_H
#define SIMX_DETECTION_DETECTIONSYSTEM_H

#include "simX/daq/SimpleDAQ.h"
#include <memory>
#include <vector>

namespace simX {
    class Detector;
    namespace detection {
        class Scaler;

        class DetectionSystem {
        public:
            using DetectorPtr = std::shared_ptr<Detector>;
            using ScalerPtr = std::shared_ptr<Scaler>;

            DetectionSystem();

            DetectionSystem(const std::vector<DetectorPtr>& detectors, const std::vector<ScalerPtr>& scalers);

            /**
             * Get a vector with the detectors in this detection system.
             */
            const std::vector<DetectorPtr>& getDetectors() const;


            const std::vector<ScalerPtr>& getScalers() const;

            /**
             * Number of detectors.
             */
            size_t size() const;


            /**
             * Get the DAQ trigger function for this detection system.
             * This will default to a global OR, if none is set.
             */
            const daq::SimpleDAQ::TriggerFunction& getTriggerFunction() const;
            void setTriggerFunction(const daq::SimpleDAQ::TriggerFunction& triggerFunction);

            const daq::SimpleDAQ::ModuleTrigger& getADCTrigger() const;
            void setADCTrigger(const daq::SimpleDAQ::ModuleTrigger&);

            const daq::SimpleDAQ::ModuleTrigger& getTDCTrigger() const;
            void setTDCTrigger(const daq::SimpleDAQ::ModuleTrigger&);

        private:
            std::vector<DetectorPtr> detectors;
            std::vector<ScalerPtr> scalers;
            daq::SimpleDAQ::TriggerFunction triggerFunction;
            daq::SimpleDAQ::ModuleTrigger adc, tdc;

            void init();
        };


    }
}

#endif //SIMX_DETECTION_DETECTIONSYSTEM_H
