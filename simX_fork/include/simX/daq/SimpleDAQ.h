
#ifndef SIMPLEDAQ_h
#define SIMPLEDAQ_h 1

#include "simX/daq/DAQ.h"

// C++ and ROOT header files
#include <vector>
#include <memory>
#include <boost/dynamic_bitset.hpp>


namespace simX {
    namespace daq {
        /**
          * Class for simplified data acquisition system with no time structure
          */
        class SimpleDAQ : public DAQ {
        public:
            /**
             * Contains the trigger status of the individual detectors.
             */
            using TriggerStatus = boost::dynamic_bitset<size_t>;

            /**
              * Allows for user-specified trigger condition.
              * The default is OR on all channels.
              */
            using TriggerFunction = std::function<bool(const TriggerStatus&)>;


            using ModuleTrigger = std::function<bool(int detector, int channel, const DAQ::BufferItem& b)>;

            SimpleDAQ( std::vector<std::shared_ptr<Detector>> detectors,
                       ModuleTrigger tdcTrigger, ModuleTrigger adcTrigger );
            
            virtual ~SimpleDAQ();

            virtual void feed( size_t n, const simX::Detector::DetectorOutput& dO ) override;
            virtual void feed(const simX::DetectionSimulator::Detections& detections) override;

            virtual bool getData( Output& out ) override;
            virtual void clear() override;
            
            void setCommonTDCThreshold( double thres );
            void setCommonADCThreshold( double thres );


            void setTrigger(const TriggerFunction& trigger);

        private:
            void clearBuffer();
            double addNoise();

            DAQ::Buffer buffer;
            std::vector<std::pair<int,int>> hits;
            TriggerStatus status;

            TriggerFunction trigger;
            ModuleTrigger tdcTrigger;
            ModuleTrigger adcTrigger;

            std::vector<std::shared_ptr<Detector>> detectors;
        };

    }
}

#endif	/* SIMPLEDAQ_H */
