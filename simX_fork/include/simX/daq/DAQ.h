#ifndef SIMX_DAQ_H
#define SIMX_DAQ_H

#include "simX/Detection/Detector.h"

// C++ and ROOT header files
#include <vector>
#include <simX/Detection/DetectionSimulator.h>


namespace simX {
    namespace daq {

        /**
        * Base class for DAQs
        */
        class DAQ {

        public:

            /**
             * Basic data format for DAQ buffer
             */
            struct BufferItem {
                unsigned int energy;
                bool trigger;
            };
            using Buffer = std::vector<std::vector<BufferItem>>;

            /**
             * Basic data format for DAQ output
             */
            struct OutputItem {
                int detector;
                int channel;
                unsigned int energy;
                double time;
            };            
            using Output = std::vector<OutputItem>;
            
            DAQ() = default;
            virtual ~DAQ() = default;

            /**
            * Feed package of data to the DAQ
            * @param n Index used to identify detector
            * @param dO Data package
            */
            virtual void feed( size_t n, const simX::Detector::DetectorOutput& dO ) = 0;

            virtual void feed(const simX::DetectionSimulator::Detections& detections) = 0;

            /**
             * Returns true if the DAQ has received a trigger signal.
             * @param out Contains the entries above the adc threshold.
             */
            virtual bool getData( Output& out ) = 0;

            /**
             * This will clear the readout buffer.
             *
             * @post After calling this size of the inner level be 0.
             */
            virtual void clear() = 0;
        };
    }
}

#endif	/* SIMX_DAQ_H */
