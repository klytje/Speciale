#ifndef SIMX_DETECTION_DETECTOR_H
#define SIMX_DETECTION_DETECTOR_H

// C++ and ROOT header files
#include "simX/Layer.h"
#include <vector>
#include <string>
#include <ausa/geometry/Volume.h>
#include <ausa/calibration/DetectorCalibration.h>
#include <functional>

namespace simX {

    // Forward declarations
    class Particle;

    /**
    * Base class for simX detectors
    */
    class Detector {

        public:
            using FluctuationFunction = std::function<double(int chan, double E)>;
            using Volume = AUSA::Geometry::Volume;

            struct DataTuple {
                int channel;
                double energy, time;
            };

            /**
             * Output from a detector.
             * Will contain multiple entries if there is multiple hits in the detector.
             */
            using DetectorOutput = std::vector<DataTuple>;

            Detector() = default;
            virtual ~Detector() = default;
            
            /**
             * Propagates particle in detector and determines 
             * energy losses in active volumes.
             * @param part Particle.
             */
            virtual const DetectorOutput& detect( Particle& part ) = 0;

            /**
             * Returns number of channels
             */
            virtual int getNumberOfChannels() = 0;
            
            /**
             * Returns union of layers (material not defined)
             */            
            virtual const Volume& getDetectorVolume() const = 0;

            virtual const std::string& getName() const = 0;

            virtual std::string description();

            virtual FluctuationFunction getFluctuationFunction() const = 0;

            virtual void setFluctuationFunction(FluctuationFunction f) = 0;

            /**
             * Do an inverse calibration ie. energy to ADC channel.
             * @param channel Which detector channel fired.
             * @param e The deposited energy
             * @pre channel is indexed from 0!
             * @return The channel number corresponding to deposited energy.
             */
            virtual unsigned int energyToChannel(int channel, double e) = 0;

            virtual void setCalibration(AUSA::Calibration::DetectorCalibration calibration) = 0;

            friend std::ostream& operator <<(std::ostream&, Detector& det);

        private:
   };
}

#endif	/* SIMX_DETECTION_DETECTOR_H */
