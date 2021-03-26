
#ifndef DSSD_h
#define DSSD_h 1

#include "simX/Detection/Detector.h"
#include "simX/Particle.h"
#include "simX/Layer.h"
#include "ausa/eloss/Material.h"
#include "ausa/calibration/DetectorCalibration.h"

// C++ and ROOT header files
#include <vector>
#include <array>
#include <memory>
#include <TVector3.h>

namespace simX {

    /**
    * Class for simX detectors
    */
    class SegmentedDetector : public Detector {

        public:
            using LayerPtr = std::unique_ptr<Layer>;


            virtual ~SegmentedDetector() = default;
            
            virtual const DetectorOutput& detect( Particle& part ) override;
            virtual int getNumberOfChannels() override;
            virtual const Volume& getDetectorVolume() const override;


            bool isSharingEnabled() const {
                return doSharing;
            }

/*
             * Methods to enable/disable sharing
             */            
            inline void enableSharing( bool b ) { doSharing=b; }

            /*
             * Set position
             * @param c Coordinates to center of DSSD (surface)
             */            
            void setPosition( TVector3 c );

            /*
             * Get position
             */            
            TVector3 getPosition() const;

            /*
            * Get direction
            */
            TVector3 getDirection() const;

            /*
            * Get up vector
            */
            TVector3 getUp() const;

            /*
             * Reverse ordering of the i'th segmentation
             * @param i Segmentation
             * @param b Normal order if 0 (false), reverse ordering if 1 (true)
             */
            virtual void reverseStripOrdering(size_t i, bool b );

            /*
             * The detector measures only the energy deposited via ionizing processes. 
             * @param set True/false.
             */
            virtual void detectOnlyIonizingEnergy(bool set) {doIonizing=set;}

            bool isDetectOnlyIonizingEnergyEnabled() const {return doIonizing;};

            virtual const std::string& getName() const override;


            virtual FluctuationFunction getFluctuationFunction() const override;

            virtual void setFluctuationFunction(FluctuationFunction fluctuationFunction) override;

            virtual unsigned int energyToChannel(int channel, double e) override;

            virtual void setCalibration(AUSA::Calibration::DetectorCalibration calibration) override;

            size_t nSegmentations() const {return segmentation.size();}

            UInt_t nSegments(size_t i) const {return segmentation[i].nSegments;}

    protected:
            virtual bool strikesGrid(double u) = 0;
            virtual int gridIndex() = 0;

            /**
             * Represent a group of segments of the detector.
             * Example: A DSSD has front and back strips. Each of these are considered a SegmenGroup.
             * Example: A Pad detector has a SingleGroup with only 1 pixel.
             * Example: A SSD has N strips each belonging to one SegmentGroup.
             */
            struct SegmentGroup {
                /**
                 * Width of each segment.
                 */
                double segmentWidth;

                /**
                 * Start coordinate of strip 1.
                 */
                double coordinateStart;

                /**
                 * Number of segments in this group.
                 */
                UInt_t nSegments;

                /**
                 * Whether the strip numbers should be calculated in reverse order.
                 */
                bool reversed;

                /**
                 * Function that maps the projection and the cross-orientation plane to a coordinate.
                 */
                using CoordinateFunction = std::function<double(const TVector2&)>;
                CoordinateFunction mapToCoordinate;

                /**
                 * Modulate strip numbers.
                 * This is relevant for sharing event with round detectors.
                 * In these case strip -1 should be mapped to strip N.
                 */
                using ModulationFunction = std::function<int(int)>;
                ModulationFunction modulate;

                /**
                 * Function that maps the projection on the cross-orientation plane to the between each segment.
                 */
                using GapFunction = std::function<double(const TVector2&)>;
                GapFunction calculateGap;
            };

            /**
             * Factory class for subclasses of SegmentedDetector.
             */
            struct Args {
                /**
                 * Build the individual layers of the detector.
                 * @return A vector with n layers. One of which must be active.
                 */
                virtual std::vector<LayerPtr> buildLayers() = 0;

                /**
                 * Build outer shell ("bounding box") volume of detector.
                 * @return A single volume equal to the sum of layers.
                 */
                virtual std::unique_ptr<Volume> buildDetectionVolume() = 0;

                /**
                 * @return The detector centroid
                 */
                virtual TVector3 buildCenter() = 0;

                /**
                 * @return The detector normal
                 */
                virtual TVector3 buildNormal() = 0;
                /**
                 * @return The detector orientation
                 */
                virtual TVector3 buildOrientation() = 0;

                /**
                 * @return A vector with the segmentation group of the detector.
                 */
                virtual std::vector<SegmentGroup> buildSegmentation() = 0;

                /**
                 * @return True if sharing should be enabled.
                 */
                virtual bool doSharing() = 0;
                /**
                 * @return True if only the ionizing energyloss should be collected.
                 */
                virtual bool doIonizing() = 0;

                /**
                 * @return Total number of channels in the detector.
                 */
                virtual UInt_t nChannels() = 0;
            };

            /**
             * Constructor.
             */
            SegmentedDetector(std::string name, std::unique_ptr<Args> args);

        private:
            bool doIonizing;
            bool doSharing;
            DetectorOutput detOut;
            int numberOfChannels;
            std::unique_ptr<Volume> detVol;
            std::vector<std::unique_ptr<Layer>> layers;
            TVector3 normal, up, right;
            TVector3 center0;
            std::string name;
            FluctuationFunction fluctuationFunction;
            std::vector<SegmentGroup> segmentation;

            AUSA::Calibration::DetectorCalibration calibration;

            void chargeCollection( TVector2& ru, double& eloss, DetectorOutput& detOut );
            double sharing( double x );
            double fluctuation( int chan, double& e );
   };
}

#endif	/* DSSD_H */
