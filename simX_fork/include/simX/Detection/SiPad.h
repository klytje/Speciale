
#ifndef SIPAD_h
#define SIPAD_h 1

#include "simX/Particle.h"
#include "simX/Layer.h"
#include "simX/Detection/SegmentedDetector.h"
#include "ausa/eloss/Material.h"
#include "ausa/calibration/LinearCalibration.h"

// C++ and ROOT header files
#include <vector>
#include <array>
#include <memory>
#include <TVector3.h>

namespace simX {

    /**
    * Class for simX detectors
    */
    class SiPad : public SegmentedDetector {

        public:
            SiPad(  std::string name, TVector3 center, TVector3 normal, TVector3 up,
                    double activeVolumeThickness, double deadlayerFrontThickness, double deadlayerBackThickness,
                    double transverseLength, double upLength, bool doIonizing = false);
                     
            virtual ~SiPad() = default;

        virtual std::string description() override;

    protected:
        bool strikesGrid(double u) override;

        int gridIndex() override;

        struct Args;
   };
}

#endif	/* SIPAD_H */
