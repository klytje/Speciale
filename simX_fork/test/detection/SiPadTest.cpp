
#include <unittest++/UnitTest++.h>
#include <ausa/eloss/Ion.h>
#include <TVector3.h>

#include "simX/Particle.h"
#include "simX/Detection/Detector.h"
#include "simX/Detection/SegmentedDetector.h"
#include <simX/Detection/SiPad.h>
#include "simX/Layer.h"
#include "simX/Random.h"
#include <simX/propagator/ParticlePropagator.h>
#include <simX/propagator/IonizingPropagator.h>

#include <ausa/geometry/Box.h>
#include <ausa/eloss/SRIMTabulation.h>
#include <ausa/eloss/RangeInterpolator.h>
#include <ausa/eloss/EnergyLossRangeInverter.h>
#include <ausa/util/Resource.h>
#include <ausa/util/memory>

#include <math.h>
#include <TMath.h>
#include <tuple>
#include <iostream>

using namespace simX;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Geometry;
using namespace simX::propagator;

SUITE(SiPadTest) {

    // use SRIM stopping powers for alpha in Silicon
    IonizingPropagator::LossCalc factory(const Layer& layer, const Particle& p) {
        std::string suffix = "/SRIM08/He4_Si.dat";
        SRIMTabulation tab(AUSA::getResourceDirectory() + suffix);
        return std::make_unique<EnergyLossRangeInverter>(std::make_unique<RangeInterpolator>(tab));
    }

    class SetupFixture {
        public:
            SetupFixture() 
            : proton(Ion(1,1)),
              alpha(Ion(2,4)),
              gamma(Ion(0,0)),
              ionProp( std::make_shared<IonizingPropagator>(IonizingPropagator(factory)) ),
              pad("Pad", {0,0,100}, {0,0, -1}, {0,1,0}, 1, 0, 0, 40, 50, false)
            {
                alpha.setPropagator(ionProp);
                proton.setPropagator(ionProp);
            }
            Particle proton, alpha, gamma;
            std::shared_ptr<IonizingPropagator> ionProp;
            SiPad pad;
    };

    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1,1);
    }

    TEST_FIXTURE(SetupFixture, DetectoVolumeHasCorrectPosition) {
        auto& vol = pad.getDetectorVolume();
        CHECK_CLOSE( 0, vol.getCenter().X(), 1E-3 );
        CHECK_CLOSE( 0, vol.getCenter().Y(), 1E-3 );
        CHECK_CLOSE( 100.5, vol.getCenter().Z(), 1E-1 );
    }

    TEST_FIXTURE(SetupFixture, DetectoVolumeHasCorrectDimensions) {
        auto vol = &pad.getDetectorVolume();
        auto box = dynamic_cast<const Box*>(vol);

        CHECK_CLOSE( 40, box->getXDim(), 1E-3 );
        CHECK_CLOSE( 50, box->getUpDim(), 1E-3 );
        CHECK_CLOSE( 1, box->getNDim(), 1E-3 );
    }

    TEST_FIXTURE(SetupFixture, SiPadOnlyOutputSingleReading) {
        alpha.setFourMomentumLab(1E3, {0,0,1});

        auto res = pad.detect(alpha);

        CHECK_EQUAL(1, res.size());
    }

    TEST_FIXTURE(SetupFixture, AlphaWith1MeVIsDetectedInStrip0) {
        alpha.setFourMomentumLab(1E3, {0,0,1});

        auto res = pad.detect(alpha);

        CHECK_EQUAL(0, res[0].channel);
    }

    TEST_FIXTURE(SetupFixture, AlphaWith1MeVIsDeposit1MeV) {
        alpha.setFourMomentumLab(1E3, {0,0,1});

        auto res = pad.detect(alpha);

        CHECK_CLOSE(1E3, res[0].energy, 0.1);
    }

    TEST_FIXTURE(SetupFixture, AlphaAtEdgeIsDetected) {
        alpha.setFourMomentumLab(1E3, {0,0,1});
        alpha.setPosition({19, 24, 0});

        auto res = pad.detect(alpha);

        CHECK_EQUAL(0, res[0].channel);
        CHECK_CLOSE(1E3, res[0].energy, 0.1);
    }
}
