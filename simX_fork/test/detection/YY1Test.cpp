
#include <unittest++/UnitTest++.h>
#include <ausa/eloss/Ion.h>
#include <TVector3.h>

#include "simX/Particle.h"
#include "simX/Detection/Detector.h"
#include "simX/Detection/SegmentedDetector.h"
#include "simX/Detection/YY1.h"
#include <simX/Detection/W1.h>
#include "simX/Layer.h"
#include "simX/Random.h"
#include <simX/propagator/ParticlePropagator.h>
#include <simX/propagator/IonizingPropagator.h>
#include <simX/propagator/MCStragglingPropagator.h>

#include <ausa/eloss/SRIMTabulation.h>
#include "ausa/eloss/SRIM13Tabulation.h"
#include <ausa/eloss/RangeInterpolator.h>
#include <ausa/eloss/EnergyLossRangeInverter.h>
#include "ausa/eloss/StoppingPowerInterpolator.h"
#include "ausa/eloss/RangeSplineFitter.h"
#include "ausa/eloss/EnergyLossIntegrator.h"
#include <ausa/util/Resource.h>
#include <ausa/util/memory>

#include <math.h>
#include <TMath.h>
#include <tuple>
#include <iostream>
#include <simX/propagator/NonIonizingPropagator.h>

using namespace simX;
using namespace simX::Detection;
using namespace AUSA::EnergyLoss;
using namespace simX::propagator;

SUITE(YY1Test) {

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
              det("det", TVector3(0, 90, 10), TVector3(0,0,-1), TVector3(0,1,0), 1, 0, 0),
              ionProp( std::make_shared<IonizingPropagator>(IonizingPropagator(factory)) )
            {
                alpha.setPropagator(ionProp);
                proton.setPropagator(ionProp);

                alpha.setPosition({0, 90, 0}); // Default to center
                proton.setPosition({0, 90, 0}); // Default to center
            }
            Particle proton, alpha, gamma;
            YY1 det;
            std::shared_ptr<IonizingPropagator> ionProp;
    };


    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(SetupFixture, NumberOfChannels) {
        CHECK_EQUAL( 16, det.getNumberOfChannels() );
    }

    TEST_FIXTURE(SetupFixture, YY1IsOnceSegmented) {
        CHECK_EQUAL( 1, det.nSegmentations() );
    }

    TEST_FIXTURE(SetupFixture, OneStripFiresWhenHitByAlpha) {
        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        auto o = det.detect(alpha);
        CHECK_EQUAL( 1, o.size() );
    }

    TEST_FIXTURE(SetupFixture, AlphaWith1MeVDepositsAllEnergy) {
        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        auto o = det.detect(alpha);
        CHECK_CLOSE(1E3, o[0].energy, 1E-3);
    }

    TEST_FIXTURE(SetupFixture, AlphaPropagtingInDetectorWillStopAtZ11) {
        alpha.setPropagator(std::make_shared<NonIonizingPropagator>());

        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        auto o = det.detect(alpha);

        CHECK_ARRAY_CLOSE(TVector3(0, 90, 11), alpha.getPosition(), 3, 1E-3);
    }

    TEST_FIXTURE(SetupFixture, AlphaAtLowerEdgeHitsStrip1) {

        alpha.setPosition({0, 50, 0});
        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        auto o = det.detect(alpha);
        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(0, o[0].channel);
    }

    TEST_FIXTURE(SetupFixture, AlphaAtFirstStripBoundaryHitsStrip2) {
        alpha.setPosition({0, 55, 0});
        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        auto o = det.detect(alpha);
        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(1, o[0].channel);
    }

    TEST_FIXTURE(SetupFixture, AlphaAtSecondStripBoundaryHitsStrip3) {
        alpha.setPosition({0, 60, 0});
        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        auto o = det.detect(alpha);
        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(2, o[0].channel);
    }

    TEST_FIXTURE(SetupFixture, AlphaAtUpperBoundaryHitsStrip16) {
        alpha.setPosition({0, 130, 0});
        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        auto o = det.detect(alpha);
        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(15, o[0].channel);
    }


    TEST_FIXTURE(SetupFixture, Strip13Is40DegreesWide) {

        auto i = 16-13+0.5;
        auto r = 130-i*5;
        auto deg = 20;
        auto x = r*sin(deg*TMath::DegToRad());
        auto y = r*cos(deg*TMath::DegToRad());


        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        alpha.setPosition({x, y, 0});
        auto o = det.detect(alpha);

        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(12, o[0].channel);

        alpha.setPosition({-x, y, 0});
        auto o2 = det.detect(alpha);

        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(12, o2[0].channel);
    }


    TEST_FIXTURE(SetupFixture, Strip14Is35DegreesWide) {

        auto i = 16-14+0.5;
        auto r = 130-i*5;
        auto deg = 24/2.;
        auto x = r*sin(deg*TMath::DegToRad());
        auto y = r*cos(deg*TMath::DegToRad());


        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        alpha.setPosition({x, y, 0});
        auto o = det.detect(alpha);

        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(13, o[0].channel);

        alpha.setPosition({-x, y, 0});
        auto o2 = det.detect(alpha);

        CHECK_EQUAL(1, o2.size());
        CHECK_EQUAL(13, o2[0].channel);
    }

    TEST_FIXTURE(SetupFixture, Outside36DegreesStrip14GivesNothing) {

        auto i = 16-14+0.5;
        auto r = 130-i*5;
        auto deg = 36./2;
        auto x = r*sin(deg*TMath::DegToRad());
        auto y = r*cos(deg*TMath::DegToRad());


        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        alpha.setPosition({x, y, 0});

        auto o = det.detect(alpha);
        CHECK_EQUAL(0, o.size());
        auto o2 = det.detect(alpha);
        CHECK_EQUAL(0, o2.size());
    }

    TEST_FIXTURE(SetupFixture, Outside29DegreesStrip2GivesNothing) {

        auto r = 130-1.5*5;
        auto deg = 29./2;
        auto x = r*sin(deg*TMath::DegToRad());
        auto y = r*cos(deg*TMath::DegToRad());


        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        alpha.setPosition({x, y, 0});

        auto o = det.detect(alpha);
        CHECK_EQUAL(0, o.size());
        auto o2 = det.detect(alpha);
        CHECK_EQUAL(0, o2.size());
    }

    TEST_FIXTURE(SetupFixture, Outside19DegreesStrip1GivesNothing) {

        auto r = 130-0.5*5;
        auto deg = 19./2;
        auto x = r*sin(deg*TMath::DegToRad());
        auto y = r*cos(deg*TMath::DegToRad());


        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        alpha.setPosition({x, y, 0});

        auto o = det.detect(alpha);
        CHECK_EQUAL(0, o.size());
        auto o2 = det.detect(alpha);
        CHECK_EQUAL(0, o2.size());
    }

    TEST_FIXTURE(SetupFixture, IfReversedLowerEdgeIsStrip16) {

        det.reverseStripOrdering(0, true);

        alpha.setPosition({0, 50, 0});
        alpha.setFourMomentumLab( 1E3, {0,0,1} );
        auto o = det.detect(alpha);
        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(15, o[0].channel);
    }
}
