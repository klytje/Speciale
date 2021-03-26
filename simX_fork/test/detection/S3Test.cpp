
#include <unittest++/UnitTest++.h>
#include <ausa/eloss/Ion.h>
#include <TVector3.h>

#include "simX/Particle.h"
#include "simX/Detection/Detector.h"
#include "simX/Detection/SegmentedDetector.h"
#include <simX/Detection/S3.h>
#include "simX/Layer.h"
#include "simX/Random.h"
#include <simX/propagator/ParticlePropagator.h>
#include <simX/propagator/IonizingPropagator.h>

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
using namespace simX::propagator;

SUITE(S3Test) {

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
              ionProp( std::make_shared<IonizingPropagator>(IonizingPropagator(factory)) )
            {
                alpha.setPropagator(ionProp);
            }
            Particle proton, alpha, gamma;
            std::shared_ptr<IonizingPropagator> ionProp;
    };


    TEST_FIXTURE(SetupFixture, S3IsTwiceSegmented) {

        S3 det( "S3", TVector3(0.,0.,100), TVector3(0,0,-1), TVector3(0, 1, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, false, false );
        CHECK_EQUAL( 2, det.nSegmentations() );
    }

    TEST_FIXTURE(SetupFixture, DetectAlphaInS3) {
        // detector 
        S3 d0( "S3", TVector3(0.,0.,100), TVector3(0,0,-1), TVector3(0, 1, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, false, false );
        d0.enableSharing(false);
        // 5 MeV alpha
        double ekin = 5E3;
        double th = 10.*TMath::Pi()/180.;        // theta=10  (=> ring #8)
        double phi = 93.*TMath::Pi()/180.;       // phi=93    (=> spoke #8)
        TVector3 dir( std::sin(th)*std::cos(phi), std::sin(th)*std::sin(phi), std::cos(th) );
        alpha.setFourMomentumLab( ekin, dir );
        alpha.setPropagator( ionProp );
        // detect
        auto o = d0.detect(alpha);
        CHECK_EQUAL( 2, o.size() );
        CHECK_EQUAL( 0, o[0].channel );   // spoke #1 = channel 0
        CHECK_EQUAL( 31+8, o[1].channel );   // ring #8 = channel 31+8=39
        CHECK_CLOSE( 4913., o[0].energy, 1E0 );
    }

    TEST_FIXTURE(SetupFixture, TestParticleCanPassThroughCenterOfS3) {
        // detector 
        S3 d0( "S3", TVector3(0.,0.,100), TVector3(0, 0, -1), TVector3(0, 1, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, false, false);
        d0.enableSharing(false);
        // 5 MeV alpha
        double ekin = 5E3;
        TVector3 dir( 0,0,1 );
        alpha.setFourMomentumLab( ekin, dir );
        alpha.setPropagator( ionProp );
        // detect
        CHECK_EQUAL( ekin, alpha.getKineticEnergyLab() );
        auto o = d0.detect(alpha);
        CHECK_EQUAL( 0, o.size() );
        CHECK_EQUAL( ekin, alpha.getKineticEnergyLab() );
    }


    TEST_FIXTURE(SetupFixture, HitRightSpoke) {
        S3 d0("S3", TVector3(0.,0.,40), TVector3(0,0,-1), TVector3(0, 1, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, false, false);
        d0.enableSharing(false);
        // 5 MeV alpha
        double ekin = 5E3;

        TVector3 dir = TVector3(0, 0, 40) + (TVector3(1, 0.00000001, 0).Unit() * (11 + (24. - 0.5) * 0.886 / 2.));

        alpha.setFourMomentumLab( ekin, dir.Unit() );
        alpha.setPropagator( ionProp );
        // detect
        auto o = d0.detect(alpha);

        CHECK_EQUAL( 2, o.size() );
        CHECK_EQUAL( 8, o[0].channel );
        CHECK_EQUAL( 31+12, o[1].channel );
    }


    TEST_FIXTURE(SetupFixture, BugReport38_TwoParticlesGives3HitsInFrontAnd1InBack) {
        // detector
        S3 d0( "S3", TVector3(0.,0.,42), TVector3(0,0,-1), TVector3(-1, 0, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, false, false );
        d0.enableSharing(false);

        alpha.setPosition({-8.87275, -13.3876, 42});
        alpha.setFourMomentumLab(2510.08, {-0.197321, -0.297725, 0.934037});

        auto det0 = d0.detect(alpha);
        CHECK_EQUAL(2, det0.size());

        alpha.setPosition({0.484187, -5.58057, 42});
        alpha.setFourMomentumLab(2271.86, {0.011427, -0.131704, 0.991223});

        auto det1 = d0.detect(alpha);
        CHECK_EQUAL(0, det1.size());
    }

    TEST_FIXTURE(SetupFixture, TestThatInnerRingIs1) {
        // detector
        S3 d0("S3", TVector3(0.,0.,40), TVector3(0,0,-1), TVector3(0, 1, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, false, false);
        d0.enableSharing(false);
        // 5 MeV alpha
        double ekin = 5E3;

        TVector3 dir = TVector3(0, 0, 40) + (TVector3(0, 1, 0).Unit() * (11.0 + 0.886 *0.5));

        alpha.setFourMomentumLab( ekin, dir.Unit() );
        alpha.setPropagator( ionProp );

        auto o = d0.detect(alpha);
        CHECK_EQUAL( 0, o[0].channel );
        CHECK_EQUAL( 31+1, o[1].channel );
    }

    TEST_FIXTURE(SetupFixture, TestThatOuterRingIs24) {
        // detector
        S3 d0("S3", TVector3(0.,0.,40), TVector3(0,0,-1), TVector3(0, 1, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, false, false);
        d0.enableSharing(false);
        // 5 MeV alpha
        double ekin = 5E3;

        TVector3 dir = TVector3(0, 0, 40) + (TVector3(0, 1, 0).Unit() * (11. + (24. - 0.5) * 0.886));

        alpha.setFourMomentumLab( ekin, dir.Unit() );
        alpha.setPropagator( ionProp );

        auto o = d0.detect(alpha);
        CHECK_EQUAL( 0, o[0].channel );
        CHECK_EQUAL( 31+24, o[1].channel );
    }


    TEST_FIXTURE(SetupFixture, TestThatSharingHappensBetweenSpokes) {
        // detector
        S3 d0( "S3", TVector3(0.,0.,100), TVector3(0,0,-1), TVector3(0, 1, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, false, false );
        d0.enableSharing(true);
        // 5 MeV alpha
        double ekin = 5E3;

        double th = 10.*TMath::Pi()/180.;
        double phi = 2.*TMath::Pi()/32./2.;
        TVector3 dir( std::sin(th)*std::cos(phi), std::sin(th)*std::sin(phi), std::cos(th) );

        alpha.setFourMomentumLab( ekin, dir.Unit() );
        alpha.setPropagator( ionProp );

        auto o = d0.detect(alpha);
        CHECK_EQUAL( 3, o.size());

        CHECK_EQUAL( 8, o[0].channel );
        CHECK_EQUAL( 7, o[1].channel );
        CHECK_EQUAL( 39, o[2].channel );
    }

    TEST_FIXTURE(SetupFixture, IfRingsAreReversedThenOuterRingIs1) {
        // detector
        S3 d0("S3", TVector3(0.,0.,40), TVector3(0,0,-1), TVector3(0, 1, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, false, true);
        d0.enableSharing(false);
        // 5 MeV alpha
        double ekin = 5E3;

        TVector3 dir = TVector3(0, 0, 40) + (TVector3(0, 1, 0).Unit() * (11. + (24. - 0.5) * 0.886));

        alpha.setFourMomentumLab( ekin, dir.Unit() );
        alpha.setPropagator( ionProp );

        auto o = d0.detect(alpha);
        CHECK_EQUAL( 0, o[0].channel );
        CHECK_EQUAL( 31+1, o[1].channel );
    }

    TEST_FIXTURE(SetupFixture, IfSpokesAreReversedThenFirstSpokesOuterRingIs32) {
        // detector
        S3 d0("S3", TVector3(0.,0.,40), TVector3(0,0,-1), TVector3(0, 1, 0), 32, 24, 11, 886E-3, 100E-3, 100E-3, true, 300E-3, 0.6E-3, 1.0E-3, true, false);
        d0.enableSharing(false);
        // 5 MeV alpha
        double ekin = 5E3;

        TVector3 dir = TVector3(0, 0, 40) + (TVector3(0, 1, 0).Unit() * (11. + (24. - 0.5) * 0.886));

        alpha.setFourMomentumLab( ekin, dir.Unit() );
        alpha.setPropagator( ionProp );

        auto o = d0.detect(alpha);
        CHECK_EQUAL( 31, o[0].channel );
        CHECK_EQUAL( 31+24, o[1].channel );
    }
}
