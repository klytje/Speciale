
#include <unittest++/UnitTest++.h>
#include <ausa/eloss/Ion.h>
#include <TVector3.h>

#include "simX/Particle.h"
#include "simX/Detection/Detector.h"
#include "simX/Detection/SegmentedDetector.h"
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

using namespace simX;
using namespace AUSA::EnergyLoss;
using namespace simX::propagator;

SUITE(W1Test) {

    // use SRIM stopping powers for alpha in Silicon
    IonizingPropagator::LossCalc factory(const Layer& layer, const Particle& p) {
        std::string suffix = "/SRIM08/He4_Si.dat";
        SRIMTabulation tab(AUSA::getResourceDirectory() + suffix);
        return std::make_unique<EnergyLossRangeInverter>(std::make_unique<RangeInterpolator>(tab));
    }

    // Factory for creating SRIM13 stopping-power integrators
    MCStragglingPropagator::LossCalc lossFactory(const Layer& layer, const Particle& p) 
    {
        SRIM13Tabulation tabulation{p.getZ(), p.getA(), layer.getMaterial()};
        auto interpolator = std::make_unique<StoppingPowerInterpolator>(tabulation);
        return std::make_unique<EnergyLossIntegrator>(move(interpolator));
    }

    // Factory for creating SRIM13 range calculators
    MCStragglingPropagator::RangeCalc rangeFactory(const Layer& layer, const Particle& p) 
    {
        SRIM13Tabulation tabulation{p.getZ(), p.getA(), layer.getMaterial()};
        return std::make_unique<RangeSplineFitter>(tabulation);
    }

    class SetupFixture {
        public:
            SetupFixture() 
            : proton(Ion(1,1)),
              alpha(Ion(2,4)),
              gamma(Ion(0,0)),
              det("det", 16, 16, TVector3(1.5,1.5,100),TVector3(0,0,-1), TVector3(0,1,0), true, 60E-3, 0.1E-3, 0.4E-3, 0.5E-3, 0.2E-3, 3120E-3, 3120E-3, 100E-3, 100E-3,90E-3), // W1 10 cm downstream of origo
              ionProp( std::make_shared<IonizingPropagator>(IonizingPropagator(factory)) )
            {}
            Particle proton, alpha, gamma;
            W1 det;
            std::shared_ptr<IonizingPropagator> ionProp;
    };


    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(SetupFixture, NumberOfChannels) {
        CHECK_EQUAL( 32, det.getNumberOfChannels() );
    }

    TEST_FIXTURE(SetupFixture, W1IsTwiceSegmented) {
        CHECK_EQUAL( 2, det.nSegmentations() );
    }

    TEST_FIXTURE(SetupFixture, DetectoVolumeHasCorrectPosition) {
        auto& vol = det.getDetectorVolume();
        CHECK_EQUAL( 1.5, vol.getCenter().X() );
        CHECK_EQUAL( 1.5, vol.getCenter().Y() );
        CHECK_CLOSE( 100, vol.getCenter().Z(), 1E-1 );
    }

    TEST_FIXTURE(SetupFixture, DetectOneMeVAlphaInW1) {
        double ekin = 1E3;
        TVector3 dir(0,0,1);
        alpha.setFourMomentumLab( ekin, dir );
        alpha.setPropagator( ionProp );
        auto o = det.detect(alpha);
        CHECK_EQUAL( 2, o.size() );
        CHECK_EQUAL( 9-1, o[0].channel );      // vertical strip no. 9 = channel 8
        CHECK_EQUAL( 16+9-1, o[1].channel );   // horizontal strip no. 9 = channel 24
        CHECK_CLOSE( 969., o[0].energy, 1E0 );
    }
    
    TEST_FIXTURE(SetupFixture, ReverseStripOrderingForW1) {
        double ekin = 1E3;
        TVector3 dir(0,0,1);
        alpha.setFourMomentumLab( ekin, dir );
        alpha.setPropagator( ionProp );
        W1 d0("W1", 16,16, TVector3(1.5,1.5,100),TVector3(0,0,-1), TVector3(0,1,0), true, 60E-3, 0.1E-3, 0.4E-3, 0.5E-3, 0.2E-3, 3120E-3, 3120E-3, 100E-3, 100E-3,90E-3);
        d0.reverseStripOrdering(0, true);
        auto o1 = d0.detect(alpha);
        CHECK_EQUAL( 7, o1[0].channel );
        CHECK_EQUAL( 24, o1[1].channel );
        // reset alpha properties
        alpha.setPosition(TVector3(0,0,0));
        alpha.setFourMomentumLab( ekin, dir );
        d0.reverseStripOrdering(1, true);
        auto o2 = d0.detect(alpha);
        CHECK_EQUAL( 23, o2[1].channel );
    }

    TEST_FIXTURE(SetupFixture, ReverseStripOrderingFrontInTheCtorWorksForW1) {
        double ekin = 1E3;
        TVector3 dir(0,0,1);
        alpha.setFourMomentumLab( ekin, dir );
        alpha.setPropagator( ionProp );
        W1 d0("W1", 16,16, TVector3(1.5,1.5,100),TVector3(0,0,-1), TVector3(0,1,0), true, 60E-3, 0.1E-3, 0.4E-3, 0.5E-3, 0.2E-3, 3120E-3, 3120E-3, 100E-3, 100E-3,90E-3, true);
//        d0.reverseStripOrdering(0, true);
        auto o1 = d0.detect(alpha);
        CHECK_EQUAL( 7, o1[0].channel );
        CHECK_EQUAL( 24, o1[1].channel );
    }

    TEST_FIXTURE(SetupFixture, ReverseStripOrderingBackInTheCtorWorksForW1) {
        double ekin = 1E3;
        TVector3 dir(0,0,1);
        alpha.setFourMomentumLab( ekin, dir );
        alpha.setPropagator( ionProp );
        W1 d0("W1", 16,16, TVector3(1.5,1.5,100),TVector3(0,0,-1), TVector3(0,1,0), true, 60E-3, 0.1E-3, 0.4E-3, 0.5E-3, 0.2E-3, 3120E-3, 3120E-3, 100E-3, 100E-3,90E-3, true, true);
//        d0.reverseStripOrdering(0, true);
        auto o1 = d0.detect(alpha);
        CHECK_EQUAL( 23, o1[1].channel );
    }

    TEST_FIXTURE(SetupFixture, TestSetAndGetMethodsForPosition) {
        det.setPosition(TVector3(1.5,-3.0,100.));
        CHECK_EQUAL( -3.0, det.getPosition().Y() );
    }

    TEST_FIXTURE(SetupFixture, TestSharing) {
        double ekin = 1E3;
        TVector3 dir(0,0,1);
        alpha.setFourMomentumLab( ekin, dir );
        alpha.setPropagator( ionProp );
        det.setPosition(TVector3(1.5,1.E-6,100.)); // DSSD (almost) centered at y=0, so we expect sharing between horizontal strips 8 and 9
        auto o = det.detect(alpha);
        CHECK_EQUAL( 3, o.size() );
        // check that correct strips are firing
        CHECK_EQUAL( 24, o[1].channel );
        CHECK_EQUAL( 23, o[2].channel );
        // check that energy is split in the correct ratio (roughly 50-50)
        CHECK_CLOSE( 0.5*969., o[1].energy, 1. );        
        CHECK_CLOSE( 0.5*969., o[2].energy, 1. );        
    }

    TEST_FIXTURE(SetupFixture, TestDetectOnlyIonizingEnergy) {
		// alpha particle
        double ekin = 1E3;
        TVector3 dir(0,0,1);
        alpha.setFourMomentumLab(ekin, dir);
		// assign propagator
        auto strag = std::make_shared<MCStragglingPropagator> (lossFactory, rangeFactory);
        alpha.setPropagator(strag);
		// detector settings 
        det.setPosition(TVector3(0, 0, 100));
        det.enableSharing(false); // disable sharing
        det.detectOnlyIonizingEnergy(true); // detect only ionizing energy loss
        // detect particle 100 times
        TH1F *h = new TH1F("h_W1_test_01", "detected energy", 1500, 0., 1500.);
		for (size_t i=0; i<100; i++) {
	        alpha.setFourMomentumLab(ekin, dir);
			alpha.setPosition(TVector3(0,0,0));		    
			auto o = det.detect(alpha);
		    auto e = o[0].energy;
			h->Fill(e);
		}
		CHECK_CLOSE(969-8, h->GetMean(), 1);        
    }

    TEST_FIXTURE(SetupFixture, TestSharingCanBeDisabled) {
        double ekin = 1E3;
        TVector3 dir(0,0,1);
        alpha.setFourMomentumLab( ekin, dir );
        alpha.setPropagator( ionProp );
        det.setPosition(TVector3(1.5,1.E-6,100.));
        det.enableSharing(false);
        auto o = det.detect(alpha);
        CHECK_EQUAL( 2, o.size() );
        // check that correct strips is firing
        CHECK_EQUAL( 24, o[1].channel );
        // check that energy is correct
        CHECK_CLOSE( 969., o[1].energy, 1. );        
    }

    TEST_FIXTURE(SetupFixture, TestSharingProbabilityIsAsExpectedForW1) {
        double ekin = 1E3;
        alpha.setPropagator( ionProp );
        det.setPosition(TVector3(0.,0.,500.));
        det.enableSharing(true);
        int n0=1e4;
        int ns[5] = {0,0,0,0,0};
        for (int i=0; i<n0; i++) {
            double costh = 1. - 2.5e-3*rnd();//0.07*rnd(); // 0-21.6 deg
            double sinth = std::sqrt(1.-std::pow(costh,2));
            double phi = 2.*TMath::Pi()*rnd();
            TVector3 dir( sinth*std::cos(phi), sinth*std::sin(phi), costh );
            alpha.setPosition(TVector3(0,0,0));
            alpha.setFourMomentumLab( ekin, dir );
            auto o = det.detect(alpha);
            if (o.size()<=4) ns[o.size()] += 1;
        }
        double hits = ns[1]+ns[2]+ns[3]+ns[4];
        double sharing = ns[3]+ns[4];
        double sharingProb = sharing/hits;
        double threeSigma = 3.*std::sqrt(sharing)/hits;
        CHECK_CLOSE( 0.059, sharingProb, threeSigma );
    }

    TEST_FIXTURE(SetupFixture, TestAlphaDetectedWithCorrectEnergyWhenEnteringFromBackSideOfW1) {
        double ekin = 1E3;
        W1 d0("W1", 16,16, TVector3(0.,0.,100), TVector3(0,0,1), TVector3(0,1,0), true, 60E-3, 0.1E-3, 0.4E-3, 0.5E-3, 0.2E-3, 3120E-3, 3120E-3, 100E-3, 100E-3,90E-3 ); // backside facing source
        alpha.setPropagator( ionProp );
        d0.enableSharing(false);
        TVector3 dir(0,0,1);
        alpha.setPosition(TVector3(0,0,0));
        alpha.setFourMomentumLab( ekin, dir );
        auto o = d0.detect(alpha);
        CHECK_CLOSE( 809., o[0].energy, 1. );
    }
    TEST_FIXTURE(SetupFixture, TestAlphaHasZeroEnergyAfterBeingStopped) {
        double ekin = 1E3;
        alpha.setPropagator( ionProp );
        TVector3 dir(0,0,1);
        alpha.setPosition(TVector3(0,0,0));
        alpha.setFourMomentumLab( ekin, dir );
        det.detect(alpha);
        CHECK_CLOSE( 0., alpha.getKineticEnergyLab(), 1E-9 );
    }

    TEST_FIXTURE(SetupFixture, TestW1CanHandleMultipleParticles) {
        // alpha particle
        Particle alpha1(Ion(2,4));
        double ekin1 = 1E3;
        alpha1.setPropagator( ionProp );
        double th1 = 10.*TMath::Pi()/180.;
        double phi1 = 0.;
        TVector3 dir1( std::sin(th1)*std::cos(phi1), std::sin(th1)*std::sin(phi1), std::cos(th1) );
        alpha1.setPosition(TVector3(0,0,0));
        alpha1.setFourMomentumLab( ekin1, dir1 );
        // another alpha particle
        Particle alpha2(Ion(2,4));
        double ekin2 = 2E3;
        alpha2.setPropagator( ionProp );
        double th2 = 10.*TMath::Pi()/180.;
        double phi2 = TMath::Pi();
        TVector3 dir2( std::sin(th2)*std::cos(phi2), std::sin(th2)*std::sin(phi2), std::cos(th2) );
        alpha2.setPosition(TVector3(0,0,0));
        alpha2.setFourMomentumLab( ekin2, dir2 );
        // detector 
        W1 d0("W1", 16,16 ,TVector3(0.,0.,100), TVector3(0,0,-1), TVector3(0,1,0), true, 60E-3, 0.1E-3, 0.4E-3, 0.5E-3, 0.2E-3, 3120E-3, 3120E-3, 100E-3, 100E-3,90E-3 );
        d0.enableSharing(false);
        // detect particles
        auto o1 = d0.detect(alpha1);
        auto o2 = d0.detect(alpha2);
        int channels = o1.size() + o2.size();
        CHECK_EQUAL( channels, 4 );
        CHECK_EQUAL( 2, o1[0].channel);
        CHECK_EQUAL( 23, o1[1].channel);
        CHECK_EQUAL( 13, o2[0].channel);
        CHECK_EQUAL( 23, o2[1].channel);
    }

    TEST_FIXTURE(SetupFixture, TestThatUpperHorizontalStripIsStrip1) {
        auto p = det.getPosition() + det.getUp()*7.5*3120E-3;
        p-= TVector3(0,0,100);

        alpha.setPosition(p);
        alpha.setFourMomentumLab(1000, {0,0,1});

        alpha.setPropagator(ionProp);

        CHECK(det.getDetectorVolume().isIntersecting(alpha.getPosition(), alpha.getDirectionLab()));

        det.enableSharing(false);
        auto d = det.detect(alpha);

        CHECK_EQUAL(2, d.size());

        CHECK_EQUAL( 9-1, d[0].channel );
        CHECK_EQUAL( 16+1-1, d[1].channel );
    }

    TEST_FIXTURE(SetupFixture, TestThatLowerHorizontalStripIsStrip16) {
        auto p = det.getPosition() - det.getUp()*7.5*3120E-3;
        p-= TVector3(0,0,100);

        alpha.setPosition(p);
        alpha.setFourMomentumLab(1000, {0,0,1});

        alpha.setPropagator(ionProp);

        CHECK(det.getDetectorVolume().isIntersecting(alpha.getPosition(), alpha.getDirectionLab()));

        det.enableSharing(false);
        auto d = det.detect(alpha);

        CHECK_EQUAL(2, d.size());

        CHECK_EQUAL( 9-1, d[0].channel );
        CHECK_EQUAL( 16+16-1, d[1].channel );
    }
}
