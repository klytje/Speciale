
#include <unittest++/UnitTest++.h>

#include <ausa/eloss/Ion.h>
#include <ausa/constants/Mass.h>
#include <simX/propagator/IonizingPropagator.h>
#include <simX/propagator/GaussianStragglingPropagator.h>

#include "simX/Particle.h"
#include "simX/Layer.h"

#include <ausa/eloss/SRIMTabulation.h>
#include <ausa/eloss/RangeInterpolator.h>
#include <ausa/eloss/EnergyLossRangeInverter.h>
#include <ausa/util/Resource.h>
#include <ausa/util/memory>
#include <ausa/geometry/Box.h>

#include <TH1.h>

using namespace std;
using namespace simX;
using namespace simX::propagator;
using namespace AUSA::EnergyLoss;

using AUSA::Geometry::Box;

SUITE(StragglingTest) {

    IonizingPropagator::LossCalc factory(const Layer& layer, const Particle& p) {
        std::string suffix = "/SRIM08/";
        if (layer.getMaterial().getName() == "Silicon")
            suffix += "He4_Si.dat";
        else
            suffix += "He4_Fe.dat";

        SRIMTabulation tab(AUSA::getResourceDirectory() + suffix);

        return std::make_unique<EnergyLossRangeInverter>(std::make_unique<RangeInterpolator>(tab));
    }
    
    class SetupFixture {
    public:
        SetupFixture()
                : silicon("Silicon", 2.3290, 14, AUSA::Constants::atomicMass(14)),
                  iron("Iron", 7.874, 26, AUSA::Constants::atomicMass(26)),
                  b0(silicon,std::make_unique<Box>(100.,100.,1.E-3,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false),
                  b1(iron,std::make_unique<Box>(100.,100.,1.E-3,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false),
                  p(Ion(2,4))
        {
            prop = std::make_shared<IonizingPropagator> (factory);
            strag = std::make_unique<GaussianStragglingPropagator> (move(prop));
        }
        Material silicon, iron;
        Layer b0, b1;
        Particle p;
        std::shared_ptr<IonizingPropagator> prop;
        std::unique_ptr<GaussianStragglingPropagator> strag;
    };

    TEST_FIXTURE(SetupFixture, Sanity) {
        factory(b0, p); // Check this does not crash

        CHECK_EQUAL(1,1);
    }
    
    TEST_FIXTURE(SetupFixture, CheckDeflectionFor1MeVAlphaIn1MicronSi) {
        double z0 = -50;
        p.setPosition({0, 0, z0});
        p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1000keV, dir=(0,0,1)

        strag -> propagate(b0, p);

        // check half-scattering angle of 1 MeV alpha in 1 um Si
        auto a = strag -> getHalfScatteringAngle();
        CHECK_CLOSE(2.54086, a, 1E-2);     

        // check deflection
        auto& dir = p.getDirectionLab();
        double theta = dir.Theta() * 180. / TMath::Pi();
        CHECK_EQUAL(true, theta < 2.1385 * 3.); // sampled angle is less that 3 times half-scattering angle        
    }

    TEST_FIXTURE(SetupFixture, CheckDeflectionForIonTravellingInXDirection) {
        double x0 = -50;
        p.setPosition({x0, 0, 0});
        p.setFourMomentumLab(1E3, {1,0,0}); // Ekin=1000keV, dir=(1,0,0)

        Layer b(silicon,std::make_unique<Box>(1.E-3,100.,100.,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)),false);
        strag -> propagate(b, p);

        // check deflection
        auto& dir = p.getDirectionLab();
        double theta = TMath::ACos(dir.Dot(TVector3(1,0,0))) * 180. / TMath::Pi();
        CHECK_EQUAL(true, theta < 2.1385 * 3.); // sampled angle is less that 3 times half-scattering angle        
    }

    TEST_FIXTURE(SetupFixture, CheckDeflectionFor2MeVAlphaIn1MicronFe) {
        double z0 = -50;
        p.setPosition({0, 0, z0});
        p.setFourMomentumLab(2E3, {0,0,1}); // Ekin=2000keV, dir=(0,0,1)

        strag -> propagate(b1, p);

        // check half-scattering angle of 1 MeV alpha in 1 um Si
        auto a = strag -> getHalfScatteringAngle();
        CHECK_CLOSE(3.11305, a, 1E-2);     
    }

    TEST_FIXTURE(SetupFixture, CheckDeflectionForFullyStoppedIon) {
        double z0 = -50;
        p.setPosition({0, 0, z0});
        p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1000keV, dir=(0,0,1)

        // 1 mm silicon
        Layer b(silicon,std::make_unique<Box>(100.,100.,1.,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)),false);
        
        strag -> propagate(b, p);

        // check that deflection is zero
        auto& dir = p.getDirectionLab();
        double theta = dir.Theta() * 180. / TMath::Pi();
        CHECK_EQUAL(0., theta);       
    }

    TEST_FIXTURE(SetupFixture, CheckMultipleElementsThickLayer) {
        double z0 = -50;
        int N = 10000;
        TH1F *h1 = new TH1F("h_test_01", "h1 title", 100, 0., 10.);

        Material fake("fake", 7.874, {26, 26}, {1, 2});

        // iron
        Layer bIron(iron,std::make_unique<Box>(100.,100.,1.E-3,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false);
        h1 -> Reset();
        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition({0, 0, z0});
            p.setFourMomentumLab(2E3, {0,0,1}); // Ekin=2000keV, dir=(0,0,1)
            strag -> propagate(bIron, p); // 1 micron Fe
            auto& dir = p.getDirectionLab();
            double theta = dir.Theta() * 180. / TMath::Pi();
            h1 -> Fill(theta);
        }
        double meanIron = h1->GetMean();

        // fake
        Layer bFake(fake,std::make_unique<Box>(100.,100.,1.E-3,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false);
        h1 -> Reset();
        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition({0, 0, z0});
            p.setFourMomentumLab(2E3, {0,0,1}); // Ekin=2000keV, dir=(0,0,1)
            strag -> propagate(bFake, p); // 1 micron Fe
            auto& dir = p.getDirectionLab();
            double theta = dir.Theta() * 180. / TMath::Pi();
            h1 -> Fill(theta);
        }
        double meanFake = h1->GetMean();

        double dev = (meanFake - meanIron) / meanIron;
        CHECK_EQUAL(true, abs(dev) < 5E-2); // deviate by less than 5%
   }

    TEST_FIXTURE(SetupFixture, CheckMultipleElementsThinLayer) {
        double z0 = -50;
        int N = 10000;
        TH1F *h1 = new TH1F("h_test_02", "h1 title", 100, 0., 10.);

        Material fake("fake", 7.874, {26, 26}, {1, 2});

        // iron
        Layer bIron(iron,std::make_unique<Box>(100.,100.,1.E-4,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false);
        h1 -> Reset();
        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition({0, 0, z0});
            p.setFourMomentumLab(2E3, {0,0,1}); // Ekin=2000keV, dir=(0,0,1)
            strag -> propagate(bIron, p); // 0.1 micron Fe
            auto& dir = p.getDirectionLab();
            double theta = dir.Theta() * 180. / TMath::Pi();
            h1 -> Fill(theta);
        }
        double meanIron = h1->GetMean();

        // fake
        Layer bFake(fake,std::make_unique<Box>(100.,100.,1.E-4,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false);
        h1 -> Reset();
        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition({0, 0, z0});
            p.setFourMomentumLab(2E3, {0,0,1}); // Ekin=2000keV, dir=(0,0,1)
            strag -> propagate(bFake, p); // 0.1 micron Fe
            auto& dir = p.getDirectionLab();
            double theta = dir.Theta() * 180. / TMath::Pi();
            h1 -> Fill(theta);
        }
        double meanFake = h1->GetMean();

        double dev = (meanFake - meanIron) / meanIron;
        CHECK_EQUAL(true, abs(dev) < 5E-2); // deviate by less than 5%
   }
   
   TEST_FIXTURE(SetupFixture, CheckEnergyStragglingForPenetratingParticle) {
        auto enProp  = std::make_shared<IonizingPropagator> (factory);
        auto enStrag = std::make_unique<GaussianStragglingPropagator> (move(enProp), false, true);
        
        double z0 = -50;
        int N = 10000;
        TH1F *h1 = new TH1F("h_test_03", "h1 title", 100, 0., 10.);
        TH1F *h2 = new TH1F("h_test_04", "h1 title", 3000, 0., 3000.);

        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition({0, 0, z0});
            p.setFourMomentumLab(2E3, {0,0,1}); // Ekin=2000keV, dir=(0,0,1)
            enStrag -> propagate(b1, p); // 1 micron Fe
            auto& dir = p.getDirectionLab();
            double theta = dir.Theta() * 180. / TMath::Pi();
            h1 -> Fill(theta);
            double energy = p.getKineticEnergyLab();
            h2 -> Fill(energy);
        }
        // check that deflection is zero
        double meanTheta = h1->GetMean();
        CHECK_EQUAL(0., meanTheta);     
        
        // check energy straggling 
        double stdE  = h2->GetStdDev();
        CHECK_CLOSE(14.02, stdE, 1.0);        
   }

   TEST_FIXTURE(SetupFixture, CheckTransmittedParticlesHaveNonNegativeEnergy) {
        auto enProp  = std::make_shared<IonizingPropagator> (factory);
        auto enStrag = std::make_unique<GaussianStragglingPropagator> (move(enProp), false, true);

        // 3.3 micron iron
        Layer b3um(iron,std::make_unique<Box>(100.,100.,3.1E-3,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false);
        
        double z0 = -50;
        int N = 10000;
        TH1F *h2 = new TH1F("h_test_05", "h1 title", 3000, 0., 3000.);

        double emin = 2000;
        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition({0, 0, z0});
            p.setFourMomentumLab(2E3, {0,0,1}); // Ekin=2000keV, dir=(0,0,1)
            enStrag -> propagate(b3um, p); // 3.3 micron Fe
            double energy = p.getKineticEnergyLab();
            h2 -> Fill(energy);
            if (energy < emin) emin = energy;
        }

        // check min energy 
        CHECK_EQUAL(true, emin >= 0);        
   }
}
