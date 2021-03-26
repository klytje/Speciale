
#include <unittest++/UnitTest++.h>

#include <ausa/eloss/Ion.h>
#include <ausa/constants/Mass.h>
#include <simX/propagator/IonizingPropagator.h>
#include <simX/propagator/MCStragglingPropagator.h>

#include "simX/Particle.h"
#include "simX/Layer.h"

#include "ausa/eloss/Material.h"
#include "ausa/eloss/EnergyLossIntegrator.h"
#include "ausa/eloss/SRIM13Tabulation.h"
#include "ausa/eloss/StoppingPowerInterpolator.h"
#include "ausa/eloss/RangeSplineFitter.h"
#include <ausa/util/Resource.h>
#include <ausa/util/memory>
#include <ausa/geometry/Box.h>

#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>

using namespace std;
using namespace simX;
using namespace simX::propagator;
using namespace AUSA::EnergyLoss;

using AUSA::Geometry::Box;

SUITE(MCStragglingTest) {

    MCStragglingPropagator::LossCalc lossFactory(const Layer& layer, const Particle& p) 
    {
        SRIM13Tabulation tabulation{p.getZ(), p.getA(), layer.getMaterial()};
        auto interpolator = make_unique<StoppingPowerInterpolator>(tabulation);
        return make_unique<EnergyLossIntegrator>(move(interpolator));
    }

    MCStragglingPropagator::RangeCalc rangeFactory(const Layer& layer, const Particle& p) 
    {
        SRIM13Tabulation tabulation{p.getZ(), p.getA(), layer.getMaterial()};
        return make_unique<RangeSplineFitter>(tabulation);
    }
    
    class SetupFixture {
    public:
        SetupFixture()
                : silicon("Silicon", 2.33, 14, AUSA::Constants::atomicMass(14)),
                  b(silicon,std::make_unique<Box>(100.,100.,1E-3,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false), // 1 um Si
                  p(Ion(2,4))
        {
            strag = std::make_unique<MCStragglingPropagator> (lossFactory, rangeFactory);
            strag -> setCrossSectionMultiplier(0.04);
        }
        Material silicon;
        Layer b;
        Particle p;
        std::unique_ptr<MCStragglingPropagator> strag;
    };

    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1,1);
    }
    
    TEST_FIXTURE(SetupFixture, CheckDeflectionFor1MeVAlphaIn1MicronSi) {
        double z0 = -50;
        int N = 1E3;

        // determine intersection point
        p.setPosition({0, 0, z0});
        p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1MeV, dir=(0,0,1)
        TVector3 intersect;
        double dist, thick;
        b.getIntersection(p.getPosition(), p.getDirectionLab(), intersect, thick, dist);

        TH1F *h1 = new TH1F("h_MC_test_01", "h1 title", 250, 0., 50.);
        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition(intersect);
            p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1MeV, dir=(0,0,1)
            strag -> propagate(b, p); // 1 micron Si
            auto& dir = p.getDirectionLab();
            double theta = dir.Theta() * 180. / TMath::Pi();
            h1 -> Fill(theta);
        }        
        // check deflection
        double defl = h1->GetMean();
        CHECK_CLOSE(2.70, defl, 0.15); 
    }

    TEST_FIXTURE(SetupFixture, CheckIonizingEnergyLossFor1MeVAlphaIn100MicronSi) {
        double z0 = -50;
        int N = 1E3;
        
        Layer b0(silicon,std::make_unique<Box>(100.,100.,100E-3,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false); // 100 um Si

        // determine intersection point
        p.setPosition({0, 0, z0});
        p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1MeV, dir=(0,0,1)
        TVector3 intersect;
        double dist, thick;
        b0.getIntersection(p.getPosition(), p.getDirectionLab(), intersect, thick, dist);

        TH1F *h2 = new TH1F("h_MC_test_02", "dE_ioni", 1500, 0., 1500.);
        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition(intersect);
            p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1MeV, dir=(0,0,1)
            strag -> propagate(b0, p); // 100 micron Si
            h2 -> Fill(strag->getIonizingEnergyLoss());
        }
        
        // check ionizing energy loss
        double ioni = h2->GetMean();
        CHECK_CLOSE(991., ioni, 2.); 
    }

    TEST_FIXTURE(SetupFixture, CheckMultipleElements) {
        double z0 = -50;
        int N = 1E3;

        Material mostlySi("mostlySi", 2.33, {26, 14}, {1, 9});

        // iron
        Layer b0(mostlySi,std::make_unique<Box>(100.,100.,100E-3,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false);

        // determine intersection point
        p.setPosition({0, 0, z0});
        p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1MeV, dir=(0,0,1)
        TVector3 intersect;
        double dist, thick;
        b0.getIntersection(p.getPosition(), p.getDirectionLab(), intersect, thick, dist);

        TH1F *h3 = new TH1F("h_MC_test_03", "dE_ioni", 1500, 0., 1500.);
        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition(intersect);
            p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1MeV, dir=(0,0,1)
            strag -> propagate(b0, p); // 100 micron Si
            h3 -> Fill(strag->getIonizingEnergyLoss());
        }
        
        // check ionizing energy loss
        double ioni = h3->GetMean();
        CHECK_CLOSE(991., ioni, 2.); 
    }

    TEST_FIXTURE(SetupFixture, CheckNonIonizingEnergyLossPassedByPointer) {
        double z0 = -50;
        int N = 1E3;
        
        Layer b0(silicon,std::make_unique<Box>(100.,100.,100E-3,TVector3(0,0,0),TVector3(0,0,-1),TVector3(0,1,0)), false); // 100 um Si

        // determine intersection point
        p.setPosition({0, 0, z0});
        p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1MeV, dir=(0,0,1)
        TVector3 intersect;
        double dist, thick;
        b0.getIntersection(p.getPosition(), p.getDirectionLab(), intersect, thick, dist);

        TH1F *h4 = new TH1F("h_MC_test_04", "dE_ioni", 1500, 0., 1500.);
        for (size_t i=0; i<N; i++) { // repeat N times 
            p.setPosition(intersect);
            p.setFourMomentumLab(1E3, {0,0,1}); // Ekin=1MeV, dir=(0,0,1)
            double * en = new double(0);
            strag -> propagate(b0, p, -1, en); // 100 micron Si
            h4 -> Fill(*en);
        }
        
        // check non-ionizing energy loss
        double en = h4->GetMean();
        CHECK_CLOSE(9., en, 1.); 
    }
}
