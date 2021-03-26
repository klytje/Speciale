
#include <unittest++/UnitTest++.h>

#include "simX/Beam.h"
#include <simX/Particle.h>
#include "ausa/constants/Mass.h"
#include <ausa/eloss/Material.h>
#include <ausa/AUSA.h>
#include <ausa/util/memory>

#include "TF1.h"
#include "Math/WrappedMultiTF1.h"


using namespace std;
using namespace simX;
using namespace AUSA::EnergyLoss;

SUITE(BeamTest) {

    class SetupFixture {
        public:
            SetupFixture() 
            : prot(Ion(1,1))
            {
            }
            Particle prot;
    };

    TEST_FIXTURE(SetupFixture, SimpleBeamReturnsCorrectEnergy) {
        Beam simpleBeam(1500.,0.5,-10, -0.5,0., 0);

        double e = simpleBeam.sampleEnergy(1);
        CHECK_EQUAL(1500.,e);
    }

    TEST_FIXTURE(SetupFixture, SimpleBeamReturnsCorrectNomialEnergyIfNOTPerNucleon) {
        Beam simpleBeam(1500.,0.5,-10, -0.5,0., 0);

        double e = simpleBeam.getNominalEnergy(10);
        CHECK_EQUAL(1500.,e);
    }

    TEST_FIXTURE(SetupFixture, SimpleBeamReturnsCorrectNomialEnergyPerNucleon) {
        Beam simpleBeam(1500.,0.5,-10, -0.5,0., 0, nullptr, nullptr, nullptr, nullptr, true);

        double e = simpleBeam.getNominalEnergy(2);
        CHECK_EQUAL(1500.*2,e);
    }

    TEST_FIXTURE(SetupFixture, SimpleBeamReturnsCorrectSampledEnergyIfNOTPerNucleon) {
        Beam simpleBeam(1500.,0.5,-10, -0.5,0., 0);

        double e = simpleBeam.sampleEnergy(10);
        CHECK_EQUAL(1500.,e);
    }

    TEST_FIXTURE(SetupFixture, SimpleBeamReturnsCorrectSampledEnergyPerNucleon) {
        Beam simpleBeam(1500.,0.5,-10, -0.5,0., 0, nullptr, nullptr, nullptr, nullptr, true);

        double e = simpleBeam.sampleEnergy(2);
        CHECK_EQUAL(1500.*2,e);
    }

    TEST_FIXTURE(SetupFixture, GaussianEnergySpread) {
        double sigma = 6.; // keV
        auto f = make_unique<TF1>("f1","gaus(0)",-1000.,1000.);
        f -> SetParameter(0,1.); // normalization
        f -> SetParameter(1,0.); // centroid
        f -> SetParameter(2,sigma); // sigma
//        ROOT::Math::WrappedMultiTF1 wf1(*f);
        Beam b(1500.,0.5,-10, -0.5,0., 0,move(f));//&wf1);
        double e = b.sampleEnergy(1);
        CHECK_CLOSE(1500.,e,5.*sigma);
    }

    TEST_FIXTURE(SetupFixture, GaussianAngularDistribution) {
        double sigma = 3.*TMath::Pi()/180.; // deg->rad
        auto f = make_unique<TF1>("f1","gaus(0)*sin(x)",0.,180.);
        f -> SetParameter(0,1.); // normalization
        f -> SetParameter(1,0.); // centroid
        f -> SetParameter(2,sigma); // sigma
        Beam b(1500.,0.5,-0.5,-10, 0., 0,nullptr,nullptr,move(f));
        TVector3 dir = b.sampleDirection();
        CHECK(dir.Theta() < 5*sigma);
    }

    TEST_FIXTURE(SetupFixture, UniformCircularBeamSpot) {
        double R0 = 2.; // mm
        auto f = make_unique<TF1>("Sin Function", "1",0.,R0);
        Beam b(1500.,0,0,-10,0., 0,nullptr,move(f),nullptr);
        double x, y;
        b.sampleXY(x,y);
        double r = sqrt(pow(x,2)+pow(y,2));
        CHECK(r <= R0);
    }

}
