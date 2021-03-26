
#include <unittest++/UnitTest++.h>
#include <ausa/eloss/Ion.h>
#include <TVector3.h>

#include "simX/Particle.h"

using namespace simX;
using namespace AUSA::EnergyLoss;

SUITE(ParticleTest) {

    class SetupFixture {
        public:
            SetupFixture() 
            : proton(Ion(1,1)),
              gamma(Ion(0,0)) {}
            
            Particle proton, gamma;
    };


    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(SetupFixture, TestParticleNameProton) {
        CHECK_EQUAL( "H1", proton.getName() );
    }

    TEST_FIXTURE(SetupFixture, TestParticleNameGamma) {
        CHECK_EQUAL( "gamma", gamma.getName() );
    }

    TEST_FIXTURE(SetupFixture, TestGetKineticEnergyLab) {
        TVector3 p(0.,0.,500.e3);
        proton.setFourMomentumLab(p);
        double mp = 938272.032662;
        double T0 = sqrt( pow(p.Mag(),2.) + pow(mp,2.) ) - mp;
        double T = proton.getKineticEnergyLab();
        CHECK_CLOSE(T0,T,1E-3);
    }

    TEST_FIXTURE(SetupFixture, TestGetDirectionLab) {
        TVector3 p(0.,0.,500.e3);
        proton.setFourMomentumLab(p);
        TVector3 dir(0.,0.,1.);
        CHECK_EQUAL( dir.X(), proton.getDirectionLab().X() );
        CHECK_EQUAL( dir.Y(), proton.getDirectionLab().Y() );
        CHECK_EQUAL( dir.Z(), proton.getDirectionLab().Z() );
    }

    TEST_FIXTURE(SetupFixture, TestGetMomentumLab) {
        TVector3 p(100.e3,-200.e3,500.e3);
        proton.setFourMomentumLab(p);
        CHECK_EQUAL( p*p, proton.getMomentumLab()*proton.getMomentumLab() );
    }

    TEST_FIXTURE(SetupFixture, TestGetPosition) {
        TVector3 pos(-77.,20.,0.1);
        proton.setPosition(pos);
        CHECK_EQUAL( -77., proton.getPosition().X() );
        CHECK_EQUAL(  20., proton.getPosition().Y() );
        CHECK_EQUAL(  0.1, proton.getPosition().Z() );
    }


}
