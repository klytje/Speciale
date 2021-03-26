
#include "simX/Particle.h"

#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include <math.h>
#include "ausa/util/memory"
#include <ausa/eloss/Ion.h>
#include <simX/generator/SingleBodyKinematics.h>


using namespace std;
using namespace simX;
using namespace AUSA::EnergyLoss;

SUITE(SingleBodyKinematicsTest) {
    class SetupFixture {
        public:
            SetupFixture() 
            : alpha(Particle(Ion(2,4)))
            {}
            Particle alpha;
    };

    TEST_FIXTURE(SetupFixture, TestTripleAlphaPositiveQValue) {

        TLorentzVector vec(TVector3(0, 0, 3e5), 10000);
        alpha.setFourMomentumLab(vec);
        auto dMomentum = SingleBodyKinematics(alpha).getFourMomenta();
        CHECK_EQUAL(1, dMomentum.size());
        CHECK(vec == dMomentum[0]);
    }

}
