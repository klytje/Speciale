
#include "simX/Particle.h"
#include "simX/generator/NBodyPhaseSpace.h"

#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include <math.h>
#include "ausa/util/memory"
#include <ausa/eloss/Ion.h>


using namespace std;
using namespace simX;
using namespace AUSA::EnergyLoss;

SUITE(NBodyPhaseSpaceTest) {
    class SetupFixture {
        public:
            SetupFixture() 
            : c12(Particle(Ion(6,12))),
              daughters({Ion(2,4),Ion(2,4),Ion(2,4)}),
              mc12(c12.getMass()),
              malpha(Ion(2,4).getMass())              
            {}
            Particle c12;
            vector<Particle> daughters;
            double mc12, malpha;
    };

    TEST_FIXTURE(SetupFixture, TestTripleAlphaPositiveQValue) {
        double ex=12710;
        c12.setNominalExcitationEnergy(ex);
        Particle * a1 = new Particle(Ion(2,4));
        Particle * a2 = new Particle(Ion(2,4));
        Particle * a3 = new Particle(Ion(2,4));
        vector<Particle*> daughtersPtr = {a1,a2,a3};
        NBodyPhaseSpace ps(c12,daughtersPtr);
        auto& PD = ps.getFourMomenta();
        double ksum = 0;
        for (auto& p : PD) {
            ksum += p.E()-p.Mag();
        }
        double e3alpha = 3.*malpha-mc12;
        CHECK_CLOSE(ex-e3alpha,ksum,1e-3);
    }

    TEST_FIXTURE(SetupFixture, TestTripleAlphaNegativeQValue) {
        double ex=6000.;
        c12.setNominalExcitationEnergy(ex);
        Particle * a1 = new Particle(Ion(2,4));
        Particle * a2 = new Particle(Ion(2,4));
        Particle * a3 = new Particle(Ion(2,4));
        vector<Particle*> daughtersPtr = {a1,a2,a3};
        NBodyPhaseSpace ps(c12,daughtersPtr);
        CHECK_THROW( auto& PD = ps.getFourMomenta(), std::invalid_argument );
    }

}
