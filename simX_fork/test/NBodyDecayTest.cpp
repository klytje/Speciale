
#include "simX/NBodyDecay.h"
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

SUITE(NBodyDecayTest) {
    class SetupFixture {
        public:
            SetupFixture() 
            : c12(Particle(Ion(6,12),12710.)),
              daughters({Ion(2,4),Ion(2,4),Ion(2,4)}),
              ps(NBodyDecay::withDefaultGenerator(c12,daughters))
            {}
            Particle c12;
            vector<Particle> daughters;
            NBodyDecay ps;
    };

    TEST_FIXTURE(SetupFixture, ReturnsCorrectMultiplicity) {
        int m = ps.getMultiplicity();
        CHECK_EQUAL(3,m);
    }

    TEST_FIXTURE(SetupFixture, CheckRunProcessGivesSensibleOutput) {
        // Isotope masses
        double mc12 = c12.getMass();
        double malpha = Ion(2,4).getMass();
        // run process
        ps.runProcess();
        // loop over daughters
        double esum = 0;
        TVector3 psum(0,0,0);
        auto& daughtersPtr = ps.getDaughters();
        for (auto& d : daughtersPtr) {
            esum += d->getKineticEnergyLab();
            psum += d->getMomentumLab();
        }
        double e3alpha = 3.*malpha-mc12;
        CHECK_CLOSE(12710.-e3alpha,esum,1e-3);
        CHECK_CLOSE(0.,psum.Mag(),1e-3);
    }

       
}
