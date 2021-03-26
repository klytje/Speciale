
#include <unittest++/UnitTest++.h>
#include <memory>

#include <ausa/eloss/Ion.h>
#include <ausa/constants/Mass.h>
#include "ausa/util/memory"
#include <simX/angular/IsotropicAngularCorrelation.h>
#include <simX/angular/FixedAngularCorrelation.h>
#include <simX/NBodyDecay.h>
#include <simX/parser/ReactionParser.h>

#include "simX/Particle.h"
#include "simX/TwoBodyDecay.h"
#include "simX/generator/TwoBodyKinematics.h"
#include "simX/angular/AngularCorrelationTF1.h"

using namespace simX;
using namespace simX::angular;
using namespace AUSA::EnergyLoss;

SUITE(TwoBodyDecayTest) {

    class SetupFixture {
        public:
            SetupFixture() 
            : C12x(Ion(6,12),12710.),
              alphaDecay(C12x,{Ion(2,4),Ion(4,8)}, nullptr, std::make_unique<FixedAngularCorrelation>(90,0) ) {
                auto& d = alphaDecay.getDaughters();
                d1 = d[0];
                d2 = d[1];

                alphaDecay.setFinalStateGenerator(std::make_unique<TwoBodyKinematics>(C12x, *d1, *d2));
            }
            
            Particle C12x;
            Particle *d1, *d2;
            NBodyDecay alphaDecay;
    };

    TEST_FIXTURE(SetupFixture, TestEnergySumsToParentEnergy) {
        alphaDecay.runProcess();

        auto k1 = d1 -> getKineticEnergyLab();
        auto k2 = d2 -> getKineticEnergyLab();

        CHECK_CLOSE(AUSA::Constants::isotopeMass("C12")+12710., k1+k2+d1->getMass()+d2->getMass(), 1E-5);
    }

    TEST_FIXTURE(SetupFixture, TestParticlesAreBackToBack) {
        alphaDecay.runProcess();

        auto v = d1->getMomentumLab() + d2->getMomentumLab();
        CHECK_CLOSE(0., v.X(), 1E-3);
        CHECK_CLOSE(0., v.Y(), 1E-3);
        CHECK_CLOSE(0., v.Z(), 1E-3);
    }

    TEST_FIXTURE(SetupFixture, TestMomentumMagnitudeIsCorrect) {
        alphaDecay.runProcess();

        auto E = C12x.getMass() + C12x.getExcitationEnergy();
        auto mx1 = d1->getMass() + d1->getExcitationEnergy();
        auto mx2 = d2->getMass() + d2->getExcitationEnergy();

        auto E1 = (E*E+mx1*mx1-mx2*mx2)/(2*E);
        auto p = std::sqrt(E1*E1-mx1*mx1);


        CHECK_CLOSE(p, d1->getMomentumLab().Mag(),1E-3);
    }

    TEST_FIXTURE(SetupFixture, TestDaughtersAreTransformedToLAB) {
        auto K = 1E3; // Kinetic energy
        C12x.setFourMomentumLab(K, {0,0,1});

        alphaDecay.runProcess();

        auto v = C12x.getFourMomentumLab();
        auto g = v.Gamma();
        auto b = v.Beta();

        auto E = v.E()-K; // Energy in rest frame
        auto mx1 = d1->getMass() + d1->getExcitationEnergy();
        auto mx2 = d2->getMass() + d2->getExcitationEnergy();

        auto E1 = (E*E+mx1*mx1-mx2*mx2)/(2*E);
        auto E2 = E - E1;

        // Transverse momentum
        auto p1T = sqrt(E1*E1 - mx1*mx1);
        auto p2T = sqrt(E2*E2 - mx2*mx2);

        // 1D Lorentz transform https://en.wikipedia.org/wiki/Lorentz_transformation#Boost_in_the_x-direction
        auto E1d = g*E1;        auto p1d = g*b*E1;
        auto E2d = g*E2;        auto p2d = g*b*E2;


        // Check particle 1
        auto p1 = d1->getFourMomentumLab();
        CHECK_CLOSE(p1T, p1.X(), 1E-3);
        CHECK_CLOSE(0., p1.Y(), 1E-3);
        CHECK_CLOSE(p1d, p1.Z(), 1E-3);
        CHECK_CLOSE(E1d, p1.E(), 1E-3);

        auto p2 = d2->getFourMomentumLab();
        CHECK_CLOSE(-p2T, p2.X(), 1E-3);
        CHECK_CLOSE(0., p2.Y(), 1E-3);
        CHECK_CLOSE(p2d, p2.Z(), 1E-3);
        CHECK_CLOSE(E2d, p2.E(), 1E-3);
    }

    TEST_FIXTURE(SetupFixture, TestThatHandMadeReactionIsTheSameAsParsed) {
        alphaDecay.runProcess();

        auto k1 = d1 -> getKineticEnergyLab();
        auto k2 = d2 -> getKineticEnergyLab();


        parser::ReactionParser p;

        auto c = p.parseString(R"(
            beam: C12 Ex: 12710keV
            -> {
                AD: FIXED(theta="90" phi="90")
                He4
                Be8
            }
            )");

        for (auto& i :*c) {
            i.runProcess();
        }

        auto& proc = *(begin(*c));

        auto& d = proc.getDaughters();
        CHECK_CLOSE(k1, d[0]->getKineticEnergyLab(), 1E-5);
        CHECK_CLOSE(k2, d[1]->getKineticEnergyLab(), 1E-5);
    }
}
