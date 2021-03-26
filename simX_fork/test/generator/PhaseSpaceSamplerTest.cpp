
#include <unittest++/UnitTest++.h>
#include <simX/generator/PhaseSpaceSampler.h>
#include <simX/samplers/SampleSet.h>

#include "../samplers/DummySampler.h"

using namespace std;
using namespace simX;


SUITE(PhaseSpaceSamplerTest) {



    initializer_list<initializer_list<double>> init(initializer_list<initializer_list<double>> list) {
        return list;
    }

    TEST(Sanity) {
        CHECK_EQUAL(1,1);
    }

    TEST(ZeroMomentumIsCorrectInput) {
        auto sampler = std::make_unique<DummySampler>(init({{0,0,0}}));

        auto alpha = AUSA::EnergyLoss::Ion("He4");
        Particle p0(alpha);
        vector<Particle*> particles = {&p0};

        PhaseSpaceSampler gen(move(sampler), particles, "keV");
        
        auto momenta = gen.getFourMomenta();
        CHECK_EQUAL(1, momenta.size());
        
        CHECK_CLOSE(alpha.getMass(), momenta[0].M(), 1E-8);
    }

    TEST(ZeroMomentumForTwoParticlesIsCorrectInput) {
        auto sampler = std::make_unique<DummySampler>(init({{0,0,0, 0,0,0}}));

        auto alpha = AUSA::EnergyLoss::Ion("He4");
        auto he3 = AUSA::EnergyLoss::Ion("He3");
        Particle p0(alpha);
        Particle p1(he3);
        vector<Particle*> particles = {&p0, &p1};

        PhaseSpaceSampler gen(move(sampler), particles, "keV");

        auto momenta = gen.getFourMomenta();
        CHECK_EQUAL(2, momenta.size());

        CHECK_CLOSE(alpha.getMass(), momenta[0].M(), 0.1);
        CHECK_CLOSE(he3.getMass(),   momenta[1].M(), 0.1);
    }

    TEST(Momentum1000keVMeansInvariantMassIsNotRestMass) {
        auto sampler = std::make_unique<DummySampler>(init({{1E3,0,0, 0,1E3,0}}));

        auto alpha = AUSA::EnergyLoss::Ion("He4");
        auto he3 = AUSA::EnergyLoss::Ion("He3");
        Particle p0(alpha);
        Particle p1(he3);
        vector<Particle*> particles = {&p0, &p1};

        PhaseSpaceSampler gen(move(sampler), particles, "keV");

        auto momenta = gen.getFourMomenta();
        CHECK_EQUAL(2, momenta.size());

        CHECK_CLOSE(sqrt(pow(alpha.getMass(), 2) + pow(1E3,2)), momenta[0].M(), 0.5);
        CHECK_ARRAY_CLOSE(TVector3(1E3,0,0), momenta[0].Vect(), 3, 1E-6);

        CHECK_CLOSE(sqrt(pow(he3.getMass(), 2) + pow(1E3,2)), momenta[1].M(), 1);
        CHECK_ARRAY_CLOSE(TVector3(0,1E3,0), momenta[1].Vect(), 3, 1E-6);
    }

    TEST(Momentum1MeVIsEquivalentWith1000keV) {
        auto sampler = std::make_unique<DummySampler>(init({{1,0,0, 0,1,0}}));

        auto alpha = AUSA::EnergyLoss::Ion("He4");
        auto he3 = AUSA::EnergyLoss::Ion("He3");
        Particle p0(alpha);
        Particle p1(he3);
        vector<Particle*> particles = {&p0, &p1};

        PhaseSpaceSampler gen(move(sampler), particles, "MeV");

        auto momenta = gen.getFourMomenta();
        CHECK_EQUAL(2, momenta.size());

        CHECK_CLOSE(sqrt(pow(alpha.getMass(), 2) + pow(1E3,2)), momenta[0].M(), 0.5);
        CHECK_ARRAY_CLOSE(TVector3(1E3,0,0), momenta[0].Vect(), 3, 1E-6);

        CHECK_CLOSE(sqrt(pow(he3.getMass(), 2) + pow(1E3,2)), momenta[1].M(), 1);
        CHECK_ARRAY_CLOSE(TVector3(0,1E3,0), momenta[1].Vect(), 3, 1E-6);
    }

    TEST(GeneratorWrapsAroundIfInputIsExhausted) {
        auto sampler = std::make_unique<DummySampler>(init({{1E3,0,0, 0,1E3,0}, {0,0,0, 0,0,0}}));

        auto alpha = AUSA::EnergyLoss::Ion("He4");
        auto he3 = AUSA::EnergyLoss::Ion("He3");
        Particle p0(alpha);
        Particle p1(he3);
        vector<Particle*> particles = {&p0, &p1};

        PhaseSpaceSampler gen(move(sampler), particles, "keV");

        // Event 0
        auto momenta = gen.getFourMomenta();
        CHECK_CLOSE(sqrt(pow(alpha.getMass(), 2) + pow(1E3,2)), momenta[0].M(), 0.5);
        CHECK_ARRAY_CLOSE(TVector3(1E3,0,0), momenta[0].Vect(), 3, 1E-6);

        CHECK_CLOSE(sqrt(pow(he3.getMass(), 2) + pow(1E3,2)), momenta[1].M(), 1);
        CHECK_ARRAY_CLOSE(TVector3(0,1E3,0), momenta[1].Vect(), 3, 1E-6);

        // Event 1
        momenta = gen.getFourMomenta();
        CHECK_CLOSE(alpha.getMass(), momenta[0].M(), 1E-8);
        CHECK_CLOSE(he3.getMass(),   momenta[1].M(), 1E-8);

        // Event 2 == 0
        momenta = gen.getFourMomenta();
        CHECK_CLOSE(sqrt(pow(alpha.getMass(), 2) + pow(1E3,2)), momenta[0].M(), 1);
        CHECK_ARRAY_CLOSE(TVector3(1E3,0,0), momenta[0].Vect(), 3, 1E-6);

        CHECK_CLOSE(sqrt(pow(he3.getMass(), 2) + pow(1E3,2)), momenta[1].M(), 1);
        CHECK_ARRAY_CLOSE(TVector3(0,1E3,0), momenta[1].Vect(), 3, 1E-6);
    }


}
