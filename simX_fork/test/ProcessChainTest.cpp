#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include <simX/angular/IsotropicAngularCorrelation.h>
#include <simX/parser/ReactionParser.h>

#include "simX/angular/AngularCorrelationTF1.h"
#include "simX/TwoBodyDecay.h"
#include "simX/Particle.h"
#include "simX/ProcessChain.h"

#include "ausa/util/memory"

using std::string;
using AUSA::EnergyLoss::Ion;
using namespace simX;
using namespace simX::parser;

SUITE(ProcessChainTest) {

    TEST(ChainIteratorVisitsAll) {
        ReactionParser p;
        auto c1 = p.parseString("beam: He3 target: B11 -> {d C12 -> {a a a}}");

        int count = 0;
        for (auto& n : *c1) {
            n.getDaughters(); // To shut up warning...
            count++;
        }
        CHECK_EQUAL(3, count);
    }

    TEST(FindDecayFindsNode) {
        ReactionParser p;
        auto c = p.parseString("beam: He3 target: Be9 -> {Be8 a}");

        auto ptr = c->findDecayOf("C12");

        CHECK(ptr != nullptr);

        auto& d = ptr -> getDaughters();
        CHECK_EQUAL(2, d.size());
        
        CHECK_EQUAL("Be8", d[0]->getName());
        CHECK_EQUAL("He4", d[1]->getName());
    }

    TEST(C12CanDecayTo_d_B10) {
        ReactionParser p;
        p.parseString("beam: H1 target: B11 -> {d B10}");
    }

    TEST(DaugtersMustUseAllNucleonsInBinary) {
        ReactionParser p;
        CHECK_THROW(p.parseString("beam: H1 target: B11 -> {p B10}"), std::invalid_argument);
    }

    TEST(ParticleCanDecayToItselfPlusGamma) {
        ReactionParser p;
        // first try valid input
        p.parseString("beam: C12 -> {g C12}");
    }

    TEST(ParticleCannotDecayToOtherPlusGamma) {
        ReactionParser p;
        CHECK_THROW(p.parseString("beam: Ne20 -> {g Ni56}"), std::invalid_argument);
    }    

    TEST(MultipleChainedValidBinaryReactionsAreValid) {
        ReactionParser p;
        // first try valid input
        p.parseString("beam: He3 target: B11 -> {d C12 -> {a a a}}");
    }

    TEST(LeafReactionUsingMoreNucleonsThanParentIsInvalid) {
        ReactionParser p;
        CHECK_THROW(p.parseString("beam: He3 target: B11 -> {d C12 -> {a a He6}}"), std::invalid_argument);
    }

     TEST(AlphaIsFirstFinalState) {
        ReactionParser p;
        auto c1 = p.parseString("beam: He3 target: Be9 -> {a Be8 Ex: 16MeV -> {g Be8 -> {He3 He5}}}");

        auto final = c1->findFinalStates();

        CHECK_EQUAL("He4", final[0]->getName());
    }

    TEST(GammaIsSecondFinalState) {
        ReactionParser p;
        auto c1 = p.parseString("beam: He3 target: Be9 -> {a Be8 Ex: 16MeV -> {g Be8 -> {He3 He5}}}");

        auto final = c1->findFinalStates();

        CHECK_EQUAL("gamma", final[1]->getName());
    }

}
