//
// Created by munk on 28-06-15.
//


#include <unittest++/UnitTest++.h>
#include "simX/CompoundFormation.h"
#include "simX/angular/AngularCorrelation.h"

#include "simX/parser/ReactionParser.h"

#include <typeinfo>
#include "simX/TwoBodyDecay.h"

using namespace simX::parser;
using namespace simX::angular;
using namespace simX;
using std::string;

SUITE(ReactionParserTest) {
    struct Fixture {
        ReactionParser parser;
    };

    TEST_FIXTURE(Fixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(Fixture, BeamParticleIsSetAsFirstParticleOfChain) {
        auto chain = parser.parseString(R"(
  beam: He3
)");
        auto f = chain -> getBeam();
        CHECK_EQUAL("He3", f.getName());
    }

    TEST_FIXTURE(Fixture, BeamParticleHasEmptyChain) {
        auto chain = parser.parseString(R"(
  beam: He3
)");
        auto b = begin(*chain);
        auto e = end(*chain);

        CHECK(b == e);
    }

    TEST_FIXTURE(Fixture, IfTargetIsSpecifiedFirstDecayIsFusionReaction) {
        auto chain = parser.parseString(R"(
  beam: He3
  target: Be9
)");
        auto iter = begin(*chain);
        NuclearProcess & d = *iter;

        CHECK(typeid(simX::CompoundFormation) == typeid(d));
    }

    TEST_FIXTURE(Fixture, IfTwoDecayProdoctsSpecifiedANodeWithTwoDecayIsGiven) {
        auto chain = parser.parseString(R"(
  beam: He3
  target: Be9
  -> {Li6 Li6}
)");
        auto iter = begin(*chain);
        ++iter;
        NuclearProcess & d = *iter;

        CHECK(typeid(simX::NBodyDecay) == typeid(d));
    }

    TEST_FIXTURE(Fixture, LvalueIsForwardedToDecay) {
        auto chain = parser.parseString(R"(
  beam: He3
  target: Be9
  -> {L:6 Li6 Li6}
)");
        auto iter = begin(*chain);
        ++iter;
        NuclearProcess & d = *iter;

        auto& t = static_cast<NBodyDecay &>(d);
        CHECK_EQUAL(6, t.getL());
    }

    TEST_FIXTURE(Fixture, ReactionParserParsesEntireChain) {
        auto chain = parser.parseString(R"(
  beam: He3
  target: Be9
  -> {L:6
      Li6 -> {L:9 a d}
      Li6}
)");
        auto iter = begin(*chain);
        ++iter;
        ++iter;
        NuclearProcess & d = *iter;

        auto& t = static_cast<NBodyDecay &>(d);
        CHECK_EQUAL(9, t.getL());
    }

    TEST_FIXTURE(Fixture, EmptyReactionBlockThrows) {
        auto input = R"(
  beam: He3
  target: Be9
  -> {}
)";
        CHECK_THROW(parser.parseString(input), std::invalid_argument);
    }

    TEST_FIXTURE(Fixture, ParentOfFirstBreakupIsFirstParticle) {
        auto chain = parser.parseString(R"(
  beam: Be8
  -> {a a}
)");
        auto iter = begin(*chain);
        NuclearProcess & d = *iter;

        auto& first = chain -> getBeam();
        auto& prod = d.getDaughters();

        CHECK_EQUAL(&first, prod[0] -> getParent());
        CHECK_EQUAL(&first, prod[1] -> getParent());
    }

    TEST_FIXTURE(Fixture, ReactionParserUsesCustomAD) {
        int n = 0;
        parser.registerFactory("HEST", [&](const ReactionParser::Limits& phi, const ReactionParser::Limits& theta, const::ReactionParser::Options& m) {
            n++;
            return std::unique_ptr<AngularCorrelation>{};
        });


        parser.parseString(R"(
  beam: He3
  target: Be9
  -> {
     AD: HEST
     Li6
     Li6
     }
)");
        CHECK_EQUAL(1, n);
    }

    TEST_FIXTURE(Fixture, ReactionParserPassesAdOptions) {
        int n = 0;
        parser.registerFactory("HEST", [&](const ReactionParser::Limits& phi, const ReactionParser::Limits& theta, const::ReactionParser::Options& m) {
            n+=m.count("phi");
            return std::unique_ptr<AngularCorrelation>{};
        });


        parser.parseString(R"(
  beam: He3
  target: Be9
  -> {
     AD: HEST(phi="2")
     Li6
     Li6
     }
)");
        CHECK_EQUAL(1, n);
    }

    TEST_FIXTURE(Fixture, CheckTreeIsIntact) {
        auto chain = parser.parseString(R"(
  beam: He3
  target: Be9
  -> {
     Li6
     Li6
     }
)");
        auto iter = chain->getTree().begin();
        NuclearProcess & d = **iter;

        CHECK(typeid(simX::CompoundFormation) == typeid(d));

        simX::CompoundFormation& c = static_cast<simX::CompoundFormation&>(d);
        auto children = c.getDaughters();
        CHECK_EQUAL("C12", children[0] ->getName());

        ++iter;
        auto& decay = **iter;
        CHECK(typeid(simX::NBodyDecay) == typeid(decay));
    }

    TEST_FIXTURE(Fixture, NoTrackTagIsSetCorretly) {
        auto chain = parser.parseString(R"(
  beam: He3 notrack
  target: Be9
  -> {
     Li6
     Li6 notrack
     }
)");
        auto b = chain -> getBeam();

        CHECK(!b.getTracking());

        auto i = ++(begin(*chain));

        auto& process = *i;
        auto& v = process.getDaughters();

        CHECK(v[0]->getTracking());
        CHECK(!v[1]->getTracking());
    }

    TEST_FIXTURE(Fixture, ParseMultiple) {
        std::vector<ReactionParser::ReactionBranch> chains = parser.parseMultiple("test/_res/reactions/multipleReactions.json");
        CHECK_EQUAL(2, chains.size());
        CHECK_EQUAL("reaction_Li7p_aa.simX", chains[0].name);
        CHECK_EQUAL("reaction_Li7p_Li7p.simX", chains[1].name);
        CHECK_EQUAL(0.6, chains[0].ratio);
        CHECK_EQUAL(0.4, chains[1].ratio);
    }

    TEST_FIXTURE(Fixture, TestGraphIsPossibleAD) {
        parser.parseString(R"(
  beam: He3
  target: Be9
  -> {
     theta : [10 160]
     AD : GRAPH(file="test/_res/dwba_output.dat")
     a
     Be8
  }
)");
    }
}
