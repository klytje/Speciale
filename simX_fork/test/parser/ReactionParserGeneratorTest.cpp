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
using namespace simX;
using std::string;

SUITE(ReactionParserGeneratorTest) {
    struct Fixture {
        ReactionParser parser;
    };

    TEST_FIXTURE(Fixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(Fixture, IfMoreThan2ParticlesNBodyDecayIsReturned) {
        auto chain = parser.parseString(R"(
  beam: C12 Ex: 15000keV
  -> {a a a}
)");
        auto iter = begin(*chain);
        NuclearProcess & d = *iter;

        CHECK(typeid(simX::NBodyDecay) == typeid(d));
    }

    TEST_FIXTURE(Fixture, AllSpecifiedParticlesAreIncluded) {
        auto chain = parser.parseString(R"(
  beam: C12 Ex: 15000keV
  -> {a a a}
)");
        auto iter = begin(*chain);
        NuclearProcess & d = *iter;

        CHECK_EQUAL(3, d.getDaughters().size());
    }

    TEST_FIXTURE(Fixture, ThatDefaultPhaseSpaceIsSet) {
        auto chain = parser.parseString(R"(
  beam: C12 Ex: 15000keV
  -> {a a a}
)");
        auto iter = begin(*chain);
        NuclearProcess & d = *iter;

        d.runProcess();
    }

    TEST_FIXTURE(Fixture, ThatSpecialPurposeGeneratorFactoryFunctionsAreCalled) {
        int n = 0;

        parser.registerFactory("HEST", [&](const NuclearProcess& p, ReactionParser::Options& o) {
            n++;
            return std::unique_ptr<FinalStateGenerator>();
        });

        auto chain = parser.parseString(R"(
  beam: He3
  target: Be9
  -> {
generator: HEST
a a a
}
)");
        CHECK_EQUAL(1, n);
    }

    TEST_FIXTURE(Fixture, ThatOptionsMapIsParsed) {
        int n = 0;

        parser.registerFactory("HEST", [&](const NuclearProcess& p, ReactionParser::Options& o) {
            n+= o.count("phi");
            return std::unique_ptr<FinalStateGenerator>();
        });

        auto chain = parser.parseString(R"(
  beam: He3
  target: Be9
  -> {
generator: HEST(phi="a")
a a a
}
)");
        CHECK_EQUAL(1, n);
    }

    TEST_FIXTURE(Fixture, IfKeywordIsPREGREN_SamplerIsUsed) {
        auto chain = parser.parseString(R"(
  beam: Be8
  -> {
generator: PREGEN(file="test/_res/samplers/input_with_weight.txt" unit="keV")
a a}
)");
        auto iter = begin(*chain);
        NuclearProcess & d = *iter;

        CHECK(typeid(simX::NBodyDecay) == typeid(d));
    }

    TEST_FIXTURE(Fixture, IfPregenFileEndsWithRootRootInputIsUsed) {
        auto chain = parser.parseString(R"(
  beam: Be8
  -> {
generator: PREGEN(file="test/_res/samplers/sample_without_weight.root" unit="MeV")
a a}
)");
        auto iter = begin(*chain);
        NuclearProcess & d = *iter;

        CHECK(typeid(simX::NBodyDecay) == typeid(d));
    }
}
