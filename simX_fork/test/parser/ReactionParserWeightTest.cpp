//
// Created by munk on 28-06-15.
//


#include <unittest++/UnitTest++.h>
#include "simX/CompoundFormation.h"
#include "simX/angular/AngularCorrelation.h"

#include "simX/parser/ReactionParser.h"
#include "simX/weight/GammaPhaseSpace.h"
#include "simX/weight/BreitWignerWithPenetrability.h"
#include "simX/weight/SimpleBreitWigner.h"
#include "simX/weight/TF1WeightCalculator.h"

#include <typeinfo>
#include <simX/weight/ProductWeightCalculator.h>
#include "simX/TwoBodyDecay.h"

using namespace simX::parser;
using namespace simX;
using std::string;

SUITE(ReactionParserWeightTest) {
    struct Fixture {
        ReactionParser parser;
    };

    TEST_FIXTURE(Fixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }


    TEST_FIXTURE(Fixture, TwoBodyDecayWithGammaWillBeAssignedSimpleBreitWignerByDefault) {
        auto c = parser.parseString(R"(beam: He4 Ex: 2MeV -> {L:3 He4 g})");

        auto& process = *begin(*c);

        NBodyDecay& d = static_cast<NBodyDecay&>(process);

        auto calc = d.getWeightCalculator();

        auto prodCalc = dynamic_cast<const ProductWeightCalculator *>(calc);
        CHECK(prodCalc != nullptr);

        auto& gCalc = prodCalc->getFunction2();

        CHECK(dynamic_cast<const GammaPhaseSpace*>(&gCalc) != nullptr);
    }

    TEST_FIXTURE(Fixture, TwoBodyParticleDecayWillHaveBreitWignerPen) {
        auto c = parser.parseString(R"(beam: Be8 -> {a a})");

        auto& process = *begin(*c);

        NBodyDecay& d = static_cast<NBodyDecay&>(process);

        auto calc = d.getWeightCalculator();

        CHECK(dynamic_cast<const BreitWignerWithPenetrability*>(calc) != nullptr);
    }

    TEST_FIXTURE(Fixture, NBodyWillBeAssignedSimpleBreitWigner) {
        auto c = parser.parseString(R"(beam: C12 Ex: 20MeV -> {a a a})");

        auto& process = *begin(*c);

        NBodyDecay& d = static_cast<NBodyDecay&>(process);

        auto calc = d.getWeightCalculator();

        CHECK(dynamic_cast<const SimpleBreitWigner*>(calc) != nullptr);
    }

    TEST_FIXTURE(Fixture, UserInputGammaOverridesDefault) {
        auto c = parser.parseString(R"(beam: Be8 -> {weight: GAMMA a a})");

        auto& process = *begin(*c);

        NBodyDecay& d = static_cast<NBodyDecay&>(process);

        auto calc = d.getWeightCalculator();

        CHECK(dynamic_cast<const GammaPhaseSpace*>(calc) != nullptr);
    }

    TEST_FIXTURE(Fixture, UserInputBWGivesSimpleBreitWigner) {
        auto c = parser.parseString(R"(beam: Be8 -> {weight: BW a a})");

        auto& process = *begin(*c);

        NBodyDecay& d = static_cast<NBodyDecay&>(process);

        auto calc = d.getWeightCalculator();

        CHECK(dynamic_cast<const SimpleBreitWigner*>(calc) != nullptr);
    }

    TEST_FIXTURE(Fixture, UserInputPENGivesBreitWignerWithPenetrability) {
        auto c = parser.parseString(R"(beam: Be8 -> {weight: PEN a a})");

        auto& process = *begin(*c);

        NBodyDecay& d = static_cast<NBodyDecay&>(process);

        auto calc = d.getWeightCalculator();

        CHECK(dynamic_cast<const BreitWignerWithPenetrability*>(calc) != nullptr);
    }

    TEST_FIXTURE(Fixture, UserInput_GAMMA_BW_GivesBreitWignerWithGamma) {
        auto c = parser.parseString(R"(beam: Be8 -> {weight: PEN*BW a a})");

        auto& process = *begin(*c);

        NBodyDecay& d = static_cast<NBodyDecay&>(process);

        auto calc = d.getWeightCalculator();

        auto pCalc = dynamic_cast<const ProductWeightCalculator*>(calc);
        CHECK(pCalc != nullptr);

        CHECK(dynamic_cast<const BreitWignerWithPenetrability*>(&pCalc->getFunction1()) != nullptr);
        CHECK(dynamic_cast<const SimpleBreitWigner*>(&pCalc->getFunction2()) != nullptr);
    }

    TEST_FIXTURE(Fixture, TwoBodyParticleDecayWillHaveCustomWeight) {
        auto c = parser.parseString("beam: Li11 -> {weight: TF1(f=\"sin(x)\") n Li10}");

        auto& process = *begin(*c);

        NBodyDecay& d = static_cast<NBodyDecay&>(process);

        auto calc = d.getWeightCalculator();
        CHECK(calc != nullptr);

        CHECK(dynamic_cast<const TF1WeightCalculator*>(calc) != nullptr);
    }
}
