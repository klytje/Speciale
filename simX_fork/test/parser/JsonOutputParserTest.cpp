//
// Created by munk on 28-06-15.
//


#include <unittest++/UnitTest++.h>

#include "simX/parser/JsonOutputParser.h"

#include <ausa/util/FileUtil.h>

using namespace simX::parser;
using std::string;

SUITE(JsonOutputParserTest) {
    struct Fixture {
        JsonOutputParser parser;

        Fixture() : parser(AUSA::asciiFileToString("test/_res/setup/setup.json"))
        {

        }
    };

    TEST_FIXTURE(Fixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(Fixture, SUHaveMapping) {
        parser.getDoubleMapping("SU");
    }

    TEST_FIXTURE(Fixture, SUHaveMultiplicityBranch_SU_S) {
        auto m = parser.getDoubleMapping("SU");

        CHECK_EQUAL("SU_S", m.front.mul);
    }

    TEST_FIXTURE(Fixture, SDHaveMultiplicityBranch_SD_S) {
        auto m = parser.getDoubleMapping("SD");

        CHECK_EQUAL("SD_S", m.front.mul);
    }

    TEST_FIXTURE(Fixture, Pad2HaveEnergyBranch_PAD2E) {
        auto m = parser.getMapping("Pad2");

        CHECK_EQUAL("PAD2E", m.adc);
    }

    TEST_FIXTURE(Fixture, SUFrontTDCIsDead) {
        auto m = parser.getDoubleMapping("SU");

        CHECK_EQUAL("-DEAD-", m.front.tdc);
    }
}
