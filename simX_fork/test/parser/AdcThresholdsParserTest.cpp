//
// Created by munk on 28-06-15.
//


#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include <simX/parser/AdcThresholdParser.h>

using namespace simX::parser;
using namespace simX;
using std::string;

SUITE(AdcThresholdsParserTest) {

    TEST(ThresholdsWithKeVUnitIsNotTransformed) {

        auto res = parseAdcThresholds(R"(unit: keV
8
9)", 2);
        double expected[] = {8, 9};

        CHECK_ARRAY_CLOSE(expected, res.getThresholds(), 2, 1E-5);
    }

    TEST(ThresholdsWith_eVUnitIsDividedBy1000) {

        auto res = parseAdcThresholds(R"(unit: eV
8
9)", 2);
        double expected[] = {8/1000., 9/1000.};

        CHECK_ARRAY_CLOSE(expected, res.getThresholds(), 2, 1E-5);
    }

    TEST(ParserThrowIfTwoFewEntries) {
        string input = R"(unit: eV
8
9)";
        CHECK_THROW(parseAdcThresholds(input, 1), std::invalid_argument);
    }
}
