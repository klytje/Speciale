//
// Created by munk on 28-06-15.
//


#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include <simX/parser/TdcThresholdParser.h>

using namespace simX::parser;
using namespace simX;
using std::string;

SUITE(TdcThresholdsParserTest) {

//    TEST(ThresholdsWithKeVUnitIsNotTransformed) {
//
//        parseTdcThresholds(R"(
//{
//  "file": "test/_res/TDC.root",
//  "histograms": [
//    "eff-spec/SU/F/effSU-1",
//    "eff-spec/SU/F/effSU-2"
//  ]
//}
//)", 2);


//        CHECK_EQUAL(1, 1);
//    }

//    TEST(ThresholdsWith_eVUnitIsDividedBy1000) {
//
//        auto res = parseAdcThresholds(R"(unit: eV
//8
//9)", 2);
//        double expected[] = {8/1000., 9/1000.};
//
//        CHECK_ARRAY_CLOSE(expected, res.getThresholds(), 2, 1E-5);
//    }
//
//    TEST(ParserThrowIfTwoFewEntries) {
//        string input = R"(unit: eV
//8
//9)";
//        CHECK_THROW(parseAdcThresholds(input, 1), std::invalid_argument);
//    }
}
