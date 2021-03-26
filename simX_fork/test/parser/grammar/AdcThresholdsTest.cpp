#include <unittest++/UnitTest++.h>
#include <memory>

#include "simX/parser/grammar/AdcThresholdsGrammar.h"

using namespace simX::parser;
using std::string;

SUITE(AdcThresholdsGrammarTest) {

    using Val = std::pair<int, std::vector<double>>;

    Val parseString(const string& s) {
        using boost::spirit::qi::phrase_parse;
        using boost::spirit::ascii::blank;

        Val v;
        bool r = phrase_parse(s.begin(), s.end(), AdcThresholdsGrammar{}, blank, v);
        if (!r) throw std::invalid_argument("failed to parse '" + s + "'");
        return v;
    }


    TEST(TwoNumbersSeparatedBySpaceShouldGivesVectorWithTwo) {
        auto res = parseString(R"(unit: keV
8
9)");
        double expected[] = {8, 9};

        CHECK_ARRAY_CLOSE(expected, res.second, 2, 1E-5);
    }

    TEST(HandlesSpaceCorrectly) {
        auto res = parseString(R"(unit: keV
  8
9)");
        double expected[] = {8, 9};

        CHECK_ARRAY_CLOSE(expected, res.second, 2, 1E-5);
    }


    TEST(KeVGivesInt3) {
        auto res = parseString(R"(unit: keV
8
9)");

        CHECK_EQUAL(3, res.first);
    }
}
