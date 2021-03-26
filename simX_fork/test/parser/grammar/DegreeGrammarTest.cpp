#include <unittest++/UnitTest++.h>
#include <memory>
#include <TMath.h>

#include "simX/parser/grammar/DegreeGrammar.h"

using namespace simX::parser;
using std::string;

SUITE(DegreeGrammarTest) {

    double parseString(const string& s) {
        using boost::spirit::qi::phrase_parse;
        using boost::spirit::ascii::space;

        double v;
        bool r = phrase_parse(s.begin(), s.end(), DegreeGrammar{}, space, v);
        if (!r) throw std::invalid_argument("failed to parse '" + s + "'");
        return v;
    }

    TEST(Radian1IsConvertedTo1) {
        auto v = parseString("1rad");
        CHECK_CLOSE(1, v, 1E-2);
    }

    TEST(Radian2IsConvertedTo2) {
        auto v = parseString("2rad");
        CHECK_CLOSE(2, v, 1E-2);
    }

    TEST(Degree1IsMultipliedByPiOver180) {
        auto v = parseString("1deg");
        CHECK_CLOSE(TMath::Pi()/180, v, 1E-7);
    }

    TEST(Degree180IsMultipliedByPiOver180) {
        auto v = parseString("180deg");
        CHECK_CLOSE(TMath::Pi(), v, 1E-7);
    }
}
