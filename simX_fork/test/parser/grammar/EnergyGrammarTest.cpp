#include <unittest++/UnitTest++.h>
#include <memory>

#include "simX/parser/grammar/EnergyGrammar.h"

using namespace simX::parser;
using std::string;

SUITE(EnergyParserTest) {

    double parseString(const string& s) {
        using boost::spirit::qi::phrase_parse;
        using boost::spirit::ascii::space;

        double v;
        bool r = phrase_parse(s.begin(), s.end(), EnergyGrammar{}, space, v);
        if (!r) throw std::invalid_argument("failed to parse '" + s + "'");
        return v;
    }

    TEST(KeVIsBaseUnit) {
        auto v = parseString("2keV");
        CHECK_CLOSE(2, v, 1E-2);
    }

    TEST(MeVIsFactor3Above) {
        auto v = parseString("2MeV");
        CHECK_CLOSE(2E3, v, 1E-2);
    }

    TEST(meVIsFactor3Below) {
        auto v = parseString("2meV");
        CHECK_CLOSE(2E-6, v, 1E-7);
    }

    TEST(AllowedToInsertWhitespace) {
        auto v = parseString("2 meV");
        CHECK_CLOSE(2E-6, v, 1E-7);
    }

    TEST(UnitWithoutSIPrefixIsParsed) {
        auto v = parseString("2 eV");
        CHECK_CLOSE(2E-3, v, 1E-7);
    }
}
