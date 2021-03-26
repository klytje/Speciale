#include <unittest++/UnitTest++.h>
#include <memory>

#include "simX/parser/grammar/WidthGrammar.h"

#include <ausa/constants/Constants.h>

using namespace simX::parser;
using std::string;

SUITE(WidthParserTest) {

        double parseString(const string& s) {
            using boost::spirit::qi::phrase_parse;
            using boost::spirit::ascii::space;

            double v;
            bool r = phrase_parse(s.begin(), s.end(), WidthGrammar{}, space, v);
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
            CHECK_CLOSE(2E-6, v, 1E-10);
        }

        TEST(AllowedToInsertWhitespace) {
            auto v = parseString("2 meV");
            CHECK_CLOSE(2E-6, v, 1E-10);
        }

	TEST(AllowedToInsert_eV) {
	  auto v = parseString("2 eV");
	  CHECK_CLOSE(2E-3, v, 1E-10);
        }

    TEST(FsParses) {
        auto v = parseString("2 fs");

        CHECK(v == v);
    }

    TEST(WidthInTimeIsCorrectlyConverted) {
        auto expect = AUSA::HBAR / 2E-3;
        auto v = parseString("2 ms");

        CHECK_CLOSE(expect, v, 1E-6);
    }
}
