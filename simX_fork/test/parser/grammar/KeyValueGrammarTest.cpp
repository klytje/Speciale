#include <unittest++/UnitTest++.h>
#include <memory>

#include "simX/parser/grammar/KeyValueGrammar.h"

using namespace simX::parser;
using std::string;

SUITE(KeyValueGrammarTest) {

    KeyValueMap parseString(const string& s) {
        using boost::spirit::qi::phrase_parse;
        using boost::spirit::ascii::space;

        KeyValueMap map;
        bool r = phrase_parse(s.begin(), s.end(), KeyValueGrammar{}, space, map);
        if (!r) throw std::invalid_argument("failed to parse '" + s + "'");
        return map;
    }

    TEST(Sanity) {
        CHECK(true);
    }

    TEST(TestIs2) {
        auto m = parseString("Test=\"2\"");

        CHECK_EQUAL("2", m["Test"]);
    }

    TEST(HestIs6) {
        auto m = parseString("Test=\"2\" Hest=\"6\"");

        CHECK_EQUAL("2", m["Test"]);
        CHECK_EQUAL("6", m["Hest"]);
    }

    TEST(NoKeyIsTakenAsEmptyString) {
        auto m = parseString("Test");

        CHECK_EQUAL(1, m.count("Test"));
    }
}
