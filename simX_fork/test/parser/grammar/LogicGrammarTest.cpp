#include <unittest++/UnitTest++.h>
#include <memory>
#include <boost/type_traits/remove_cv.hpp>
#include <boost/variant/get.hpp>

#include "simX/parser/grammar/LogicGrammar.h"

using namespace simX::parser;
using namespace simX::parser::logic;
using std::string;

SUITE(LogicGrammarTest) {

    using expr = LogicExpression<std::string>;

    expr parse(const std::string& s) {
        auto f(std::begin(s)), l(std::end(s));
        LogicGrammar<std::string, decltype(f)> p;

        expr result;
        bool ok = qi::phrase_parse(f,l,p,qi::space,result);

        if (!ok)
                throw std::runtime_error("h");
        return result;
    }


    TEST(SingleLetterParses) {
        parse("A");
    }

    TEST(NamesCanHaveMultipleLetters) {
        auto r = parse("AB");
        CHECK_EQUAL(0, r.which()); // We only have a named var

        auto str = boost::get<std::string>(r);
        CHECK_EQUAL("AB", str);
    }

    TEST(NamesCanHaveNumbers) {
        auto r = parse("A1");
        CHECK_EQUAL(0, r.which()); // We only have a named var

        auto str = boost::get<std::string>(r);
        CHECK_EQUAL("A1", str);
    }

}
