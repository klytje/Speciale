#include <unittest++/UnitTest++.h>
#include <memory>
#include <boost/type_traits/remove_cv.hpp>
#include <boost/variant/get.hpp>

#include "simX/parser/grammar/LogicGrammar.h"
#include <boost/fusion/adapted/struct/adapt_struct.hpp>


using namespace simX::parser;
using namespace simX::parser::logic;
using std::string;

namespace {
    struct Indexed {
        std::string name;
        size_t reduction;
        size_t index;
        size_t count;
    };
}

BOOST_FUSION_ADAPT_STRUCT(
        Indexed,
(std::string, name)
(std::size_t, reduction)
)


SUITE(DownscaledGrammarTest) {

    using expr = LogicExpression<Indexed>;

    expr parse(const std::string& s) {
        auto f(std::begin(s)), l(std::end(s));
        DownscaledGrammar<Indexed, decltype(f)> p;

        expr result;
        bool ok = qi::phrase_parse(f,l,p,qi::space,result);

        if (!ok)
                throw std::runtime_error("h");
        return result;
    }


    TEST(DefaultReductionIs1) {
        auto r = parse("A");
        CHECK_EQUAL(0, r.which()); // We only have a named var

        auto i = boost::get<Indexed>(r);
        CHECK_EQUAL("A", i.name);
        CHECK_EQUAL(1, i.reduction);
    }

    TEST(If2InParenthesisThenReductionIs2) {
        auto r = parse("A(2)");
        CHECK_EQUAL(0, r.which()); // We only have a named var

        auto i = boost::get<Indexed>(r);
        CHECK_EQUAL("A", i.name);
        CHECK_EQUAL(2, i.reduction);
    }


    TEST(If200InParenthesisThenReductionIs200) {
        auto r = parse("A(200)");
        CHECK_EQUAL(0, r.which()); // We only have a named var

        auto i = boost::get<Indexed>(r);
        CHECK_EQUAL("A", i.name);
        CHECK_EQUAL(200, i.reduction);
    }


    TEST(OrOperatorIs) {
        auto r = parse("A(7) | B(999)");
        CHECK_EQUAL(4, r.which()); // We have OR op

        auto or_ = boost::get<binop<op_or, Indexed>>(r);
        // Two simple vars
        CHECK_EQUAL(0, or_.oper1.which());
        CHECK_EQUAL(0, or_.oper2.which());

        // Left is A
        auto left = boost::get<Indexed>(or_.oper1);
        CHECK_EQUAL("A", left.name);
        CHECK_EQUAL(7, left.reduction);

        // right is B
        auto right = boost::get<Indexed>(or_.oper2);
        CHECK_EQUAL("B", right.name);
        CHECK_EQUAL(999, right.reduction);
    }


    TEST(DoubleAndParsesToAnd) {
        auto r = parse("A(7) && B(999)");
        CHECK_EQUAL(2, r.which()); // We only have a named var

        auto and_ = boost::get<binop<op_and, Indexed>>(r);
        // Two simple vars
        CHECK_EQUAL(0, and_.oper1.which());
        CHECK_EQUAL(0, and_.oper2.which());

        // Left is A
        auto left = boost::get<Indexed>(and_.oper1);
        CHECK_EQUAL("A", left.name);
        CHECK_EQUAL(7, left.reduction);

        // right is B
        auto right = boost::get<Indexed>(and_.oper2);
        CHECK_EQUAL("B", right.name);
        CHECK_EQUAL(999, right.reduction);
    }
}
