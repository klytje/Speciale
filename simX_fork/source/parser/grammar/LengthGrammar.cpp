//
// Created by munk on 25-06-15.
//

#include "simX/parser/grammar/LengthGrammar.h"

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

namespace {
    struct lazy_pow_
    {
        template <typename X, typename Y>
        struct result { typedef X type; };

        template <typename X, typename Y>
        X operator()(X x, Y y) const
        {
            return std::pow(x, y);
        }
    };
}

simX::parser::LengthGrammar::LengthGrammar() : LengthGrammar::base_type(start) {
    using qi::double_;
    using qi::_1;
    using qi::_val;
    using qi::lit;
    using ascii::no_case;

    boost::phoenix::function<lazy_pow_>   lazy_pow;

    start = double_[_val = _1*std::pow(10, 3)] >> -prefix[_val *= lazy_pow(10., _1)] >> no_case[lit("m")];
}
