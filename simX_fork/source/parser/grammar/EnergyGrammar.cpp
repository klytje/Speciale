//
// Created by munk on 25-06-15.
//

#include "simX/parser/grammar/EnergyGrammar.h"

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <string>

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

simX::parser::EnergyGrammar::EnergyGrammar() : EnergyGrammar::base_type(start, "energy") {
    using qi::double_;
    using qi::_1;
    using qi::_val;
    using qi::lit;
    using ascii::no_case;

    using boost::phoenix::construct;
    using boost::phoenix::val;
    using qi::on_error;
    using qi::fail;

    boost::phoenix::function<lazy_pow_>   lazy_pow;

    start = double_[_val = _1*std::pow(10,-3)] > -prefix[_val *= lazy_pow(10., _1)] > no_case[lit("eV")];

    on_error<fail>
      (
       start
       , std::cerr
       << val("Error! Expecting ")
       << qi::_4                               // what failed?
       << val(" here: \"")
       << construct<std::string>(qi::_3, qi::_2)   // iterators to error-pos, end
       << val("\"")
       << std::endl
       );
}
