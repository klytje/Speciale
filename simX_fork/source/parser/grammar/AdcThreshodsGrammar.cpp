//
// Created by munk on 25-06-15.
//

#include "simX/parser/grammar/AdcThresholdsGrammar.h"

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <utility>
#include <vector>
#include <string>

simX::parser::AdcThresholdsGrammar::AdcThresholdsGrammar() : AdcThresholdsGrammar::base_type(start, "AdcThreshold") {
    using qi::double_;
    using qi::_1;
    using qi::_a;
    using qi::_val;
    using qi::lit;
    using ascii::no_case;

    using boost::phoenix::construct;
    using boost::phoenix::val;
    using qi::on_error;
    using qi::fail;


    start %= qi::lexeme["unit:"] > (prefix | qi::attr(0)) > qi::lexeme["eV"] > '\n' > sub;
    sub %= double_ % '\n';

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
