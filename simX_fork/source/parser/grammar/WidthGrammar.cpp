//
// Created by munk on 26-06-15.
//

#include "simX/parser/grammar/WidthGrammar.h"

#include <ausa/constants/Constants.h>

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/foreach.hpp>

namespace {
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

        struct lazy_hbar_
        {
            template <typename X>
            struct result { typedef X type; };

            template <typename X>
            X operator()(X x) const
            {
                return AUSA::HBAR / x;
            }
        };
    }
}

simX::parser::WidthGrammar::WidthGrammar() : WidthGrammar::base_type(start) {
    using qi::double_;
    using qi::_1;
    using qi::_val;
    using qi::lit;
    using qi::eps;
    using ascii::no_case;

    using boost::phoenix::construct;
    using boost::phoenix::val;
    using qi::on_error;
    using qi::fail;


    boost::phoenix::function<lazy_pow_>   lazy_pow;
    boost::phoenix::function<lazy_hbar_>  lazy_hbar;

//    width %= double_[_val = _1] >> -prefix[_val *= lazy_pow(10., _1)] >> no_case[lit("s")] >> eps[_val = lazy_hbar(_val)];
//    start %= (energy | width);

    start %= double_[_val = _1] > -prefix[_val *= lazy_pow(10., _1)] > (no_case["s"] >> eps[_val = lazy_hbar(_val)] | no_case["eV"] >> eps[_val *= lazy_pow(10., -3)]);

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
