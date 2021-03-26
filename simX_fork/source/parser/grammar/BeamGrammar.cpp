//
// Created by munk on 18-11-15.
//

#include "simX/parser/grammar/BeamGrammar.h"

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/boost_tuple.hpp>

#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <string>

using __Pair = boost::tuple<double, double, double>;
using __EPair = std::pair<double, bool>;
using __Map = boost::optional<std::pair<std::string, boost::optional<simX::parser::KeyValueMap>>>;

BOOST_FUSION_ADAPT_STRUCT(
        simX::parser::BeamPrototype,
        (__Map, thetaDist)
        (__Map, phiDist)
	    (__Map, energyDist)
	    (__Map, xyDist)
        (boost::optional<__EPair>, energy)
        (boost::optional<__Pair>, center)
        (boost::optional<double>, theta)
        (boost::optional<double>, phi)
)

simX::parser::BeamGrammar::BeamGrammar() : BeamGrammar::base_type(start, "beam") {
    using qi::double_;
    using qi::lit;
    using qi::attr;
    using qi::matches;
    using boost::spirit::ascii::no_case;
    using boost::spirit::eoi;

using namespace qi::labels;

using boost::phoenix::construct;
using boost::phoenix::val;
using qi::on_error;
using qi::fail;

    start %=
            (
                    (no_case["dist_theta"]             > ":" >> id >> -("(" > keyValue > ")")) ^
                    (no_case["dist_phi"]             > ":" >> id >> -("(" > keyValue > ")")) ^
		            (no_case["dist_energy"]             > ":" >> id >> -("(" > keyValue > ")")) ^
		            (no_case["dist_xy"]             > ":" >> id >> -("(" > keyValue > ")")) ^
                    (no_case["energy"] > ':' > energy >> matches["/A"]) ^
                    (no_case["center"] > ':' > '(' > length > ',' > length > ',' > length > ')') ^
                    (no_case["theta"] > ':' > degree) ^
                    (no_case["phi"] > ':' > degree)
            ) >> eoi;

    id %= +qi::char_("a-zA-Z_0-9");

    on_error<fail>
      (
       start
       , std::cerr
       << val("Error! Expecting ")
       << qi::labels::_4                               // what failed?
       << val(" here: \"")
       << construct<std::string>(qi::labels::_3, qi::labels::_2)   // iterators to error-pos, end
       << val("\"")
       << std::endl
       );
}
