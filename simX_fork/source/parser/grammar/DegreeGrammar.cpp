//
// Created by munk on 25-06-15.
//

#include "simX/parser/grammar/DegreeGrammar.h"
#include <TMath.h>

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

simX::parser::DegreeMultiplier::DegreeMultiplier() {
    add
            ("deg", TMath::Pi()/180)
            ("rad", 1)
            ;
}

simX::parser::DegreeGrammar::DegreeGrammar() : DegreeGrammar::base_type(start) {
    using qi::double_;
    using qi::_1;
    using qi::_val;
    using qi::lit;
    using ascii::no_case;

    start = double_[_val = _1] >> multiplier[_val *= _1];
}
