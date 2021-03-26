//
// Created by munk on 29-06-15.
//

#include "simX/parser/grammar/KeyValueGrammar.h"

#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/adapted/std_pair.hpp>

simX::parser::KeyValueGrammar::KeyValueGrammar() : KeyValueGrammar::base_type(start) {
    start = +pair;
    pair = key >> -('=' >> value);
    key = +qi::char_("a-zA-Z_0-9*");
    value %= qi::lexeme['"' > +(qi::char_ - '"') > '"'];
}