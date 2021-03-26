//
// Created by munk on 20-11-15.
//

#include "simX/parser/IonParser.h"

namespace simX {
    namespace parser {
        Ion parseIon(const std::string& s) {
            using boost::spirit::qi::phrase_parse;
            using boost::spirit::ascii::space;

            Ion result;
            bool r = phrase_parse(s.begin(), s.end(), IonGrammar{}, space, result);
            if (!r) throw std::invalid_argument("failed to parse '" + s + "'");
            return result;
        }
    }
}