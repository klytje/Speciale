//
// Created by munk on 25-06-15.
//

#ifndef SIMX_ENERGYPARSER_H
#define SIMX_ENERGYPARSER_H

#include "SiPrefix.h"
#include <boost/spirit/include/qi.hpp>

namespace simX {
    namespace parser {
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        typedef std::string::const_iterator Iterator;

        struct EnergyGrammar : qi::grammar<Iterator, double(), ascii::space_type>
        {
            EnergyGrammar();

            qi::rule<Iterator, double(), ascii::space_type> start;
            SiPrefix prefix;
        };
    }
}
#endif //SIMX_ENERGYPARSER_H
