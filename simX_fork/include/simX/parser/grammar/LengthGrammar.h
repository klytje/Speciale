//
// Created by munk on 25-06-15.
//

#ifndef SIMX_LENGTHPARSER_H
#define SIMX_LENGTHPARSER_H

#include "SiPrefix.h"
#include <boost/spirit/include/qi.hpp>

namespace simX {
    namespace parser {
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        typedef std::string::const_iterator Iterator;

        struct LengthGrammar : qi::grammar<Iterator, double(), ascii::space_type>
        {
            LengthGrammar();

            qi::rule<Iterator, double(), ascii::space_type> start;
            SiPrefix prefix;
        };
    }
}
#endif //SIMX_LENGTHPARSER_H
