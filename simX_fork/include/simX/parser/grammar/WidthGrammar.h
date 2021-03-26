//
// Created by munk on 25-06-15.
//

#ifndef SIMX_WIDTHGRAMMAR_H
#define SIMX_WIDTHGRAMMAR_H

#include "SiPrefix.h"
#include "EnergyGrammar.h"
#include <boost/spirit/include/qi.hpp>

namespace simX {
    namespace parser {
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        typedef std::string::const_iterator Iterator;

        struct WidthGrammar : qi::grammar<Iterator, double(), ascii::space_type>
        {
            WidthGrammar();

            qi::rule<Iterator, double(), ascii::space_type> start;
            qi::rule<Iterator, double(), ascii::space_type> width;
            SiPrefix prefix;
            EnergyGrammar energy;
        };
    }
}
#endif //SIMX_WIDTHGRAMMAR_H
