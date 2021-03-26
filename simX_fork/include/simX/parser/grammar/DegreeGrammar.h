//
// Created by munk on 25-06-15.
//

#ifndef SIMX_DEGREE_MULTIPLIER_H
#define SIMX_DEGREE_MULTIPLIER_H

#include <boost/spirit/include/qi.hpp>

namespace simX {
    namespace parser {
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        struct DegreeMultiplier : qi::symbols<char, double> {
            DegreeMultiplier();
        };


        typedef std::string::const_iterator Iterator;

        struct DegreeGrammar : qi::grammar<Iterator, double(), ascii::space_type>
        {
            DegreeGrammar();

            qi::rule<Iterator, double(), ascii::space_type> start;
            DegreeMultiplier multiplier;
        };
    }
}
#endif //SIMX_DEGREE_MULTIPLIER_H
