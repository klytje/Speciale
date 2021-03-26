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

        struct AdcThresholdsGrammar : qi::grammar<Iterator, std::pair<int, std::vector<double>>(), ascii::blank_type>
        {
            AdcThresholdsGrammar();

            qi::rule<Iterator, std::pair<int, std::vector<double>>(), ascii::blank_type> start;
            qi::rule<Iterator, std::vector<double>(), ascii::blank_type> sub;

            SiPrefix prefix;
        };
    }
}
#endif //SIMX_ENERGYPARSER_H
