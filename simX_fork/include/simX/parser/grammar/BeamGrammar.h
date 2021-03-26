//
// Created by munk on 25-06-15.
//

#ifndef SIMX_PARSER_GRAMMAR_BEAM_GRAMMAR_H
#define SIMX_PARSER_GRAMMAR_BEAM_GRAMMAR_H

#include "EnergyGrammar.h"
#include "LengthGrammar.h"
#include "DegreeGrammar.h"
#include "KeyValueGrammar.h"

#include <boost/optional.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>
#include <vector>

namespace simX {
    namespace parser {
        
        struct BeamPrototype {
            boost::optional<std::pair<double, bool>> energy;
            boost::optional<boost::tuple<double, double, double>> center;
            boost::optional<double> theta, phi;

            boost::optional<std::pair<std::string, boost::optional<KeyValueMap>>> thetaDist, energyDist, xyDist, phiDist;
        };
        
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        typedef std::string::const_iterator Iterator;

        struct BeamGrammar : qi::grammar<Iterator, BeamPrototype(), ascii::space_type> {
            BeamGrammar();

        private:
            qi::rule<Iterator, BeamPrototype(), ascii::space_type> start;

            EnergyGrammar energy;
            LengthGrammar length;
            DegreeGrammar degree;
            KeyValueGrammar keyValue;
            qi::rule<Iterator, std::string()> id;
        };
    }
}

#endif //SIMX_PARSER_GRAMMAR_BEAM_GRAMMAR_H
