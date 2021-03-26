//
// Created by munk on 25-06-15.
//

#ifndef SIMX_PARSER_GRAMMAR_ION_GRAMMAR_H
#define SIMX_PARSER_GRAMMAR_ION_GRAMMAR_H

#include <ausa/eloss/Ion.h>

#include <boost/spirit/include/qi.hpp>

namespace simX {
    namespace parser {
        namespace qi = boost::spirit::qi;
        typedef std::string::const_iterator Iterator;
        typedef AUSA::EnergyLoss::Ion Ion;

        struct ElementSymbolTable : qi::symbols<char, unsigned> {
            ElementSymbolTable();
        };

        struct PredefinedIons : qi::symbols<char, Ion> {
            PredefinedIons();
        };


        struct AzeIsotopeGrammar : qi::grammar<Iterator, Ion()>
        {
            AzeIsotopeGrammar();

            qi::rule<Iterator, AUSA::EnergyLoss::Ion()> start;
            ElementSymbolTable element_name;
        };

        struct IonGrammar : qi::grammar<Iterator, Ion()>
        {
            IonGrammar();

            qi::rule<Iterator, AUSA::EnergyLoss::Ion()> start;
            PredefinedIons predefined;
            AzeIsotopeGrammar aze;
        };
    }
}
#endif //SIMX_PARSER_GRAMMAR_ION_GRAMMAR_H
