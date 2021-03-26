//
// Created by munk on 25-06-15.
//

#ifndef SIMX_PARSER_GRAMMAR_REACTION_GRAMMAR_H
#define SIMX_PARSER_GRAMMAR_REACTION_GRAMMAR_H

#include "simX/Particle.h"
#include "ParticleGrammar.h"
#include "IonGrammar.h"
#include "KeyValueGrammar.h"

#include <boost/optional.hpp>
#include <boost/spirit/include/qi.hpp>
#include <string>
#include <vector>

namespace simX {
    namespace parser {
        using Limits = std::pair<double, double>;

        struct ReactionSettings {
            boost::optional<Limits> phi;
            boost::optional<Limits> theta;
            boost::optional<std::pair<std::string, boost::optional<KeyValueMap>>> ad, gen, weight;
            boost::optional<unsigned> L;
        };

        struct ReactionDescription {
            typedef std::pair<Particle, boost::optional<ReactionDescription>> SubReaction;

            boost::optional<ReactionSettings> settings;
            std::vector<SubReaction> products;
        };

        struct ReactionChainDescription {

            Particle beam;
            boost::optional<Particle> target;
            boost::optional<ReactionDescription> chain;
        };


        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        typedef std::string::const_iterator Iterator;

        struct ReactionChainGrammar : qi::grammar<Iterator, ReactionChainDescription(), ascii::space_type> {
            ReactionChainGrammar();

        private:
            qi::rule<Iterator, ReactionChainDescription(), ascii::space_type> start;
            qi::rule<Iterator, Particle(), ascii::space_type> beam;
            qi::rule<Iterator, Particle(), ascii::space_type> target;
            qi::rule<Iterator, ReactionDescription(), ascii::space_type> chain;
            qi::rule<Iterator, ReactionSettings(), ascii::space_type> settings;
            qi::rule<Iterator, Limits(), ascii::space_type> pair_;
            qi::rule<Iterator, std::string()> id;

            IonGrammar ion;
            ParticleGrammar particle;
            KeyValueGrammar keyValue;
        };
    }
}

#endif //SIMX_PARSER_GRAMMAR_REACTION_GRAMMAR_H
