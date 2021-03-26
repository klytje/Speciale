//
// Created by munk on 25-06-15.
//

#include "simX/parser/grammar/ParticleGrammar.h"

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(
        simX::parser::ParticleSettings,
        (boost::optional<double>, ex)
        (boost::optional<double>, G0)
        (bool, noTrack)
)

BOOST_FUSION_ADAPT_STRUCT(
        simX::parser::Particle,
        (AUSA::EnergyLoss::Ion, ion)
        (boost::optional<simX::parser::ParticleSettings>, settings)
        (bool, noTrack)
)


simX::parser::ParticleGrammar::ParticleGrammar() : ParticleGrammar::base_type(start, "particle") {
    using qi::double_;
    using qi::lit;
    using qi::attr;
    using qi::matches;
    using boost::spirit::ascii::no_case;


    noTrack %= matches[no_case["notrack"]];
    settings %= (no_case[lit("Ex")] >> lit(":") >> energy) ^
                (no_case[lit("G0")] >> lit(":") >> width) ^
                noTrack
                ;

    start = ion >> -(-lit("(") >> settings >> -lit(")")) >> noTrack;
}
