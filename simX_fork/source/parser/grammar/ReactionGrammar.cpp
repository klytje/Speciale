//
// Created by munk on 25-06-15.
//

#include "simX/parser/grammar/ReactionGrammar.h"

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/adapted/std_pair.hpp>

#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <string>

BOOST_FUSION_ADAPT_STRUCT(
        simX::parser::ReactionChainDescription,
        (simX::parser::Particle, beam)
        (boost::optional<simX::parser::Particle>, target)
        (boost::optional<simX::parser::ReactionDescription>, chain)
)

typedef boost::optional<std::pair<double, double>> Limits;
typedef std::pair<simX::parser::Particle, boost::optional<simX::parser::ReactionDescription>> SubReaction;
typedef boost::optional<std::pair<std::string, boost::optional<simX::parser::KeyValueMap>>> IdWithMap_t;

BOOST_FUSION_ADAPT_STRUCT(
        simX::parser::ReactionSettings,
        (Limits, phi)
        (Limits, theta)
        (IdWithMap_t, ad)
        (IdWithMap_t, gen)
        (IdWithMap_t, weight)
        (boost::optional<unsigned>, L)
)

BOOST_FUSION_ADAPT_STRUCT(
        simX::parser::ReactionDescription,
        (boost::optional<simX::parser::ReactionSettings>, settings)
                (std::vector<SubReaction>, products)
)


simX::parser::ReactionChainGrammar::ReactionChainGrammar() : ReactionChainGrammar::base_type(start){
    using qi::lit;
    using qi::double_;
    using qi::uint_;
    using qi::char_;
    using qi::lexeme;
    using boost::spirit::ascii::no_case;

using boost::phoenix::construct;
using boost::phoenix::val;
using qi::on_error;
using qi::fail;

    pair_ %= (-lit("[") >> double_ >> -lit(",") >> double_ >> -lit("]")); 
    beam %= no_case[lit("beam")] > ":" > particle;
    target %= no_case[lit("target")] > ":" > particle;
    id %= +qi::char_("a-zA-Z_0-9*");
    settings %= (no_case[lit("phi")]            > lit(":") > pair_) ^
                (no_case[lit("theta")]          > lit(":") > pair_) ^
                (no_case[lit("AD")]             > lit(":") > id >> -("(" > keyValue > ")")) ^
                (no_case[lit("generator")]      > lit(":") > id >> -("(" > keyValue > ")")) ^
                (no_case[lit("weight")]      > lit(":") > id >> -("(" > keyValue > ")")) ^
                (no_case[lit("L")]              >> lit(":") >> uint_);

    chain %= lit("->") > lit("{") >>
                                 -settings >> -(*(particle >> -chain)) >>
            lit("}");

    start %= beam ^ target ^ (-no_case[lit("chain")] >> -lit(":") >> chain);

    pair_.name("pair");
    beam.name("beam_particle");
    target.name("target_particle");
    id.name("id");
    particle.name("particle");

    
    on_error<fail>
      (
       start
       , std::cerr
       << val("Error! Expecting ")
       << qi::_4                               // what failed?
       << val(" here: \"")
       << construct<std::string>(qi::_3, qi::_2)   // iterators to error-pos, end
       << val("\"")
       << std::endl
       );
}
