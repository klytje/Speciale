//
// Created by munk on 25-06-15.
//

#ifndef SIMX_PARSER_PARTICLEPARSER_H
#define SIMX_PARSER_PARTICLEPARSER_H

#include "IonGrammar.h"
#include "EnergyGrammar.h"
#include "WidthGrammar.h"

#include <boost/spirit/include/qi.hpp>

namespace simX {
    namespace parser {
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        typedef std::string::const_iterator Iterator;

        struct ParticleSettings {
            boost::optional<double> ex;
            boost::optional<double> G0;
            bool noTrack;
        };

        struct Particle {
            AUSA::EnergyLoss::Ion ion;
            boost::optional<ParticleSettings> settings;
            bool noTrack;

            double ex() const {
                return settings.is_initialized() ? settings -> ex.get_value_or(0.) : 0.;
            }

            double G0() const {
                return settings.is_initialized() ? settings -> G0.get_value_or(0.) : 0.;
            }

            bool doTracking() const {
                return !noTrack && !(settings.is_initialized() && settings->noTrack);
            }
        };

        struct ParticleGrammar : qi::grammar<Iterator, Particle(), ascii::space_type> {
            ParticleGrammar();

        private:
            qi::rule<Iterator, Particle(), ascii::space_type> start;
            qi::rule<Iterator, ParticleSettings(), ascii::space_type> settings;
            qi::rule<Iterator, bool(), ascii::space_type> noTrack;

            EnergyGrammar energy;
            WidthGrammar width;
            IonGrammar ion;
        };
    }
}
#endif //SIMX_PARSER_PARTICLEPARSER_H
