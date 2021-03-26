//
// Created by munk on 13-12-15.
//

#ifndef SIMX_PARSER_BEAMPARSER_H
#define SIMX_PARSER_BEAMPARSER_H

#include "simX/Beam.h"

#include <string>
#include <map>
#include <memory>
#include <ausa/util/UnitParser.h>

class TF1;

namespace simX {
    namespace parser {
        class BeamParser {
        public:
            BeamParser();

            /**
             * A string -> string map with additional user options.
             */
            using Options = std::map<std::string, std::string>;


            using DistributionFactory = std::function<std::unique_ptr<TF1>(Options&, AUSA::UnitParser& parser)>;

            /**
             * Register a DistributionFactory for 1D distributions.
             */
            void registerFactory(std::string key, DistributionFactory f);

            /**
             * Register a DistributionFactory for 2D XY distributions.
             */
            void registerXYFactory(std::string key, DistributionFactory f);

            Beam parse(const std::string& input);
            Beam parseFile(const std::string& file);

        private:
            std::map<std::string, DistributionFactory> distFactories, xyFactories;
        };
    }
}
#endif //SIMX_PARSER_BEAMPARSER_H
