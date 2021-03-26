//
// Created by munk on 07-08-15.
//

#ifndef SIMX_CONFIGPARSER_H
#define SIMX_CONFIGPARSER_H

#include "simX/Config.h"
#include "ausa/util/UnitParser.h"

#include <string>
#include <map>
#include <rapidjson/document.h>

class TF1;

namespace simX {
    namespace parser {
        class ConfigParser {
        public:
            /**
             * A string -> string map with additional user options.
             */
            using Options = rapidjson::Value;

            /**
             * Factory function that will create a Config::PropagatorFactory based on user string input.
             */
            using PropagatorFactory = std::function<Config::PropagatorFactory(Options&)>;

            ConfigParser();

            Config parse(const std::string& input);
            Config parseFile(const std::string& file);

            void registerFactory(std::string key, PropagatorFactory f);

        private:
            std::map<std::string, PropagatorFactory> propFactories;
        };
    }
}

#endif //SIMX_CONFIGPARSER_H
