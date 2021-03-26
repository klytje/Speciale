//
// Created by munk on 23-09-15.
//

#ifndef SIMX_JSONOUTPUTPARSER_H
#define SIMX_JSONOUTPUTPARSER_H

#include <unordered_map>
#include <string>

namespace simX {
    namespace parser {
        class JsonOutputParser {
        public:
            struct Mapping {
                std::string mul, seg, adc, tdc;
            };

            struct DoubleMapping {
                Mapping front, back;
            };

            JsonOutputParser(const std::string& input);

            Mapping& getMapping(const std::string& s);
            bool hasMapping(const std::string& s);

            DoubleMapping& getDoubleMapping(const std::string& s);
            bool hasDoubleMapping(const std::string& s);

        private:
            std::unordered_map<std::string, Mapping> singleMap;
            std::unordered_map<std::string, DoubleMapping> doubleMap;
        };
    }
}
#endif //SIMX_JSONOUTPUTPARSER_H
