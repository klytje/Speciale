//
// Created by munk on 22-09-15.
//

#ifndef SIMX_PARSER_DETECTIONSYSTEMPARSER_H
#define SIMX_PARSER_DETECTIONSYSTEMPARSER_H

#include "simX/Detection/DetectionSystem.h"
#include "simX/Detection/Detector.h"

#include <string>
#include <typeindex>
#include <unordered_map>
#include <ausa/setup/Detector.h>
#include <ausa/setup/Setup.h>

#include <rapidjson/document.h>

namespace simX {
    namespace parser {
        class DetectionSystemParser {
        public:
            struct Result {
                std::string ausalibFile;
                detection::DetectionSystem system;
            };


            using DetectorBuilder = std::function<std::shared_ptr<simX::Detector>(AUSA::Detector&, rapidjson::Value&)>;

            DetectionSystemParser(bool calibrate);
            ~DetectionSystemParser();

            Result parseString(const std::string& input, std::string prefix = "");
            Result parseFile(const std::string& input);

            detection::DetectionSystem buildSimXSetup(AUSA::Setup& s);

            template <class Detector>
            void registerBuilder(DetectorBuilder b) {
                builders[std::type_index(typeid(Detector))] = b;
            }

        private:
            struct Impl;
            std::unique_ptr<Impl> pimpl;

            std::unordered_map<std::type_index, DetectorBuilder> builders;
        };
    }
}
#endif //SIMX_PARSER_DETECTIONSYSTEMPARSER_H
