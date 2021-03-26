//
// Created by munk on 23-09-15.
//

#include "simX/parser/JsonOutputParser.h"
#include "ausa/json/JSONUtil.h"

#include <iostream>

using namespace std;
using namespace simX;
using namespace simX::parser;
using namespace AUSA::JSON;

namespace {
    using Mapping = JsonOutputParser::Mapping;

    Mapping parseMapping(const rapidjson::Value& v) {
        auto prefix = v.HasMember("prefix") ? readString(v, "prefix") : "";

        std::string b[] = {"multiplicity", "segment", "adc", "tdc"};
        std::string out[4];
        for (size_t i = 0; i < 4; ++i) {
            if (v.HasMember(b[i].c_str())) {
                auto s = readString(v, b[i].c_str());
                if (s != "-DEAD-") out[i] = prefix + s;
                else out[i] = "-DEAD-";
            }
            else out[i] = "-DEAD-";
        }

        return Mapping{out[0], out[1], out[2], out[3]};
    }
}

JsonOutputParser::JsonOutputParser(const std::string& input) {
    rapidjson::Document doc;
    tryParse(input, doc);

    auto& detectors = readArray(doc, "detectors");
    
    for (auto i = 0U, n = detectors.Size(); i < n; i++) {
        auto& d = detectors[i];
        auto name = readString(d, "name");

        if (d.HasMember("frontMapping")) doubleMap[name] = {parseMapping(d["frontMapping"]), parseMapping(d["backMapping"])};
        else if (d.HasMember("mapping")) singleMap[name] = parseMapping(d["mapping"]);
        else
            cerr << "No mapping entry found for " << name << endl;
    }
    
}

JsonOutputParser::Mapping& JsonOutputParser::getMapping(const std::string& s) {
    return singleMap[s];
}

JsonOutputParser::DoubleMapping& JsonOutputParser::getDoubleMapping(const std::string& s) {
    return doubleMap[s];
}

bool JsonOutputParser::hasMapping(const std::string& s) {
    return singleMap.count(s) != 0;
}

bool JsonOutputParser::hasDoubleMapping(const std::string& s) {
    return doubleMap.count(s) != 0;
}
