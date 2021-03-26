//
// Created by munk on 09-11-16.
//

#include <boost/algorithm/string.hpp>
#include <simX/daq/AdcThresholds.h>
#include "simX/parser/AdcThresholdParser.h"
#include "simX/parser/grammar/AdcThresholdsGrammar.h"
#include "ausa/util/FileUtil.h"

using namespace simX::daq;
using namespace simX::parser;

AdcThresholds simX::parser::parseAdcThresholds(const std::string &input, size_t required) {
    using boost::spirit::qi::phrase_parse;
    using boost::spirit::ascii::blank;

    std::pair<int, std::vector<double>> tmp;
    bool r = phrase_parse(input.begin(), input.end(), AdcThresholdsGrammar{}, blank, tmp);
    if (!r) throw std::invalid_argument("failed to parse '" + input + "'");

    std::vector<double> &thresholds = tmp.second;
    if (thresholds.size() != required)
        throw std::invalid_argument("Should have " + std::to_string(required) +
                                            " adc thresholds, but only have " + std::to_string(thresholds.size()));
    
    auto mul = pow(10, tmp.first-3);
    for (auto& entry : thresholds) entry *= mul;

    return AdcThresholds{thresholds};
}

AdcThresholds simX::parser::parseAdcThresholdsFromFile(const std::string &file, size_t required) {
    auto str = AUSA::asciiFileToString(file);
    return parseAdcThresholds(str, required);
}
