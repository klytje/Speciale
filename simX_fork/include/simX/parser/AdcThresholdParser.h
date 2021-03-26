//
// Created by munk on 09-11-16.
//

#ifndef SIMX_ADCTHRESHOLDPARSER_H
#define SIMX_ADCTHRESHOLDPARSER_H

#include "simX/daq/AdcThresholds.h"

namespace simX {
    namespace parser {
        daq::AdcThresholds parseAdcThresholds(const std::string& input, size_t required);
        daq::AdcThresholds parseAdcThresholdsFromFile(const std::string& file, size_t required);
    }
}
#endif //SIMX_ADCTHRESHOLDPARSER_H
