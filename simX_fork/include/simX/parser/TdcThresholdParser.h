//
// Created by munk on 09-11-16.
//

#ifndef SIMX_TDC_THRESHOLD_PARSER_H
#define SIMX_TDC_THRESHOLD_PARSER_H

#include "simX/daq/TdcThresholds.h"

namespace simX {
    namespace parser {
        daq::TdcThresholds parseTdcThresholds(const std::string& input, size_t required);
        daq::TdcThresholds parseTdcThresholdsFromFile(const std::string& file, size_t required);
    }
}


#endif //SIMX_TDC_THRESHOLD_PARSER_H
