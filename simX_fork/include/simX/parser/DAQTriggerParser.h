//
// Created by munk on 13-12-15.
//

#ifndef SIMX_DAQTRIGGERPARSER_H
#define SIMX_DAQTRIGGERPARSER_H

#include "simX/daq/SimpleDAQ.h"
#include "simX/Detection/DetectionSystem.h"

namespace simX {
    namespace parser {
        daq::SimpleDAQ::TriggerFunction parseTriggerFunction(const std::string& s,
                                                             const detection::DetectionSystem& detectionSystem);

        daq::SimpleDAQ::TriggerFunction parseTriggerFunction(const std::string& s,
                                                             const std::vector<std::string>& names);
    }
}
#endif //SIMX_DAQTRIGGERPARSER_H
