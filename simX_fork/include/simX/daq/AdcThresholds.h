//
// Created by munk on 09-11-16.
//

#ifndef SIMX_ADCTHRESHOLDS_H
#define SIMX_ADCTHRESHOLDS_H

#include <vector>

namespace simX {
    namespace daq {
        class AdcThresholds {
        public:
            AdcThresholds(const std::vector<double> &thresholds) : thresholds(thresholds) {}

            bool isTriggering(int channel, double E) const {
                return thresholds[channel] < E;
            }

            const std::vector<double>& getThresholds() const {
                return thresholds;
            }

        private:
            std::vector<double> thresholds;
        };
    }
}
#endif //SIMX_ADCTHRESHOLDS_H
