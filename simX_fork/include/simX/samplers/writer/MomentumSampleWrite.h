//
// Created by munk on 6/22/16.
//

#ifndef SIMX_MOMENTUMSAMPLEWRITE_H
#define SIMX_MOMENTUMSAMPLEWRITE_H

#include <vector>
#include <TVector3.h>

namespace simX {
    namespace samplers {
        namespace writer {
            class MomentumSampleWriter {
            public:
                virtual void addSample(const std::vector<TVector3>& sample, double w = 1) = 0;
            };
        }
    }
}
#endif //SIMX_MOMENTUMSAMPLEWRITE_H
