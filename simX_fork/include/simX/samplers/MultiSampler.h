//
// Created by munk on 09-11-17.
//

#ifndef SIMX_MULTISAMPLER_H
#define SIMX_MULTISAMPLER_H

#include "SampleSet.h"
#include <cstdint>
#include <memory>

namespace simX {
    namespace samplers {
        class MultiSampler : public SampleSet {
        public:

            using SamplePtr = std::unique_ptr<SampleSet>;
            using SampleCollection = std::vector<SamplePtr>;

            MultiSampler(SampleCollection&& samples);

            const std::vector<double>& getSample() override;

            const double getWeight() override;

            bool isSampled() const override;

            void next() override;

            bool hasNext() override;

            void reset() override;

        private:
            const SampleCollection samples;
            size_t current, _next;
            bool sampled;

        };
    }
}



#endif //SIMX_MULTISAMPLER_H
