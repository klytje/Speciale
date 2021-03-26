//
// Created by munk on 6/19/16.
//

#ifndef SIMX_SAMPLER_INPUT_SAMPLER_H
#define SIMX_SAMPLER_INPUT_SAMPLER_H

#include "SampleSet.h"

#include <memory>

namespace simX {
    namespace samplers {

        /**
         * A specilization of SampleInput, which samples another SampleInput using AUSA::Stat::InverseTransformSampler.
         *
         * This class samples another finite SampleInput providing an infitite set of samples.
         */
        class SampleSampler : public SampleSet {
        public:
            using Inner = std::unique_ptr<SampleSet>;

            /**
             * Create of SamplerInputSampler from another sampler.
             *
             * @param inner Some SampleInput
             * @pre inner must be finite ie. inner.hasNext() must return false at some point.
             * @pre inner must have at least 1 element.
             */
            SampleSampler(Inner inner);

            virtual ~SampleSampler();

            virtual const std::vector<double> &getSample() override;

            virtual const double getWeight() override;

            virtual bool isSampled() const override;

            virtual void next() override;

            virtual bool hasNext() override;

            /**
             * Currently not possible
             */
            virtual void reset() override;

        private:
            struct Impl;
            std::unique_ptr<Impl> pimpl;
        };
    }
}
#endif //SIMX_SAMPLER_INPUT_SAMPLER_H
