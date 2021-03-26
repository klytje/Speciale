//
// Created by munk on 6/16/16.
//

#ifndef SIMX_SAMPLEINPUT_H
#define SIMX_SAMPLEINPUT_H

#include <vector>

namespace simX {
    namespace samplers {

        /**
         * Abstract class that represent either a discrete set of samples or a set of samples with weights.
         */
        class SampleSet {
        public:

            /**
             * Get a const ref to the current sample.
             *
             * @pre This is undefined before calling next().
             * @return Const ref to the current sample
             */
            virtual const std::vector<double>& getSample() = 0;

            /**
             * Get the weight of the current sample.
             * @pre This is undefined before calling next().
             * @pre This is undefined is isSampled() returns true.
             * @return Weight of the current sample
             */
            virtual const double getWeight() = 0;

            /**
             * Return whether the dataset have been sampled or is provided with weights.
             * @return False if the datasets has weights, true else.
             */
            virtual bool isSampled() const = 0;

            /**
             * Increment to the next element.
             *
             * @pre hasNext() must return true.
             */
            virtual void next() = 0;

            /**
             * Return whether the input contains more samples.
             * @return True if the input has more samples.
             */
            virtual bool hasNext() = 0;

            /**
             * Reset the input to its initial state.
             *
             * @post getWeight() and getSample() is undefined.
             */
            virtual void reset() = 0;
        };
    }
}
#endif //SIMX_SAMPLEINPUT_H
