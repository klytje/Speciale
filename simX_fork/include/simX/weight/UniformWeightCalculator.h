//
// Created by munk on 16-12-16.
//

#ifndef SIMX_UNIFORMWEIGHTCALCULATOR_H
#define SIMX_UNIFORMWEIGHTCALCULATOR_H

#include "WeightCalculator.h"

namespace simX {

    /**
    * A uniform weight of 1
    */
    class UniformWeightCalculator : public WeightCalculator {
    public:
        UniformWeightCalculator() = default;
        virtual ~UniformWeightCalculator() = default;

        /**
        * Returns weight
        */
        virtual double getWeight() const override {
            return 1;
        }
    };
}
#endif //SIMX_UNIFORMWEIGHTCALCULATOR_H
