//
// Created by munk on 01-12-15.
//

#ifndef SIMX_PRODUCTWEIGHT_FUNCTION_H
#define SIMX_PRODUCTWEIGHT_FUNCTION_H

#include "WeightCalculator.h"

#include <memory>

namespace simX {
    class ProductWeightCalculator : public WeightCalculator {
    public:
        ProductWeightCalculator(std::unique_ptr<WeightCalculator> w1, std::unique_ptr<WeightCalculator> w2);

        virtual double getWeight() const;

        const WeightCalculator& getFunction1() const;
        const WeightCalculator& getFunction2() const;

    private:
        std::unique_ptr<WeightCalculator> w1, w2;
    };
}
#endif //SIMX_PRODUCTWEIGHT_FUNCTION_H
