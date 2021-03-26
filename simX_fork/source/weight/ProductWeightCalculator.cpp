//
// Created by munk on 01-12-15.
//

#include "simX/weight/ProductWeightCalculator.h"

using namespace std;
using namespace simX;

ProductWeightCalculator::ProductWeightCalculator(std::unique_ptr<WeightCalculator> w1,
                                                 std::unique_ptr<WeightCalculator> w2)
    : w1(move(w1)), w2(move(w2))
{

}

double ProductWeightCalculator::getWeight() const {
    return w1->getWeight()*w2->getWeight();
}

const WeightCalculator& ProductWeightCalculator::getFunction1() const {
    return *w1;
}

const WeightCalculator& ProductWeightCalculator::getFunction2() const {
    return *w2;
}
