//
// Created by munk on 6/19/16.
//

#ifndef SIMX_DUMMYSAMPLER_H
#define SIMX_DUMMYSAMPLER_H

#include "simX/samplers/SampleSet.h"

class DummySampler : public simX::samplers::SampleSet {
public:
    DummySampler(std::initializer_list<std::initializer_list<double>> list) {
        for (auto& l : list) {
            input.emplace_back(l);
        }
        reset();
    }

    virtual const std::vector<double> &getSample() override {
        return input[current];
    }

    virtual const double getWeight() override {
        return !weight.empty() ? weight[current] : 0;
    }

    virtual void next() override {
        current++;
    }

    virtual bool hasNext() override {
        return current+1 < input.size();
    }

    virtual bool isSampled() const override {
        return weight.empty();
    }

    virtual void reset() override {
        current = -1;
    }

    size_t current;
    std::vector<std::vector<double>> input;
    std::vector<double> weight;
};
#endif //SIMX_DUMMYSAMPLER_H
