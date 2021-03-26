//
// Created by munk on 6/19/16.
//

#include <limits>
#include <iostream>
#include <ausa/stat/InverseTransformSampler.h>
#include <ausa/util/memory>
#include "simX/samplers/SampleSampler.h"

using namespace simX::samplers;
using namespace AUSA::Stat;
using namespace std;

namespace {
    struct Element {
        std::vector<double> value;
    };

    struct WeightCalc {
        double weight(const Element& t, size_t i) {
            return w[i];
        }

        const std::vector<double> w;
    };
}

using Sampler = InverseTransformSampler<Element, TRandom3Generator, WeightCalc>;


struct SampleSampler::Impl {

    Impl(SampleSampler::Inner inner) {
        std::vector<Element> elements;
        std::vector<double> weights;

        while (inner->hasNext()) {
            inner->next();
            weights.emplace_back(inner->getWeight());
            elements.emplace_back(Element{inner->getSample()});
        }

        sampler = std::make_unique<Sampler>(move(elements), TRandom3Generator{}, WeightCalc{weights});
    }

    void next() {
        current = &(sampler->sample());
    }

    std::unique_ptr<Sampler> sampler;
    const Element* current;
};

SampleSampler::SampleSampler(SampleSampler::Inner inner)
    : pimpl(std::make_unique<Impl>(move(inner)))
{
    // See pimpl
}


const std::vector<double> &SampleSampler::getSample() {
    return pimpl->current->value;
}

const double SampleSampler::getWeight() {
    return std::numeric_limits<double>::quiet_NaN();
}

bool SampleSampler::isSampled() const {
    return true;
}

void SampleSampler::next() {
    pimpl->next();
}

bool SampleSampler::hasNext() {
    return true;
}

void SampleSampler::reset() {

}

SampleSampler::~SampleSampler() {
    // Needs to be here. Using PIMPL.
}

