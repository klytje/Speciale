//
// Created by munk on 09-11-17.
//

#include "simX/samplers/MultiSampler.h"
#include <cassert>
#include <stdexcept>

using namespace std;
using namespace simX;
using namespace simX::samplers;

const vector<double>& MultiSampler::getSample() {
    assert(current < samples.size());
    return samples[current]->getSample();
}

const double MultiSampler::getWeight() {
    assert(current < samples.size());
    return samples[current]->getWeight();
}

bool MultiSampler::isSampled() const {
    return sampled;
}

void MultiSampler::next() {
    for (size_t i = _next; i < samples.size(); ++i) {
        if (samples[i]->hasNext()) {
            _next = current = i;
            samples[i]->next();
            break;
        }
    }
    while (_next < samples.size() && !samples[_next]->hasNext()) _next++;
}

bool MultiSampler::hasNext() {
    return _next < samples.size() ? samples[_next]->hasNext() : false;
}

void MultiSampler::reset() {
    current = _next = 0;
    for (auto& s : samples) s->reset();

    while (_next < samples.size() && !samples[_next]->hasNext()) _next++;

    assert(_next < samples.size());
}

MultiSampler::MultiSampler(vector<unique_ptr<SampleSet>>&& _samples)
        : samples(move(_samples)) {
    assert(!samples.empty());

    reset();

    sampled = samples[0]->isSampled();
    for (auto& s : samples) {
        if (s->isSampled() != sampled)
            throw std::invalid_argument("All SampleSets must either be sampled or not sampled!");
    }
}
