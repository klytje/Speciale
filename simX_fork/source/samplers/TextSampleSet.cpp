//
// Created by munk on 6/16/16.
//

#include "simX/samplers/TextSampleSet.h"

#include <stdexcept>
#include <cstring>
#include <sstream>

using namespace std;

simX::samplers::TextSampleSet::TextSampleSet(std::string path_, size_t elements)
    : stream(path_), elements(elements), path(path_)
{
    if (stream.bad()) throw std::runtime_error("Failed to open file for text sampler due to: " + std::string(strerror(errno)));

    sample.resize(elements);
    reset();
    next();
    reset();
}

const std::vector<double> &simX::samplers::TextSampleSet::getSample() {
    return sample;
}

const double simX::samplers::TextSampleSet::getWeight() {
    return weight;
}

void simX::samplers::TextSampleSet::next() {
    if (!read) readSample();

    current++;
    bool flag = false;

    stringstream ss(line);
    for (size_t i = 0; i < elements; ++i) {
        flag |= !(ss >> sample[i]);
    }

    if (flag)
        throw std::runtime_error("Failed to parse line " + to_string(current) + " of " + path);

    sampled = !(ss >> weight);
    read = false;
}

bool simX::samplers::TextSampleSet::isSampled() const {
    return sampled;
}

void simX::samplers::TextSampleSet::readSample() {
    if (!getline(stream, line)) line = "";
}

bool simX::samplers::TextSampleSet::hasNext() {
    if (!read) {
        readSample();
        read = true;
    }

    return !line.empty();
}

void simX::samplers::TextSampleSet::reset() {
    stream.seekg(0, std::ios::beg);
    read = false;
    current = -1;
}


