//
// Created by munk on 6/22/16.
//

#include "simX/samplers/writer/TextMomentumWriter.h"

using namespace simX::samplers::writer;

TextMomentumWriter::TextMomentumWriter(std::string path, bool weight)
    : stream(path), weight(weight)
{

}

void TextMomentumWriter::addSample(const std::vector<TVector3> &sample, double w) {
    for (auto& v : sample) {
        for (int i = 0; i < 3; ++i) {
            stream << v[i] << " ";
        }
    }
    if (weight) stream << w;
    stream << std::endl;
}


