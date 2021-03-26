//
// Created by munk on 08-11-17.
//

#ifndef JONAS_TEXTWRITER_H
#define JONAS_TEXTWRITER_H

#include "Decay.h"
#include <fstream>
#include <iostream>

struct TextWriter {

    explicit TextWriter(std::string output);

    virtual ~TextWriter();

    void fill(const TriAlphaDecay& decay);

    std::ofstream f;
};

#endif //JONAS_TEXTWRITER_H
