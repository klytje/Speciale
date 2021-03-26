//
// Created by munk on 08-11-17.
//

#include "TextWriter.h"

TextWriter::TextWriter(std::string output) {
    f.open(output);
    if (f.bad()) {
        std::cerr << "Could not open " << output << std::endl;
        exit(2);
    }
}

TextWriter::~TextWriter() {
    f.flush();
    f.close();
}

void TextWriter::fill(const TriAlphaDecay& decay) {
    for (auto ptr : decay.p) {
        auto& p = *ptr;
        for (size_t j = 0; j < 3; ++j) {
            f << p[j] << " ";
        }
    }
    f << decay.w << std::endl;
}