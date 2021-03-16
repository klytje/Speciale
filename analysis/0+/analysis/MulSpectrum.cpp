//
// Created by munk on 11-10-16.
//

#include "MulSpectrum.h"

using namespace AUSA::Event;

std::vector<bool, std::allocator<bool>> MulSpectrum::cut(const IdVector &ids) const {
    size_t n = 0;
    for (auto& i : ids) {
        auto& ion = i->getParticleType()->ion;
        n += ion.getA() == 4 && ion.getZ() == 2;
    }
    mul.Fill(n);

    return std::vector<bool>{ids.size(), true, std::vector<bool>::allocator_type()};
}

MulSpectrum::MulSpectrum() : AbstractIdCutter("MulSpecturm") {
    mul = TH1I("mul", "mul", 16, -0.5, 15.5);
}
