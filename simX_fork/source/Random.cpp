//
// Created by munk on 29-06-15.
//

#include "simX/Random.h"
#include <TRandom3.h>

namespace {
    TRandom& getRandom() {
        static TRandom3 rand;
        return rand;
    }
}

double ::simX::rnd() {
    return getRandom().Rndm();
}

double ::simX::gaus() {
    return getRandom().Gaus(0,1);
}

void ::simX::setSeed(unsigned seed) {
    getRandom().SetSeed(seed);
}
