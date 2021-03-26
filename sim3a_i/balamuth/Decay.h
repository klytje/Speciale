//
// Created by munk on 08-11-17.
//

#ifndef JONAS_DECAY_H
#define JONAS_DECAY_H

#include <TLorentzVector.h>

struct TriAlphaDecay {
    /** Total weight */
    double w;

    /** Phase space weight */
    double wU;

    /** Balamuth weight */
    double wB;

    /** CM energies */
    double eCM[3];

    /** CM Lorentz vectors */
    const TLorentzVector* p[3];

    /** Sort order */
    int order[3];

    /** Weight factors for interference */
    std::vector<std::vector<double>> factors;
};

#endif //JONAS_DECAY_H
