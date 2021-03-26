//
// Created by munk on 08-11-17.
//

#ifndef JONAS_ROOTWRITER_H
#define JONAS_ROOTWRITER_H

#include "Decay.h"
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <memory>

struct RootWriter {

    explicit RootWriter(std::string output);

    virtual ~RootWriter();

    void fill(const TriAlphaDecay& decay);

    std::unique_ptr<TFile> f;
    TTree* t;

    double w, x, y, wB, wU;
    TClonesArray* pOut;
    double ecm[3];
};

#endif //JONAS_ROOTWRITER_H
