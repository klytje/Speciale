//
// Created by kuhlwein on 1/13/20.
//

#ifndef JONAS_ROOTINTERFERENCEWRITER_H
#define JONAS_ROOTINTERFERENCEWRITER_H

#include "Decay.h"
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <memory>

struct RootInterferenceWriter {

    explicit RootInterferenceWriter(std::string output);

    virtual ~RootInterferenceWriter();

    void fill(const TriAlphaDecay& decay);

    TFile* f;
    TTree* t;

    double x, y, wU, p[9];
    TClonesArray* pOut;
    double ecm[3];
    std::vector<std::vector<double>> b;

};

#endif //JONAS_ROOTINTERFERENCEWRITER_H
