//
// Created by munk on 08-11-17.
//

#include "RootWriter.h"
#include <ausa/util/FileUtil.h>
#include <TRandom.h>


RootWriter::RootWriter(std::string output) {
    f = AUSA::openFileEnsureDirectory(output, "RECREATE");
    t = new TTree("sim","Tree with simulated decays");

    pOut = new TClonesArray("TLorentzVector",3);

    t->Branch("wB",&wB,"wB/D");
    t->Branch("wU",&wU,"wU/D");
    t->Branch("w",&w,"w/D");
    t->Branch("x",&x,"x/D");
    t->Branch("y",&y,"y/D");
    t->Branch("ecm",ecm,"ecm[3]/D");
    t->Branch("pAlpha","TClonesArray",&pOut,32000,1);
}

RootWriter::~RootWriter() {
    f->WriteTObject(t);
    f->Flush();
    f->Close();
}

void RootWriter::fill(const TriAlphaDecay& decay) {
    wB = decay.wB;
    wU = decay.wU;
    w = decay.w;

    auto& p = *pOut;

    for (size_t i = 0; i < 3; ++i) {
        ecm[i] = decay.eCM[decay.order[i]];

        new(p[i]) TLorentzVector(*decay.p[decay.order[i]]);
    }

    auto i = gRandom->Integer(3);
    auto j = i;
    while (j == i) j = gRandom->Integer(3);

    const auto tot = ecm[0] + ecm[1] + ecm[2];
    double e1 = ecm[i] / tot;
    double e2 = ecm[j] / tot;
    x = (e1 + 2. * e2 - 1.) / sqrt(3.);
    y = e1 - 1./3.;

    t->Fill();
}