//
// Created by kuhlwein on 1/13/20.
//

#include "RootInterferenceWriter.h"
#include <ausa/util/FileUtil.h>
#include <TRandom.h>
#include <TParameter.h>
#include <TF1.h>


RootInterferenceWriter::RootInterferenceWriter(std::string output) {
    f = new TFile(output.c_str(), "RECREATE");
    t = new TTree("sample","sample");
    auto para = new TParameter<int>("N_ELM",9);
    t->GetUserInfo()->Add(para);

    pOut = new TClonesArray("TLorentzVector",3);

    //t->Branch("wB",&wB,"wB/D");
    t->Branch("wU",&wU,"wU/D");
    //t->Branch("w",&w,"w/D");
    t->Branch("x",&x,"x/D");
    t->Branch("y",&y,"y/D");
    //t->Branch("ecm",ecm,"ecm[3]/D");
    //t->Branch("pAlpha","TClonesArray",&pOut,32000,1);
    t->Branch("p", p,std::string("p[9]/D").c_str());
    t->Branch("f",&b);
}

RootInterferenceWriter::~RootInterferenceWriter() {
    f->WriteTObject(t);

    std::vector<int> J = {1,2,3};

    f->WriteObject(&J,"J");
    f->WriteObject(&J,"P");
    f->WriteObject(&J,"L1");
    f->WriteObject(&J,"L2");
    f->Flush();
    f->Close();
}

void RootInterferenceWriter::fill(const TriAlphaDecay& decay) {
    wU = decay.wU;

    for (size_t i = 0; i < 3; ++i) {
        ecm[i] = decay.eCM[decay.order[i]];
    }

    int counter = 0;
    for (auto &v : decay.p) {
        auto vect = v->Vect();
        for (int i=0; i<3; i++) {
            p[counter] = vect[i];
            counter++;
        }
    }

    auto i = gRandom->Integer(3);
    auto j = i;
    while (j == i) j = gRandom->Integer(3);

    const auto tot = ecm[0] + ecm[1] + ecm[2];
    double e1 = ecm[i] / tot;
    double e2 = ecm[j] / tot;
    x = (e1 + 2. * e2 - 1.) / sqrt(3.);
    y = e1 - 1./3.;

    b.clear();
    b = decay.factors;

    t->Fill();
}