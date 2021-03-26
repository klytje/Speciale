//
// Created by munk on 6/22/16.
//

#include "simX/samplers/writer/RootMomentumWriter.h"

#include <TFile.h>
#include <TTree.h>

using namespace simX::samplers::writer;

void RootMomentumWriter::addSample(const std::vector<TVector3> &s, double w) {
    size_t j = 0;
    for (auto& v : s) {
        for (int i = 0; i < 3; ++i) {
            sample[j++] = v[i];
        }
    }
    weight = w;

    tree->Fill();
}

RootMomentumWriter::RootMomentumWriter(std::string path, bool hasWeight, size_t nElements) {
    file = new TFile(path.c_str(), "RECREATE");
    tree = new TTree("sample", "sample");

    sample.resize(nElements);

    tree->Branch("p", sample.data(), std::string("p[" + std::to_string(nElements) + "]/D").c_str());
    if (hasWeight) tree->Branch("w", &weight);
}

RootMomentumWriter::~RootMomentumWriter() {
    if (tree && file) {
        file->WriteTObject(tree);
    }
    if (file) {
        file->Flush();
        file->Close();

        delete file;
    }
}
