//
// Created by munk on 6/19/16.
//

#include "simX/samplers/RootSampleSet.h"
#include <TTree.h>
#include <TFile.h>
#include <TParameter.h>
#include "simX/Context.h"

simX::samplers::RootSampleSet::RootSampleSet(std::string path, size_t nElements)
    :  result(nElements), sampled(true)
{
    SIMX_TCONTEXT;
    file = new TFile(path.c_str(), "READ");

    if (!file->IsOpen()) throw std::runtime_error("Failed to open: " + path);

    tree = dynamic_cast<TTree*>(file->Get("sample"));
    if (!tree) throw std::invalid_argument("Tree named sample not found in ROOT file " + path);

    auto param = dynamic_cast<TParameter<int>*>(tree->GetUserInfo()->FindObject("N_ELM"));
    if (!param) throw std::invalid_argument("Tree named sample does not contain user parameter named N_ELM!");

    if (param->GetVal() != nElements)
        throw std::runtime_error("Supplied ROOT file is for " + std::to_string(param->GetVal())
                                 + "elements. Not " + std::to_string(nElements) + " as expected");

    tree->SetBranchAddress("p", result.data());

    if (tree->GetBranch("w")) {
        sampled = false;
        tree->SetBranchAddress("w", &w);
    }

    reset();
}

simX::samplers::RootSampleSet::~RootSampleSet() {
    if (file) {
        file->Close();
        delete file;
    }
}

const std::vector<double> &simX::samplers::RootSampleSet::getSample() {
    return result;
}

const double simX::samplers::RootSampleSet::getWeight() {
    return w;
}

bool simX::samplers::RootSampleSet::isSampled() const {
    return sampled;
}

void simX::samplers::RootSampleSet::next() {
    tree->GetEntry(++current);
}

bool simX::samplers::RootSampleSet::hasNext() {
    return current+1 < tree->GetEntries();
}

void simX::samplers::RootSampleSet::reset() {
    current = -1;
}




