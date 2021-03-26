//
// Created by munk on 6/19/16.
//

#include "simX/samplers/JonasLibSampleSet.h"
#include "simX/Logger.h"
#include "simX/Context.h"
#include <TChain.h>
#include <TParameter.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <ausa/util/memory>

simX::samplers::JonasLibSampleSet::JonasLibSampleSet(std::string path, size_t nParticles)
    :  result(300*nParticles), arr(new TClonesArray("TLorentzVector"))
{
    static auto LOGGER = simX::log::getLogger("JonasLibSampleSet");
    SIMX_TCONTEXT;
    tree = std::make_unique<TChain>("sim");

    LOGGER->debug("path: {}, nParticles {}", path, nParticles);
    auto n = tree->Add(path.c_str());
    if (n == 0)
        throw std::invalid_argument("Jonas lib files not found");

    if (!tree -> GetBranch("pAlpha"))
        throw std::invalid_argument("Missing branch pAlpha!");

    if (!tree -> GetBranch("pAlpha.fP"))
        throw std::invalid_argument("Missing branch pAlpha.fP!");

    if (!tree -> GetBranch("w"))
        throw std::invalid_argument("Missing branch w!");

    tree->SetBranchAddress("w", &w);
    tree->SetBranchAddress("pAlpha", &arr);

    tree->GetEntry(0);
    if (arr->GetEntries() != nParticles)
        throw std::runtime_error("Supplied ROOT file is for " + std::to_string(n)
                                 + "particles. Not " + std::to_string(nParticles) + " as expected");


    tree->SetBranchStatus("*", false);
    tree->SetBranchStatus("pAlpha*", true);
    tree->SetBranchStatus("w", true);

    reset();
}

simX::samplers::JonasLibSampleSet::~JonasLibSampleSet() {

}

const std::vector<double> &simX::samplers::JonasLibSampleSet::getSample() {
    return result;
}

const double simX::samplers::JonasLibSampleSet::getWeight() {
    return w;
}

bool simX::samplers::JonasLibSampleSet::isSampled() const {
    return false;
}

void simX::samplers::JonasLibSampleSet::next() {
    tree->GetEntry(++current);

    for (int i = 0; i < arr->GetEntries(); ++i) {
        for (int j = 0; j < 3; ++j) {
            result[3*i+j] = ((TLorentzVector&)*arr->At(i))[j];
        }
    }
}

bool simX::samplers::JonasLibSampleSet::hasNext() {
    return current+1 < tree->GetEntries();
}

void simX::samplers::JonasLibSampleSet::reset() {
    current = -1;
}




