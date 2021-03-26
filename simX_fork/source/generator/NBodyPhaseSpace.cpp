
#include "simX/generator/NBodyPhaseSpace.h"
#include "simX/Random.h"
#include "simX/NuclearProcess.h"

#include "ausa/util/memory"

using namespace std;
using namespace simX;

namespace {
	double calcTotalEnergyGeV(const Particle& parent) {
		return (parent.getMass() + parent.getNominalExcitationEnergy())*1E-6;
	}

	vector<double> extractMassesGeV(const vector<Particle*>& d) {
		vector<double> masses;
		masses.reserve(d.size());

		for (auto& p : d) masses.push_back((p->getMass() + p->getNominalExcitationEnergy())*1E-6);

		return masses;
	}
}

NBodyPhaseSpace::NBodyPhaseSpace(const Particle& p, const vector<Particle*>& d)
		: parent(p), 
		  daughtersPtr(d),
		  isInitialized(false),
		  mult(d.size()),
		  pDaughters(d.size())
{}


std::vector<TLorentzVector>& NBodyPhaseSpace::getFourMomenta() {
    // Initialize sampler if needed
    // (nominal excitation energies are being used)
    if (!isInitialized) init();

    // Generate final state from phase-space distribution
	double w=0, r=1;
	while (r>w) { // discard events if r>w (Von Neumann sampling)
		w = phaseSpaceGenerator.Generate() / Wmax; // normalized weight
		r = rnd(); // random number between 0 and 1
	}

    // Access generated four-momenta
	for ( int i=0; i<mult; i++ ) {
	    pDaughters[i] = *(phaseSpaceGenerator.GetDecay(i));
	    pDaughters[i] *= 1e6;    // GeV->keV
    }

    return pDaughters;
}


void NBodyPhaseSpace::init() { 
    // Total energy
    double etotGeV = calcTotalEnergyGeV(parent);
    // Total mass of daughters
    auto massesGeV = extractMassesGeV(daughtersPtr);
    // Check if decay is energetically allowed
    double msumGeV = 0;
   	for (auto& m : massesGeV) msumGeV += m;
    if (msumGeV>etotGeV)
        throw invalid_argument("Phase-space generator cannot be initialized because decay is not energetically allowed - exit.");
	TLorentzVector P(0.,0.,0.,etotGeV);
	phaseSpaceGenerator.SetDecay(P, static_cast<int>(massesGeV.size()), massesGeV.data() );

	Wmax = 0;
	// thePhaseSpaceGenerator->GetWtMax() overestimates the max weight (ROOT BUG)
	// So we must determine maximum weight ourselves:

	for (Int_t i=0; i< N_SAMPLES; i++) {
		double w = phaseSpaceGenerator.Generate();
		Wmax = std::max(w, Wmax);
	}

	isInitialized = true;
}



