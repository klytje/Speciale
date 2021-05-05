#include "StandardAnalyzer.h"
#include "TFile.h"
#include "Math/ProbFunc.h"

#include <ausa/constants/Mass.h> 
#include <iostream>
#include <TTree.h>

namespace {
    const auto C12_MASS = AUSA::Constants::isotopeMass("C12");
    const auto B11_MASS = AUSA::Constants::isotopeMass("B11");
    const auto BE8_MASS = AUSA::Constants::isotopeMass("Be8");
    const auto PROTON_MASS = AUSA::Constants::isotopeMass("H1");
    const auto ALPHA = AUSA::EnergyLoss::Ion::predefined("He4");

    auto Q = PROTON_MASS + B11_MASS - 3*ALPHA.getMass();
    auto E_SUM = B11_MASS/C12_MASS * 161 + Q;

    auto Z_AXIS = TVector3{0, 0, 1};
}

using namespace std;
using namespace AUSA;
using namespace AUSA::Event;

StandardAnalyzer::StandardAnalyzer() {}
StandardAnalyzer::~StandardAnalyzer() {}

void StandardAnalyzer::setup(std::shared_ptr<Setup> setup) {
    tree = make_unique<TTree>("tree", "analyzed tree");

    // define the private variables as the leaf nodes
    tree->Branch("mul", &mul);
    tree->Branch("exC12", &exC12);
    tree->Branch("p_tot", &ptot);
    tree->Branch("deltaE", &dE);
    tree->Branch("N", &N);
    tree->Branch("mi", &mi, "mi[3]/i"); // ids of the alphas in mul
    tree->Branch("exBe8", &exBe8, "exBe8[3]/D");
    tree->Branch("E_cm", &E_cm, "E_cm[3]/D");
    tree->Branch("E_dep", &E_dep, "E_dep[3]/D");
    tree->Branch("E_lab", &E_lab, "E_lab[3]/D");
    tree->Branch("theta_lab", &theta_lab, "theta_lab[3]/D");
    tree->Branch("phi_lab", &phi_lab, "phi_lab[3]/D");
    tree->Branch("px", &px, "px[3]/D"); // momenta of the alphas in the cm frame
    tree->Branch("py", &py, "py[3]/D");
    tree->Branch("pz", &pz, "pz[3]/D");

    // I've disabled a few since I don't use them, and so they only waste space. They are fully functional however, and should work if you uncomment them.
    // tree->Branch("theta_cm", theta_cm, "theta_cm[3]/D");
    // tree->Branch("phi_cm", phi_cm, "phi_cm[3]/D");
//    tree->Branch("z_angle", &zangle);
//    tree->Branch("vZ", vZ, "vZ[3]/F"); 
//    tree->Branch("vZL", vZL, "vZL[3]/F");
//    tree->Branch("prob", prob, "prob[3]/F");
}

void StandardAnalyzer::setScalers(const ScalerOutput &s) {
    scaler = s;
}

void StandardAnalyzer::analyze(const std::vector<PhysicsEvent> &events) {
    // iterate through all events
    for (auto& event : events) {
        N = scaler.getOutput("___N___").getValue();
        size_t m = event.getDetectionMultiplicity(); // event multiplicity
        if (m < 2) continue; // we need at least two particle interpretations to continue

        TLorentzVector p_beam = event.getInitialLorentzVector();
        double be = p_beam.E(); // beam energy
        TVector3 bboost = p_beam.BoostVector(); // beam boost vector

        mul = 0; // number of alpha particles identified
        TLorentzVector p[3]; // momenta of the alphas
        TLorentzVector p_tot; // total momentum of the alphas
        for (int i = 0; i < m; i++) {
            if (mul == 3) break; // stop when we've found all three
            
            auto ion = event.getIon(i); 
            if (*ion != ALPHA) continue; // we are looking for alphas
            
            TLorentzVector pi = event.getLorentzVector(i);
            mi[mul] = i; // store the index of the particle
            p[mul] = pi;
            p_tot += pi; // update the total momentum
            mul++;
        }

        // if we only found two alpha particles, we can reconstruct the third
        if (mul == 2) {
            mul = -1; // mul is set to something distinctive, so we can easily filter all reconstructed events later
            p[2] = p_beam - p_tot; // pbeam - (p0 + p1 + p2) = 0
            p_tot += p[2]; // ptot = pbeam
            mi[2] = -1; // indicate that this particle is reconstructed
        }
        
        // at this stage we work only with the three alphas, possibly reconstructed, and any other event is discarded
        if (!(mul == 3 || mul == -1)) continue; 

        // calculate the momentum in the cm frame
        TLorentzVector p_tot_cm = p_tot;
        p_tot_cm.Boost(-bboost);

        zangle = p[0].Vect().Cross(p[1].Vect()).Angle(Z_AXIS)*TMath::RadToDeg();
        dE = p_beam.E() - p_tot.E(); // energy difference between beam and alphas
        exC12 = p_tot_cm.E() - C12_MASS; // excitation energy of C12
        ptot = p_tot_cm.P(); 

        // loop over the three alpha particles
        for (int i = 0; i < 3; i++) {
            TLorentzVector pi = p[i]; // copy the momentum so we can boost it
            vz[i] = pi.Angle(Z_AXIS)*TMath::RadToDeg(); 
            E_lab[i] = pi.E() - Constants::ALPHA_MASS;
            theta_lab[i] = pi.Theta()*TMath::RadToDeg();
            phi_lab[i] = pi.Phi()*TMath::RadToDeg();
            pi.Boost(-bboost);

            // from this point on, pi is boosted to the cm frame
            px[i] = pi.X();
            py[i] = pi.Y();
            pz[i] = pi.Z();

            E_cm[i] = pi.E() - Constants::ALPHA_MASS;
            theta_cm[i] = pi.Theta()*TMath::RadToDeg();
            phi_cm[i] = pi.Phi()*TMath::RadToDeg();
            TLorentzVector pBe8 = p_tot - pi;
            pBe8.Boost(-pBe8.BoostVector());
            exBe8[i] = pBe8.E() - BE8_MASS;
            prob[i] = ROOT::Math::breitwigner_cdf_c(abs(this->exBe8[i]-3030),1513,0)*2;
        }
        tree->Fill();
    }
}

void StandardAnalyzer::terminate() {
    gDirectory->WriteTObject(tree.release());
}

void StandardAnalyzer::saveToRootFile(TFileWrapper &file) {}
void StandardAnalyzer::reset() {}
