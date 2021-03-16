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
    tree = make_unique<TTree>("a", "a");

    // define the private variables as the leaf nodes
    tree->Branch("mul", &mul);
    tree->Branch("exC12", &exC12);
    tree->Branch("exBe8", &exBe8, "exBe8[3]/D");
    tree->Branch("pt", &ptot);
    tree->Branch("px", &px);
    tree->Branch("py", &py);
    tree->Branch("pz", &pz);
//    tree->Branch("sV", &sumAng);
//    tree->Branch("dV", &devAng);
    tree->Branch("deltaE", &dE);
    tree->Branch("z_angle", &zangle);
    tree->Branch("N", &N);
//    tree->Branch("d", d , "d[3]/I");
    tree->Branch("E_cm", E_cm, "ecm[3]/D");
    tree->Branch("E_dep", E_dep, "edep[3]/D");
    tree->Branch("E_lab", E_lab, "elab[3]/D");
    tree->Branch("theta_cm", theta_cm, "thetacm[3]/D");
    tree->Branch("theta_lab", theta_lab, "thetalab[3]/D");
    tree->Branch("phi_cm", phi_cm, "phicm[3]/D");
    tree->Branch("phi_lab", phi_lab, "philab[3]/D");
//    tree->Branch("vZ", vZ, "vZ[3]/D");
//    tree->Branch("vZL", vZL, "vZL[3]/D");
    tree->Branch("prob", prob, "prob/D");
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
        int indices[3]; // indices of the alphas
        TLorentzVector p[3]; // momenta of the alphas
        TLorentzVector p_tot; // total momentum of the alphas
        for (int i = 0; i < m; i++) {
            if (mul == 3) break; // stop when we've found all three
            
            auto ion = event.getIon(i); 
            if (*ion != ALPHA) continue; // we are looking for alphas
            
            auto& pi = event.getLorentzVector(i);
            indices[mul] = i; // store the index of the particle
            p[mul] = pi;
            p_tot += pi; // update the total momentum
            mul++;
        }

        // if we only found two alpha particles, we can reconstruct the third
        if (m == 2) {
            m++; // increment m for the next if check
            p[2] = p_beam - p_tot; // pbeam - (p0 + p1 + p2) = 0
            p_tot += p[2]; // ptot = pbeam
            indices[2] = -1; // indicate that this particle is reconstructed
        }
        
        // at this stage we work only with the three alphas, so any other event is discarded
        if (m != 3) continue; 

        // calculate the momentum in the cm frame
        TLorentzVector p_tot_cm = p_tot;
        p_tot_cm.Boost(-bboost);

        zangle = p[0].Vect().Cross(p[1].Vect()).Angle(Z_AXIS)*TMath::RadToDeg();
        dE = p_beam.E() - p_tot.E(); // energy difference between beam and alphas
        exC12 = p_tot_cm.E() - C12_MASS; // excitation energy of C12
        ptot = p_tot_cm.P(); 
        px = p_tot.X(); // note: not cm frame
        py = p_tot.Y(); // note: not cm frame
        pz = p_tot.Z(); // note: not cm frame

        // loop over the three alpha particles
        // sumAng = 0; devAng = 0;
        for (int i = 0; i < 3; i++) {
            TLorentzVector pi = p[i]; // copy the momentum so we can boost it
            vz[i] = pi.Angle(Z_AXIS)*TMath::RadToDeg(); // idk why we have two
            vzl[i] = pi.Angle(Z_AXIS)*TMath::RadToDeg();
            E_lab[i] = pi.E() - Constants::ALPHA_MASS;
            theta_lab[i] = pi.Theta()*TMath::RadToDeg();
            phi_lab[i] = pi.Phi()*TMath::RadToDeg();
            pi.Boost(-bboost);

        //    if (i > 0) 
        //        sumAng += abs(pi.Angle(p[0].Vect()));
        //    if (i == 2)
        //        devAng = abs(pi.Angle(p[0].Vect().Cross(p[1].Vect())));

            E_cm[i] = pi.E() - Constants::ALPHA_MASS;
            theta_cm[i] = pi.Theta()*TMath::RadToDeg();
            phi_cm[i] = pi.Phi()*TMath::RadToDeg();
            TLorentzVector pBe8 = p_tot - pi;
            pBe8.Boost(-pBe8.BoostVector());
            exBe8[i] = pBe8.E() - BE8_MASS;
            prob[i] = ROOT::Math::breitwigner_cdf_c(abs(this->exBe8[i]-3030),1513,0)*2;
        }
        // sumAng += abs(p[1].Angle(p[2].Vect()));
        tree->Fill();
    }
}

void StandardAnalyzer::terminate() {
    gDirectory->WriteTObject(tree.release());
}

void StandardAnalyzer::saveToRootFile(TFileWrapper &file) {}
void StandardAnalyzer::reset() {}
