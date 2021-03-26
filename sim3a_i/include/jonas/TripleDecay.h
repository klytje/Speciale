#ifndef TRIPLE_DECAY_H
#define TRIPLE_DECAY_H
#include <vector>
#include <memory>
#include <TLorentzVector.h>
#include "Recoil.h"
#include "DecayWeight.h"

using namespace std;

class TripleDecay {

  private:
    double Q;
    //Decay products before beta-nu recoil
    //Products after boost
    TLorentzVector beta;
    TLorentzVector neutrino;
    vector<double> cmEnergies;
    std::unique_ptr<Recoil> recoilGenerator;

protected:
    vector<TLorentzVector> products;
    double wGenerator;
    bool doGroundState;

    double GenerateGroundState();

    double wCalculator;

    double GenerateUniform();

    unique_ptr<DecayWeight> weight;
    vector<TLorentzVector> rawProducts;

    void RecoilBoost();

    bool doRecoil;
public:
    /**
    * The TripleDecay will per default generate triple-alpha decays and return a weight
    * which describes a uniform phase-space decay.
    */
    TripleDecay();
    ~TripleDecay();

    /**
    * Set the available energy for the triple-alpha breakup. Units: keV.
    */
    void SetQ(double);

    /**
    * You can enable beta-neutrino recoil calculation here. The Q-value of the beta-
    * decay is calculated from the ground state of 12N. The TLorentzVectors of the
    * decay products are boosted with the recoil-momentum immediately after generation.
    */
    void DoRecoil(bool input = true);

    /**
    * If you need to generate decays through the 8Be ground state.
    */
    void DoGroundState(bool input = true);

    /**
    * Generates a decay and returns the weight of that decay. The result of the decay
    * can be retrieved through the GetProduct() and GetCmEnergy() methods.
    */
    double Generate();
  
    const TLorentzVector & GetProduct(int);

    const TLorentzVector & GetRawProduct(int);

    const TLorentzVector & GetBeta();

    const TLorentzVector & GetNeutrino();

    /**
    * Returns the center of mass energies of the three decay products. This
    * corrects for the beta-neutrino recoil effect.
    */
    double GetCmEnergy(int);

    double GetMaxWeight();

    double GetGeneratorWeight();
    double GetCalculatorWeight();

    /**
    * 
    */
    void SetDecayWeight(unique_ptr<DecayWeight>);
};
#endif
