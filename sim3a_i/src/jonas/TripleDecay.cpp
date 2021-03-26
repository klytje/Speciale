#include <TGenPhaseSpace.h>
#include <ausa/constants/Mass.h>
#include <assert.h>
#include <ausa/util/memory>
#include <ausa/constants/Mass.h>
#include <jonas/Constants.h>
#include "jonas/Nucleus.h"
#include "jonas/DecayWeight.h"
#include "jonas/TripleDecay.h"

using namespace std;

TripleDecay::TripleDecay() : weight(new DecayWeight())
{
  Q = 1000.;
  doRecoil = true;
  doGroundState = false;
}

TripleDecay::~TripleDecay() {}

void TripleDecay::SetQ(double q)
{
  Q = q;
}

void TripleDecay::DoRecoil(bool input)
{
  doRecoil = input;
  if(doRecoil){
    //We initialise the recoil generator for th 12N decay.
    double Ei = AUSA::Constants::isotopeMass(7,12) - 510.9989461;
    double Ef = 3 * AUSA::Constants::isotopeMass(2,4);
    double QBetaMax = Ei - Ef;
    recoilGenerator = std::make_unique<Recoil>(QBetaMax,1000);
  }
}

void TripleDecay::DoGroundState(bool input)
{
  doGroundState = input;
}
 
double TripleDecay::Generate()
{
  double w = 1.;
  if(doGroundState){
    wGenerator = GenerateGroundState();
    wCalculator = 1.;
  }
  else{
    wGenerator = GenerateUniform();
    wCalculator = weight->Calculate(products);
    //factors = {weight->CalculateInterference(products)};
  }

  w *= wGenerator;
  w *= wCalculator;

  rawProducts = products; //Store raw decay.
  if(doRecoil) RecoilBoost();

  return w;
}

double TripleDecay::GenerateUniform()
{
  //We'll need the mass of the alpha-particle.
  double initialMass = (3 * ALPHA_MASS + Q)/1.e6; //Remember, GeV.
  TLorentzVector mother(0.0,0.0,0.0,initialMass);
  double masses[3] = {ALPHA_MASS/1e6, ALPHA_MASS/1e6, ALPHA_MASS/1e6}; //Remember, GeV.

  //Generate the decay.
  TGenPhaseSpace event;
  event.SetDecay(mother, 3, masses);
  double w = event.Generate();

  //Then we store the results in a safe place.
  products.clear();
  products.push_back((*(event.GetDecay(0)))*1e6); //We like the results in keV.
  products.push_back((*(event.GetDecay(1)))*1e6);
  products.push_back((*(event.GetDecay(2)))*1e6);

  cmEnergies.clear();
  cmEnergies.push_back(products.at(0).Energy() - ALPHA_MASS);
  cmEnergies.push_back(products.at(1).Energy() - ALPHA_MASS);
  cmEnergies.push_back(products.at(2).Energy() - ALPHA_MASS);

  return w;
}

double TripleDecay::GenerateGroundState()
{
  double initialMass = (3 * ALPHA_MASS + Q)/1.e6; //Remember, GeV.
  TLorentzVector mother(0.0,0.0,0.0,initialMass);
  double masses[2] = {ALPHA_MASS/1e6, BE8_MASS/1e6}; //Remember, GeV.
  //Generate the decay.
  TGenPhaseSpace event;
  event.SetDecay(mother, 2, masses);
  double w = event.Generate();

   //Then we store the results in a safe place.
  products.clear();
  products.push_back((*(event.GetDecay(0)))*1e6); //We like the results in keV.

  //Next, the 8Be breaks up: New mother and daughters.
  mother = *(event.GetDecay(1));
  masses[0] = ALPHA_MASS/1e6;
  masses[1] = ALPHA_MASS/1e6;
  event.SetDecay(mother, 2, masses); 
  w *= event.Generate();
  products.push_back((*(event.GetDecay(0)))*1e6);
  products.push_back((*(event.GetDecay(1)))*1e6);

  cmEnergies.clear();
  cmEnergies.push_back(products.at(0).Energy() - ALPHA_MASS);
  cmEnergies.push_back(products.at(1).Energy() - ALPHA_MASS);
  cmEnergies.push_back(products.at(2).Energy() - ALPHA_MASS);

  return w;
}

void TripleDecay::RecoilBoost()
{
  Nucleus He4(4,2,0,1);
  double Ei = AUSA::Constants::isotopeMass(7,12) - 510.9989461;
  double Ef = 3 * ALPHA_MASS + Q;
  double Qbeta = Ei - Ef;
  TVector3 p = recoilGenerator->Generate(Qbeta,&beta,&neutrino);
  double Erecoil = TMath::Sqrt(TMath::Power(Ef,2) + p.Mag2());
  TLorentzVector recoil(p,Erecoil);
  for(auto & element : products) {
    element.Boost(-recoil.BoostVector());
  }
}

const TLorentzVector & TripleDecay::GetProduct(int i)
{
  return products.at(i);
}

const TLorentzVector & TripleDecay::GetRawProduct(int i)
{
  return rawProducts.at(i);
}

const TLorentzVector & TripleDecay::GetBeta()
{
  return beta;
}

const TLorentzVector & TripleDecay::GetNeutrino()
{
  return neutrino;
}
 
double TripleDecay::GetCmEnergy(int i)
{
  return cmEnergies.at(i);
}

double TripleDecay::GetMaxWeight()
{
  return 0.5;
}

void TripleDecay::SetDecayWeight(unique_ptr<DecayWeight> w)
{
  weight = move(w);
}

double TripleDecay::GetGeneratorWeight()
{
  return wGenerator;
}

double TripleDecay::GetCalculatorWeight()
{
  return wCalculator;
}
