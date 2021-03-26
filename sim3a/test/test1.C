#include <memory>
#include <iostream>
#include <TH1I.h>
#include <TH2I.h>
#include <TRandom3.h>
#include <TMath.h>
#include "jonas/TripleDecay.h"
#include "jonas/Pair.h"
#include "jonas/Nucleus.h"
#include "jonas/Level.h"
#include "jonas/Channel.h"
#include "jonas/BalamuthWeight.h"

using namespace std;

TH2I * test(double q)
{
  TripleDecay decay;
  //decay.DoRecoil();
  //double q = 1000.; //5434.2;
  decay.SetQ(q);
  
  //We define the first channel
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  Pair pPair(Be8,He4);
  Channel pChannel(pPair,2,35.,6.7);  //reduced width not really important here, arbitrary.
  Level pLevel(12.71,0,1);            //energy not relevant.
  
  //The second channel
  Pair sPair(He4,He4);
  Channel sChannel(sPair,2,32.787,4.5);
  Level sLevel(3037.,2,1);            //Energy above ground state.

  //We create the weight calculator for the triple-alpha decay.
  unique_ptr<BalamuthWeight> bal(new BalamuthWeight());
  bal->SetPrimary(pLevel, pChannel);
  bal->SetSecondary(sLevel, sChannel);
  //bal->SetScale(1.8217e4);     //Appropriate for 1+ through 2+ at Q = 5435.2keV.
  //bal->SetScale(9.4e5);          //Appropriate for 1+ through 2+ at Q = 3000keV.
  //bal->SetScale(3.09e5);     //Appropriate for 0+ through L=2 at Q = 3000keV.
  //bal->SetScale(6.57e4);     //Appropriate for 2+ through L=0 at Q = 3000keV.
  //bal->SetScale(5.39e4);       //Appropriate for 2+ through L=2 at Q = 3000keV.
  //bal->SetScale(3.55e6);        //Appropriate for 2+ through L=4 at Q = 3000keV.

  decay.SetDecayWeight(move(bal));

  TH1I *h = new TH1I("h","",1000,0,1e-3);
  TH2I *dalitz = new TH2I("dalitz","",140,0,1.4,140,-0.7,0.7);

  TRandom3 rGen(0);

  int nDec = 0;
  double maxWeight = 0;
  while(nDec < 2e5){
    double weight = decay.Generate();
    h->Fill(weight);
    if(weight > maxWeight) maxWeight = weight;
    //if(rGen.Uniform() > weight) continue;
    nDec++;

    /*
    double alphaEnergies[3];
    alphaEnergies[0] = decay.GetProduct(0).Energy() * 1e6 - He4.M();
    alphaEnergies[1] = decay.GetProduct(1).Energy() * 1e6 - He4.M();
    alphaEnergies[2] = decay.GetProduct(2).Energy() * 1e6 - He4.M();
    int alphaOrder[3];
    TMath::Sort(3,alphaEnergies,alphaOrder,true);
    double  e1 = alphaEnergies[alphaOrder[0]];
    double  e2 = alphaEnergies[alphaOrder[1]];
    double  e3 = alphaEnergies[alphaOrder[2]];
    double x = Sqrt(3) * (e1 - e3) / (e1 + e2 + e3);
    double y = (2*e2 - e1 - e3) / (e1 + e2 + e3);
    dalitz->Fill(x,y);
    */
    /*
    double EA1 = decay.GetCmEnergy(0);
    double EA2 = decay.GetCmEnergy(1);
    double EA3 = decay.GetCmEnergy(2);

    double ea1 = decay.GetProduct(0).Energy()*1e6 - He4.M();
    double ea2 = decay.GetProduct(1).Energy()*1e6 - He4.M();
    double ea3 = decay.GetProduct(2).Energy()*1e6 - He4.M();

    TLorentzVector total = decay.GetProduct(0) + decay.GetProduct(1) + decay.GetProduct(2);
    double Etot = total.Energy()*1e6 - 3*He4.M();

    cout << "weight = " << weight << endl;
    cout << "E1 = " << EA1 << endl;
    cout << "E2 = " << EA2 << endl;
    cout << "E3 = " << EA3 << endl;
    cout << "Etot = " << EA1 + EA2 + EA3 << endl;
    cout << "In lab system: " << endl;
    cout << "E1 = " << ea1 << endl;
    cout << "E2 = " << ea2 << endl;
    cout << "E3 = " << ea3 << endl;
    cout << "Etot = " << Etot << endl;
    cout << endl;
    */
  }

  cout << maxWeight << endl;
  h->Draw();
  return dalitz;
}
