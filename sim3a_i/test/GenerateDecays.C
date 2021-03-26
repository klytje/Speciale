#include <memory>
#include <iostream>
#include <string>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TH1I.h>
#include <TDatime.h>
#include "jonas/TripleDecay.h"
#include "jonas/Pair.h"
#include "jonas/Nucleus.h"
#include "jonas/Level.h"
#include "jonas/Channel.h"
#include "jonas/BalamuthWeight.h"

using namespace std;

void GenerateDecays(double q, int N = 1e6, string file = "test.root")
{
  TripleDecay decay;
  decay.DoRecoil();
  decay.SetQ(q);

  //XXX --------------------------------------- XXX  
  //We define the first channel
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  Pair pPair(Be8,He4);
  Channel pChannel(pPair,2,35.,6.7);  //reduced width not really important here, arbitrary.
  Level pLevel(12.71,1,1);            //energy not relevant.
  
  //The second channel
  Pair sPair(He4,He4);
  Channel sChannel(sPair,2,32.787,4.5);
  Level sLevel(3037.,2,1);            //Energy above ground state.

  //We create the weight calculator for the triple-alpha decay.
  unique_ptr<BalamuthWeight> bal(new BalamuthWeight());
  bal->SetPrimary(pLevel, pChannel);
  bal->SetSecondary(sLevel, sChannel);
  decay.SetDecayWeight(move(bal));
  //XXX --------------------------------------- XXX

  TFile *outFile = new TFile(file.c_str(),"recreate");
  TTree *outTree = new TTree("sim","Tree with simulated decays");

  double w, x, y;
  double ecm[3];
  TClonesArray *pAlpha = new TClonesArray("TLorentzVector",3);
  TClonesArray &p = *pAlpha;
  TLorentzVector pb;
  TLorentzVector *pBeta = &pb;
  TLorentzVector pn;
  TLorentzVector *pNeutrino = &pn;

  //Now we add the branches.
  outTree->Branch("w",&w,"w/D");
  outTree->Branch("x",&x,"x/D");
  outTree->Branch("y",&y,"y/D");
  outTree->Branch("ecm",ecm,"ecm[3]/D");
  outTree->Branch("pAlpha","TClonesArray",&pAlpha,32000,1); //Split level = 1 -> accessible through TTreeReaderArray (bug in ROOT).
  outTree->Branch("pBeta","TLorentzVector",&pBeta);
  outTree->Branch("pNeutrino","TLorentzVector",&pNeutrino);

  TDatime tBegin, tNow;
  tBegin.Set(); printf("*==* ---------- Begin of Job ---------- "); 
  tBegin.Print();

  int nDec = 0;
  while(nDec < N){
    p.Clear();
    w = decay.Generate();

    double alphaEnergies[3];
    alphaEnergies[0] = decay.GetCmEnergy(0);
    alphaEnergies[1] = decay.GetCmEnergy(1);
    alphaEnergies[2] = decay.GetCmEnergy(2);
    int alphaOrder[3];
    TMath::Sort(3,alphaEnergies,alphaOrder,true);

    ecm[0] = alphaEnergies[alphaOrder[0]];
    ecm[1] = alphaEnergies[alphaOrder[1]];
    ecm[2] = alphaEnergies[alphaOrder[2]];
    x = Sqrt(3) * (ecm[0] - ecm[2]) / (ecm[0] + ecm[1] + ecm[2]);
    y = (2*ecm[1] - ecm[0] - ecm[2]) / (ecm[0] + ecm[1] + ecm[2]);

    TLorentzVector p1 = decay.GetProduct(alphaOrder[0]);
    TLorentzVector p2 = decay.GetProduct(alphaOrder[1]);
    TLorentzVector p3 = decay.GetProduct(alphaOrder[2]);

    new(p[0]) TLorentzVector(p1.Px(),p1.Py(),p1.Pz(),p1.E());
    new(p[1]) TLorentzVector(p2.Px(),p2.Py(),p2.Pz(),p2.E());
    new(p[2]) TLorentzVector(p3.Px(),p3.Py(),p3.Pz(),p3.E());

    TLorentzVector pBetaTmp = decay.GetBeta();
    pBeta = &pBetaTmp;
    TLorentzVector pNeutrinoTmp = decay.GetNeutrino();
    pNeutrino = &pNeutrinoTmp;

    outTree->Fill();
    nDec++;
  }

  tNow.Set(); printf("*==* ---------- End of Job ------------ "); 
  tNow.Print();

  outTree->Write("",TObject::kOverwrite);
  outFile->Close();
}
