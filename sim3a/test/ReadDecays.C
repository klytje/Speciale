#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <string>
#include <iostream>

using namespace std;

void ReadDecays()
{
  string simName = "/home/jonas/data/i161/simulated/decays/3000keV/J1L2.root";
  TFile *simFile = new TFile(simName.c_str(),"READ");
  TTree *simTree = (TTree*)simFile->Get("sim");
  double maxWeight = simTree->GetMaximum("w");
  TTreeReader reader("sim",simFile);
  TTreeReaderValue<double> w(reader,"w");
  TTreeReaderValue<double> x(reader,"x");
  TTreeReaderValue<double> y(reader,"y");
  TTreeReaderArray<double> ecm(reader,"ecm");
  TTreeReaderArray<TLorentzVector> pAlpha(reader,"pAlpha");

  TH2I *h = new TH2I("h","",140,0,1.4,140,-0.7,0.7);
  TH1I *h2 = new TH1I("h2","",400,0,2000);
  TH1I *h3 = new TH1I("h3","",400,0,2000);
  TRandom3 rGen;
  rGen.SetSeed(0);

  while(reader.Next()) {
    //Just access the data as if myPx and myPy were iterators (note the '*'
    //in front of them):
    if(rGen.Uniform(maxWeight) < *w){
      h->Fill(*x,*y);
      for(double E : ecm){
        h2->Fill(E);
      }
      for(TLorentzVector alpha : pAlpha){
        h3->Fill(alpha.E() - alpha.M());
      }
    }
  }

  //h->Draw("colz");
  //h2->Draw();
  h3->Draw();
}
