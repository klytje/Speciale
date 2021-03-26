#include <TMath.h>
#include <iostream>
#include <assert.h>
#include "jonas/Recoil.h"

using namespace std;

Recoil::Recoil(double q, int n)
{
  Initialise(q,n);
}

Recoil::~Recoil() {}

void Recoil::Initialise(double q, int n)
{
  fTbeta.clear();
  fCos.clear();
  rGen.SetSeed(0);

  Qmax = q;
  Wmax = (Qmax + mBeta) / mBeta;
  int nBins = n;

  for(int i=0; i<nBins; i++){
    double Q = Qmax/nBins * (i + 0.5); //'Bin' center.
    char fName[16];
    sprintf(fName,"fTbeta_%.1f",Q);
    TF1 fT(fName,"TMath::Sqrt(x*x+2*x*[1])*TMath::Power([0]-x,2)*(x+[1])",0,Q);
    fT.SetParameter(0,Q);
    fT.SetParameter(1,mBeta);
    fTbeta.push_back(fT);

    double w = (Wmax - 1)/nBins * (i + 0.5) + 1.0;
    sprintf(fName,"fW_%.3f",w);   
    TF1 fC(fName,"0.5*(1.-TMath::Sqrt([0]*[0]-1.)/(3.*[0])*x)",-1,1);
    fC.SetParameter(0,w);
    fCos.push_back(fC);
  }
}



int Recoil::FindQbin(double Q)
{
  int bin = Floor(Q/(Qmax/NBins()));
  return bin;
}

int Recoil::FindWbin(double W)
{
  int bin = Floor((W - 1) / (Wmax - 1) * NBins());
  return bin;
}

TVector3 Recoil::Generate(double Q, TLorentzVector *pb, TLorentzVector *pn)
{
  assert(Q>0 && Q<Qmax);
  int binQ = FindQbin(Q);
  double Tbeta = fTbeta[binQ].GetRandom();    //Find the kinetic energy of the electron.
  double Ebeta = Tbeta + mBeta;
  double w = Ebeta / mBeta;
  int binW = FindWbin(w);
  double x = fCos[binW].GetRandom();          //Find cos(angle between electron and neutrino momenta).
  double theta = ACos(x);
  double phi = rGen.Uniform(0,2*Pi());        //Uniform distribution of phi.
  double Eneutrino = Q - Tbeta;
  double pBetaMag = Sqrt(Power(Ebeta,2)-Power(mBeta,2));
  double pNeutrinoMag = Eneutrino;
 
  //We choose the direction of the beta particle randomly.   
  double px, py, pz;
  rGen.Sphere(px,py,pz,pBetaMag);
  TVector3 pBeta(px,py,pz);                      //Units: keV/c,keV.
  TVector3 pNeutrino(pBeta.Unit()*pNeutrinoMag); //Initially in same direction as pBeta.

  //The direction of the neutrino must follow the correct angular distribution, so we rotate it.
  pNeutrino.Rotate(theta,pBeta.Cross(TVector3(0,0,1)));
  pNeutrino.Rotate(phi,pBeta);

  TVector3 p = pBeta + pNeutrino;  //Total recoil momentum.

  if(pb){
    pb->SetVect(pBeta);
    pb->SetE(Ebeta);
  }
  if(pn){
    pn->SetVect(pNeutrino);
    pn->SetE(Eneutrino);
  }
  return p;
}

int Recoil::NBins()
{
  return fCos.size();
}
