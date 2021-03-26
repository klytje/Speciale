#include "jonas/DetectorResponse.h"
#include "jonas/LineShapes.h"
#include <iostream>
#include <TMath.h>

using namespace std;
using namespace TMath;

DetectorResponse::DetectorResponse(double sigPhys, double tau1, double tau2, double eta, double sigF, double sigB, double Eref)
{
  //We want a physical response function for each 10keV bin in the pulse height defect.
  nBins = 1000;
  Emax = 10000.;
  stripSort = true;

  TF1 fPHD("fPHD","x/([0]+x)*([1]+[2]*x)",0,10000);   //Approximation to pulse height defect for alphas in silicon.
  fPHD.SetParameters(68,10.56,6.2e-4);

  for(int i=0; i<nBins; i++){
    double E = Emax/nBins * (i + 0.5); //'Bin' center.
    char fName[32];
    sprintf(fName,"fPhys_%.3f",E);
    double scale = 1.0;
    if(Eref > 0) scale = fPHD.Eval(E) / fPHD.Eval(Eref);
    double tau1Scaled = tau1 * scale;
    double tau2Scaled = tau2 * scale; 
    double lowerLim = -20 * tau2Scaled;        //This is when the distribution has dropped to ~1e-9 of its maximum value.
    double upperLim = 6.5 * sigPhys;           //Same.
    int nPoints = Floor(upperLim - lowerLim);  //We want one point pr. keV.
    TF1 temp(fName,MonoAlpha,lowerLim,upperLim,6);
    temp.SetParameters(1,0,sigPhys,tau1Scaled,tau2Scaled,eta);
    temp.SetNpx(nPoints);
    fPhys.push_back(temp);
  }

  fFront = new TF1("fFront","TMath::Gaus(x,0,[0],1)",-6.5*sigF,6.5*sigF);
  fFront->SetParameter(0,sigF);

  fBack = new TF1("fBack","TMath::Gaus(x,0,[0],1)",-6.5*sigB,6.5*sigB);
  fBack->SetParameter(0,sigB);

  fThres = new TF1("fThres","1",0,10000);

  rGen = new TRandom3(0);
}


DetectorResponse::~DetectorResponse()
{
  delete fFront;
  delete fBack;
  delete fThres;
  delete rGen;
}

void DetectorResponse::Generate(int mul, int *fi, int *bi, double *energy)
{
  frontMul = 0;
  backMul = 0;

  if(mul == 0) return;

  double fe[16] = {0};
  double be[16] = {0};

  //Calculate physical response and take care of summing.
  for(int i=0; i<mul; i++){
    int bin = FindEbin(energy[i]);
    double depositedEnergy;
    depositedEnergy = energy[i] + fPhys[bin].GetRandom();
    fe[fi[i]-1] += depositedEnergy > 0 ? depositedEnergy : -1;
    be[bi[i]-1] += depositedEnergy > 0 ? depositedEnergy : -1;
  }

  //Calculate electronic signal, check for triggering and store values.
  for(int i=0; i<16; i++){
    if(fe[i] > 0){
      double signal = fe[i] + fFront->GetRandom();
      signal = signal > 0 ? signal : 0;
      frontSignals[frontMul] = signal;
      frontStrips[frontMul] = i + 1;
      double eff = fThres->Eval(signal);
      double randomNum = rGen->Uniform();
      if(randomNum > eff) frontTriggers[frontMul] = false;
      else frontTriggers[frontMul] = true;
      frontMul++;
    }
    else{ fe[i] = 0.0;}
    if(be[i] > 0){
      double signal = be[i] + fBack->GetRandom();
      signal = signal > 0 ? signal : 0;
      backSignals[backMul] = signal;
      backStrips[backMul] = i + 1;
      double eff = fThres->Eval(signal);
      double randomNum = rGen->Uniform();
      if(randomNum > eff) backTriggers[backMul] = false;
      else backTriggers[backMul] = true;
      backMul++;
    }
    else{ be[i] = 0.0;}
  }

  if(stripSort){
    Sort(frontMul,frontStrips,frontOrder,false);
    Sort(backMul,backStrips,backOrder,false);
  }
}

int DetectorResponse::GetFrontMultiplicity()
{
  return frontMul;
}

double DetectorResponse::GetFrontSignal(int hit)
{
  return frontSignals[frontOrder[hit]];
}

int DetectorResponse::GetBackMultiplicity()
{
  return backMul;
}

double DetectorResponse::GetBackSignal(int hit)
{
  return backSignals[backOrder[hit]];
}

void DetectorResponse::SetThreshold(double level, double width)
{
  delete fThres;
  fThres = new TF1("fThres","0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1])))",0,10000);
  fThres->SetParameters(level,width);
}

bool DetectorResponse::HasFrontTrigger(int hit)
{
  return frontTriggers[frontOrder[hit]];
}

bool DetectorResponse::HasBackTrigger(int hit)
{
  return backTriggers[backOrder[hit]];
}

int DetectorResponse::GetFrontStrip(int hit)
{
  return frontStrips[frontOrder[hit]];
}

int DetectorResponse::GetBackStrip(int hit)
{
  return backStrips[backOrder[hit]];
}

int DetectorResponse::FindEbin(double E)
{
  int bin = Floor(E/(Emax/nBins));
  bin = bin < 1000 ? bin : 999;
  bin = bin >= 0 ? bin : 0;
  return bin;
}

void DetectorResponse::SortByStrip(bool input = true)
{
  stripSort = input;
  if(input == false){
    for(int i=0; i<16; i++){
      frontOrder[i] = i;
      backOrder[i] = i;
    }
  }
}
