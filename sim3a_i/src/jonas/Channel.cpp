#include "jonas/Channel.h"
#include "jonas/RMatrixUtils.h"
#include <TMath.h>
#include <TF1.h>
#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace TMath;

Channel::Channel(Pair pc, int lc, double gc, double rc, bool interpc)
{
  pair = pc;
  l = lc;
  g = gc;
  r = rc * pair.Radius();  //rc is the 'nucleon radius'.
  UseInterpolation(interpc);
}

Channel::~Channel(){}

void Channel::SetPair(Pair& p)
{
  pair = p;
}

double Channel::Penetrability(double E)
{
  /*
  * Calculates the penetrability. Input should be given in units of keV.  
  */
  if(E <= 0) return 0.0;

  if(interpolate) return FastPenetrability(E);
  else return ExactPenetrability(E);
}

double Channel::ShiftFunction(double E)
{
  if(interpolate) return FastShiftFunction(E);
  else return ExactShiftFunction(E);
}

double Channel::ShiftFunctionDeriv(double E)
{
  if(interpolate) return FastShiftDeriv(E);
  else return ExactShiftDeriv(E);
}

double Channel::ExactShiftFunction(double E)
{
  /*
  * Calculates the R-matrix Shift function. Input should be given in keV.  
  */

  //Calculate the Sommerfeld-parameter and k*r
  double eta = pair.Eta(E);
  double rho = Rho(E);

  //Is the energy positive or negative?
  double S = 0;
  double epsilon = 5.;  //Small energy for which we can still evaluate the Coulomb functions.
  if(E > epsilon){
    S = ShiftPos(l,eta,rho);
  }
  else if(E < -epsilon){
    S = ShiftNeg(l,eta,rho);
  }
  else{
    //We assume the shift function to be continuous and make a linear approximation around zero.
    double kappa = (E + epsilon) / (2. * epsilon); //Variable between zero and one.
    double kSmall = Sqrt(2 * epsilon * pair.RedMass()) / hbarc; //We take a small step on either side of zero.
    double etaSmall = pair.RedMass() / (kSmall * hbarc) * alpha * pair.QProduct();
    double rhoSmall = kSmall * r;
    S = kappa * ShiftPos(l,etaSmall,rhoSmall) + (1. - kappa) * ShiftNeg(l,etaSmall,rhoSmall);
  }

  return S;
}

double Channel::CoulombShift(double E)
{
  if(E <= 0){
    if(l == 0) return 0.;
    else if(l > 0) return Pi()/2.;
  }

  double eta = pair.Eta(E);
  double omega = 0.0;
  for(int i=1; i<=l; i++){ omega += ATan2(eta,i);}

  return omega;
}

double Channel::HardSphereShift(double E)
{
  if(E <= 0) return 0.0;

  if(interpolate) return FastHardSphere(E);
  else return ExactHardSphere(E);
}

double Channel::ExactHardSphere(double E)
{
  if(E <= 0) return 0.0;

  double eta = pair.Eta(E);
  double rho = Rho(E); 

  double F;
  double G;
  CoulombFunctions(l,eta,rho,F,G);
  
  //double phi = - ATan2(F,G);   //Not entirely sure about the minus-sign here.
  double phi = ATan2(F,G);       //29/6/2017: This is probably the correct sign.
  //double phi = F / G;
  return phi;
}

double Channel::PartialWidth(double E)
{
  return 2. * Penetrability(E) * Power(g,2);
}

void Channel::SetReducedWidth(double gc)
{
  g = gc;
}

void Channel::SetL(int lc)
{
  l = lc;
  if(interpolate) MakeInterpolation(emax,estep);
}

void Channel::SetRadius(double rc)
{
  r = rc;
  if(interpolate) MakeInterpolation(emax,estep);
}

int Channel::L()
{
  return l;
}

double Channel::Radius()
{
  return r;
}

const Pair & Channel::GetPair()
{
  return pair;
}

double Channel::ReducedWidth()
{
  return g;
}

double Channel::BarrierHeight()
{
  return pair.QProduct() * hbarc * alpha / r;
}

void Channel::MakeInterpolation(double Emax, double Estep)
{
  //cout << "Setting up interpolation..." << endl;
  if(Emax < 0. || Estep < 0.){
    cout << "  Channel::MakeInterpolation(): Bad interpolation parameters." << endl;
    exit (EXIT_FAILURE);
  }
  emax = Emax;
  estep = Estep;
  emin = -0.5*emax;
 
  /*
  * First step is to calculate the penetrability and interpolate the result. The penetrability
  * is interpolated in the region 0 <= E <= Emax, while the Shift-function is interpolated in
  * the region -0.5*Emax <= E <= Emax.
  */

  //XXX:First, we interpolate the penetrability
  vector<double> ev, fv;
  double e = 0;
  double delta = 1.;  //Value at zero energy is approximated with the value at 1keV.
  double f = ExactPenetrability(e + delta) / pair.GamowFactor(e + delta);  //XXX
  do{
    ev.push_back(e);
    fv.push_back(f);
    e += Estep;
    f = ExactPenetrability(e) / pair.GamowFactor(e); //XXX
    //cout << e << "  " << f << endl;
  }
  while(e < emax + 3*estep);  //We go a bit beyond the necessary limit.
  int N = ev.size();

  Pl = unique_ptr<ROOT::Math::Interpolator>(new ROOT::Math::Interpolator(N,ROOT::Math::Interpolation::kCSPLINE));
  Pl->SetData(ev, fv);


  //XXX:We interpolate the Shift-function
  ev.clear();
  fv.clear();
  e = emin - 3*estep;
  f = ExactShiftFunction(e);
  do{
    ev.push_back(e);
    fv.push_back(f);
    e += Estep;
    f = ExactShiftFunction(e);
    //cout << e << "  " << f << endl;
  }
  while(e < emax + 3*estep);  //We go a bit beyond the necessary limit.
  N = ev.size();

  S = unique_ptr<ROOT::Math::Interpolator>(new ROOT::Math::Interpolator(N,ROOT::Math::Interpolation::kCSPLINE));
  S->SetData(ev, fv);

  //XXX:We interpolate the derivative of the Shift-function
  ev.clear();
  fv.clear();
  e = emin - 3*estep;
  f = ExactShiftDeriv(e);
  do{
    ev.push_back(e);
    fv.push_back(f);
    e += Estep;
    f = ExactShiftDeriv(e);
    //cout << e << "  " << f << endl;
  }
  while(e < emax + 3*estep);  //We go a bit beyond the necessary limit.
  N = ev.size();

  dS = unique_ptr<ROOT::Math::Interpolator>(new ROOT::Math::Interpolator(N,ROOT::Math::Interpolation::kCSPLINE));
  dS->SetData(ev, fv);

  //XXX:We interpolate the hard-sphere phase shift. This is a bit tricky, since
  //    this function is discontinuous, but we just add 2*pi to the function whenever
  //    it changes sign. Being a phase shift, this shouldn't affect the final results.
  ev.clear();
  fv.clear();
  e = 0.;
  f = 0.;
  double fOld = f;
  int jumps = 0;
  do{
    ev.push_back(e);
    fv.push_back(f + jumps*2*Pi());
    e += Estep;
    f = ExactHardSphere(e);
    if(f < fOld) jumps++;
    fOld = f;
  }
  while(e < emax + 3*estep);  //We go a bit beyond the necessary limit.
  N = ev.size();

  Phi = unique_ptr<ROOT::Math::Interpolator>(new ROOT::Math::Interpolator(N,ROOT::Math::Interpolation::kCSPLINE));
  Phi->SetData(ev, fv);
}

double Channel::Rho(double E)
{
  double k = pair.WaveNumber(E);
  return k * r;
}

double Channel::ExactPenetrability(double E)
{
  //Calculate the Sommerfeld-parameter and k*r
  double eta = pair.Eta(E);
  double rho = Rho(E);
  return P(l,eta,rho);
}

double Channel::ExactShiftDeriv(double E)
{
  TF1 shiftFcn("ShiftFunction",[&](double* x, double* p) { return ExactShiftFunction(*x);},0,2*E,0);
  return shiftFcn.Derivative(E,0,0.0001);
}

double Channel::FastPenetrability(double E)
{
  CheckAndExpand(E);
  double f = Pl->Eval(E);
  return f * pair.GamowFactor(E);
}

double Channel::FastShiftFunction(double E)
{
  CheckAndExpand(E);
  return S->Eval(E);
}

double Channel::FastShiftDeriv(double E)
{
  CheckAndExpand(E);
  return dS->Eval(E);
}

double Channel::FastHardSphere(double E)
{
  CheckAndExpand(E);
  return Phi->Eval(E);
}

void Channel::UseInterpolation(bool interp, double Emax, double Estep)
{
  interpolate = interp;
  if(interpolate) MakeInterpolation(Emax, Estep);
}

void Channel::CheckAndExpand(double E)
{
  //Check whether energy is within bounds, otherwise expand range.
  if(E > emax) MakeInterpolation(1.5*E,estep);
  if(E < emin) MakeInterpolation(-3.*E,estep);
}
