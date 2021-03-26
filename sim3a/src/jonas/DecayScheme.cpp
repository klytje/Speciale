#include "jonas/DecayScheme.h"
#include "jonas/Nucleus.h"
#include <algorithm>  //std::sort
#include <iostream>
#include <stdlib.h>  //exit(), EXIT_FAILURE
#include <float.h>
#include <TMath.h>
#include <TF1.h>

using namespace std;

DecayScheme::DecayScheme(vector<int> channels, vector<SecondaryLevel> levels, double r1, double r2)
{
//  double r0 = 1.42; //We construct with default channel radius, but allow for later modification.
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  Pair p1(Be8,He4);
  Pair p2(He4,He4);

  //We construct the primary channels
  for(int l : channels){
    unique_ptr<Channel> cl = unique_ptr<Channel>(new Channel(p1,l,1.,r1,true));
    primChannels.push_back(move(cl));
  }

  //Then we count how many secondary channels we need.
  vector<int> ls; 
  for(SecondaryLevel level : levels){
    int li = level.J;
    bool foundL = false;
    for(int l : ls) if(l == li) foundL = true;
    if(!foundL) ls.push_back(li);
  }

  //The channels are sorted in order of ascending l.
  sort(ls.begin(), ls.end());

  //We check that the user provided secondary levels with the proper number of
  //widths (somewhat stupid).
  for(SecondaryLevel level : levels){
    if(level.widths.size() != ls.size()){
      cout << "DecayScheme: Number of widths doesn't match number of secondary channels! N(widths) = " << level.widths.size() << ",  N(channels) = " << ls.size() << endl;
      exit(EXIT_FAILURE);
    }
  }

  //Everything seems ok, we adopt the secondary levels.
  secLevels = levels;

  //We construct the secondary channels. 
  for(int l : ls){
    unique_ptr<Channel> cl = unique_ptr<Channel>(new Channel(p2,l,1.,r2,true));
    secChannels.push_back(move(cl));  
  }

  //Finally, we calculate the boundary conditions, and we're ready to go.
  CalculateBoundaries();
}

DecayScheme::~DecayScheme() {}

void DecayScheme::AddPrimaryLevel(PrimaryLevel level)
{
  //We check that the user provided secondary levels with the proper width matrix.
  if(level.GetWidths().n_rows != primChannels.size() || level.GetWidths().n_cols != secLevels.size()){
    cout << "DecayScheme: Dimension of width matrix is wrong!" ;
    cout << "  N(rows) = " << level.GetWidths().n_rows << ",  N(channels) = " << primChannels.size();
    cout << ",  N(cols) = " << level.GetWidths().n_cols << ",  N(levels) = " << secLevels.size() << endl;
    exit(EXIT_FAILURE);
  }
  //Seems ok, we adopt the level in our scheme.
  primLevels.push_back(level);
  CalculateBoundaries();
}

void DecayScheme::SetPrimaryRadius(double r0)
{  
  for(int c=0; c<primChannels.size(); c++){
    double rp = primChannels.at(c)->GetPair().Radius();
    primChannels.at(c)->SetRadius(r0 * rp);
  }
  CalculateBoundaries();
}

void DecayScheme::SetSecondaryRadius(double r0)
{
  for(int c=0; c<secChannels.size(); c++){
    double rp = secChannels.at(c)->GetPair().Radius();
    secChannels.at(c)->SetRadius(r0 * rp);
  }
  CalculateBoundaries();
}

void DecayScheme::CalculateBoundaries()
{
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  Nucleus C12(12,6,0,1);

  primBoundaries = Mat<double>(primChannels.size(),secLevels.size());
  secBoundaries = vector<double>(secChannels.size());

  //Set thresholds in order to get the 'channel energies'.
  double primThreshold = He4.M() + Be8.M() - C12.M();  //Positive number ~7366keV
  double secThreshold = 2. * He4.M() - Be8.M();  //Negative number (Be8 is unbound) ~ -92keV

  //We establish boundary conditions for primary channels distinguished by l and mu.
  for(int il=0; il<primChannels.size(); il++){
    for(int mu=0; mu<secLevels.size(); mu++){
      //For each channel we identify the lowest-energy level with a
      //non-zero widdth for decay through the respective channel.
      int lambda0 = -1;
      double Emin = DBL_MAX;
      for(int i=0; i<primLevels.size(); i++){
        double gammac = primLevels.at(i).GetWidth(il,mu);
        if(!IsAlmostZero(gammac)){
          double El = primLevels.at(i).El;
          if(El < Emin){
            Emin = El;
            lambda0 = i;
          }
        }
      }
      if(lambda0 >= 0){
        //Ok, lambda0 is the relevant level.
        //cout << "Lowest-energy level for the (" << il << "," << mu << ")-channel is " ;
        //cout << "primLevel(" << lambda0 << ") at " << primLevels.at(lambda0).El << "keV, ";
        double Bc = primChannels.at(il)->ShiftFunction(primLevels.at(lambda0).El - primThreshold);
        primBoundaries(il,mu) = Bc;
      }
      else{
        primBoundaries(il,mu) = 0.;
      }
      //cout << "Bc(" << il << "," << mu << ") = " << primBoundaries(il,mu) << endl;
    }
  }

  //Next, we calculate the boundary conditions for the secondary system.
  for(int il=0; il<secChannels.size(); il++){
    //For each channel we identify the lowest-energy level with a
    //non-zero widdth for decay through the respective channel.
    int lambda0 = -1;
    double Emin = DBL_MAX;
    for(int i=0; i<secLevels.size(); i++){
      double gammac = secLevels.at(i).widths.at(il);
      if(!IsAlmostZero(gammac)){
        double El = secLevels.at(i).El;
        if(El < Emin){
          Emin = El;
          lambda0 = i;
        }
      }
    }
    if(lambda0 >= 0){
      //Ok, lambda0 is the relevant level.
      //cout << "Lowest-energy level for the " << il <<"-channel is " ;
      //cout << "secLevel(" << lambda0 << ") at " << secLevels.at(lambda0).El << "keV, ";
      double Bc = secChannels.at(il)->ShiftFunction(secLevels.at(lambda0).El - secThreshold);
      secBoundaries.at(il) = Bc;
      //cout << ",  Bc = " << Bc << endl;
    }
    else{
      secBoundaries.at(il) = 0.;
    }
    //cout << "Bc(" << il << ") = " << secBoundaries.at(il) << endl;
  }  
}

vector<PrimaryLevel> & DecayScheme::GetPrimaryLevels()
{
  return primLevels;
}

vector<unique_ptr<Channel>> & DecayScheme::GetPrimaryChannels()
{
  return primChannels;
}

vector<SecondaryLevel> & DecayScheme::GetSecondaryLevels()
{
  return secLevels;
}

vector<unique_ptr<Channel>> & DecayScheme::GetSecondaryChannels()
{
  return secChannels;
}

double DecayScheme::GetPrimaryBoundary(int il, int mu)
{
  return primBoundaries(il,mu);
}

double DecayScheme::GetSecondaryBoundary(int il)
{
  return secBoundaries.at(il);
}

bool DecayScheme::IsAlmostZero(double num)
{
  int smudgeFactor = 10;
  if(abs(num) < smudgeFactor * DBL_EPSILON) return true;
  else return false;
}

//A dirty hack, which makes one-level, one-channel approximation for the secondary system.
double DecayScheme::DensityFunction(int ib, double E23)
{
  int il = ib;  //Only works if the first level decays through the first channel, etc.
  double El = secLevels.at(ib).El + 91.84;
  double gamma = secLevels.at(ib).widths.at(il);
  double gamma2 = TMath::Power(gamma,2);
  double P = secChannels.at(il)->Penetrability(E23);
  double S = secChannels.at(il)->ShiftFunction(E23);
  double B = GetSecondaryBoundary(il);

  //cout << "El = " << El << ",  gamma2 = " << gamma2 << ",  P = " << P << ",  S = " << S << ",  B = " << B << endl;

  return 1./TMath::Pi() * P * gamma2 / (TMath::Power(El-E23-gamma2*(S - B),2) + TMath::Power(P * gamma2,2));
}

double DecayScheme::Norm(int ib, double E3a)
{
  auto integrand = [&] (double e){
    return DensityFunction(ib,e);   
  };
  TF1 f("dos",[&](double* x, double* p) { return integrand(*x);},0,E3a,0);
  double N = 0.;
  if(E3a < 100.){ N = f.Integral(10.,E3a);}
  else{
    N = f.Integral(10.,100.); 
    N += f.Integral(100.,E3a);
  }
  return N;
}

double DecayScheme::AvgPenetrability(int il, int ib, double E3a)
{
  auto integrand = [&] (double e){
    double P = primChannels.at(il)->Penetrability(E3a-e);
    double rho = DensityFunction(ib,e);
    return P * rho;
  };
  TF1 f("i1",[&](double* x, double* p) { return integrand(*x);},0,E3a,0);
  double Pm = 0.;
  if(E3a < 100.){ Pm = f.Integral(10.,E3a);}
  else{
    Pm = f.Integral(10.,100.); 
    Pm += f.Integral(100.,E3a);
  }
  return Pm;// / Norm(ib,E3a); 
}

double DecayScheme::AvgShiftFunction(int il, int ib, double E3a)
{
  auto integrand = [&] (double e){
    double S = primChannels.at(il)->ShiftFunction(E3a-e);
    double rho = DensityFunction(ib,e);
    return S * rho;
  };
  TF1 f("i2",[&](double* x, double* p) { return integrand(*x);},0,E3a,0);
  double Sm = 0.;
  if(E3a < 100.){ Sm = f.Integral(10.,E3a);}
  else{
    Sm = f.Integral(10.,100.); 
    Sm += f.Integral(100.,E3a);
  }
  return Sm / Norm(ib,E3a);  
}
