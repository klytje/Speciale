#include <complex>
#include <iostream>
#include <cmath>
#include <memory>
#include <gsl/gsl_sf_coupling.h>
#include <TMath.h>
#include <TF1.h>
#include <ausa/util/memory>
#include "jonas/InterferenceWeight.h"
#include "jonas/Nucleus.h"
#include "jonas/Pair.h"
#include "jonas/SphericalHarmonic.h"

using namespace std;

InterferenceWeight::InterferenceWeight()
{
  DoCorrection(false);
}

InterferenceWeight::~InterferenceWeight() {}

double InterferenceWeight::Calculate(vector<TLorentzVector> &p)
{
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  Pair pPair(Be8,He4);
  Pair sPair(He4,He4);
  double alphaMass = He4.M();

  double Q = 0;
  for(auto & alpha : p){ Q += (alpha.Energy() - alphaMass);}       //Q-value of the alpha breakup.

  //This only applies to 12C, but we can only have Ja = 0, 1 or 2.
  complex<double> amplitude0(0,0);
  complex<double> amplitude1(0,0);
  complex<double> amplitude2(0,0);

  for(int k=0; k<weights.size(); k++){
    int Ja = primLevels.at(k).J();
    int L = primChannels.at(k)->L();
    int Jb = secLevels.at(k).J();
    double w = weights[k];

    //cout << "k = " << k << ",  w = " << w << ",  Ja = " << Ja << ",  L = " << L << ",  Jb = " << Jb ;

    double Eb = (secLevels.at(k).E() + Be8.M()) - 2 * He4.M();              //Energy of intermediate state above threshold.
    //cout << ",  Eb = " << Eb << endl;

    //Initial spin direction to be averaged over.
    for(int ma = -Ja; ma <= Ja ;ma++){
      //Symmetrisation in the alpha particle labels.
      for(int j=0; j<3; j++){    
        TLorentzVector alpha1 = p.at(j);
        TLorentzVector alpha2 = p.at((j+1)%3);
        TLorentzVector alpha3 = p.at((j+2)%3);

        double E1 = (alpha1.Energy() - alphaMass);// * 1e6;
        TLorentzVector rcm = alpha2 + alpha3;     //Recoil center of mass system.

        double theta1 = alpha1.Theta();
        double phi1 = alpha1.Phi();
        alpha2.Boost(-rcm.BoostVector());
        double theta2 = alpha2.Theta();
        double phi2 = alpha2.Phi();
        alpha2.Boost(rcm.BoostVector());

        rcm.Boost(-rcm.BoostVector());
        double E23 = (rcm.Energy() - 2 * alphaMass);// * 1e6;  //Relative energy of alpha2 and alpha3.
        double E = Q - E23;                                 //Relative energy of alpha1 and recoil.

        double Gamma1 = primChannels.at(k)->PartialWidth(E);

        //XXX:Coulomb correction.
        if(doCorrection){
          //For the Coulomb correction we need the relative energies of alpha 1-2 and alpha 1-3.
          TLorentzVector cm12 = alpha1 + alpha2;
          cm12.Boost(-cm12.BoostVector());
          double E12 = cm12.Energy() - 2 * alphaMass;
          TLorentzVector cm13 = alpha1 + alpha3;
          cm13.Boost(-cm13.BoostVector());
          double E13 = cm13.Energy() - 2 * alphaMass;
          
          double rtilde = CorrectionRadius(primChannels.at(k),secChannels.at(k),secLevels.at(k),E,E23);
          c1.at(k)->SetRadius(rtilde);
          c2.at(k)->SetRadius(rtilde);

          //cout << "E12 = " << E12 << endl;
          //cout << "E13 = " << E13 << endl;
 
          //Then we multiply by the correction factors.
          Gamma1 *= c1.at(k)->Rho(E) / (c2.at(k)->Rho(E12) * c2.at(k)->Rho(E13));
          Gamma1 *= c2.at(k)->Penetrability(E12) * c2.at(k)->Penetrability(E13) / c1.at(k)->Penetrability(E);
        }
        //XXX

        double Gamma2 = secChannels.at(k)->PartialWidth(E23);
        //cout << "Gamma1 = " << Gamma1 << endl;
        //cout << "Gamma2 = " << Gamma2 << endl;
        double phase1 = primChannels.at(k)->CoulombShift(E) - primChannels.at(k)->HardSphereShift(E);
        double phase2 = secChannels.at(k)->CoulombShift(E23) - secChannels.at(k)->HardSphereShift(E23);
 
        //We start calculating...
        complex<double> i(0,1);
        complex<double> a(0,0);
 
        //XXX: We have modified the expression to include rho's instead of sqrt(E)'s.
        a += sqrt(Gamma1 * Gamma2 / (primChannels.at(k)->Rho(E) * secChannels.at(k)->Rho(E23))) * exp(i * phase1) * exp(i * phase2);
        a /= Eb - pow(secChannels.at(k)->ReducedWidth(),2) * (secChannels.at(k)->ShiftFunction(E23) - secChannels.at(k)->ShiftFunction(Eb)) - E23 - 0.5 * i * Gamma2;
        complex<double> sum(0,0);
        for(int mb = -Jb ; mb <= Jb; mb ++){
          complex<double> Cmmj(ClebschGordan(L,Jb,Ja,ma-mb,mb,ma),0);
          complex<double> Ylm1(spherical_harmonic(L,ma-mb,theta1,phi1));
          complex<double> Ylm2(spherical_harmonic(secChannels.at(k)->L(),mb,theta2,phi2));
          sum += Cmmj * Ylm1 * Ylm2;
        }
        //cout << "contribution = " << w * a * sum << endl;
        if(Ja == 0){ amplitude0 += w * a * sum;}// cout << " a0 = " << w * a * sum << endl;}
        else if(Ja == 1){ amplitude1 += w * a * sum;}
        else if(Ja == 2){ amplitude2 += w * a * sum;}// cout << " a2 = " << w * a * sum << endl;}
        //amplitude += w * a * sum;
        //cout << "amplitude = " << amplitude << endl;
      }
    }
  }

  double weight = norm(amplitude0) + norm(amplitude1) + norm(amplitude2);
  //cout << "norm(0) = " << norm(amplitude0) << ",  norm(2) = " << norm(amplitude2) << endl;
  return weight;
}

void InterferenceWeight::SetContributions(vector<vector<int>> models, vector<double> w, double r0)
{
  if(models.size() != w.size()) return;
//  contributions = models;
  weights = w;

  //cout << "weights.size() = " << weights.size() << endl;

  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  Pair pPair(Be8,He4);
  Pair sPair(He4,He4);
  double primRadius = r0;
  double secRadius = r0;
  double primEnergy = 1000.; //Not relevant.
  double secEnergy0 = 0.0;   //Be8(gs).
  double secEnergy2 = 3037.; //Be8(ex).
  double primWidth = 1.;     //Not relevant.
  double secWidth0 = 28.81;   //Sqrt of reduced width of Be8(0+). //XXX: Is in fact dependent on channel radius
  double secWidth2 = 32.787; //Sqrt of reduced width of Be8(2+).  //XXX: and should be calculated from the partial width.

  for(int i=0; i<models.size(); i++){
    int Ja = models.at(i).at(0);
    int L = models.at(i).at(1);
    int Jb = models.at(i).at(2);

    Level primLevel(primEnergy,Ja,1);                  //energy not relevant.
    unique_ptr<Channel> pci = std::make_unique<Channel>(pPair,L,primWidth,primRadius); //reduced width not really important here, arbitrary.
    unique_ptr<Channel> sci;
    Level secLevel(0.0,Jb,1);
    if(Jb == 0){
      secLevel.SetEnergy(secEnergy0);
      sci = std::make_unique<Channel>(sPair,Jb,secWidth0,secRadius);
    }
    else if(Jb == 2){
      secLevel.SetEnergy(secEnergy2);
      sci = std::make_unique<Channel>(sPair,Jb,secWidth2,secRadius);
    }

    primLevels.push_back(primLevel);
    primChannels.push_back(std::move(pci));
    secLevels.push_back(secLevel);
    secChannels.push_back(std::move(sci));
  }
}

double InterferenceWeight::ClebschGordan(double j1, double j2, double j3, double m1, double m2, double m3)
{
  m3=-m3;
  int j1x2=(int)lrint(2*j1);
  int j2x2=(int)lrint(2*j2);
  int j3x2=(int)lrint(2*j3);
  int m1x2=(int)lrint(2*m1);
  int m2x2=(int)lrint(2*m2);
  int m3x2=(int)lrint(2*m3);

  double w3j=gsl_sf_coupling_3j(j1x2,j2x2,j3x2,m1x2,m2x2,m3x2);

  return pow(-1.0,lrint(j1-j2-m3)) * sqrt(2.0*j3+1.) * w3j;
}

void InterferenceWeight::DoCorrection(bool input, double radius)
{
  doCorrection = input;
  if(!doCorrection) return;

  //Cannot do correction if the channels are not defined.
  if(primChannels.size() == 0){
    cout << "  InterferenceWeight::DoCorrection(): Channels need to be defined!" << endl;
    exit(EXIT_FAILURE);
  }

  //Correction for Coulomb repulsion.
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  Pair p1(He4,Be8);
  Pair p2(He4,He4);
  for(int i=0; i<primChannels.size(); i++){
    int L1 = primChannels.at(i)->L();
    int L2 = secLevels.at(i).J();
    unique_ptr<Channel> c1i = unique_ptr<Channel>(new Channel(p1,L1,1.,radius/p1.Radius()));
    unique_ptr<Channel> c2i = unique_ptr<Channel>(new Channel(p2,L2,1.,radius/p2.Radius())); //We guess that the L's are the same as the 8Be spin.
    c1.push_back(move(c1i));
    c2.push_back(move(c2i));
  }
}

void InterferenceWeight::UseInterpolation(bool input)
{
  if(primChannels.size() == 0){
    cout << "  InterferenceWeight::UseInterpolation(): Define channels first!" << endl;
    exit(EXIT_FAILURE);   
  }

  for(int i=0; i<primChannels.size(); i++) primChannels.at(i)->UseInterpolation();
  for(int i=0; i<secChannels.size(); i++) secChannels.at(i)->UseInterpolation();
  for(int i=0; i<c1.size(); i++) c1.at(i)->UseInterpolation();
  for(int i=0; i<c2.size(); i++) c2.at(i)->UseInterpolation();
}

double InterferenceWeight::CorrectionRadius(unique_ptr<Channel> & cp, unique_ptr<Channel> & cs, Level & l, double E, double E23)
{
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);

  auto phase = [&] (double e){
    double Eb = (l.E() + Be8.M()) - 2 * He4.M(); //Energy of intermediate state above threshold.
    double phi = cs->HardSphereShift(e);
    double Gamma = cs->PartialWidth(e);
    double S = cs->ShiftFunction(e);
    double B = cs->ShiftFunction(Eb);
    double gamma2 = pow(cs->ReducedWidth(),2);

    double arg = (Gamma/2.) / (Eb - e - (S - B) * gamma2);
    double y = abs(arg);
    double x = TMath::Sign(1.,arg);

    double delta = atan2(y,x) - phi;
    return delta;   
  };

  const Pair& p1 = cp->GetPair();
  const Pair& p2 = cs->GetPair();
  double r1 = cp->Radius();
  double m1 = p1.RedMass();
  double r2 = cs->Radius();
  double m2 = p2.RedMass();

  TF1 delta("delta",[&](double* x, double* p) { return phase(*x);},0,2*E23,0);
  double deriv = delta.Derivative(E23,0,0.000001);
 
  double rtilde = r1 + hbarc * sqrt(2.*E/m1) * deriv + r2 * sqrt(E*m2/(E23*m1));

  return rtilde;
}
