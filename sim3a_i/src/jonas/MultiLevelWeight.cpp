#include <vector>
#include <complex>
#include <gsl/gsl_sf_coupling.h>
#include <armadillo>

#include "jonas/Nucleus.h"
#include "jonas/Pair.h"
#include "jonas/Channel.h"
#include "jonas/MultiLevelWeight.h"
#include "jonas/SphericalHarmonic.h"

#include <TLorentzVector.h>
#include <TF2.h>
#include <TH2D.h>
#include <TMath.h>

using namespace std;
using namespace arma;
using namespace TMath;


/*
Special class calculating weight for 1+ decay through several states in 8Be.
*/
MultiLevelWeight::MultiLevelWeight() {}

MultiLevelWeight::~MultiLevelWeight() {}

void MultiLevelWeight::SetPrimary(Level& l, unique_ptr<Channel> c)
{
  primLevel = l;
  primChannel = std::move(c);
}

void MultiLevelWeight::SetSecondary(Level& l, unique_ptr<Channel> c)
{
  secLevel = l;
  secChannel = std::move(c);
}

double MultiLevelWeight::Calculate(vector<TLorentzVector> &p)
{
  int Ja = 1;
  int Jb = 2;
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  double alphaMass = He4.M();// / 1e6;  //In GeV.
  double Egs = 91.84;
  double energies[3] = {3037., 16626., 16922.}; //energies above ground state.
  //double widths[3] = {1075., 10.89, 7.46};
  double widths[3] = {1075., 10.89, 7.46};
  int signs[3] = {1,1,1};

  //Primary channel
  Pair pPair(Be8,He4);
  Channel primChannel(pPair,2,35.,5.1);  //reduced width not really important here, arbitrary.
  
  //The second channel
  Pair sPair(He4,He4);
  Channel secChannel(sPair,2,32.787,4.5);

  double Q = 0;
  for(auto & alpha : p){ Q += (alpha.Energy() - alphaMass);}// * 1e6;}  //Q-value of the alpha breakup.

  //XXX double Eb = (secLevel.E() + Be8.M()) - 2 * He4.M();                 //Energy of intermediate state above threshold.


  double prob = 0;

  //Initial spin direction to be averaged over.
  for(int ma = -Ja; ma <= Ja ; ma++){
    //cout << "ma = " << ma << endl;
    complex<double> amplitude(0,0);
    for(int j=0; j<3; j++){
      //cout << "j = " << j << endl;
      //First some kinematics...
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

      //cout << "E = " << E << endl;
      //cout << "E23 = " << E23 << endl;

      //We start calculating the weight function.
      double Gamma1 = primChannel.PartialWidth(E);
      //cout << "Gamma1 = " << Gamma1 << endl;
      //double Gamma2 = secChannel.PartialWidth(E23);

      double phase1 = primChannel.CoulombShift(E) - primChannel.HardSphereShift(E);
      double phase2 = secChannel.CoulombShift(E23) - secChannel.HardSphereShift(E23);
      //cout << "Phase1 = " << phase1 << endl;
      //cout << "Phase2 = " << phase2 << endl;

      complex<double> i(0,1);
      //complex<double> a(0,0);
 
      //a += sqrt(Gamma1 * Gamma2 / sqrt(E1 * E23)) * exp(i * phase1) * exp(i * phase2);
      //a /= Eb - pow(secChannel.ReducedWidth(),2) * (secChannel.ShiftFunction(E23) - secChannel.ShiftFunction(Eb)) - E23 - 0.5 * i * Gamma2;
   
      //We calculate some standard R-matrix functions.
      double Sc = secChannel.ShiftFunction(E23);
      double Bc = secChannel.ShiftFunction(energies[0] + Egs);  //Boundary = shift at first level.
      double Pc = secChannel.Penetrability(E23);

      //We construct the inverse level matrix.
      Mat<complex<double>> Ainv(3,3);
      for(int lambda=0; lambda<3; lambda++){
        for(int mu=0; mu<3; mu++){
          complex<double> val(0,0);
          if(lambda == mu) val += energies[lambda] - (E23 - Egs);      //level energies are relative to ground state.
          val -= (Sc - Bc + i*Pc) * sqrt(widths[lambda] * widths[mu]);  //sqrt since widths are actually gamma^2.
          Ainv(lambda,mu) = val;
        }
      }

      //cout << "Ainv = " << Ainv << endl;
      //Then we find the level matrix proper.
      Mat<complex<double>> A = inv(Ainv);
      //cout << "A = " << A << endl;

      //Calculate the resonance amplitude
      complex<double> a(0,0);
      for(int mu = 0; mu<3; mu++){
        double Gamma2 = 2 * Pc * widths[mu];
        //cout << "Gamma2 = " << Gamma2 << endl;
        a += sqrt(Gamma2) * A(0,mu);
        //cout << "amu = " << a << endl;
      }
      //cout << "amu = " << a << endl;
      //a *= sqrt(Gamma1 / sqrt(E1 * E23)) * exp(i * phase1) * exp(i * phase2);
      complex<double> factor(0,0);
      //cout << "Gamma1 = " << Gamma1 << endl;
      //cout << "E = " << E << endl;
      //cout << "E23 = " << E23 << endl;
      //cout << "phase1 = " << phase1 << endl;
      //cout << "phase2 = " << phase2 << endl;
      factor += sqrt(Gamma1 / sqrt(E * E23)) * exp(i * phase1) * exp(i * phase2);
      //cout << "factor = " << factor << endl;
      //a *= sqrt(Gamma1 / sqrt(E * E23)) * exp(i * phase1) * exp(i * phase2);
      a *= factor;
      
      //cout << "a = " << a << endl;

      //Multiply by geometrical part of the amplitude.
      complex<double> sum(0,0);
      for(int mb = -Jb ; mb <= Jb; mb ++){
        complex<double> Cmmj(ClebschGordan(primChannel.L(),Jb,Ja,ma-mb,mb,ma),0);
        complex<double> Ylm1(spherical_harmonic(primChannel.L(),ma-mb,theta1,phi1));
        complex<double> Ylm2(spherical_harmonic(secChannel.L(),mb,theta2,phi2));
        sum += Cmmj * Ylm1 * Ylm2;
      }
      amplitude += a * sum;
    }
    prob += norm(amplitude);
  }

  return prob;
}


double MultiLevelWeight::ClebschGordan(double j1, double j2, double j3, double m1, double m2, double m3)
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
