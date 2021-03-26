#include <iostream>
#include <cmath>
#include <float.h>
#include <gsl/gsl_sf_coupling.h>
#include <TMath.h>
#include <TF1.h>
#include "jonas/GrandWeight.h"
#include "jonas/SphericalHarmonic.h"
#include <logft/phase_space.h>

using namespace std;
using namespace logft;

GrandWeight::GrandWeight(unique_ptr<DecayScheme> input)
{
  scheme = std::move(input);

  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  Nucleus C12(12,6,0,1);

  //Set thresholds in order to get the 'channel energies'.
  primThreshold = He4.M() + Be8.M() - C12.M();  //Positive number ~7366.59keV
  secThreshold = 2. * He4.M() - Be8.M();  //Negative number (Be8 is unbound) ~ -91.84keV
  tripleThreshold = 3. * He4.M() - C12.M();   //Positive number ~7274.75keV

  doCorrection = false;
}

GrandWeight::~GrandWeight() {}

double GrandWeight::ClebschGordan(double j1, double j2, double j3, double m1, double m2, double m3)
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

int GrandWeight::delta(int i, int j)
{
  if(i == j) return 1;
  else return 0;
}

bool GrandWeight::IsAlmostZero(double num)
{
  int smudgeFactor = 10;
  if(abs(num) < smudgeFactor * DBL_EPSILON) return true;
  else return false;
}

vector<int> GrandWeight::GetPrimaryIndices(int Ja)
{
  vector<int> indices;
  vector<PrimaryLevel> &primLevels = scheme->GetPrimaryLevels();
  for(int i=0; i<primLevels.size(); i++){
    if(primLevels.at(i).J == Ja){
      indices.push_back(i);
    }
  }

  return indices;
}

vector<int> GrandWeight::GetSecondaryIndices(int Jb)
{
  vector<int> indices;
  vector<SecondaryLevel> &secLevels = scheme->GetSecondaryLevels();
  for(int i=0; i<secLevels.size(); i++){
    if(secLevels.at(i).J == Jb){
      indices.push_back(i);
    }
  }

  return indices;
}

Mat<complex<double>> GrandWeight::ConstructPrimaryMatrix(int Ja, double Q3a)
{
  //We calculate with energies relative to the gs.
  double E = Q3a + tripleThreshold;
  //Define the complex unit.
  complex<double> i(0.,1.);

  vector<PrimaryLevel> & primLevels = scheme->GetPrimaryLevels();

  //First, we figure out which levels to include in the matrix.
  vector<int> indices = GetPrimaryIndices(Ja);

  //Then we construct the (inverse) level matrix.
  int N = indices.size();
  Mat<complex<double>> A = Mat<complex<double>>(N,N);
  
  for(int lambda = 0; lambda<N; lambda++){
    int ilambda = indices.at(lambda);
    double Elambda = primLevels.at(ilambda).El;  //Elambda is energy above gs.
    for(int mu = 0; mu<N; mu++){
      int imu = indices.at(mu);
      A(lambda,mu) = (Elambda - E) * delta(lambda,mu);
      //We need a loop over all channels.
      for(int il=0; il<scheme->GetPrimaryChannels().size(); il++){
        for(int mup=0; mup<scheme->GetSecondaryLevels().size(); mup++){
          double Ec = E - primThreshold;  //Channel energy.
          double gammalc = primLevels.at(ilambda).GetWidth(il,mup);
          double gammamc = primLevels.at(imu).GetWidth(il,mup);
          double Sc = scheme->GetPrimaryChannels().at(il)->ShiftFunction(Ec);
          double Bc = scheme->GetPrimaryBoundary(il,mup);
          double Pc = scheme->GetPrimaryChannels().at(il)->Penetrability(Ec);
          A(lambda,mu) -= gammalc * gammamc * (Sc - Bc + i * Pc);
        }
      }   
    }
  }

  //We invert to get get the proper level matrix
  A = inv(A);

  return A;
}


Mat<complex<double>> GrandWeight::ConstructSecondaryMatrix(int Jb, double Q2a)
{
  //We calculate energies relative to the g.s.
  double E = Q2a + secThreshold;
  //Define the complex unit
  complex<double> i(0.,1.);

  vector<SecondaryLevel> & secLevels = scheme->GetSecondaryLevels();

  //First, we figure out which levels to include in the matrix.
  vector<int> indices = GetSecondaryIndices(Jb);

  /*
  cout << "Indices to include in level matrix: ";
  for(int i : indices) cout << i << ", " ;
  cout << endl;
  */

  //Then we construct the (inverse) level matrix.
  int N = indices.size();
  Mat<complex<double>> A = Mat<complex<double>>(N,N);
  
  for(int lambda = 0; lambda<N; lambda++){
    int ilambda = indices.at(lambda);
    double Elambda = secLevels.at(ilambda).El;  //Elambda is energy above gs.
    for(int mu = 0; mu<N; mu++){
      int imu = indices.at(mu);
      A(lambda,mu) = (Elambda - E) * delta(lambda,mu);
      //cout << "(Elambda - E)*delta = " << (Elambda - E) * delta(lambda,mu) ;
      //We loop over channels.
      for(int il=0; il<scheme->GetSecondaryChannels().size(); il++){
        double Ec = E - secThreshold;  //Channel energy.
        double gammalc = secLevels.at(ilambda).widths.at(il);
        double gammamc = secLevels.at(imu).widths.at(il);
        double Sc = scheme->GetSecondaryChannels().at(il)->ShiftFunction(Ec);
        double Bc = scheme->GetSecondaryBoundary(il);
        double Pc = scheme->GetSecondaryChannels().at(il)->Penetrability(Ec);
        //cout << ", Ec = " << Ec << ", Sc = " << Sc << ", Bc = " << Bc << ", Pc = " << Pc << endl;
        A(lambda,mu) -= gammalc * gammamc * (Sc - Bc + i * Pc);
      }   
    }
  }

  //cout << "Inverse level matrix: " << endl << A << endl;

  //We invert to get get the proper level matrix
  A = inv(A);

  //cout << "Proper level matrix: " << endl << A << endl;

  return A;
}

void GrandWeight::DoCorrection(double rt)
{
  doCorrection = true;
  correctionRadius = rt;

  //We construct the channels with variable channel radius.
  for(int il=0; il<scheme->GetPrimaryChannels().size(); il++){
    //Nucleus n1 = scheme->GetPrimaryChannels().at(il)->GetPair().Particles().at(0);
    //Nucleus n2 = scheme->GetPrimaryChannels().at(il)->GetPair().Particles().at(1);
    //Nucleus nn1(n1.A(),n1.Q(),n1.J(),n1.Pi());
    //Nucleus nn2(n2.A(),n2.Q(),n2.J(),n2.Pi());
    //Nucleus nn1(8,4,0,1);
    //Nucleus nn2(4,2,0,1);
    //Pair p(nn1,nn2);
    Pair p = scheme->GetPrimaryChannels().at(il)->GetPair();
    int l = scheme->GetPrimaryChannels().at(il)->L();
    unique_ptr<Channel> c1i = unique_ptr<Channel>(new Channel(p,l));
    c1i->SetRadius(rt);
    c1.push_back(std::move(c1i));
  }
  //Also for the secondary channels.
  for(int ilp=0; ilp<scheme->GetSecondaryChannels().size(); ilp++){
    //Nucleus n1 = scheme->GetSecondaryChannels().at(ilp)->GetPair().Particles().at(0);
    //Nucleus n2 = scheme->GetSecondaryChannels().at(ilp)->GetPair().Particles().at(1);
    //Nucleus nn1(n1.A(),n1.Q(),n1.J(),n1.Pi());
    //Nucleus nn2(n2.A(),n2.Q(),n2.J(),n2.Pi());
    //Nucleus nn1(4,2,0,1);
    //Nucleus nn2(4,2,0,1);
    //Pair p(nn1,nn2);
    Pair p = scheme->GetSecondaryChannels().at(ilp)->GetPair();
    int l = scheme->GetSecondaryChannels().at(ilp)->L();
    unique_ptr<Channel> c2i = unique_ptr<Channel>(new Channel(p,l));
    c2i->SetRadius(rt);
    c2.push_back(std::move(c2i));
  }
}


//This implementation prevents the interpolation to be used in the correction channels.
double GrandWeight::CoulombCorrection(double E, double E23, double E12, double E13, int il, int ilp)
{
  if(!doCorrection) return 1.;

  double rtilde;
  //We figure out if we need to do the final-state Coulomb correction with variable radius.
  if(correctionRadius < 0.){
    rtilde = CorrectionRadius(E,E23,il,ilp);
  }
  else{
    rtilde = correctionRadius;
  }

  c1.at(il)->SetRadius(rtilde);  //Primary decay channel with channel radius at rtilde.
  c2.at(ilp)->SetRadius(rtilde);  //Secondary decay channel with channel radius at rtilde.
 
  double FSCI = c1.at(il)->Rho(E) / (c2.at(ilp)->Rho(E12) * c2.at(ilp)->Rho(E13))
              * c2.at(ilp)->Penetrability(E12) * c2.at(ilp)->Penetrability(E13) / c1.at(il)->Penetrability(E);

  //cout << "E23 = " << E23 << ",  E12 = " << E12 << ",  E13 = " << E13 << endl;
  return FSCI;
}

double GrandWeight::CorrectionRadius(double E, double E23, int il, int ilp)//(unique_ptr<Channel> & c, Level & l, double E, double E23)
{
  int L = scheme->GetSecondaryChannels().at(ilp)->L();
  vector<int> sIndices = GetSecondaryIndices(L);

  auto phase = [&] (double e){
    double Pc = scheme->GetSecondaryChannels().at(ilp)->Penetrability(e);
    double Sc = scheme->GetSecondaryChannels().at(ilp)->ShiftFunction(e);
    double Bc = scheme->GetSecondaryBoundary(ilp);
    double phic = scheme->GetSecondaryChannels().at(ilp)->HardSphereShift(e);
    double sum = 0.;
    for(int m=0; m<sIndices.size(); m++){    //Level index with respect to other levels of the same spin.
      int lambda = sIndices.at(m);           //Global level index;
      double gammalc = scheme->GetSecondaryLevels().at(lambda).widths.at(ilp);
      double El = scheme->GetSecondaryLevels().at(lambda).El; //El is energy above gs.
      El -= secThreshold;  //Now with respect to threshold;
      sum += TMath::Power(gammalc,2) / (El - e);
      //cout  << "El = " << El << "Gamma = " << 2. * Pc * TMath::Power(gammalc,2) << ",  gamma2 = " << TMath::Power(gammalc,2);
    }

    //cout << ",  Sc = " << Sc << ",  Bc = " << Bc << ",  phic = " << phic;

    double arg = Pc / (1./sum - Sc + Bc);
    double y = abs(arg);
    double x = TMath::Sign(1.,arg);

    double del = atan2(y,x) - phic;
    //cout << ",  d = " << del << endl;
    return del;   
  };
  TF1 delta("delta",[&](double* x, double* p) { return phase(*x);},0,2*E23,0);
  double deriv = delta.Derivative(E23,0,0.000001);

  //cout << endl << ",  d(" << L << ") = " << delta.Eval(E23) << ",  deriv = " << deriv ;

  double hbarc = 197326.9788;

  const Pair& p1 = scheme->GetPrimaryChannels().at(il)->GetPair();
  const Pair& p2 = scheme->GetSecondaryChannels().at(ilp)->GetPair();
  double r1 = scheme->GetPrimaryChannels().at(il)->Radius();
  double m1 = p1.RedMass();
  double r2 = scheme->GetSecondaryChannels().at(ilp)->Radius();
  double m2 = p2.RedMass();

  double rtilde = r1 + hbarc * sqrt(2.*E/m1) * deriv + r2 * sqrt(E*m2/(E23*m1));
  rtilde = rtilde > r1 ? rtilde : r1;  //Check that rtilde is never smaller than r1.
  //rtilde = rtilde > 40. ? rtilde : 40.;

  //cout << ",  rtilde = " << rtilde << endl;

  return rtilde;
}

double GrandWeight::Calculate(vector<TLorentzVector> &p)
{
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);
  Pair pPair(Be8,He4);
  Pair sPair(He4,He4);
  double alphaMass = He4.M();

  double Q = 0; //Q-value of the alpha breakup.
  for(auto & alpha : p){ Q += (alpha.Energy() - alphaMass);} 

  //We determine the maximum Ja we need to include in the sum.
  int JaMax = -1;
  for(PrimaryLevel pi : scheme->GetPrimaryLevels()) JaMax = pi.J > JaMax ? pi.J : JaMax;

  //We construct a matrix to hold the complex amplitudes from the giant summation below.
  Mat<complex<double>> amplitudes(JaMax+1,2*JaMax+1);
  complex<double> zero(0.,0.);
  complex<double> i(0.,1.);
  amplitudes.fill(zero);

  //Same for Jb
  int JbMax = -1;
  for(SecondaryLevel si : scheme->GetSecondaryLevels()) JbMax = si.J > JbMax ? si.J : JbMax; 

  //Then we start the proper calculation.
  for(int Ja = 0; Ja <= JaMax; Ja++){
    Mat<complex<double>> A1 = ConstructPrimaryMatrix(Ja, Q);  //Primary level matrix.
    if(A1.size() == 0) continue; //Skip if A1 is empty (for instance if no levels of spin Ja exist).

    //The indices of the levels with the proper spin.
    vector<int> pIndices = GetPrimaryIndices(Ja);
      
    //Symmetrisation in the alpha particle labels.
    for(int j=0; j<3; j++){
      TLorentzVector alpha1 = p.at(j);
      TLorentzVector alpha2 = p.at((j+1)%3);
      TLorentzVector alpha3 = p.at((j+2)%3);

       //We find angles and energies.
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

      //For the Coulomb correction we need the relative energies of alpha 1-2 and alpha 1-3.
      TLorentzVector cm12 = alpha1 + alpha2;
      cm12.Boost(-cm12.BoostVector());
      double E12 = cm12.Energy() - 2 * alphaMass;
      TLorentzVector cm13 = alpha1 + alpha3;
      cm13.Boost(-cm13.BoostVector());
      double E13 = cm13.Energy() - 2 * alphaMass;

      //double Gamma1 = primChannels.at(k)->PartialWidth(E);
      for(int Jb=0; Jb <= JbMax; Jb++){
        Mat<complex<double>> A2 = ConstructSecondaryMatrix(Jb,E23);
        if(A2.size() == 0) continue; //Skip if matrix is empty.

        vector<int> sIndices = GetSecondaryIndices(Jb);
        
        //Sum over l1
        for(int il=0; il<scheme->GetPrimaryChannels().size(); il++){
          int L1 = scheme->GetPrimaryChannels().at(il)->L();
          double P1 = scheme->GetPrimaryChannels().at(il)->Penetrability(E);
          double rho1 = scheme->GetPrimaryChannels().at(il)->Rho(E);
          double phi1 = scheme->GetPrimaryChannels().at(il)->HardSphereShift(E);
          double omega1 = scheme->GetPrimaryChannels().at(il)->CoulombShift(E);
          complex<double> Omega1 = exp(i * (omega1 - phi1));

          //Sum over L2
          for(int ilp=0; ilp<scheme->GetSecondaryChannels().size(); ilp++){
            double FSCI = CoulombCorrection(E,E23,E12,E13,il,ilp);
            int L2 = scheme->GetSecondaryChannels().at(ilp)->L();
            double P23 = scheme->GetSecondaryChannels().at(ilp)->Penetrability(E23);
            double rho23 = scheme->GetSecondaryChannels().at(ilp)->Rho(E23);
            double phi23 = scheme->GetSecondaryChannels().at(ilp)->HardSphereShift(E23);
            double omega23 = scheme->GetSecondaryChannels().at(ilp)->CoulombShift(E23);
            complex<double> Omega23 = exp(i * (omega23 - phi23));

            //We use the index vectors to only loop over levels with the proper J.
            for(int m=0; m<pIndices.size(); m++){    //Level index with respect to other levels of the same spin.
              int lambda = pIndices.at(m);           //Global level index
              for(int n=0; n<pIndices.size(); n++){
                int mu = pIndices.at(n);        
                double Gmu = scheme->GetPrimaryLevels().at(mu).Bg;
                if(IsAlmostZero(Gmu)) continue;
                for(int mp=0; mp<sIndices.size(); mp++){
                  int lambdap = sIndices.at(mp);
                  double gammaLambdap = scheme->GetSecondaryLevels().at(lambdap).widths.at(ilp);
                  if(IsAlmostZero(gammaLambdap)) continue;
                  for(int np=0; np<sIndices.size(); np++){  
                    int mup = sIndices.at(np);
                    double gammaLambda = scheme->GetPrimaryLevels().at(lambda).GetWidth(il,mup);
                    if(IsAlmostZero(gammaLambda)) continue;
                    /*
                    cout << "Ja = " << Ja ;
                    cout << ",  j = " << j;   
                    cout << ",  Jb = " << Jb ;
                    cout << ",  L1 = " << L1 ;
                    cout << ",  L2 = " << L2 ;
                    cout << ",  lambda = " << lambda ;
                    cout << ",  mu = " << mu ;
                    cout << ",  lambdap = " << lambdap ;
                    cout << ",  mup = " << mup ;     
                    cout << ",  Gmu = " << Gmu << ",  gammaLambda = " << gammaLambda << ",  gammaLambdap = " << gammaLambdap ;
                    */
                    for(int ma=-Ja; ma<=Ja; ma++){
                      //cout << ",  ma = " << ma ;
                      for(int mb = -Jb; mb <= Jb; mb++){
                        //cout << FSCI << endl;
                        //cout << ",  mb = " << mb ;
                        complex<double> Cmmj(ClebschGordan(L1,Jb,Ja,ma-mb,mb,ma),0);
                        complex<double> Ylm1(spherical_harmonic(L1,ma-mb,theta1,phi1));
                        Ylm1 *= pow(i,L1);  //Supposedly good for time-reversal symmetry.
                        complex<double> Ylm2(spherical_harmonic(L2,mb,theta2,phi2));
                        Ylm2 *= pow(i,L2);

                        complex<double> f = Cmmj * Ylm1 * Ylm2 * Omega1 * Omega23  * A1(m,n) * gammaLambda*sqrt(2.*P1/rho1*FSCI)
                                            * Gmu * A2(mp,np) * gammaLambdap*sqrt(2.*P23/rho23);
                        //cout << ",  accessing element (" << Ja << "," << ma+JaMax << ")." << endl;
                        amplitudes(Ja,ma+JaMax) += f;
                        //cout << f << endl;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  double weight = sum(sum(square(abs(amplitudes))));

  //We calculate beta-decay phase space.
  double Ex = Q + tripleThreshold;
  //cout << Ex << endl;
  if(Ex > 16305) return 0.;
  double fBeta = calculatePhaseSpace(7,6,12,Ex);

  return fBeta * weight;

  /*
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
  */

  //double weight = norm(amplitude0) + norm(amplitude1) + norm(amplitude2);
  //cout << "norm(0) = " << norm(amplitude0) << ",  norm(2) = " << norm(amplitude2) << endl;
  //return weight;
}
