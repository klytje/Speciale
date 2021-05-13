#include <complex>
#include <iostream>
#include <gsl/gsl_sf_coupling.h>
#include <TF1.h>
#include <TMath.h>
#include <ausa/util/memory>
#include "jonas/BalamuthWeight.h"
#include "jonas/Nucleus.h"
#include "jonas/Pair.h"
#include "jonas/Constants.h"

using namespace std;

BalamuthWeight::BalamuthWeight()
{
    DoCorrection(false);
}

BalamuthWeight::~BalamuthWeight() {}

void BalamuthWeight::SetPrimary(Level& l, unique_ptr<Channel> c)
{
    primLevel = l;
    primChannel = std::move(c);
    spinCouple();
    primSH = std::make_unique<SphericalHarmonic>(primChannel->L());
    primSH->interpolate(true, 100);
}

void BalamuthWeight::SetSecondary(Level& l, unique_ptr<Channel> c)
{
    secLevel = l;
    secChannel = std::move(c);
    spinCouple();
    secSH = std::make_unique<SphericalHarmonic>(secChannel->L());
    secSH->interpolate(true, 100);
}

vector<complex<double>> BalamuthWeight::CalculateAmplitudes(vector<TLorentzVector> &p)
{
    double Q = 0;
    for(auto & alpha : p){ Q += (alpha.Energy() - ALPHA_MASS);}// * 1e6;}  //Q-value of the alpha breakup.
    double Eb = (secLevel.E() + BE8_MASS) - 2 * ALPHA_MASS;                 //Energy of intermediate state above threshold.
	// std::cout << secLevel.E() << "\n";


    vector<complex<double>> amplitudes;
    amplitudes.resize(static_cast<size_t>(2*primLevel.J()+1));

    for(int j=0; j<3; j++){
        TLorentzVector alpha1 = p.at(j);
        TLorentzVector alpha2 = p.at((j+1)%3);
        TLorentzVector alpha3 = p.at((j+2)%3);

        double E1 = (alpha1.Energy() - ALPHA_MASS);// * 1e6;
        TLorentzVector rcm = alpha2 + alpha3;     //Recoil center of mass system.

        double theta1 = alpha1.Theta();
        double phi1 = alpha1.Phi();
        alpha2.Boost(-rcm.BoostVector());
        double theta2 = alpha2.Theta();
        double phi2 = alpha2.Phi();
        alpha2.Boost(rcm.BoostVector());

        rcm.Boost(-rcm.BoostVector());
        double E23 = (rcm.Energy() - 2 * ALPHA_MASS);// * 1e6;  //Relative energy of alpha2 and alpha3.
        double E = Q - E23;                                 //Relative energy of alpha1 and recoil.

        double Gamma1 = primChannel->PartialWidth(E);
        if(doCorrection){
            //For the Coulomb correction we need the relative energies of alpha 1-2 and alpha 1-3.
            TLorentzVector cm12 = alpha1 + alpha2;
            cm12.Boost(-cm12.BoostVector());
            double E12 = cm12.Energy() - 2 * ALPHA_MASS;
            TLorentzVector cm13 = alpha1 + alpha3;
            cm13.Boost(-cm13.BoostVector());
            double E13 = cm13.Energy() - 2 * ALPHA_MASS;

            //Then we multiply by the correction factors.
            Gamma1 *= c1->Rho(E) / (c2->Rho(E12) * c2->Rho(E13));
            Gamma1 *= c2->Penetrability(E12) * c2->Penetrability(E13) / c1->Penetrability(E);
        }
        double Gamma2 = secChannel->PartialWidth(E23);
        double phase1 = primChannel->CoulombShift(E) - primChannel->HardSphereShift(E);
        double phase2 = secChannel->CoulombShift(E23) - secChannel->HardSphereShift(E23);

        //We start calculating...
        complex<double> i(0,1);
        complex<double> a(0,0);

        a += sqrt(Gamma2 / secChannel->Rho(E23)); //Notice sqrt(E) -> rho(E).
        a /= Eb - pow(secChannel->ReducedWidth(),2) * (secChannel->ShiftFunction(E23) - secChannel->ShiftFunction(Eb)) - E23 - 0.5 * i * Gamma2;
        a *= sqrt(Gamma1 / primChannel->Rho(E)) * exp(i * phase1) * exp(i * phase2);

        for(int ma = (int)-primLevel.J(), ima = 0; ma <= primLevel.J() ;ma++, ima++) {
            complex<double> sum(0, 0);

            for (auto mb = (int)-secLevel.J(); mb <= secLevel.J(); mb++) {
                auto Cmmj = CG->Coefficient(ma - mb, mb, primLevel.J());
                if (Cmmj == 0) continue;

                auto Ylm1 = primSH->eval(ma - mb, theta1, phi1);
                auto Ylm2 = secSH->eval(mb, theta2, phi2);

                sum += Cmmj * Ylm1 * Ylm2;
            }
            amplitudes[ima] += a*sum;
        }
    }
    return amplitudes;
}

double BalamuthWeight::Calculate(vector<TLorentzVector> &p)
{
    double prob = 0.;
    for(auto a : CalculateAmplitudes(p)) {
        prob += norm(a);
    }
    prob *= scale;   //Apparent maximum probability;
    return prob;
}


void BalamuthWeight::DoCorrection(bool input, double radius)
{
    doCorrection = input;
    if(!doCorrection) return;

    //Cannot do correction if the channels are not defined.
    if(!primChannel || !secChannel){
        cout << "  BalamuthWeight::DoCorrection(): Channels need to be defined!" << endl;
        exit(EXIT_FAILURE);
    }

    //Correction for Coulomb repulsion.
    Nucleus He4(4,2,0,1);
    Nucleus Be8(8,4,0,1);
    Pair p1(He4,Be8);
    Pair p2(He4,He4);
    c1 = std::make_unique<Channel>(p1,primChannel->L(),1.,radius/p1.Radius()); //Width arbitrary, only radius is important.
    c2 = std::make_unique<Channel>(p2,secLevel.J(),1.,radius/p2.Radius());     //We guess that the L's are the same as the 8Be spin.

    //XXX: Interpolation impossible for varying correction radius!
    c1->UseInterpolation(true,10000.,10.);
    c2->UseInterpolation(true,10000.,10.);
}


double BalamuthWeight::GetClebschGordan(double j1, double j2, double j3, double m1, double m2, double m3)
{
    m3=-m3;
    auto j1x2=(int)lrint(2*j1);
    auto j2x2=(int)lrint(2*j2);
    auto j3x2=(int)lrint(2*j3);
    auto m1x2=(int)lrint(2*m1);
    auto m2x2=(int)lrint(2*m2);
    auto m3x2=(int)lrint(2*m3);

    double w3j=gsl_sf_coupling_3j(j1x2,j2x2,j3x2,m1x2,m2x2,m3x2);

    return pow(-1.0,lrint(j1-j2-m3)) * sqrt(2.0*j3+1.) * w3j;
}

void BalamuthWeight::spinCouple() {
    if (!primChannel || !secChannel) return;

    CG = std::make_unique<ClebschGordan>(primChannel->L(), secLevel.J());
}

/*
double BalamuthWeight::CorrectionRadius(unique_ptr<Channel> & c, Level & l, double E, double E23)
{
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);

  auto phase = [&] (double e){
    double Eb = (l.E() + Be8.M()) - 2 * ALPHA_MASS; //Energy of intermediate state above threshold.
    double phi = c->HardSphereShift(e);
    double Gamma = c->PartialWidth(e);
    double S = c->ShiftFunction(e);
    double B = c->ShiftFunction(Eb);
    double gamma2 = pow(c->ReducedWidth(),2);

    double arg = (Gamma/2.) / (Eb - e - (S - B) * gamma2);
    double y = abs(arg);
    double x = TMath::Sign(1.,arg);

    //cout << "El = " << Eb << ",  Gamma = " << Gamma << ",  S = " << S << ",  B = " << B << ", gamma2 = " << gamma2 << ", phi = " << phi;

    double del = atan2(y,x) - phi;
    return del;   
  };

  const Pair& p1 = primChannel->GetPair();
  const Pair& p2 = c->GetPair();
  double r1 = primChannel->Radius();
  double m1 = p1.RedMass();
  double r2 = c->Radius();
  double m2 = p2.RedMass();

  TF1 delta("delta",[&](double* x, double* p) { return phase(*x);},0,2*E23,0);
  double deriv = delta.Derivative(E23,0,0.000001);
 
  //double deriv = phase(E23);

  double rtilde = r1 + hbarc * sqrt(2.*E/m1) * deriv + r2 * sqrt(E*m2/(E23*m1));
  rtilde = rtilde > r1 ? rtilde : r1;  //Check that rtilde is never smaller than r1.

  //cout << endl << ",  d(" << c->L() << ") = " << delta.Eval(E23) << ",  deriv = " << deriv ;
  //cout << ",  rtilde = " << rtilde << ",  r1 = " << r1 << ",  r2 = " << r2 << endl;

  return rtilde;
}


void BalamuthWeight::Hello()
{
  cout << "Hello!" << endl;
}
*/

/*
double BalamuthWeight::SecondaryPS(double E23)
{
  Nucleus He4(4,2,0,1);
  Nucleus Be8(8,4,0,1);

  double Eb = (secLevel.E() + Be8.M()) - 2 * ALPHA_MASS; //Energy of intermediate state above threshold.
  double phi = secChannel->HardSphereShift(E23);
  double Gamma = secChannel->PartialWidth(E23);
  double S = secChannel->ShiftFunction(E23);
  double B = secChannel->ShiftFunction(Eb);
  double gamma2 = pow(secChannel->ReducedWidth(),2);

  double arg = (Gamma/2.) / (Eb - E23 - (S - B) * gamma2);
  double y = abs(arg);
  double x = TMath::Sign(1.,arg);

  double delta = atan2(y,x) - phi;
  return delta;
}

//SecondaryPS modified to satisfy ROOT::TF1.
double BalamuthWeight::ROOTPS(double *var, double *par)
{
  double E23 = var[0];
  return SecondaryPS(E23);
}

//The arguments are the relative energies of the first and secondary breakup
double BalamuthWeight::CorrectionRadius(double E, double E23)
{
  const Pair& p1 = primChannel->GetPair();
  const Pair& p2 = secChannel->GetPair();
  double r1 = primChannel->Radius();
  double m1 = p1.RedMass();
  double r2 = secChannel->Radius();
  double m2 = p2.RedMass();

  TF1 delta("delta",[&](double* x, double* p) { return SecondaryPS(*x);},0,2*E23,0);
  double deriv = delta.Derivative(E23,0,0.000001);

//  cout << "E = " << E << ",  E23 = " << E23 << ",  r1 = " << r1 << ",  rd = " << hbarc * sqrt(2.*E/m1) * deriv << ",  rt = " << 2 * sqrt(E*m2/(E23*m1)) << endl;

  double rtilde = r1 + hbarc * sqrt(2.*E/m1) * deriv + r2 * sqrt(E*m2/(E23*m1));

  return rtilde;
}
*/

