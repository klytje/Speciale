// ROOT
#include <TROOT.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TF1.h>
#include <TAxis.h>
#include <TH1.h>

// cwfcomplex setup
#include <complex>

using namespace std;
const double precision = 1E-10,sqrt_precision = 1E-5;
#define SIGN(a) (((a) < 0) ? (-1) : (1))

#include "cwfcomplex/complex_functions.H"
#include "cwfcomplex/cwfcomp.cpp"

// other stuff
#include <boost/format.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
using boost::format;

// constants
//const auto hbar = boost::units::si::constants::codata::hbar.value();
//const auto c = boost::units::si::constants::codata::c.value();

const double hbar = 197.329;
const double c = 1;
const double fm = 1;
const double MeV = 1;
const double m_p = 931.494*MeV/pow(c, 2);
const double alpha = 0.007297;

// setup
const double Z1 = 2, Z2 = 2;
const double A1 = 2, A2 = 2;
const double E_res = 3.03*MeV;
const int l = 2;

const double mu = A1*A2/(A1+A2)*m_p;
const double a = 1.45*fm*(pow(A1, 1./3) + pow(A2, 1./3));

double sommerfeld(double rho) {
	return alpha*Z1*Z2*c*mu*a/(rho*hbar);
}

double penetrability(int l, double rho) {
    double eta = sommerfeld(rho);

    class Coulomb_wave_functions cwf(true, l, eta);
    complex<double> F, dF;
    cwf.F_dF(rho, F, dF);
    
    complex<double> G, dG;
    cwf.G_dG(rho, G, dG);
    return rho/(norm(F) + norm(G));
}

double shift_factor(int l, double rho) {
    double eta = sommerfeld(rho);

    class Coulomb_wave_functions cwf(true, l, eta);
    complex<double> F, dF;
    cwf.F_dF(rho, F, dF);

    complex<double> G, dG;
    cwf.G_dG(rho, G, dG);
    return real(rho*(F*dF + G*dG)/(norm(F) + norm(G)));
}

double penetrability(double* x, double* par) {
    double rho = x[0];
    int l = par[0];
    double eta = sommerfeld(rho);

    class Coulomb_wave_functions cwf(true, l, eta);
    complex<double> F, dF;
    cwf.F_dF(rho, F, dF);
    
    complex<double> G, dG;
    cwf.G_dG(rho, G, dG);
    return rho/(pow((double)F.real(),2) + pow((double)G.real(),2));
}

double shift_factor(double* x, double* par) {
    double rho = x[0];
    int l = par[0];
    double eta = sommerfeld(rho);

    class Coulomb_wave_functions cwf(true, l, eta);
    complex<double> F, dF;
    cwf.F_dF(rho, F, dF);

    complex<double> G, dG;
    cwf.G_dG(rho, G, dG);
    return rho*(F.real()*dF.real() + G.real()*dG.real())/(pow((double)F.real(),2) + pow((double)G.real(),2));
}

double k(double E) {
	return sqrt(2*mu*E)/hbar;
}

int main(int argc, char *argv[]) {
    double gamma = atof(argv[1])*sqrt(MeV/1000);

    double rho = a*k(E_res);
    double P = penetrability(l, rho);
    double Sp = 0;
    double Gamma = 2*pow(gamma, 2)*P/(1 + pow(gamma, 2)*Sp*hbar)/MeV;

    cout << "The width is " << Gamma << " MeV" << endl;

    // penetrability & shift plot for various l
    if (true) {
        TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
        TF1* pen = new TF1("pen", penetrability, 0.01, 2, 1);
        pen->SetParameter(0, 0);
        pen->Draw();

        TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
        TF1* shift = new TF1("shift", shift_factor, 0.01, 2, 1);
        shift->SetParameter(0, 0);
        shift->GetHistogram()->SetMinimum(-3.1);
        shift->Draw();

        string name;
        for (int l = 1; l < 4; l++) {
            c1->cd();
            TF1* pen = new TF1((format("pen %1%") % l).str().c_str(), penetrability, 0.01, 2, 1);
            pen->SetParameter(0, l);
            pen->DrawClone("SAME");

            c2->cd();
            TF1* shift = new TF1((format("shift %1%") % l).str().c_str(), shift_factor, 0.01, 2, 1);
            shift->SetParameter(0, l);
            shift->DrawClone("SAME");
        }
        c1->SaveAs("pen.png");
        c2->SaveAs("shift.png");
    }
}