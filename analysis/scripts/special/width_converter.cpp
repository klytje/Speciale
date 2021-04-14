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
using boost::format;

double penetrability(double* x, double* par) {
    double rho = x[0];
    int l = par[0];
    double eta = par[1];

    class Coulomb_wave_functions cwf(true, l, eta);
    complex<double> F, dF;
    cwf.F_dF(rho, F, dF);
    
    complex<double> G, dG;
    cwf.G_dG(rho, G, dG);
    return rho/(norm(F) + norm(G));
}

double shift_factor(double* x, double* par) {
    double rho = x[0];
    int l = par[0];
    double eta = par[1];

    class Coulomb_wave_functions cwf(true, l, eta);
    complex<double> F, dF;
    cwf.F_dF(rho, F, dF);

    complex<double> G, dG;
    cwf.G_dG(rho, G, dG);
    return real(rho*(F*dF + G*dG)/(norm(F) + norm(G)));
}

int main(int argc, char *argv[]) {
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    TF1* pen = new TF1("pen", penetrability, 0.01, 2, 2);
    pen->SetParameter(0, 0);
    pen->SetParameter(1, 0);
    pen->Draw();

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    TF1* shift = new TF1("shift", shift_factor, 0.01, 2, 2);
    shift->SetParameter(0, 0);
    shift->SetParameter(1, 0);
    shift->GetHistogram()->SetMinimum(-3.1);
    shift->Draw();

    string name;
    for (int l = 1; l < 4; l++) {
        c1->cd();
        TF1* pen = new TF1((format("pen %1%") % l).str().c_str(), penetrability, 0.01, 2, 2);
        pen->SetParameter(0, l);
        pen->SetParameter(1, 0);
        pen->DrawClone("SAME");

        c2->cd();
        TF1* shift = new TF1((format("shift %1%") % l).str().c_str(), shift_factor, 0.01, 2, 2);
        shift->SetParameter(0, l);
        shift->SetParameter(1, 0);
        shift->DrawClone("SAME");
    }
    c1->SaveAs("pen.png");
    c2->SaveAs("shift.png");
}