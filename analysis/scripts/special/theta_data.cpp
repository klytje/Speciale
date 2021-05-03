// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>
#include <TVector3.h>

// other stuff
#include <boost/format.hpp>
#include <iostream>
#include <math.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cout << "Usage: ./theta_proj <output path> <input>" << endl;
        exit(1);
    }

    setup_style();
    auto theta = [] (double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double r3x, double r3y, double r3z) {
        TVector3 p23 = {r2x-r3x, r2y-r3y, r2z-r3z};
        TVector3 p1 = {r1x, r1y, r1z};
        return cos(p1.Angle(p23));

    };

    // auto max = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
    //     TVector3 p1 = {x[0], y[0], z[0]};
    //     TVector3 p2 = {x[1], y[1], z[1]};
    //     TVector3 p3 = {x[2], y[2], z[2]};
    //     vector<double> mags = {p1.Mag(), p2.Mag(), p3.Mag()};
    //     sort(mags.begin(), mags.end());
    //     if (p1.Mag() > mags[1]) {
    //         return 0;
    //     } else if (p2.Mag() > mags[1]) {
    //         return 1;
    //     } else {
    //         return 2;
    //     }
    // };

    // auto mid = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
    //     TVector3 p1 = {x[0], y[0], z[0]};
    //     TVector3 p2 = {x[1], y[1], z[1]};
    //     TVector3 p3 = {x[2], y[2], z[2]};
    //     vector<double> mags = {p1.Mag(), p2.Mag(), p3.Mag()};
    //     sort(mags.begin(), mags.end());
    //     if (p1.Mag() > mags[1]) {
    //         if (mags[1] > p2.Mag()) {
    //             return 2;
    //         } else {
    //             return 1;
    //         }
    //     } else if (p2.Mag() > mags[1]) {
    //         if (mags[1] > p3.Mag()) {
    //             return 0;
    //         } else {
    //             return 2;
    //         }
    //     } else {
    //         if (mags[1] > p1.Mag()) {
    //             return 1;
    //         } else {
    //             return 0;
    //         }
    //     }
    // };

    auto max = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
        TVector3 p1 = {x[0], y[0], z[0]};
        TVector3 p2 = {x[1], y[1], z[1]};
        TVector3 p3 = {x[2], y[2], z[2]};
        vector<double> m = {p1.Mag(), p2.Mag(), p3.Mag()};
        double maxv = std::max({m[0], m[1], m[2]});
        if (maxv == m[0]) {
            return 0;
        } else if (maxv == m[1]) {
            return 1;
        } else {
            return 2;
        }
    };

    auto mid = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
        TVector3 p1 = {x[0], y[0], z[0]};
        TVector3 p2 = {x[1], y[1], z[1]};
        TVector3 p3 = {x[2], y[2], z[2]};
        vector<double> m = {p1.Mag(), p2.Mag(), p3.Mag()};
        if (m[0] > m[1]) {
            if (m[1] > m[2]) { // m0 > m1 > m2
                return 1;
            } else if (m[0] > m[2]) { // m0 > m2 > m1
                return 2;
            } else { // m2 > m0 > m1
                return 0;
            }
        } else if (m[1] > m[2]) { // m1 > m0, m2
            if (m[0] > m[2]) { // m1 > m0 > m2
                return 0;
            } else { // m1 > m2 > m0
                return 2;
            }
        } else { // m2 > m1 > m0
            return 1;
        }
    };

    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", argv[2]);
    df = df.Define("i_max", max, {"px", "py", "pz"})
            .Define("i_mid", mid, {"px", "py", "pz"})
            .Define("i_min", "3 - i_max - i_mid") // i_min + i_mid + i_max = 3
            .Define("px1", "px[i_max]")
            .Define("py1", "py[i_max]")
            .Define("pz1", "pz[i_max]")
            .Define("px2", "px[i_mid]")
            .Define("py2", "py[i_mid]")
            .Define("pz2", "pz[i_mid]")
            .Define("px3", "px[i_min]")
            .Define("py3", "py[i_min]")
            .Define("pz3", "pz[i_min]");

    TCanvas* c = new TCanvas("c", "c", 600, 600);
    TH1D h = df.Define("theta", theta, {"px1", "py1", "pz1", "px2", "py2", "pz2", "px3", "py3", "pz3"}).Histo1D({"h", "h", 100, 0, 1}, "theta").GetValue();
    h.SetLineColor(kBlack);
    h.SetLineWidth(2);
    h.Scale(1./h.GetMaximum());
    h.Draw("HIST L");

    // TH1D h1 = df.Define("x", "cos(phi_cm/180)").Histo1D({"h1", "h", 100, 0, 1}, "x").GetValue();
    // TH1D h2 = df.Define("x", "cos(theta_cm/180)").Histo1D({"h2", "h", 100, 0, 1}, "x").GetValue();
    // h1.SetLineColor(kOrange+1);
    // h1.SetLineWidth(2);
    // h1.Scale(1./h1.GetMaximum());
    // h2.SetLineColor(kAzure+1);
    // h2.SetLineWidth(2);
    // h2.Scale(1./h2.GetMaximum());
    // h1.Draw("HIST L SAME");
    // h2.Draw("HIST L SAME");

    auto func = [] (double* x, double* par) {
        double scale = par[0];
        double c1 = 0.368245, c3 = 0.630681;
        double beta = M_PI - x[0];
        return scale*(c1*correlation_functions.at("2- 1")(beta) + c3*correlation_functions.at("2- 3")(beta));
    };
    TF1* prediction = new TF1("mix", func, 0, 3.14, 1);
    prediction->SetParameter(0, h.GetMaximum());

    // prediction->Draw("same");
    string path = string(argv[1]) + "breakup_angle_data.pdf";
    c->SaveAs(path.c_str());
}