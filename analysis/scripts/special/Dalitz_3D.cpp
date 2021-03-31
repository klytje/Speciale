// ROOT stuff
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>

// other stuff
#include <filesystem>
#include <boost/format.hpp>
#include <iostream>

// debug stuff
#include <TApplication.h>

// my stuff
#include "../plot_style.cpp"

using namespace std;
using namespace ROOT;
using boost::format;

int main(int argc, char *argv[]) {
    ROOT::RDF::RNode df = RDataFrame("tree", argv[1]);
    df = df.Define("E_tot","E_cm[0]+E_cm[1]+E_cm[2]");

    // set the axes    
    double x_axis[] = {100, -2, 2};
    double y_axis[] = {100, -2, 2};

    // define sorting methods
    auto max = [] (double e1, double e2, double e3) {return std::max({e1, e2, e3});};
    auto min = [] (double e1, double e2, double e3) {return std::min({e1, e2, e3});};
    auto mid = [] (double e1, double e2, double e3) {
        if (e1 > e2) {
            if (e2 > e3) {
                return e2;
            } else if (e1 > e3) {
                return e3;
            } else {
                return e1;
            }
        } else {
            if (e1 > e3) {
                return e1;
            } else if (e2 > e3) {
                return e3;
            } else {
                return e2;
            }
        }
    };
    df = df.Filter("abs(deltaE)<200")
           .Define("e_cm_1", "E_cm[0]/E_tot") // normalized such that e1 + e2 + e3 = 1
           .Define("e_cm_2", "E_cm[1]/E_tot")
           .Define("e_cm_3", "E_cm[2]/E_tot")
           .Define("e_1", max, {"e_cm_1", "e_cm_2", "e_cm_3"}) // we want e1 > e2 > e3
           .Define("e_2", mid, {"e_cm_1", "e_cm_2", "e_cm_3"})
           .Define("e_3", min, {"e_cm_1", "e_cm_2", "e_cm_3"})
        //    .Define("X","(e_cm_2 - e_cm_3 + 0.5)*sqrt(3)*2")
        //    .Define("Y","(2*e_cm_1 - e_cm_2 - e_cm_3 + 0.5)*2");
           .Define("X","sqrt(3)*(e_2 - e_3)")
           .Define("Y","3*e_1 - 1");

    TApplication *app = new TApplication("ROOT window", 0, 0);
    TCanvas* canvas = new TCanvas("3D Dalitz", "3D Dalitz", 600, 600);
    TH3D hist = df.Histo3D({"h3", "3D Dalitz", 50, 0, 1, 50, 0, 1, 50, 0, 1}, "e_cm_1", "e_cm_2", "e_cm_3").GetValue();

    hist.GetXaxis()->SetTitle("E1");
    hist.GetXaxis()->CenterTitle();
    hist.GetXaxis()->SetTitleOffset(2);

    hist.GetYaxis()->SetTitle("E2");
    hist.GetYaxis()->CenterTitle();
    hist.GetYaxis()->SetTitleOffset(2);

    hist.GetZaxis()->SetTitle("E3");
    hist.GetZaxis()->CenterTitle();
    hist.GetZaxis()->SetTitleOffset(2);
    
    hist.Draw("ISO");

    canvas->Show();
    canvas->SetPhi(225);
    canvas->SetLeftMargin(0.15);
    app->Run();

    df.Foreach([] (double e1, double e2, double e3) {cout << format("(e1, e2, e3): (%1%, %2%, %3%)") % e1 % e2 % e3 << endl; usleep(10);}, {"e_1", "e_2", "e_3"});
}