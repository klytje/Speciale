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

// my stuff
#include "../plot_style.cpp"

using namespace std;
using namespace ROOT;
using boost::format;

int main(int argc, char *argv[]) {
    // all but the last argument are input files
    TChain chain("tree");
    for (int i = 1; i < argc-1; i++) {
        chain.Add(argv[i]);
    }    
    // final argument is the destination
    string dest = argv[argc-1]; 
    dest += "dalitz.pdf";

    // prepare the dataframe
    ROOT::RDF::RNode df = RDataFrame(chain);

    // set the axes    
    double x_axis[] = {200, -1.3, 1.3};
    double y_axis[] = {200, -1.3, 1.3};

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
    // define the necessary variables
    df = df.Define("E_tot","E_cm[0] + E_cm[1] + E_cm[2]")
           .Define("e_cm_1", "E_cm[0]/E_tot") // normalized such that e1 + e2 + e3 = 1
           .Define("e_cm_2", "E_cm[1]/E_tot")
           .Define("e_cm_3", "E_cm[2]/E_tot")
           .Define("e_1", max, {"e_cm_1", "e_cm_2", "e_cm_3"}) // we want e1 > e2 > e3
           .Define("e_2", mid, {"e_cm_1", "e_cm_2", "e_cm_3"})
           .Define("e_3", min, {"e_cm_1", "e_cm_2", "e_cm_3"})
           .Define("x","sqrt(3)*(e_2 - e_3)")
           .Define("y","3*e_1 - 1")
           .Filter("pow(x,2) + pow(y,2) < 1.0")
           .Filter("y < 0.93");
    
    //*** PLOT ***//
    setup_style();

    TCanvas* canvas = new TCanvas("c", "c", 600, 600);
    TH2D* hist = new TH2D("h1", "Dalitz plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);

    // we can get the other slices simply by permutating i, j, k
    int perms[] = {1, 2, 3};
    do {
        int i = perms[0];
        int j = perms[1];
        int k = perms[2];
        TH2D htemp = df.Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                       .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                       .Histo2D({"h1", "temp", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp").GetValue();
        hist->Add(&htemp);
    } while (std::next_permutation(perms, perms+3)); // repeat for each of the 3! = 6 permutations of {1, 2, 3}
    
    hist->GetXaxis()->SetTitle("x");
    hist->GetYaxis()->SetTitle("y");
    hist->Draw("colz");

    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->SaveAs(dest.c_str());
}