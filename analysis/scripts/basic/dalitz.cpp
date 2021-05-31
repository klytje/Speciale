// ROOT stuff
#include <TFile.h>
#include <TTree.h>
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
#include "../utility.cpp"

using namespace std;
using namespace ROOT;
using boost::format;

int main(int argc, char *argv[]) {
    //*** RAW ***//
    setup_style();

    TCanvas* c1 = new TCanvas("c1", "c", 600, 600);
    TH2D* h1 = dalitz(argv[1], 200, false, false);
    setup_dalitz_plot(h1);
    h1->Draw("colz");

    string path = string(argv[2]) + "dalitz_no_cuts.pdf"; 
    c1->SetLogz();
    c1->SetRightMargin(0.15);
    c1->SaveAs(path.c_str());

    //*** w/ CUTS ***//
    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", argv[1]);
    filter(&df);
    setup_dataframe(&df);
    cut_circle(&df);
    if (string(argv[1]).find("events") != string::npos) {
        cut_gs(&df);
    }
    TCanvas* c2 = new TCanvas("c2", "c", 600, 600); 
    TH2D* h2 = dalitz(&df, 200, true, false);
    setup_dalitz_plot(h2);
    h2->Draw("colz");

    path = string(argv[2]) + "dalitz.pdf"; 
    c2->SetLogz();
    c2->SetRightMargin(0.15);
    c2->SaveAs(path.c_str());
}