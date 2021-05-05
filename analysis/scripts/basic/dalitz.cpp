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
    //*** PLOT ***//
    setup_style();

    TCanvas* canvas = new TCanvas("c", "c", 600, 600);
    TH2D* hist = dalitz(argv[1]);
    
    hist->GetXaxis()->SetTitle("x");
    hist->GetYaxis()->SetTitle("y");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetNdivisions(2);
    hist->GetYaxis()->SetNdivisions(2);
    hist->Draw("colz");

    string path = string(argv[2]) + "dalitz.pdf"; 
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->SaveAs(path.c_str());
}