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
    dest += "QpT.pdf";

    // prepare the dataframe
    ROOT::RDF::RNode df = RDataFrame(chain);
    df = df.Define("E_sum","E_cm[0]+E_cm[1]+E_cm[2]");

    // set the axes    
    double x_axis[] = {400, 0, 11e3};
    double y_axis[] = {400, 0, 250e3};

    //*** PLOT ***//
    setup_style();

    TCanvas* canvas = new TCanvas("QpT", "QpT", 600, 600);
    TH2D hist = df.Histo2D({"h2", "QpT", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "E_sum", "p_tot").GetValue();
    hist.GetXaxis()->SetTitle("exC12");
    hist.GetYaxis()->SetTitle("p_tot");
    hist.Draw("colz");

    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->SaveAs(dest.c_str());
}