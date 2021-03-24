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
    dest += "phi.pdf";

    // prepare the dataframe
    ROOT::RDF::RNode df = RDataFrame(chain);
    df = df.Define("E_sum","E_cm[0]+E_cm[1]+E_cm[2]");

    // set the axes
    double x_axis[] = {200, 0, M_PI/3};

    // define the necessary variables
    df = df.Filter("p_tot<35e3")
           .Filter("abs(deltaE)<200")
           .Define("E2","min(E_cm[0],min(E_cm[1],E_cm[2]))")
           .Define("E1","max(min(E_cm[0],E_cm[1]), min(max(E_cm[0],E_cm[1]),E_cm[2]))")\
           .Define("E0","max(E_cm[0],max(E_cm[1],E_cm[2]))")
           .Define("X","((E1-E2)/E_sum)*sqrt(3)")
           .Define("Y","3*E0/E_sum-1.0")
           .Filter("pow(X,2)+pow(Y,2)<1.0")
           .Filter("Y<0.93")
           .Filter("E2>250");

    //*** PLOT ***//
    setup_style();

    TCanvas* canvas = new TCanvas("phi", "phi", 400, 400);
    TH1D hist = df.Define("x", "atan2(X,Y)")
                  .Histo1D({"h1", "phi", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();

    hist.GetXaxis()->SetTitle("\\phi");
    hist.GetYaxis()->SetTitle("Count");
    hist.Draw("colz");

    canvas->SetRightMargin(0.15);
    canvas->SaveAs(dest.c_str());
}