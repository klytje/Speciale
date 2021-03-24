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
    dest += "eCM.pdf";

    // prepare the dataframe
    ROOT::RDF::RNode df = RDataFrame(chain);
    df = df.Define("esum","E_cm[0]+E_cm[1]+E_cm[2]");

    // set the axes    
    double x_axis[] = {300, 0, 6000};

    // define the necessary variables
    df = df.Filter("p_tot<35e3")
           .Filter("abs(deltaE)<200")
           .Define("E2","min(E_cm[0],min(E_cm[1],E_cm[2]))")
           .Define("E1","max(min(E_cm[0],E_cm[1]), min(max(E_cm[0],E_cm[1]),E_cm[2]))")\
           .Define("E0","max(E_cm[0],max(E_cm[1],E_cm[2]))")
           .Define("X","((E1-E2)/esum)*sqrt(3)")
           .Define("Y","3*E0/esum-1.0")
           .Filter("pow(X,2)+pow(Y,2)<1.0")
           .Filter("Y<0.93")
           .Filter("E2>250");

    //*** PLOT ***//
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gROOT->ForceStyle();

    TCanvas* canvas = new TCanvas("eCM", "eCM", 600, 600);
    TH1D* hist = new TH1D("h1", "eCM", int(x_axis[0]), x_axis[1], x_axis[2]);
    for (int i = 0; i < 3; i++) {
        TH1D htemp = df.Histo1D({"h1", "eCM", int(x_axis[0]), x_axis[1], x_axis[2]}, (format("E_cm[%1%]") % i).str()).GetValue();
        hist->Add(&htemp);
    }

    hist->GetXaxis()->SetTitle("$\\rho$");
    hist->GetYaxis()->SetTitle("Count");
    hist->Draw("colz");

    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->SaveAs(dest.c_str());
}