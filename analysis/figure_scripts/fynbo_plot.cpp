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
    dest += "fynbo_plot.pdf";

    // prepare the dataframe
    ROOT::RDF::RNode df = RDataFrame(chain);
    df = df.Define("E_sum","E_cm[0]+E_cm[1]+E_cm[2]");

    // set the axes    
    double x_axis[] = {400, 0, 12000};
    double y_axis[] = {400, 0, 6500};

    df = df.Filter("p_tot<50e3");

    //*** PLOT ***//
    setup_style();

    TCanvas* canvas = new TCanvas("fynbo_plot", "fynbo_plot", 400, 400);
    TH2D* hist = new TH2D("h1", "fynbo_plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    for (int i = 0; i < 3; i++) {
        TH2D htemp = df.Define("y", (format("E_cm[%1%]") % i).str())
                       .Histo2D({"h1", "fynbo_plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "E_sum", "y").GetValue();
        hist->Add(&htemp);
    }
    
    hist->GetXaxis()->SetTitle("E_{tot} [keV]");
    hist->GetYaxis()->SetTitle("E_\\alpha [keV]");
    hist->Draw("colz");

    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->SaveAs(dest.c_str());
}