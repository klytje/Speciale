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
    dest += "dalitzfig.pdf";

    // prepare the dataframe
    ROOT::RDF::RNode df = RDataFrame(chain);
    df = df.Define("E_sum","E_cm[0]+E_cm[1]+E_cm[2]");

    // set the axes    
    double x_axis[] = {200, -1.3, 1.3};
    double y_axis[] = {200, -1.3, 1.3};

    // define the necessary variables
    df = df.Filter("p_tot<35e3")
        .Filter("abs(deltaE)<200")
        .Define("E2","min(E_cm[0],min(E_cm[1],E_cm[2]))")
        .Define("E1","max(min(E_cm[0],E_cm[1]), min(max(E_cm[0],E_cm[1]),E_cm[2]))")
        .Define("E0","max(E_cm[0],max(E_cm[1],E_cm[2]))")
        .Define("X","((E1-E2)/E_sum)*sqrt(3)")
        .Define("Y","3*E0/E_sum-1.0")
        .Filter("pow(X,2)+pow(Y,2)<1.0")
//        .Filter("Y<0.93")
        .Filter("E2>250");

    //*** PLOT ***//
    setup_style();

    TCanvas* canvas = new TCanvas("Dalitz_fig", "Dalitz_fig", 600, 600);
    TH2D* hist = new TH2D("h1", "kinematic_lab", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);

    int perms[] = {0, 1, 2};
    do {
        int i = perms[0];
        int j = perms[1];
        int k = perms[2];
        TH2D htemp = df.Define("x", (format("((E_cm[%1%] - E_cm[%2%])/E_sum)*sqrt(3)") % j % k).str())
                       .Define("y", (format("3*E_cm[%1%]/E_sum-1") % i).str())
                       .Histo2D({"h1", "kinematic_lab", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x", "y").GetValue();
        hist->Add(&htemp);
    } while (std::next_permutation(perms, perms+3));
    
    hist->GetXaxis()->SetTitle("X");
    hist->GetYaxis()->SetTitle("Y");
    hist->Draw("colz");

    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->SaveAs(dest.c_str());
}