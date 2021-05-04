// ROOT stuff
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TCanvas.h>

// other stuff
#include <boost/format.hpp>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

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
    dest += "kinematic.pdf";

    // prepare the dataframe
    ROOT::RDF::RNode df = RDataFrame(chain);
    filter(&df);
    df = df.Define("E_sum","E_cm[0]+E_cm[1]+E_cm[2]");

    // set the axes
    double x_axis[] = {200, 25, 170};
    double y_axis[] = {200, 0, 6000};

    // define the necessary variables
    df = df.Define("E2","min(E_cm[0],min(E_cm[1],E_cm[2]))")
           .Define("E1","max(min(E_cm[0],E_cm[1]), min(max(E_cm[0],E_cm[1]),E_cm[2]))")\
           .Define("E0","max(E_cm[0],max(E_cm[1],E_cm[2]))")
           .Define("X","((E0+2*E1)/E_sum-1)/sqrt(3)")
           .Define("Y","E0/E_sum-1.0/3.0")
           .Filter("pow(X,2)+pow(Y,2)<1.0/9.0")
           .Filter("Y<0.31");

    //*** PLOT ***//
    setup_style();

    TCanvas* canvas = new TCanvas("kinematic", "kinematic", 400, 400);
    TH2D* hist = new TH2D("h1", "kinematic", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    for (int i = 0; i < 3; i++) {
        TH2D htemp = df.Define("x", (format("theta_lab[%1%]") % i).str())
                       .Define("y", (format("E_lab[%1%]") % i).str())
                       .Histo2D({"h1", "kinematic", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x", "y").GetValue();
        hist->Add(&htemp);
    }

    hist->GetXaxis()->SetTitle("$\\theta_{lab}$");
    hist->GetYaxis()->SetTitle("$E_{lab}$");
    hist->Draw("colz");

    canvas->SetRightMargin(0.15);
    canvas->SaveAs(dest.c_str());
}