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
    dest += "dalitz_raw.pdf";

    // prepare the dataframe
    ROOT::RDF::RNode df = RDataFrame(chain);
    filter(&df);

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
    df = df.Define("E_tot","E_cm[0] + E_cm[1] + E_cm[2]")
           .Define("e_cm_1", "E_cm[0]/E_tot") // normalized such that e1 + e2 + e3 = 1
           .Define("e_cm_2", "E_cm[1]/E_tot")
           .Define("e_cm_3", "E_cm[2]/E_tot")
           .Define("x","sqrt(3)*(e_cm_2 - e_cm_3)")
           .Define("y","3*e_cm_1 - 1");

    //*** PLOT ***//
    setup_style();

    TCanvas* canvas = new TCanvas("c", "c", 600, 600);
    TH2D hist = df.Histo2D({"h2", "Dalitz raw", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x", "y").GetValue();
    hist.GetXaxis()->SetTitle("x");
    hist.GetYaxis()->SetTitle("y");
    hist.Draw("colz");

    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->SaveAs(dest.c_str());
}