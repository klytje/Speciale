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

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using namespace ROOT;
using boost::format;

int main(int argc, char *argv[]) {
    setup_style();

    if (string(argv[2]).find("_i") == string::npos) {
        cout << "\033[1;31m" << "ERROR: Input file does not appear to be from sim3a_i." << "\033[0m" << endl;
        exit(1);
    }

    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", argv[2]);
    filter(&df);
    setup_dataframe(&df);
    cut_circle(&df);

//*** mul = 2 ***//
    ROOT::RDF::RNode m2 = df.Filter("mul == -1");
    TCanvas* c1 = new TCanvas("c1", "c", 600, 600);
    TH2D* h1 = dalitz_slice(&m2, 200, true, false);
    setup_dalitz_plot(h1);
    h1->SetNdivisions(1, "X");
    h1->SetNdivisions(1, "Y");
    h1->SetNdivisions(4, "Z");
    h1->Draw("colz");

    string path = string(argv[1]) + "phase_space_mul2.pdf"; 
    // c1->SetLogz();
    c1->SetRightMargin(0.15);
    c1->SaveAs(path.c_str());


//*** mul == 3 ***//
    ROOT::RDF::RNode m3 = df.Filter("mul == 3");
    TCanvas* c2 = new TCanvas("c2", "c", 600, 600); 
    TH2D* h2 = dalitz_slice(&m3, 200, true, false);
    setup_dalitz_plot(h2);
    h2->SetNdivisions(1, "X");
    h2->SetNdivisions(1, "Y");
    h2->SetNdivisions(4, "Z");
    h2->Draw("colz");

    path = string(argv[1]) + "phase_space_mul3.pdf"; 
    // c2->SetLogz();
    c2->SetRightMargin(0.15);
    c2->SaveAs(path.c_str());
}