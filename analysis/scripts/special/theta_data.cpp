// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>
#include <TVector3.h>

// other stuff
#include <boost/format.hpp>
#include <iostream>
#include <math.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cout << "Usage: ./theta_proj <output path> <input>" << endl;
        exit(1);
    }

    setup_style();
    string path = string(argv[1]) + "angular_projection.pdf";
    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", argv[2]);

    auto theta = [] (double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double r3x, double r3y, double r3z) {
        TVector3 p23 = {r2x-r3x, r2y-r3y, r2z-r3z};
        TVector3 p1 = {r1x, r1y, r1z};
        return p1.Angle(p23);

    };

    TCanvas* c = new TCanvas("c", "c", 600, 600);
    for (int i = 0; i < 3; i++) { // flatten the momentum vectors
        df = df.Define((format("px%1%") % i).str().c_str(), (format("px[%1%]") % i).str().c_str())
                .Define((format("py%1%") % i).str().c_str(), (format("py[%1%]") % i).str().c_str())
                .Define((format("pz%1%") % i).str().c_str(), (format("pz[%1%]") % i).str().c_str());
    }
    TH1D h = df.Define("theta", theta, {"px0", "py0", "pz0", "px1", "py1", "pz1", "px2", "py2", "pz2"}).Histo1D({"h", "h", 100, 0, 3.1415}, "theta").GetValue();

    h.Draw();
    c->SaveAs(path.c_str());
}