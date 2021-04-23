// ROOT stuff
#include <TAxis.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH2D.h>

// other stuff
#include <boost/format.hpp>

// my stuff
#include "../plot_style.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[]) {
    setup_style();
    int bins = 200;
    TCanvas* c = new TCanvas("c", "c", 600, 600);
    TH2D* hist = new TH2D("h1", "Dalitz plot", bins, -1, 1, bins, -1, 1);
    
    for (double x = -1; x < 1; x+=2./bins) {
        for (double y = -1; y < 1; y+=2./bins) {
            if (pow(x, 2) + pow(y, 2) <= 1) {
                double theta = acos(x/sqrt(1-pow(y, 2)));
                int binx = hist->GetXaxis()->FindBin(x);
                int biny = hist->GetYaxis()->FindBin(y);
                hist->SetBinContent(hist->GetBin(binx, biny), theta);
            }
        }
    }

    auto circle = [] (double* x, double* par) {
        return par[0]*sqrt(1-pow(x[0], 2));
    };

    hist->SetContour(100);
    hist->GetXaxis()->SetNdivisions(-2);
    hist->GetYaxis()->SetNdivisions(-2);
    hist->GetZaxis()->SetNdivisions(-2);
    hist->GetZaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "0");
    hist->GetZaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "#pi/2");
    hist->GetZaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "#pi");
    hist->Draw("colz");

    string path = string(argv[1]) + "theta.pdf";
    c->SaveAs(path.c_str());
}