// ROOT stuff
#include <TAxis.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH2D.h>
#include <TLegend.h>

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

    //*** theta(y) PLOT ***//
    TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
    auto func = [] (double* x, double* par) {
        double y = par[0];
        return acos(x[0]/sqrt(1-pow(y, 2)));
    };

    TF1 *dummy = new TF1("dummy", "3.14", -1, 1);
    dummy->SetTitle("Angular correlation functions");
    dummy->SetMaximum(M_PI);
    dummy->SetMinimum(0);
    dummy->SetLineColor(kBlack);
    dummy->Draw();

    auto legend = new TLegend(0.2, 0.4);
    for (int i = 0; i < 8; i++) {
        double y = 0 + double(i)/8;
        double x = sqrt(1-pow(y, 2));
        string label = (format("y = %1%") % y).str();
        TF1 *corr = new TF1(label.c_str(), func, -x+0.0001, x-0.0001, 1);
        corr->SetParameter(0, y);
        corr->SetLineColor(i+1);
        corr->DrawClone("same");

        legend->AddEntry(label.c_str(), label.c_str(), "l");
    }
    legend->Draw();
    path = string(argv[1]) + "theta(y).pdf";
    c1->SaveAs(path.c_str());
}