// ROOT stuff
#include <TROOT.h>
#include <TCanvas.h>
#include <TArc.h>
#include <TLine.h>

// other stuff
#include <boost/format.hpp>
#include <math.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using namespace ROOT;
using boost::format;

int main(int argc, char *argv[]) {
    setup_style();
    double y = 0.93; // cutoff
    double x = sqrt(1-pow(y, 2));
    double theta = atan(y/x)*180./M_PI;

    TH2D hist = dalitz_slice(argv[1], 100, false);

    //*** "before" plot with the cuts shown ***//
    TCanvas* c1 = new TCanvas("c1", "c", 600, 600);
    hist.GetXaxis()->SetTitle("x");
    hist.GetXaxis()->CenterTitle();
    hist.GetXaxis()->SetNdivisions(2);
    hist.GetYaxis()->SetTitle("y");
    hist.GetYaxis()->CenterTitle();
    hist.GetYaxis()->SetNdivisions(2);
    hist.Draw("colz");

    TArc* a1 = new TArc(0, 0, 1, 0, theta);
    a1->SetLineColor(kRed);
    a1->SetLineWidth(3);
    a1->SetFillColorAlpha(0, 0);
    a1->Draw("only same");

    TArc* a2 = new TArc(0, 0, 1, theta, 90);
    a2->SetLineColor(kRed);
    a2->SetLineWidth(3);
    a2->SetLineStyle(12);
    a2->SetFillColorAlpha(0, 0);
    a2->Draw("only same");

    TLine* l1 = new TLine(0, y, x, y);
    l1->SetLineColor(kRed);
    l1->SetLineWidth(3);
    l1->Draw("same");

    c1->SetLogz();
    c1->SetLogz();
    c1->SetRightMargin(0.15);

    string path = string(argv[2]) + "dalitz_slice.pdf";
    c1->SaveAs(path.c_str());

    //*** "after" plot with the cuts shown ***//
    hist = dalitz_slice(argv[1], 100, true);

    TCanvas* c2 = new TCanvas("c2", "c", 600, 600);
    hist.GetXaxis()->SetTitle("x");
    hist.GetXaxis()->CenterTitle();
    hist.GetXaxis()->SetNdivisions(1);
    hist.GetYaxis()->SetTitle("y");
    hist.GetYaxis()->CenterTitle();
    hist.GetYaxis()->SetNdivisions(1);
    hist.Draw("colz");
    a1->Draw("only same");
    l1->Draw("same");

    c2->SetLogz();
    c2->SetLogz();
    c2->SetRightMargin(0.15);

    path = string(argv[2]) + "dalitz_slice_cut.pdf";
    c2->SaveAs(path.c_str());
}