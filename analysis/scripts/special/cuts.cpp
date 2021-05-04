// ROOT stuff
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TCanvas.h>
#include <TLine.h>

// other stuff
#include <boost/format.hpp>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using namespace ROOT;
using boost::format;

int main(int argc, char *argv[]) {
    setup_style();
    TChain chain("tree");
    for (int i = 2; i < argc-1; i++) {
        chain.Add(argv[i]);
    }    

    // prepare the dataframe
    ROOT::RDF::RNode raw = RDataFrame(chain).Define("mp", "p_tot/1000").Define("mexC12", "exC12/1000"); // scale p to MeV
    ROOT::RDF::RNode tdc = RDataFrame("tree", argv[1]).Define("mp", "p_tot/1000").Define("mexC12", "exC12/1000");

    // set the axes
    double x_axis[] = {400, 8, 20};
    double y_axis[] = {400, 0, 250};

    //*** RAW PLOT ***//
    TCanvas* c1 = new TCanvas("raw", "", 600, 600);
    TH2D hraw = raw.Histo2D({"hraw", "", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "mexC12", "mp").GetValue();
    hraw.GetXaxis()->SetTitle("E_{tot} [MeV]");
    hraw.GetXaxis()->CenterTitle();
    hraw.GetXaxis()->SetNdivisions(205);
    hraw.GetYaxis()->SetTitle("p_{tot} [MeV]");
    hraw.GetYaxis()->CenterTitle();
    hraw.GetYaxis()->SetNdivisions(205);
    hraw.Draw("colz");

    c1->SetLogz();
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.12);
    string path = string(argv[argc-1]) + "cuts_raw.pdf";
    c1->SaveAs(path.c_str());

    //*** TDC CALIBRATED PLOT ***//
    TCanvas* c2 = new TCanvas("tdc", "", 600, 600);
    TH2D htdc = tdc.Histo2D({"htdc", "", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "mexC12", "mp").GetValue();
    htdc.GetXaxis()->SetTitle("E_{tot} [MeV]");
    htdc.GetXaxis()->CenterTitle();
    htdc.GetXaxis()->SetNdivisions(205);
    htdc.GetYaxis()->SetTitle("p_{tot} [MeV]");
    htdc.GetYaxis()->CenterTitle();
    htdc.GetYaxis()->SetNdivisions(205);
    htdc.Draw("colz");

    c2->SetLogz();
    c2->SetRightMargin(0.15);
    c2->SetLeftMargin(0.12);
    path = string(argv[argc-1]) + "cuts_tdc.pdf";
    c2->SaveAs(path.c_str());

    //*** ENERGY & MOMENTUM CUT ***//
    double pmax = 100;
    double p_cut = 40;
    double deltaE_cut = 0.3;
    double resonance_E = 17.76;
    x_axis[1] = 14; x_axis[2] = 20;
    TCanvas* c3 = new TCanvas("c3", "", 600, 600);
    htdc = tdc.Histo2D({"htdc", "", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), 0, pmax}, "mexC12", "mp").GetValue();
    htdc.GetXaxis()->SetTitle("E_{tot} [MeV]");
    htdc.GetXaxis()->CenterTitle();
    htdc.GetXaxis()->SetNdivisions(205);
    htdc.GetYaxis()->SetTitle("p_{tot} [MeV]");
    htdc.GetYaxis()->CenterTitle();
    htdc.GetYaxis()->SetNdivisions(205);
    htdc.Draw("colz");

    TLine* l1 = new TLine(x_axis[1], p_cut, x_axis[2], p_cut);
    l1->SetLineColor(kRed);
    l1->SetLineWidth(3);
    l1->Draw("same");

    TLine* l2 = new TLine(resonance_E-deltaE_cut, y_axis[1], resonance_E-deltaE_cut, pmax);
    l2->SetLineColor(kRed);
    l2->SetLineWidth(3);
    l2->Draw("same");

    TLine* l3 = new TLine(resonance_E+deltaE_cut, y_axis[1], resonance_E+deltaE_cut, pmax);
    l3->SetLineColor(kRed);
    l3->SetLineWidth(3);
    l3->Draw("same");

    c3->SetLogz();
    c3->SetLeftMargin(0.12);
    c3->SetRightMargin(0.15);
    path = string(argv[argc-1]) + "cuts_Ep.pdf";
    c3->SaveAs(path.c_str());
}