// ROOT stuff
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>

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
    ROOT::RDF::RNode raw = RDataFrame(chain).Define("mp", "pt/1000").Define("mexC12", "exC12/1000"); // scale p to MeV
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

//*** TIME OF ARRIVAL ***//
    // NO CUT PERFORMED! The output file looks is a very nice-looking Gaussian peak, indicating that this is unnecessary. 
    // This is probably due to it being very similar to the already-performed TDC cuts.

    // find the largest time interval given three arrival times
    auto dt = [] (double t1, double t2, double t3, int mul) {
        if (mul == -1) { // for the reconstructed data, t3 = 0
            return t1-t2;
        }

        double dt1 = t1 - t3;
        double dt2 = t2 - t1;
        double dt3 = t3 - t2;
        if (abs(dt1) > abs(dt2) && abs(dt1) > abs(dt3)) {
            return dt1;
        } else if (abs(dt2) > abs(dt1) && abs(dt2) > abs(dt3)) {
            return dt2;
        } else {
            return dt3;
        }
    };

    TCanvas* c3 = new TCanvas("c3", "", 600, 600);
    TH1D hist = tdc.Define("BT0", "BT[0]")
                    .Define("BT1", "BT[1]")
                    .Define("BT2", "BT[2]")
                    .Define("dt", dt, {"BT0", "BT1", "BT2", "mul"})
                    .Histo1D("dt").GetValue();
    TF1* gauss = new TF1("gauss", "gaus", -50, 50);
    hist.Fit(gauss, "LQR");

    hist.Draw();
    gauss->Draw("same");
    path = string(argv[argc-1]) + "time_of_arrival_peak.pdf";
    c3->SaveAs(path.c_str());


//*** ENERGY & MOMENTUM CUT ***//
    double pmax = 100;
    double p_cut = 40;
    double deltaE_cut = 0.3;
    double resonance_E;
    if (string(argv[argc-1]).find("0+") != string::npos) {
        resonance_E = 17.76;
    } else if (string(argv[argc-1]).find("2-") != string::npos) {
        resonance_E = 16.62;
    } else if (string(argv[argc-1]).find("3-") != string::npos) {
        resonance_E = 18.35;
    } else {
        cout << "\033[1;31m" << "Could not determine Jp. Cannot generate final figure." << "\033[0m" << endl;
        exit(1);
    }
    x_axis[1] = 14; x_axis[2] = 20;
    TCanvas* c4 = new TCanvas("c4", "", 600, 600);
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

    c4->SetLogz();
    c4->SetLeftMargin(0.12);
    c4->SetRightMargin(0.15);
    path = string(argv[argc-1]) + "cuts_Ep.pdf";
    c4->SaveAs(path.c_str());
}