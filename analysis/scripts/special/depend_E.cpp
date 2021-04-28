// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TCanvas.h>
#include <TChain.h>
#include <TLegend.h>
#include <TF1.h>

// other stuff
#include <boost/format.hpp>

// my stuff
#include "../utility.cpp"
#include "../calibrate/utility.cpp"

using namespace std;
using boost::format;

int main(int argc, char *argv[]) {    
    double x_axis[] = {500, 65100, 65550}; // area of the peaks we're looking at
    int detector = 1; // filter to the SD detector

    data_container data;
    prepare_match(argc, argv, &data, "mul==3");
    vector<double>* BT = data.get_double("BT");
    vector<double>* BE = data.get_double("BE");
    vector<int>* ID = data.get_int("ID");
    vector<double> bt = detector_filter(BT, ID, detector);
    vector<double> be = detector_filter(BE, ID, detector);

    setup_style();
    vector<int> e1 = {10000, 4000, 3000, 2000, 1000}; // the energy bounds for each plot
    vector<int> e2 = {800, 600, 400, 200, 100}; // same, but for the second plot
    TH1D h1[e1.size()];
    TH1D h2[e2.size()];
    for (int i = 0; i < e1.size(); i++) {
        h1[i] = TH1D(("h1" + to_string(i)).c_str(), ("h1" + to_string(i)).c_str(), int(x_axis[0]), x_axis[1], x_axis[2]);
    } 
    for (int i = 0; i < e2.size(); i++) {
        h2[i] = TH1D(("h2" + to_string(i)).c_str(), ("h2" + to_string(i)).c_str(), int(x_axis[0]), x_axis[1], x_axis[2]);
    }

    for (int i = 0; i < bt.size(); i++) {
        double e = be[i];
        double t = bt[i];
        for (int j = 0; j < e1.size(); j++) {
            if (e <= e1[j]) {
                h1[j].Fill(t);
            } 
        }
        for (int j = 0; j < e2.size(); j++) {
            if (e <= e2[j]) {
                h2[j].Fill(t);
            } 
        }
    }    

    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 600);
    c1->Divide(e1.size(), 1, 0, 0);

    h1[0].GetYaxis()->SetTitle("Arbitrary units");
    h1[0].GetYaxis()->SetTitleSize(0.1);
    h1[0].GetYaxis()->SetTitleOffset(0.3);
    h1[0].GetYaxis()->CenterTitle();
    for (int i = 0; i < e1.size(); i++) {
        c1->cd(i+1);
        h1[i].Scale(1./h1[i].GetMaximum());
        h1[i].SetLineColor(kBlack);
        h1[i].SetLabelSize(0, "X");
        h1[i].SetTickLength(0, "X");
        h1[i].SetLabelSize(0, "Y");
        h1[i].SetTickLength(0, "Y");
        h1[i].GetXaxis()->SetTitle(to_string(e1[i]).c_str());
        h1[i].GetXaxis()->SetTitleSize(0.1);
        h1[i].GetXaxis()->SetTitleOffset(0.3);
        h1[i].GetXaxis()->CenterTitle();
        h1[i].Draw("HIST L");
    }
    string path = string(argv[argc-1]) + "large_E.pdf";
    c1->SaveAs(path.c_str());

    TCanvas* c2 = new TCanvas("c2", "c2", 1200, 600);
    c2->Divide(e2.size(), 1, 0, 0);

    h2[0].GetYaxis()->SetTitle("Arbitrary units");
    h2[0].GetYaxis()->SetTitleSize(0.1);
    h2[0].GetYaxis()->SetTitleOffset(0.3);
    h2[0].GetYaxis()->CenterTitle();
    for (int i = 0; i < e2.size(); i++) {
        c2->cd(i+1);
        h2[i].Scale(1./h2[i].GetMaximum());
        h2[i].SetLineColor(kBlack);
        h2[i].SetLabelSize(0, "X");
        h2[i].SetTickLength(0, "X");
        h2[i].SetLabelSize(0, "Y");
        h2[i].SetTickLength(0, "Y");
        h2[i].GetXaxis()->SetTitle(to_string(e1[i]).c_str());
        h2[i].GetXaxis()->SetTitleSize(0.1);
        h2[i].GetXaxis()->SetTitleOffset(0.3);
        h2[i].GetXaxis()->CenterTitle();
        h2[i].Draw("HIST L");
    }
    path = string(argv[argc-1]) + "small_E.pdf";
    c2->SaveAs(path.c_str());
}