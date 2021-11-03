// ROOT stuff
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TExec.h>
#include <TLatex.h>

// other stuff
#include <filesystem>
#include <fstream>
#include <boost/format.hpp>
#include <math.h>
#include <time.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using namespace ROOT;

int main(int argc, char *argv[]) {
    if (argc != 6) {
        cout << "Usage: ./3-_dalitz <output path> <data> <3- l = 1 sim> <3- l = 3 sim> <3- l = 5 sim>" << endl;
        exit(1);
    }
    string folder = string(argv[1]) + "3-_article/";

    // path to fit with free delta
    string path_fit = "figures/dalitz_fit/article/";
    double k1 = 0.918074;
    double k2 = 0.0374003;
    double k3 = 0.445252;

    std::cout << "Using hardcoded values for k and delta." << endl;
    filesystem::create_directories(folder);

    int bins = 200;
    labelsize = 0.06;
    titlesize = 0.07;
    xlabeloffset = 0.8;
    ylabeloffset = 0.65;
    setup_style();

    // prepare the data sets
    ROOT::RDF::RNode data = RDataFrame("tree", argv[2]);
    ROOT::RDF::RNode sim1 = RDataFrame("tree", argv[3]);
    ROOT::RDF::RNode sim3 = RDataFrame("tree", argv[4]);
    ROOT::RDF::RNode sim5 = RDataFrame("tree", argv[5]);
    filter(&data); // energy & momentum cuts
    filter(&sim1);
    filter(&sim3);
    filter(&sim5);
    setup_dataframe(&data); // define Dalitz coordinates
    setup_dataframe(&sim1);
    setup_dataframe(&sim3);
    setup_dataframe(&sim5);
    cut_circle(&data); // cut everything outside the unit circle
    cut_circle(&sim1);
    cut_circle(&sim3);
    cut_circle(&sim5);
    cut_gs(&data); // remove ground state decays
    cut_gs(&sim1);
    cut_gs(&sim3);
    cut_gs(&sim5);

//**********************************//
// FIRST PANEL COLUMN, DALITZ PLOTS //
//**********************************//
if (true) {
    TCanvas* c1 = new TCanvas("c1", "c", 600, 2000);
    c1->Divide(1, 5, 0, 0);
    double rightmargin = 0.14;
    double topmargin = 0.03;
    double bottommargin = 0.03;

    auto setup_plot_dalitz = [] (TH2D* h, string ylabel = "y") {
        h->GetXaxis()->CenterTitle();
        h->GetXaxis()->SetNdivisions(2);
        h->GetYaxis()->SetTitle(ylabel.c_str());
        h->GetYaxis()->CenterTitle();
        h->GetYaxis()->SetNdivisions(2);
        h->GetZaxis()->SetNdivisions(2);
    };

    double contours[10] = {0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7};

// first panel: l = 1 sim
    c1->cd(1);

    TH2D* dalitz_sim1 = dalitz(&sim1, bins, true);
    dalitz_sim1->Scale(1./dalitz_sim1->GetMaximum());
    setup_plot_dalitz(dalitz_sim1, "L = 1");

    // cosmetics
    dalitz_sim1->GetXaxis()->SetLabelSize(0);
    dalitz_sim1->SetContour(7, contours);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(bottommargin);

    dalitz_sim1->DrawCopy("colz");

// second panel: l = 3 sim
    c1->cd(2);

    TH2D* dalitz_sim3 = dalitz(&sim3, bins, true);
    dalitz_sim3->Scale(1/dalitz_sim3->GetMaximum());
    setup_plot_dalitz(dalitz_sim3, "L = 3");

    // cosmetics
    dalitz_sim3->GetXaxis()->SetLabelSize(0);
    dalitz_sim3->SetContour(7, contours);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(bottommargin);

    dalitz_sim3->DrawCopy("colz");

// third panel: l = 5 sim
    c1->cd(3);

    TH2D* dalitz_sim5 = dalitz(&sim5, bins, true);
    dalitz_sim5->Scale(1/dalitz_sim5->GetMaximum());
    setup_plot_dalitz(dalitz_sim5, "L = 5");

    // cosmetics
    dalitz_sim5->GetXaxis()->SetLabelSize(0);
    dalitz_sim5->SetContour(7, contours);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(bottommargin);

    dalitz_sim5->DrawCopy("colz");

// fourth panel: fit
    c1->cd(4);

    dalitz_sim1->Scale(k1);
    dalitz_sim3->Scale(k2);
    dalitz_sim5->Scale(k3);

    dalitz_sim1->Add(dalitz_sim3);
    dalitz_sim1->Add(dalitz_sim5);
    dalitz_sim1->Scale(1/dalitz_sim1->GetMaximum());

    setup_plot_dalitz(dalitz_sim1, "Fit");

    // cosmetics
    dalitz_sim1->GetXaxis()->SetLabelSize(0);
    dalitz_sim1->SetContour(7, contours);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(bottommargin);

    dalitz_sim1->Draw("colz");

// fifth panel: data
    c1->cd(5);
    TH2D* dalitz_data = dalitz(&data, bins, true);
    dalitz_data->Scale(1./dalitz_data->GetMaximum());
    setup_plot_dalitz(dalitz_data, "Data");

    // cosmetics
    dalitz_data->SetXTitle("x");
    dalitz_data->SetContour(7, contours);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(0.05);

    dalitz_data->Draw("colz");

// save
    string path = folder + "dalitz_panels.pdf";
    c1->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;
}

//*************************//
//*** BINWISE CHI2 PLOT ***//
//*************************//
if (true) {
    gStyle->SetPalette(kThermometer);
    TCanvas* c2 = new TCanvas("c2", "c", 600, 600);

    // read binwise chi2 distribution from disk
    bins = 100;
    vector<double> binwise_chi(pow(bins, 2));
    ifstream file(path_fit + "binwise_chi2.txt", ios::in);
    string line; int counter = 0;
    while (std::getline(file, line)) {
        binwise_chi[counter] = atof(line.c_str());
        counter++;
    }
    file.close();

    // perform the actual plot
    TH2D* chi2 = new TH2D("chi", "h", bins, 0, 1, bins, 0, 1);
    double minval = 0;
    int scaler = 1;
    for (int y = 0; y < bins; y++) {
        for (int x = 0; x < bins; x++) {
            int my_bin = x + y*bins; // the logical choice of bins
            int root_bin = (x + 1) + (y + 1)*(bins + 2); // index 0 is underflow bin and index bins+1 is overflow bin

            if (binwise_chi[my_bin] >= 1e6) { 
                chi2->SetBinContent(root_bin, -1e6);
            } else {
                chi2->SetBinContent(root_bin, binwise_chi[my_bin]);
            }
        }
    }
    chi2->SetMinimum(0);

    // cosmetics
    c2->SetRightMargin(0.15);
    c2->SetLeftMargin(0.13);
    chi2->SetLabelSize(0.035, "X");
    chi2->SetLabelSize(0.035, "Y");
    chi2->SetLabelSize(0.035, "Z");

    chi2->SetTitleSize(0.045, "X");
    chi2->SetTitleSize(0.045, "Y");
    chi2->SetXTitle("x");
    chi2->SetYTitle("y");
    chi2->GetXaxis()->CenterTitle();
    chi2->GetYaxis()->CenterTitle();
    chi2->GetXaxis()->SetNdivisions(2);
    chi2->GetYaxis()->SetNdivisions(2);
    chi2->GetXaxis()->SetTitleOffset(1);
    chi2->GetYaxis()->SetTitleOffset(1);
    chi2->GetZaxis()->SetTitleOffset(1.2);

    TLatex* zlabel = new TLatex(0.98, 0.5, "\\chi^{2}");
    zlabel->SetNDC();
    zlabel->SetTextAngle(270);
    zlabel->SetTextSize(0.04);
    zlabel->SetTextAlign(22);

    chi2->Draw("colz");
    zlabel->Draw();
    
    // save
    string path = folder + "chi2.pdf";
    c2->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

    chi2->SetMaximum(25);
    path = folder + "chi2_max.pdf";
    c2->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;
}

//***********************//
//*** DIFFERENCE PLOT ***//
//***********************//
if (true) {
    gStyle->SetPalette(kThermometer);
    TCanvas* c4 = new TCanvas("c4", "c", 600, 600);
    bins = 100;

    TH2D* slice_data = dalitz_slice(&data, bins, true);
    TH2D* slice_sim1 = dalitz_slice(&sim1, bins, true);
    TH2D* slice_sim3 = dalitz_slice(&sim3, bins, true);
    TH2D* slice_sim5 = dalitz_slice(&sim5, bins, true);

    slice_sim1->Scale(k1/slice_sim1->GetMaximum());
    slice_sim3->Scale(k2/slice_sim3->GetMaximum());
    slice_sim5->Scale(k3/slice_sim5->GetMaximum());

    slice_sim1->Add(slice_sim3);
    slice_sim1->Add(slice_sim5);

    TH2D* diff = new TH2D("diff", "h", bins, 0, 1, bins, 0, 1);
    slice_data->Scale(1/slice_data->GetMaximum());
    slice_sim1->Scale(1/slice_sim1->GetMaximum());
    for (int x = 1; x < bins+1; x++) {
        for (int y = 1; y < bins+1; y++) {
            double bdat = slice_data->GetBinContent(x, y);
            double bfit = slice_sim1->GetBinContent(x, y);
            if (bdat > 0 && bfit > 0) {
                diff->SetBinContent(x, y, bdat-bfit);
            } 
        }
    }
    diff->Scale(1/max(diff->GetMaximum(), abs(diff->GetMinimum())));

    // cosmetics
    c4->SetRightMargin(0.15);

    // double diffcont[] = {diff->GetMinimum(), -40, -20, -10, -5, 5, 10, 20, 40, diff->GetMaximum()}; // upscaled contours
    // double diffcont[] = {diff->GetMinimum(), -.2, -.1, -.05, -.025, .025, .05, .1, .2, diff->GetMaximum()}; // raw difference contours
    double diffcont[] = {-1, -.6, -.3, -.15, -.075, .075, .15, .3, .6, 1}; // normalized difference contours
    if (diff->GetMaximum() < 1) {
        double m = diff->GetMaximum();
        for (int i = 5; i < 10; i++) {
            diffcont[i] *= m;            
        }
    }
    diff->SetContour(10, diffcont);

    diff->GetXaxis()->SetTitle("x");
    diff->GetXaxis()->SetTitleSize(0.045);
    diff->GetXaxis()->SetLabelSize(0.035);
    diff->GetXaxis()->SetTitleOffset(1);
    diff->GetXaxis()->CenterTitle();
    diff->GetXaxis()->SetNdivisions(2);

    diff->GetYaxis()->SetTitle("y");
    diff->GetYaxis()->SetTitleSize(0.045);
    diff->GetYaxis()->SetLabelSize(0.035);
    diff->GetYaxis()->SetTitleOffset(1);
    diff->GetYaxis()->CenterTitle();
    diff->GetYaxis()->SetNdivisions(2);

    diff->GetZaxis()->SetLabelSize(0.035);
    diff->GetZaxis()->SetTitleOffset(1.2);

    TText* zlabel2 = new TText(0.98, 0.53, "Difference");
    zlabel2->SetNDC();
    zlabel2->SetTextAngle(270);
    zlabel2->SetTextSize(0.04);
    zlabel2->SetTextAlign(22);
    zlabel2->SetTextFont(42);

    diff->Draw("colz1");
    zlabel2->Draw();

    // save
    string path = folder + "difference.pdf";
    c4->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;
}

//**********************************************//
// SECOND PANEL COLUMN, ENERGY PROJECTION PLOTS //
//**********************************************//
if (true) {
    TCanvas* c3 = new TCanvas("c2", "c", 600, 600);
    // c3->Divide(1, 3, 0, 0);

    bins = 100;
    vector<double> x_axis = {0, 6000};
    TH1D* data_E = energy_projection(&data, false, "data_E", "none");
    TH1D* sim1_E = energy_projection(&sim1, false, "sim1_E", "none");
    TH1D* sim3_E = energy_projection(&sim3, false, "sim3_E", "none");
    TH1D* sim5_E = energy_projection(&sim5, false, "sim5_E", "none");

    data_E->Scale(1/data_E->Integral());
    sim1_E->Scale(k1/sim1_E->Integral());
    sim3_E->Scale(k2/sim3_E->Integral());
    sim5_E->Scale(k3/sim5_E->Integral());
    sim1_E->Add(sim3_E);
    sim1_E->Add(sim5_E);
    sim1_E->Scale(1/sim1_E->Integral());

    double scale = max(data_E->GetMaximum(), sim1_E->GetMaximum());
    data_E->Scale(1/scale);
    sim1_E->Scale(1/scale);
    data_E->SetMaximum(1);

    auto setup_plot_energy = [] (TH1D* hdat, TH1D* hsim, string xlabel, string ylabel) {
        hdat->GetXaxis()->SetTitle(xlabel.c_str());
        hdat->GetYaxis()->SetTitle(ylabel.c_str());
        hdat->GetXaxis()->CenterTitle();
        hdat->GetYaxis()->CenterTitle();
        hdat->GetXaxis()->SetNdivisions(205);
        hdat->GetYaxis()->SetNdivisions(203);
        hdat->GetYaxis()->SetTitleOffset(1);
        hdat->GetYaxis()->SetTitleSize(titlesize);
        hdat->GetXaxis()->SetTitleSize(titlesize);
        
        hsim->SetLineColor(kOrange+1);
        hdat->SetLineColor(kBlack);
        hsim->SetLineWidth(2);
        hdat->SetLineWidth(2);
    };

    data_E->GetXaxis()->SetLabelSize(0);
    double leftmargin = 0.15;
    double rightmargin = 0.055;
    double topmargin = 0.03;
    double bottommargin = 0.03;

// first panel: data + fit
    c3->cd(1);
    setup_plot_energy(data_E, sim1_E, "E_{cm}", "Arbitrary scale");
    data_E->DrawClone("hist l");
    sim1_E->Draw("hist l same");

    // cosmetics
    data_E->GetXaxis()->SetLabelSize(labelsize);
    gPad->SetTopMargin(topmargin);
    gPad->SetLeftMargin(leftmargin);
    gPad->SetRightMargin(rightmargin);
    gPad->SetBottomMargin(0.05);

// save
    string path = folder + "energy_panels.pdf";
    c3->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;
}
}