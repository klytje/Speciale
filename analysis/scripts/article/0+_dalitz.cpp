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


// weights defined by eq 42 in Morten's thesis
double calc_weight(vector<vector<double>> f, double wU, double k, double delta) { 
    return wU*(k*f[0][0]+(1-k)*f[0][1] + 2*sqrt(k*(1-k))*(f[0][2]*cos(delta) + f[0][3]*sin(delta)));
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        cout << "Usage: ./0+_dalitz <output path> <data> <0+ sim> <2- sim3a_i> <3- sim3a_i>" << endl;
        exit(1);
    }
    string folder = string(argv[1]) + "0+_article/";

    // path to fit with free delta
    string path_fit = "figures/dalitz_fit/article_free/";
    double k1 = 0.0469975;
    double delta1 = 0.750007*2*M_PI;
    double c = 0.822108;
    double k2 = 0;
    double delta2 = 0;

    // path to fit with fixed delta
    // string path_fit = "figures/dalitz_fit/article_fixed/";
    // double k = 0.154447;
    // double delta = 0;

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
    ROOT::RDF::RNode sim0 = RDataFrame("tree", argv[3]);
    ROOT::RDF::RNode sim2 = RDataFrame("tree", argv[4]);
    ROOT::RDF::RNode sim3 = RDataFrame("tree", argv[5]);
    filter(&data); // energy & momentum cuts
    filter(&sim0);
    filter(&sim2);
    filter(&sim3);
    setup_dataframe(&data); // define Dalitz coordinates
    setup_dataframe(&sim0);
    setup_dataframe(&sim2);
    setup_dataframe(&sim3);
    cut_circle(&data); // cut everything outside the unit circle
    cut_circle(&sim0);
    cut_circle(&sim2);
    cut_circle(&sim3);
    cut_gs(&data); // remove ground state decays
    cut_gs(&sim0);
    cut_gs(&sim2);
    cut_gs(&sim3);

    auto w2 = [&k1, &delta1] (vector<vector<double>> f, double wU) {return calc_weight(f, wU, k1, delta1);};
    auto w3 = [&k2, &delta2] (vector<vector<double>> f, double wU) {return calc_weight(f, wU, k2, delta2);};
    sim2 = sim2.Define("w", w2, {"f", "wU"});
    sim3 = sim3.Define("w", w3, {"f", "wU"});

//**********************************//
// FIRST PANEL COLUMN, DALITZ PLOTS //
//**********************************//
if (true) {
    TCanvas* c1 = new TCanvas("c1", "c", 600, 1400);
    c1->Divide(1, 3, 0, 0);
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

    auto fill_empty = [&bins] (TH2D* h) {
        for (int bin = 1; bin < pow(bins+1, 2); bin++) {
            int xbin, ybin, zbin;
            h->GetBinXYZ(bin, xbin, ybin, zbin);
            double x = h->GetXaxis()->GetBinCenter(xbin), y = h->GetYaxis()->GetBinCenter(ybin);
            // cout << "bin: " << bin << ", x: " << x << ", y: " << y << endl;
            if (sqrt(pow(x, 2) + pow(y, 2)) <= 1 && h->GetBinContent(bin) == 0) {
                h->SetBinContent(bin, 1e-5); // just needs to be low enough to get colored as "0"
            }
        }
    };

    double contours[10] = {0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7};

// first panel: 0+ sim
    c1->cd(1);
    TH2D* dalitz_sim0 = dalitz(&sim0, bins, true);
    dalitz_sim0->Scale(1./dalitz_sim0->GetMaximum());
    setup_plot_dalitz(dalitz_sim0, "0^{+}");

    // cosmetics
    // fill_empty(dalitz_sim2l1);
    dalitz_sim0->GetXaxis()->SetLabelSize(0);
    dalitz_sim0->SetContour(7, contours);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(bottommargin);

    dalitz_sim0->Draw("colz");

// second panel: fit
    c1->cd(2);
    TH2D* dalitz_sim2 = dalitz(&sim2, bins, true, true);
    TH2D* dalitz_sim3 = dalitz(&sim3, bins, true, true);

    dalitz_sim2->Scale(c/dalitz_sim2->GetMaximum());
    dalitz_sim3->Scale((1-c)/dalitz_sim3->GetMaximum());
    dalitz_sim2->Add(dalitz_sim3);
    setup_plot_dalitz(dalitz_sim2, "Fit");

    // cosmetics
    dalitz_sim2->GetXaxis()->SetLabelSize(0);
    dalitz_sim2->SetContour(7, contours);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(bottommargin);

    dalitz_sim2->Draw("colz");

// third panel: data
    c1->cd(3);
    TH2D* dalitz_data = dalitz(&data, bins, true);
    double data_scale = dalitz_data->GetMaximum();
    dalitz_data->Scale(1./data_scale);
    setup_plot_dalitz(dalitz_data, "Data");

    // cosmetics
    // fill_empty(dalitz_data);
    dalitz_data->SetXTitle("x");
    dalitz_data->SetContour(7, contours);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(0.05);

    dalitz_data->DrawCopy("colz");

// save
    string path = folder + "dalitz_panels.pdf";
    c1->SetLogz();
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
    TH2D* slice_sim2 = dalitz_slice(&sim2, bins, true, true);
    TH2D* slice_sim3 = dalitz_slice(&sim3, bins, true, true);

    slice_sim2->Scale(c/slice_sim2->GetMaximum());
    slice_sim3->Scale((1-c)/slice_sim3->GetMaximum());
    slice_sim2->Add(slice_sim3);

    TH2D* diff = new TH2D("diff", "h", bins, 0, 1, bins, 0, 1);
    slice_data->Scale(1/slice_data->GetMaximum());
    slice_sim2->Scale(1/slice_sim2->GetMaximum());
    for (int x = 1; x < bins+1; x++) {
        for (int y = 1; y < bins+1; y++) {
            double bdat = slice_data->GetBinContent(x, y);
            double bfit = slice_sim2->GetBinContent(x, y);
            if (bdat > 0 && bfit > 0) {
                diff->SetBinContent(x, y, bdat-bfit);
            } 
        }
    }
    diff->Scale(1/diff->GetMaximum());

    // cosmetics
    c4->SetRightMargin(0.15);

    // double diffcont[] = {diff->GetMinimum(), -40, -20, -10, -5, 5, 10, 20, 40, diff->GetMaximum()}; // upscaled contours
    // double diffcont[] = {diff->GetMinimum(), -.2, -.1, -.05, -.025, .025, .05, .1, .2, diff->GetMaximum()}; // raw difference contours
    double diffcont[] = {-1, -.6, -.3, -.15, -.075, .075, .15, .3, .6, 1}; // normalized difference contours
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
    TCanvas* c3 = new TCanvas("c2", "c", 600, 1000);
    c3->Divide(1, 2, 0, 0);

    bins = 100;
    vector<double> x_axis = {0, 6000};
    TH1D* data_E = energy_projection(&data, false, "data_E", "integral");
    TH1D* sim0_E = energy_projection(&sim0, false, "sim0_E", "integral");
    TH1D* fit_E = energy_projection(&sim2, &sim3, c, true, true, "fit_E", "integral");

    double data_scale = max({data_E->GetMaximum(), sim0_E->GetMaximum(), fit_E->GetMaximum()});
    data_E->Scale(1/data_scale);
    sim0_E->Scale(1/data_scale);
    fit_E->Scale(1/data_scale);
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

// first panel: data + 0+ sim
    c3->cd(1);
    setup_plot_energy(data_E, sim0_E, "E_{cm}", "0+ simulation");
    data_E->DrawClone("hist l");
    sim0_E->Draw("hist l same");

    // cosmetics
    gPad->SetTopMargin(topmargin);
    gPad->SetLeftMargin(leftmargin);
    gPad->SetRightMargin(rightmargin);
    gPad->SetBottomMargin(bottommargin);

// second panel: data + fit
    c3->cd(2);
    setup_plot_energy(data_E, fit_E, "E_{cm}", "2- & 3- fit");
    data_E->Draw("hist l");
    fit_E->Draw("hist l same");

    // cosmetics
    data_E->GetXaxis()->SetLabelSize(labelsize);
    gPad->SetTopMargin(topmargin);
    gPad->SetLeftMargin(leftmargin);
    gPad->SetRightMargin(rightmargin);
    gPad->SetBottomMargin(0.06);

// save
    string path = folder + "energy_panels.pdf";
    c3->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;
}
}