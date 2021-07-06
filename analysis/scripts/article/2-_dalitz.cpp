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
        cout << "Usage: ./2-_dalitz <output path> <data> <l = 1 sim> <l = 3 sim> <both sim>" << endl;
        exit(1);
    }
    bool delta_is_fixed = true;

    string folder = string(argv[1]) + "2-_article/";
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
    ROOT::RDF::RNode sim2 = RDataFrame("tree", argv[4]);
    ROOT::RDF::RNode both = RDataFrame("tree", argv[5]);
    filter(&data); // energy & momentum cuts
    filter(&sim1);
    filter(&sim2);
    filter(&both);
    setup_dataframe(&data); // define Dalitz coordinates
    setup_dataframe(&sim1);
    setup_dataframe(&sim2);
    setup_dataframe(&both);
    cut_circle(&data);
    cut_circle(&sim1);
    cut_circle(&sim2);
    cut_circle(&both);
    cut_gs(&data);
    cut_gs(&sim1);
    cut_gs(&sim2);
    cut_gs(&both);

//**********************************//
// FIRST PANEL COLUMN, DALITZ PLOTS //
//**********************************//
    TCanvas* c1 = new TCanvas("c1", "c", 600, 1800);
    c1->Divide(1, 4, 0, 0);
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

// first panel: data
    c1->cd(1);
    TH2D* dalitz_data = dalitz(&data, bins, true);
    double data_scale = dalitz_data->GetMaximum();
    dalitz_data->Scale(1./data_scale);
    setup_plot_dalitz(dalitz_data, "Data");

    // cosmetics
    fill_empty(dalitz_data);
    dalitz_data->GetXaxis()->SetLabelSize(0);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(bottommargin);

    dalitz_data->DrawCopy("colz");

// second panel: l = 1
    c1->cd(2);
    TH2D* dalitz_sim1 = dalitz(&sim1, bins, true);
    dalitz_sim1->Scale(1./dalitz_sim1->GetMaximum());
    setup_plot_dalitz(dalitz_sim1, "L = 1");

    // cosmetics
    fill_empty(dalitz_sim1);
    dalitz_sim1->GetXaxis()->SetLabelSize(0);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(bottommargin);

    dalitz_sim1->Draw("colz");

// third panel: l = 3
    c1->cd(3);
    TH2D* dalitz_sim2 = dalitz(&sim2, bins, true);
    dalitz_sim2->Scale(1./dalitz_sim2->GetMaximum());
    setup_plot_dalitz(dalitz_sim2, "L = 3");

    // cosmetics
    fill_empty(dalitz_sim2);
    dalitz_sim2->GetXaxis()->SetLabelSize(0);
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);
    gPad->SetBottomMargin(bottommargin);

    dalitz_sim2->Draw("colz");

// fourth panel: mix
    c1->cd(4);
    std::cout << "Using hardcoded values for k and delta." << endl;
    double k = 0.197965;
    double delta = 0.673937*2*M_PI;

    if (delta_is_fixed) {
        k = 0.139746;
        delta = 0;
    }

    auto calc_weights = [&k, &delta] (vector<vector<double>> f, double wU) { 
        return wU*(k*f[0][0]+(1-k)*f[0][1] + 2*sqrt(k*(1-k))*(f[0][2]*cos(delta) + f[0][3]*sin(delta)));
    };
    both = both.Define("w", calc_weights, {"f", "wU"});

    TH2D* dalitz_mix = dalitz(&both, bins, true, true);
    dalitz_mix->Scale(1./dalitz_mix->GetMaximum());
    setup_plot_dalitz(dalitz_mix, "Fit");

    // cosmetics
    dalitz_mix->SetXTitle("x");
    gPad->SetRightMargin(rightmargin);
    gPad->SetTopMargin(topmargin);

    dalitz_mix->Draw("colz");

// save
    string path = folder + "dalitz_panels.pdf";
    c1->SetLogz();
    c1->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

// chi2 slice plot
    TCanvas* c2 = new TCanvas("c2", "c", 600, 600);

    //*** DIFFERENCE ***//
    // dalitz_data->Scale(data_scale);
    // dalitz_mix->Scale(data_scale);

    // dalitz_data->Add(dalitz_mix, -1); // subtract the two plots
    // double minval = dalitz_data->GetMinimum();
    // for (int x = 1; x < bins+1; x++) {
    //     for (int y = 1; y < bins+1; y++) {
    //         if (abs(dalitz_data->GetBinContent(x, y)) < 1) {
    //             dalitz_data->SetBinContent(x, y, -1e6);
    //         }
    //     }
    // }
    // dalitz_data->SetMinimum(minval);

    //*** BINWISE CHI2 ***//
    // read binwise chi2 distribution from disk
    bins = 100;
    vector<double> binwise_chi(pow(bins, 2));
    string load_path = "figures/dalitz_fit/true_events 2-_i/binwise_chi2.txt";

    if (delta_is_fixed) {
        load_path = "figures/dalitz_fit/true_events 2-_i fixed/binwise_chi2.txt";
    }
    ifstream file(load_path, ios::in);
    string line; int counter = 0;
    while (std::getline(file, line)) {
        binwise_chi[counter] = atof(line.c_str());
        counter++;
    }
    file.close();

    // perform the actual plot
    TH2D* chi2 = new TH2D("chi", "h", bins, 0, 1, bins, 0, 1);
    double minval = 0;
    for (int y = 0; y < bins; y++) {
        for (int x = 0; x < bins; x++) {
            int my_bin = x + y*bins; // the logical choice of bins
            int root_bin = (x + 1) + (y + 1)*(bins + 2); // index 0 is underflow bin and index bins+1 is overflow bin

            if (binwise_chi[my_bin] >= 1e6) { 
                chi2->SetBinContent(root_bin, -1e6);
            } else {
                chi2->SetBinContent(root_bin, binwise_chi[my_bin]);

                // change sign depending on if data or simulation is higher
                // if (dalitz_data->GetBinContent(root_bin) < dalitz_mix->GetBinContent(root_bin)) {
                //     minval = min(minval, -binwise_chi[my_bin]);
                //     chi2->SetBinContent(root_bin, -binwise_chi[my_bin]);
                // } else {
                //     chi2->SetBinContent(root_bin, binwise_chi[my_bin]);
                // }
            } 
        }
    }
    chi2->SetMinimum(minval);

    // Mirror the plot to form the full Dalitz-plot
    // TH2D* chi2 = new TH2D("chi", "h", 2*bins, -1, 1, 2*bins, -1, 1);
    
    // const vector<double> cost = {cos(0), cos(2./3*M_PI), cos(4./3*M_PI)};
    // const vector<double> sint = {sin(0), sin(2./3*M_PI), sin(4./3*M_PI)};

    // // in this first loop, each event is mirrored across the radial line at 60 degrees
    // for (int xbin = 0; xbin < bins; xbin++) {
    //     for (int ybin = 0; ybin < bins; ybin++) {
    //         int my_bin = xbin + ybin*bins; // the logical choice of bins
    //         int root_bin = (xbin + 1) + (ybin + 1)*(bins + 2); // index 0 is underflow bin and index bins+1 is overflow bin
    //         double val = binwise_chi[my_bin];
    //         if (val >= 1e6) {
    //             continue;
    //         }

    //         double x = -1./(2*bins) + (double) xbin/bins; // cast to double to avoid integer division
    //         double y = -1./(2*bins) + (double) ybin/bins;
    //         double phi = atan2(y, x);
    //         double X = x*cos(2*phi-M_PI/3) + y*sin(2*phi-M_PI/3); // mirrored x
    //         double Y = -x*sin(2*phi-M_PI/3) + y*cos(2*phi-M_PI/3); // mirrored y

    //         // in this second loop, both the original and mirrored event are duplicated and rotated by 120 degrees
    //         // the rotation matrix is precalculated for efficiency
    //         for (int r = 0; r < 3; r++) {
    //             double x1 = x*cost[r] + y*sint[r];
    //             double x2 = X*cost[r] + Y*sint[r];
    //             double y1 = -x*sint[r] + y*cost[r];
    //             double y2 = -X*sint[r] + Y*cost[r];

    //             chi2->SetBinContent(chi2->FindBin(x1, y1), val);
    //             chi2->SetBinContent(chi2->FindBin(x2, y2), val);
    //         }
    //     }
    // }

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

    // hack solution to change palette
    TExec *ex = new TExec("ex","gStyle->SetPalette(kThermometer);");
    chi2->Draw("colz");
    ex->Draw();
    chi2->Draw("colz same");
    zlabel->Draw();
    
    // save
    path = folder + "chi2.pdf";
    c2->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

//**********************************************//
// SECOND PANEL COLUMN, ENERGY PROJECTION PLOTS //
//**********************************************//
    TCanvas* c3 = new TCanvas("c2", "c", 600, 1400);
    c3->Divide(1, 3, 0, 0);

    vector<double> x_axis = {0, 6000};
    TH1D* data_E = new TH1D("data_E", "h", bins, x_axis[0], x_axis[1]);
    TH1D* sim1_E = new TH1D("sim1_E", "h", bins, x_axis[0], x_axis[1]);
    TH1D* sim2_E = new TH1D("sim2_E", "h", bins, x_axis[0], x_axis[1]);
    TH1D* both_E = new TH1D("both_E", "h", bins, x_axis[0], x_axis[1]);

    for (int i = 0; i < 3; i++) {
        TH1D data_temp = data.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"dat_temp", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
        TH1D sim1_temp = sim1.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"sim1_temp", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
        TH1D sim2_temp = sim2.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"sim2_temp", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
        TH1D both_temp = both.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"both_temp", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
        data_E->Add(&data_temp);
        sim1_E->Add(&sim1_temp);
        sim2_E->Add(&sim2_temp);
        both_E->Add(&both_temp);
    }
    sim1_E->Scale(1/sim1_E->Integral());
    sim2_E->Scale(1/sim2_E->Integral());
    both_E->Scale(1/both_E->Integral());
    data_E->Scale(1/data_E->Integral());

    sim1_E->Scale(1/data_E->GetMaximum());
    sim2_E->Scale(1/data_E->GetMaximum());
    both_E->Scale(1/data_E->GetMaximum());
    data_E->Scale(1/data_E->GetMaximum());

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
    rightmargin = 0.05;
// first panel: data + fit l = 1
    c3->cd(1);
    setup_plot_energy(data_E, sim1_E, "E_{cm}", "L = 1");
    data_E->DrawClone("hist l");
    sim1_E->Draw("hist l same");

    // cosmetics
    gPad->SetTopMargin(topmargin);
    gPad->SetLeftMargin(leftmargin);
    gPad->SetRightMargin(rightmargin);
    gPad->SetBottomMargin(bottommargin);

// second panel: data + fit l = 3
    c3->cd(2);
    setup_plot_energy(data_E, sim2_E, "E_{cm}", "L = 3");
    data_E->DrawClone("hist l");
    sim2_E->Draw("hist l same");

    // cosmetics
    gPad->SetTopMargin(topmargin);
    gPad->SetLeftMargin(leftmargin);
    gPad->SetRightMargin(rightmargin);
    gPad->SetBottomMargin(bottommargin);

// third panel: data + fit both
    c3->cd(3);
    setup_plot_energy(data_E, both_E, "E_{cm}", "Fit");
    data_E->Draw("hist l");
    both_E->Draw("hist l same");

    // cosmetics
    data_E->GetXaxis()->SetLabelSize(labelsize);
    gPad->SetTopMargin(topmargin);
    gPad->SetLeftMargin(leftmargin);
    gPad->SetRightMargin(rightmargin);

// save
    path = folder + "energy_panels.pdf";
    c3->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;
}