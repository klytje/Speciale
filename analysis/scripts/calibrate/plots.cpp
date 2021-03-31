#pragma once
#include <iostream>
#include <cmath>
#include <TROOT.h>
#include <TF1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLine.h>

using namespace std;

// plot settings
struct plot {
    static const bool save = true; // note that only "enabled" figures are actually created and can be saved
    static inline string path = "analysis/figures/";
    static const inline string format = ".png";

    // detector plot controllers. prioritized over any individual figure (mainly used for debugging)
    static const bool W = true;
    static const bool S = true;

    // DATA CORRECTION FIGURES
    static const bool repair_peaks = true;
    static const bool align_peaks = true;
    static const bool gauss_cut = true; // FT and BT cut

    // TDC CALIBRATION FIGURES
    static const bool raw = true;
    static const bool centered = true;
    static const bool fit = true;
    static const bool gauss = true; // DT cut
    static const bool gauss_sigma = true;

    // a simple wrapper around the two plot controllers; needed by the code itself
    static const inline vector<bool> plot_det = {S, S, W, W};

    // axes for the histograms
    struct x_axes {
        static inline vector<double> energy = {100, -500, 500}; // standard limits but with fewer bins; needed for the energy histograms
        static inline vector<double> standard = {1000, -500, 500}; // standard limits; also used as bounds for the dt filters
        static inline vector<double> centered = {100, -50, 50}; // centered limits; used after the data has been centered on 0
        static inline vector<double> gauss = {500, -50, 50}; // x-axis for the gaussian plots. more bins since it's 1d
        static inline vector<double> FT = {500, 65000, 65600}; // x-axis for the gaussian plots. more bins since it's 1d
        static inline vector<double> BT = {500, 65000, 65600}; // x-axis for the gaussian plots. more bins since it's 1d
    };

    // a struct to avoid polluting this namespace with y-axis stuff. use the y_axis map to access them
    struct y_axes {
        static const inline vector<double> SU = {24, 1, 25};
        static const inline vector<double> SD = {24, 1, 25};
        static const inline vector<double> Det1 = {16, 1, 17};
        static const inline vector<double> Det2 = {16, 1, 17};
    };
    static const inline map<int, vector<double>> y_axis = {{0, y_axes::SU}, {1, y_axes::SD}, {2, y_axes::Det1}, {3, y_axes::Det2}}; // accessor for the y-axes
};

// a simple container for histogram information
class hist_info {
    public: 
        // 1D histograms - no y-axis
        hist_info(vector<double> x_axis, string filename, string title, string x_label, string y_label) {
            this->x_axis = x_axis;
            this->x_label = x_label;
            this->y_label = y_label;
            this->title = title;
            this->filename = plot::path + filename + plot::format;
        }

        // 2D histograms
        hist_info(vector<double> x_axis, vector<double> y_axis, string filename, string title, string x_label, string y_label) {
            this->x_axis = x_axis;
            this->x_label = x_label;
            this->y_axis = y_axis;
            this->y_label = y_label;
            this->title = title;
            this->filename = plot::path + filename + plot::format;
        }

        vector<double> x_axis;
        vector<double> y_axis;
        string title;
        string x_label;
        string y_label;
        string filename;
};

// plot a 2d histogram
void hist2D(const vector<double> *DT, const vector<int> *BI, const hist_info info) {
    const vector<double> &DTref = *DT;
    const vector<int> &BIref = *BI;

    TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);
    TH2D *h = new TH2D("2D Histogram", info.title.c_str(), int(info.x_axis[0]), info.x_axis[1], info.x_axis[2], int(info.y_axis[0]), info.y_axis[1], info.y_axis[2]);

    // fill the histogram
    for (int i = 0; i < DTref.size(); i++) {
        h->Fill(DTref[i], BIref[i]);
    }
    std::cout << "    Histogram filled with " << h->GetEntries() << " data points." << endl;

    // set up the figure
    h->GetXaxis()->SetTitle(info.x_label.c_str());
    h->GetYaxis()->SetTitle(info.y_label.c_str());
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->Draw("colz");
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->Modified(); canvas->Update();

    if (plot::save) {
        canvas->SaveAs(info.filename.c_str());
    }

    canvas->WaitPrimitive();
    canvas->Close();
    h->Delete();
}

// plot a 2d histogram with vertical lines from a gaussian fit
void hist2D(const vector<double> *DT, const vector<int> *BI, double sigma, int sigma_n, const hist_info info) {
    const vector<double> &DTref = *DT;
    const vector<int> &BIref = *BI;

    TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);
    TH2D *h = new TH2D("2D Histogram", info.title.c_str(), int(info.x_axis[0]), info.x_axis[1], info.x_axis[2], int(info.y_axis[0]), info.y_axis[1], info.y_axis[2]);
    TLine *vline_l = new TLine(sigma_n*sigma, info.y_axis[1], sigma_n*sigma, info.y_axis[2]);
    TLine *vline_r = new TLine(-sigma_n*sigma, info.y_axis[1], -sigma_n*sigma, info.y_axis[2]);
    vline_l->SetLineWidth(2);
    vline_r->SetLineWidth(2);
    
    // fill the histogram
    for (int i = 0; i < DTref.size(); i++) {
        h->Fill(DTref[i], BIref[i]);
    }
    std::cout << "    Histogram filled with " << h->GetEntries() << " data points." << endl;

    // set up the figure
    h->GetXaxis()->SetTitle(info.x_label.c_str());
    h->GetYaxis()->SetTitle(info.y_label.c_str());
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->Draw("colz");
    vline_l->Draw();
    vline_r->Draw();

    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->Modified(); canvas->Update();

    if (plot::save) {
        canvas->SaveAs(info.filename.c_str());
    }

    canvas->WaitPrimitive();
    canvas->Close();
    h->Delete();
}

// plot a 1d histogram with vertical lines from a gaussian fit
template <typename T>
void hist1D(const vector<T> *input, double mean, int sigma, const hist_info info) {
    const vector<T> &data = *input;

    TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);
    TH1D *h = new TH1D("1D Histogram", info.title.c_str(), int(info.x_axis[0]), info.x_axis[1], info.x_axis[2]);    
    // fill the histogram
    for (int i = 0; i < data.size(); i++) {
        h->Fill(data[i]);
    }
    std::cout << "    Histogram filled with " << h->GetEntries() << " data points." << endl;

    TLine *vline_l = new TLine(mean + sigma, 0, mean + sigma, 10000);
    TLine *vline_r = new TLine(mean - sigma, 0, mean - sigma, 10000);
    vline_l->SetLineWidth(2);
    vline_r->SetLineWidth(2);

    // set up the figure
    h->GetXaxis()->SetTitle(info.x_label.c_str());
    h->GetYaxis()->SetTitle(info.y_label.c_str());
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetMaxDigits(3);
    h->Draw();
    vline_l->Draw();
    vline_r->Draw();
    
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->Modified(); canvas->Update();

    if (plot::save) {
        canvas->SaveAs(info.filename.c_str());
    }

    canvas->WaitPrimitive();
    canvas->Close();
    h->Delete();
    return;
}

// plot a 1d histogram
template <typename T>
void hist1D(const vector<T> *input, const hist_info info) {
    const vector<T> &data = *input;

    TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);
    TH1D *h = new TH1D("1D Histogram", info.title.c_str(), int(info.x_axis[0]), info.x_axis[1], info.x_axis[2]);    

    // fill the histogram
    for (int i = 0; i < data.size(); i++) {
        h->Fill(data[i]);
    }
    std::cout << "    Histogram filled with " << h->GetEntries() << " data points." << endl;

    // set up the figure
    h->GetXaxis()->SetTitle(info.x_label.c_str());
    h->GetYaxis()->SetTitle(info.y_label.c_str());
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetMaxDigits(3);
    h->Draw();
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->Modified(); canvas->Update();

    if (plot::save) {
        canvas->SaveAs(info.filename.c_str());
    }

    canvas->WaitPrimitive();
    canvas->Close();
    h->Delete();
    return;
}
