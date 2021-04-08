#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <TChain.h>
#include <TStyle.h>

#include <filesystem>
#include <boost/format.hpp>
#include <iostream>
#include <cmath>

#include "plots.cpp"
#include "utility.cpp"

using namespace std;
using boost::format;

void repair_peaks(data_container* data) {
    print_title("### Attempting to repair broken peaks in the data ###");
    vector<double>* FT = data->get_double("FT");
    vector<double>* BT = data->get_double("BT");
    vector<int>* ID = data->get_int("ID");

    TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);

    // check FT or BT for misaligned peaks
    auto check_col = [&canvas] (vector<double> *T, vector<int> *t_hist, TH1D* t_plot, string name, string filename) {
        vector<double> &t = *T;
        vector<int> &bin = *t_hist;
        TH1D &plot = *t_plot;

        int left[2] = {0, 0}, right = 0;

        // calculate the mean so we can base step-req and background on the data instead of some arbitrary values
        double mean;
        {
            double sum = std::accumulate(std::begin(bin), std::end(bin), 0.0);
            mean =  sum / bin.size();
        }
        int step_req = mean*4; // if it jumps more than step_req in a single step, assume it is the edge we're looking for
        int background = mean/4; // background limit, previous bin must be lower than this
        for (int i = 1; i < bin.size(); i++) {
            if (left[1] == 0) { // we are still searching for the right edge of the first peak
                if (left[0] == 0) { // we are still searching for the left edge of the first peak
                    if (bin[i-1] < background) { // if the last bin was background noise
                        if (bin[i] - step_req > bin[i-1]) { // if there is a large jump in bin count
                            left[0] = plot.GetBinLowEdge(i); // we've found the left edge
                        }
                    }
                }
                else { // we have found the left edge of the first peak, but not the right
                    if (bin[i] < background) { // if this bin is background noise
                        left[1] = plot.GetBinLowEdge(i-1) + plot.GetBinWidth(i-1); // we've found the right edge
                    }
                }
            } else { // we have found the first peak, and are searching for the second
                if (right == 0) { // we are still searching for the right edge
                    if (bin[i] < background) { // if this bin is background noise
                        if (bin[i] + step_req < bin[i-1]) { // if there is a large jump in bin count
                            right = plot.GetBinLowEdge(i);
                        }
                    }
                }
            }
        }

        int offset = right - left[0]; // everything between left[0] and left[1] will be offset by this value
        if (right == 0) {
            offset = 0;
            std::cout << "    " << name << ": Couldn't find any split peaks." << endl;
        } else {
            std::cout << "    " << name << ": Found a split peak. Attempting to fix it." << endl;
        }
        // std::cout << "    mean: " << mean << endl; 
        // cout << "offset: " << offset << ", right edge: " << right << ", left edge: " << left[0] << endl;

        // show "before" plot
        if (plot::repair_peaks) {
            plot.GetXaxis()->SetTitle(name.c_str());
            plot.GetYaxis()->SetTitle("Count");
            plot.GetXaxis()->CenterTitle();
            plot.GetYaxis()->CenterTitle();
            plot.SetNdivisions(205, "X");
            plot.SetNdivisions(205, "Y");
            plot.Draw();

            canvas->SetLeftMargin(0.15);
            canvas->Modified(); canvas->Update();
            if (plot::save) {
                string file = plot::path + filename + "_before" + plot::format;
                canvas->SaveAs(file.c_str());
            }
            canvas->WaitPrimitive();
        }

        plot.Reset("ICESM"); // remove all data already in the plot

        if (offset != 0) {
            for (int i = 0; i < t.size(); i++) {
                if (left[0] < t[i] && t[i] < left[1]) {
                    t[i] += offset;
                }
                plot.Fill(t[i]);
            }

            // show "after" plot
            if (plot::repair_peaks) {
                plot.GetXaxis()->SetTitle(name.c_str());
                plot.GetYaxis()->SetTitle("Count");
                plot.GetXaxis()->CenterTitle();
                plot.GetYaxis()->CenterTitle();
                plot.SetNdivisions(205, "X");
                plot.SetNdivisions(205, "Y");
                plot.Draw();

                canvas->SetLeftMargin(0.15);
                canvas->Modified(); canvas->Update();
                if (plot::save) {
                    string file = plot::path + filename + "_after" + plot::format;
                    canvas->SaveAs(file.c_str());
                }
                canvas->WaitPrimitive();
            }
        }

        vector<int> result = {left[0], left[1], offset};
        return result;
    };

    // applies the offset in offset[2] to the peak bounded by offset[0] and offset[1] in detector det
    auto apply_offset = [] (vector<double> *T, vector<int> *ID, vector<int> offset, int det) {
        vector<double> &t = *T;
        vector<int> &id = *ID;
        int a = offset[0], b = offset[1], shift = offset[2];
        for (int i = 0; i < t.size(); i++) {
            if (id[i] == det) { // detector filter
                if (a < t[i] && t[i] < b) { // check if entry is within the given area
                    t[i] += shift;
                }
            }
        }
        return;
    };

    for (int det = 0; det < 4; det++) {
        TH1D ft_plot("", "", int(plot::x_axes::FT[0]), plot::x_axes::FT[1], plot::x_axes::FT[2]);
        TH1D bt_plot("", "", int(plot::x_axes::BT[0]), plot::x_axes::BT[1], plot::x_axes::BT[2]);
        vector<double> ft = detector_filter(FT, ID, det);
        vector<double> bt = detector_filter(BT, ID, det);

        for (int i = 0; i < ft.size(); i++) {
            ft_plot.Fill(ft[i]);
            bt_plot.Fill(bt[i]);
        }

        // extract histogram information
        vector<int> ft_data(ft_plot.GetNbinsX()-1);
        vector<int> bt_data(bt_plot.GetNbinsX()-1);
        for (int i = 1; i < ft_data.size(); i++) { // the first bin is the underflow bin, not sure if overflow bin is also included here
            ft_data[i] = ft_plot.GetBinContent(i);
            bt_data[i] = bt_plot.GetBinContent(i);
        }

        std::cout << "\033[1;34m" << "\nLooking for split peaks in detector " << detector_map.at(det) << "." << "\033[0m" << endl; 
        vector<int> offset;

        // check FT
        if (det == 2 || det == 3) { // FT is only well-defined for the W1 detectors
            offset = check_col(&ft, &ft_data, &ft_plot, "FT", (format("repair_peaks_%1%_FT") % detector_map.at(det)).str());
            if (offset[2] != 0) { // offset[2] is set to 0 if no peak was found
                apply_offset(FT, ID, offset, det);
            }
        }
        
        // check BT
        offset = check_col(&bt, &bt_data, &bt_plot, "BT", (format("repair_peaks_%1%_BT") % detector_map.at(det)).str()); // always correct BT
        if (offset[2] != 0) { // offset[2] is set to 0 if no peak was found
            apply_offset(BT, ID, offset, det);
        }
        ft_plot.Delete();
        bt_plot.Delete();
    }
    canvas->Close();
}

void align_peaks(data_container* data, vector<vector<int>> ft_peaks, vector<vector<int>> bt_peaks) {
    print_title("### Attempting to align peaks from different detectors ###");
    vector<double>& ft = *data->get_double("FT");
    vector<double>& bt = *data->get_double("BT");
    vector<int>& id = *data->get_int("ID");

    TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);
    TH1D ft_plot("", "", int(plot::x_axes::FT[0]), plot::x_axes::FT[1], plot::x_axes::FT[2]);
    TH1D bt_plot("", "", int(plot::x_axes::BT[0]), plot::x_axes::BT[1], plot::x_axes::BT[2]);

    for (int i = 0; i < ft.size(); i++) {
        ft_plot.Fill(ft[i]);
        bt_plot.Fill(bt[i]);
    }

    auto multi_gauss_fit = [&canvas] (TH1* hist, vector<vector<int>> peaks, TF1* total, string name) {
        // before we can fit multiple gauss simultaneously, we fit them indiviually within their ranges to get some nice initial parameter estimates
        vector<TF1*> gauss(peaks.size()); 
        for (int i = 0; i < peaks.size(); i++) {
            gauss[i] = new TF1((format("m%1%") % i).str().c_str(), "gaus", peaks[i][0], peaks[i][1]);
        }

        // perform the actual fitting
        double params[3*peaks.size()];
        hist->Fit(gauss[0], "QR"); // fit the first gauss only within the range it is defined on
        gauss[0]->GetParameters(&params[0]); // extract the parameters
        for (int i = 1; i < peaks.size(); i++) {
            hist->Fit(gauss[i], "QR+"); // R+ tells ROOT that we want to fit multiple functions (otherwise the previous one is simply deleted)
            gauss[i]->GetParameters(&params[3*i]); // extract the parameters
        }

        if (plot::align_peaks) {
            std::cout << "Results of the individual Gauss fits is shown on the figure." << endl;
            hist->GetXaxis()->SetTitle(name.c_str());
            hist->GetXaxis()->CenterTitle();
            hist->GetXaxis()->SetLabelOffset(0.005);
            hist->GetYaxis()->SetTitle("Count");
            hist->GetYaxis()->CenterTitle();
            hist->SetNdivisions(205, "X");
            hist->SetNdivisions(205, "Y");
            hist->Draw();

            canvas->SetLeftMargin(0.15);
            canvas->Modified(); canvas->Update();
            if (plot::save) {
                string file = plot::path + "align_peaks_" + name + "_individual" + plot::format;
                canvas->SaveAs(file.c_str());
            }
            canvas->WaitPrimitive();
            canvas->Clear();
        }

        total->SetParameters(params);
        hist->Fit(total, "Q");

        if (plot::align_peaks) {
            std::cout << "Results of the total Gauss fit is shown on the figure." << endl;
            hist->GetXaxis()->SetTitle(name.c_str());
            hist->GetXaxis()->CenterTitle();
            hist->GetYaxis()->SetTitle("Count");
            hist->GetYaxis()->CenterTitle();
            hist->SetNdivisions(205, "X");
            hist->SetNdivisions(205, "Y");
            hist->Draw();

            canvas->SetLeftMargin(0.15);
            canvas->Modified(); canvas->Update();
            if (plot::save) {
                string file = plot::path + "align_peaks_" + name + "_simultaneous" + plot::format;
                canvas->SaveAs(file.c_str());
            }
            canvas->WaitPrimitive();
            canvas->Clear();
        }
    };

    auto perform_shift = [&canvas, &data] (vector<double>* t, TH1* hist, TF1* total, int n_sigma, string name) {
        int n = total->GetNumberFreeParameters()/3; // for looping over all means
        // cout << "Number of free parameters: " << n << endl;
        double result[3*n]; // allocate space for the result vector
        vector<double> amp(n), mean(n), sigma(n);
        total->GetParameters(&result[0]); // extract the parameters
        for (int i = 0; i < n; i++) { // put them into sensible containers
            amp[i] = result[3*i];
            mean[i] = result[3*i + 1];
            sigma[i] = result[3*i + 2];
        }
        
        vector<double> &tref = *t; // dereference the data
        TH1 &href = *hist;

        // determine the edges of all peaks
        vector<vector<int>> edges(n, vector<int>(2));
        vector<int> offsets(n);
        for (int i = 0; i < n; i++) {
            edges[i][0] = mean[i] - n_sigma*sigma[i];
            edges[i][1] = mean[i] + n_sigma*sigma[i];
            offsets[i] = mean[i] - mean[0]; // we shift everything on top of the first peak
            // cout << "Borders: (" << edges[i][0] << ", " << edges[i][1] << "), offset: " << offsets[i] << endl;
        }

        // perform the actual shifts
        href.Reset("ICESM");
        for (int i = 0; i < tref.size(); i++) {
            int e = tref[i]; // element
            for (int j = 0; j < n; j++) {
                if (edges[j][0] < e && e < edges[j][1]) {
                    tref[i] -= offsets[j];
                    break;
                }
            }
            href.Fill(tref[i]);
        }

        if (plot::align_peaks) {
            hist->GetXaxis()->SetTitle(name.c_str());
            hist->GetXaxis()->CenterTitle();
            hist->GetYaxis()->SetTitle("Count");
            hist->GetYaxis()->CenterTitle();
            hist->SetNdivisions(205, "X");
            hist->SetNdivisions(205, "Y");
            href.Draw();

            canvas->SetLeftMargin(0.15);
            canvas->Modified(); canvas->Update();
            if (plot::save) {
                string file = plot::path + "align_peaks_" + name + "_aligned" + plot::format;
                canvas->SaveAs(file.c_str());
            }
            canvas->WaitPrimitive();
            canvas->Clear();
        }
    };
    //### FT ###//
    // "arg" must be a string describing how many gauss functions we expect, in a specific notation 
    // if we are fitting 5 gaussians, it should be "gaus(0) + gaus(3) + gaus(6) + gaus(9) + gaus(12)"
    string arg = "gaus(0)";
    for (int i = 1; i < ft_peaks.size(); i++) {
        arg += (format(" + gaus(%1%)") % (3*i)).str();
    }
    TF1* total = new TF1("total", arg.c_str(), plot::x_axes::FT[1], plot::x_axes::FT[2]); // the multi-gauss fit function
    multi_gauss_fit(&ft_plot, ft_peaks, total, "FT");
    perform_shift(&ft, &ft_plot, total, 3, "FT");

    //### BT ###//
    arg = "gaus(0)";
    for (int i = 1; i < 6; i++) {
        arg += (format(" + gaus(%1%)") % (3*i)).str();
    } 
    total = new TF1("total", arg.c_str(), plot::x_axes::BT[1], plot::x_axes::BT[2]); // the multi-gauss fit function
    multi_gauss_fit(&bt_plot, bt_peaks, total, "BT");
    perform_shift(&bt, &bt_plot, total, 3, "BT");
}

// perform a gaussian fit to both FT and BT and remove outliers beyond n_sigma
void gauss_filter(data_container* data, int n_sigma, vector<double> ft_peak, vector<double> bt_peak) {
    print_title("### Removing outliers ###");
    vector<double>& ft = *data->get_double("FT");
    vector<double>& bt = *data->get_double("BT");
    vector<int>& id = *data->get_int("ID");

    map<string, double> ft_result = gauss_fit(&ft, ft_peak, plot::gauss_cut);
    map<string, double> bt_result = gauss_fit(&bt, bt_peak, plot::gauss_cut);

    vector<vector<double>> peak = {{ft_result.at("mean") - n_sigma*ft_result.at("sigma"), ft_result.at("mean") + n_sigma*ft_result.at("sigma")}, 
                                   {bt_result.at("mean") - n_sigma*bt_result.at("sigma"), bt_result.at("mean") + n_sigma*bt_result.at("sigma")}};

    int n = ft.size();
    int m = n/3;
    int c = 0; // counter for output information
    vector<bool> filter(n, false);
    for (int i1 = 0; i1 < m; i1++) {
        int i2 = i1 + m; 
        int i3 = i2 + m;

        // FT check, only applicable to the W1 detectors
        if ((id[i1] == 2 || id[i1] == 3) && !(peak[0][0] < ft[i1] && ft[i1] < peak[0][1])) continue;
        if ((id[i2] == 2 || id[i2] == 3) && !(peak[0][0] < ft[i2] && ft[i2] < peak[0][1])) continue;
        if ((id[i3] == 2 || id[i3] == 3) && !(peak[0][0] < ft[i3] && ft[i3] < peak[0][1])) continue;

        // BT check
        if (!(peak[1][0] < bt[i1] && bt[i1] < peak[1][1])) continue;
        if (!(peak[1][0] < bt[i2] && bt[i2] < peak[1][1])) continue;
        if (!(peak[1][0] < bt[i3] && bt[i3] < peak[1][1])) continue;

        filter[i1] = filter[i2] = filter[i3] = true;
        c++;
    }

    if (plot::gauss_cut) {
        std::cout << format("Imposing a cut at %1% sigma, keeping %2% of %3% elements.") % n_sigma % (3*c) % n << endl; 
        hist_info info_ft(plot::x_axes::FT, "gauss_cut_ft", "title", "ft", "count");
        hist1D(&ft, ft_result.at("mean"), n_sigma*ft_result.at("sigma"), info_ft);

        hist_info info_bt(plot::x_axes::BT, "gauss_cut_bt", "title", "bt", "count");
        hist1D(&bt, bt_result.at("mean"), n_sigma*bt_result.at("sigma"), info_bt);
    }
    data->filter(&filter);
}

// perform a gaussian fit to both FT and BT and remove outliers beyond n_sigma
// note this is almost the same as the gauss_filter above, except we only check the non-zero BT/FT alphas here
void gauss_filter_custom(data_container* data, int n_sigma, vector<double> ft_peak, vector<double> bt_peak) {
    print_title("### Removing outliers ###");
    vector<double>& ft = *data->get_double("FT");
    vector<double>& bt = *data->get_double("BT");
    vector<int>& id = *data->get_int("ID");

    map<string, double> ft_result = gauss_fit(&ft, ft_peak, plot::gauss_cut);
    map<string, double> bt_result = gauss_fit(&bt, bt_peak, plot::gauss_cut);

    vector<vector<double>> peak = {{ft_result.at("mean") - n_sigma*ft_result.at("sigma"), ft_result.at("mean") + n_sigma*ft_result.at("sigma")}, 
                                   {bt_result.at("mean") - n_sigma*bt_result.at("sigma"), bt_result.at("mean") + n_sigma*bt_result.at("sigma")}};

    int n = ft.size();
    int m = n/3;
    int c = 0; // counter for output information
    vector<bool> filter(n, false);
    for (int i1 = 0; i1 < m; i1++) {
        int i2 = i1 + m; 
        int i3 = i2 + m;

        // FT check, only applicable to the W1 detectors
        if ((ft[i1] != 0) && (id[i1] == 2 || id[i1] == 3) && !(peak[0][0] < ft[i1] && ft[i1] < peak[0][1])) continue;
        if ((ft[i2] != 0) && (id[i2] == 2 || id[i2] == 3) && !(peak[0][0] < ft[i2] && ft[i2] < peak[0][1])) continue;
        if ((ft[i3] != 0) && (id[i3] == 2 || id[i3] == 3) && !(peak[0][0] < ft[i3] && ft[i3] < peak[0][1])) continue;

        // BT check
        if ((bt[i1] != 0) && !(peak[1][0] < bt[i1] && bt[i1] < peak[1][1])) continue;
        if ((bt[i2] != 0) && !(peak[1][0] < bt[i2] && bt[i2] < peak[1][1])) continue;
        if ((bt[i3] != 0) && !(peak[1][0] < bt[i3] && bt[i3] < peak[1][1])) continue;

        filter[i1] = filter[i2] = filter[i3] = true;
        c++;
    }

    if (plot::gauss_cut) {
        std::cout << format("Imposing a cut at %1% sigma, keeping %2% of %3% elements.") % n_sigma % (3*c) % n << endl; 
        hist_info info_ft(plot::x_axes::FT, "gauss_cut_ft", "title", "ft", "count");
        hist1D(&ft, ft_result.at("mean"), n_sigma*ft_result.at("sigma"), info_ft);

        hist_info info_bt(plot::x_axes::BT, "gauss_cut_bt", "title", "bt", "count");
        hist1D(&bt, bt_result.at("mean"), n_sigma*bt_result.at("sigma"), info_bt);

        auto dt = calc_dt(&ft, &bt);
        hist_info info_dt(plot::x_axes::standard, "gauss_cut_dt", "title", "dt", "count");
        hist1D(&dt, info_dt);

    }
    data->filter(&filter);
}

// int main(int argc, char *argv[]) {
//     data_container data;
//     prepare_data(argc, argv, &data, "");
//     gStyle->SetOptStat(0);

//     TApplication *app = new TApplication("ROOT window", 0, 0);
//     repair_peaks(&data, false);
//     align_peaks(&data, true);
//     debug(&data);
//     return 0;
// }