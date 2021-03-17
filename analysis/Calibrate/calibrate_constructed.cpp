// root stuff
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>

// other stuff
#include <filesystem>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>

// my own stuff
#include "plots.cpp"
#include "utility.cpp"
#include "data_fix.cpp"

using namespace std;
using boost::format;

// perform a gaussian fit to both FT and BT and remove outliers beyond n_sigma
// note this is almost the same as gauss_filter from data_fix.cpp, except we only check the non-zero BT/FT alphas here
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
    }
    data->filter(&filter);
}

// load a fit calibration file
tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> load_calibration() {
    print_title("*** LOADING CALIBRATION RESULTS ***");
    vector<vector<double>> doubleF;
    vector<vector<double>> doubleB;
    vector<vector<double>> offset;
    int c = 0; // counter in current vector

    string line;
    string currently_reading;
    ifstream file(plot::path + "../fit_info.txt");
    if (file.is_open()) {
        while(getline(file, line)) {
            if (line == "doubleF") {
                c = 0;
                currently_reading = "doubleF";
                continue;
            } else if (line == "doubleB") {
                c = 0;
                currently_reading = "doubleB";
                continue;
            } else if (line == "offset") {
                c = 0;
                currently_reading = "offset";
                continue;
            } else if (line == "") {
                continue;
            }

            vector<string> vals;
            if (currently_reading == "doubleF") {
                boost::split(vals, line, boost::is_any_of(" "));
                doubleF.push_back(vector<double>(vals.size()));
                for (int i = 0; i < vals.size(); i++) {
                    doubleF[c][i] = atof(vals[i].c_str());
                }
                c++;
            } else if (currently_reading == "doubleB") {
                boost::split(vals, line, boost::is_any_of(" "));
                doubleB.push_back(vector<double>(vals.size()));
                for (int i = 0; i < vals.size(); i++) {
                    doubleB[c][i] = atof(vals[i].c_str());
                }
                c++;
            } else if (currently_reading == "offset") {
                boost::split(vals, line, boost::is_any_of(" "));
                offset.push_back(vector<double>(vals.size()));
                for (int i = 0; i < vals.size(); i++) {
                    offset[c][i] = atof(vals[i].c_str());
                }
                c++;
            }
        }
    }

    // verify that we have read something
    if (currently_reading != "offset") {
        std::cout << "\033[1;31m" << "Something went wrong while loading the calibration file. Is it really located at " << plot::path + "../fit_info.txt ?"  << "\033[0m" << endl;
        exit(1);
    }

    std::cout << "Calibration loaded." << endl;
    return tuple(doubleF, doubleB, offset);
}

// calculate the mean, ignoring all events outside the x_axis, as well as those that are 0
vector<double> calc_offset_custom(const vector<int> *BI, const vector<double> *DT, const vector<double> x_axis, int strips) {
    const vector<int> &bi = *BI;
    const vector<double> &dt = *DT;
    vector<double> sum(strips);
    vector<int> count(strips);

    // sum over each of the strips
    int c = 0;
    for (int i = 0; i < bi.size(); i++) {
        if (dt[i] != 0 && x_axis[1] < dt[i] && dt[i] < x_axis[2]) {
            sum[bi[i]-1] += dt[i];
            count[bi[i]-1]++;
        }
    }
    
    // divide each sum by its number of elements (i.e. calculate the mean)
    for (int i = 0; i < strips; i++) {
        if (count[i] != 0) {
            sum[i] /= count[i];
        } else {
            sum[i] = 0;
        }
    }
    return sum;
}

// takes the mean value of each strip and subtracts it from dt and bt. ignores any event outside the x_axis, as well as any for which dt == 0
vector<double> center_dt_custom(vector<double>* DT, vector<double>* BT, const vector<int>* BI, const vector<double> x_axis, int strips) {
    vector<double>& dt = *DT;
    vector<double>& bt = *BT;

    vector<double> offset = calc_offset_custom(BI, DT, x_axis, strips);
    apply_offset(DT, BI, offset);
    apply_offset(BT, BI, offset);
    return offset;
}

// applies the offsets from all detectors to bt. ignores any event for which id[i] == -1 (reconstructed event)
void apply_offset_custom(data_container* data, const vector<vector<double>> offset) {
    data->add_dt();
    const vector<int> &bi = *data->get_int("BI");
    const vector<int> &id = *data->get_int("ID");
    vector<double> &bt = *data->get_double("BT");
    
    for (int i = 0; i < bt.size(); i++) {
        if (id[i] != -1) {
            bt[i] -= offset[id[i]][bi[i]-1];
        }
    }
    return;
}

// applies deltaF and deltaB to FT and BT according to the fitted equation
// this method applies the fit on data from all detectors simultaneously
// ignores any reconstructed event (id[i] == -1)
void apply_fit_custom(data_container* data, const vector<vector<double>> deltaF, const vector<vector<double>> deltaB) {
    vector<double> &ft = *data->get_double("FT");
    vector<double> &bt = *data->get_double("BT");
    const vector<int> &fi = *data->get_int("FI");
    const vector<int> &bi = *data->get_int("BI");
    const vector<int> &id = *data->get_int("ID");

    for (int i = 0; i < ft.size(); i++) {
        if (id[i] != -1) {
            ft[i] += deltaF[id[i]][fi[i]-1];
            bt[i] += deltaB[id[i]][bi[i]-1];
        }
    }
    return;
}

// apply a given tdc calibration to the data
void apply_tdc_calibration(data_container* data, vector<vector<double>> deltaF, vector<vector<double>> deltaB) {
    print_title("*** APPLYING EXISTING TDC CALIBRATION ***");
    // dereference all the necessary data from the container
    vector<double>* FT = data->get_double("FT");
    vector<double>* BT = data->get_double("BT");
    vector<int>* FI = data->get_int("FI");
    vector<int>* BI = data->get_int("BI");
    vector<int>* ID = data->get_int("ID");

    // since the peaks are shifted slightly (not all peaks could be found in align_peaks), we must calculate the offset again
    vector<vector<double>> offset(4);

    // the purpose of this loop is to calculate the offsets of each detector
    vector<int> strips = {24, 24, 16, 16}; // number of strips on each detector
    for (int det = 0; det <= 3; det++) {
        cout << "\033[1;34m" << "\nStarting calibration of detector " << detector_map.at(det) << "." << "\033[0m" << endl; 
        vector<int> fi, bi, id;
        vector<double> ft, bt;
        if (det == 0 || det == 1) { // S3
            tie(ft, bt, fi, bi, id) = find_coincidences(FT, BT, FI, BI, ID, det); 
        } else { // W1
            tie(ft, bt, fi, bi, id) = detector_filter(FT, BT, FI, BI, ID, det); 
        }

        // plot the raw data
        vector<double> dt = calc_dt(&ft, &bt);
        if (plot::raw) {
            hist_info info(plot::x_axes::standard, plot::y_axis.at(det), (format("tdc_calibration_constructed_raw_%1%") % detector_map.at(det)).str(), "Raw data", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, info);
        }

        // center the columns on 0
        offset[det] = center_dt_custom(&dt, &bt, &bi, plot::x_axes::standard, strips[det]);
        if (plot::centered) {
            hist_info info(plot::x_axes::centered, plot::y_axis.at(det), (format("tdc_calibration_constructed_centered_%1%") % detector_map.at(det)).str(), "Raw data", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, info);
        }
        if (plot::fit) {
            // apply the fit on the local data
            apply_fit(&ft, &bt, &fi, &bi, deltaF[det], deltaB[det]);
            hist_info info(plot::x_axes::centered, plot::y_axis.at(det), (format("tdc_calibration_fit_%1%") % detector_map.at(det)).str(), "Raw data", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, info);
        }
    }

    apply_offset_custom(data, offset);
    apply_fit_custom(data, deltaF, deltaB);
    return;
}

// fits a Gaussian to each dt histogram, and cuts it at a given x*sigma. this way we remove the outliers in a statistically sound manner
void gauss_cut_custom(data_container* data, int n_sigma) {
    print_title("*** GAUSS FITTING INITIALIZED ***");
    
    // the limits for the dt cut. these will be n_sigma*sigma from the Gauss fits
    vector<double> a(4);
    vector<double> b(4);
    vector<double> sigma(4); // can be deleted later

    // dereference all the necessary data from the container
    vector<double>* FT = data->get_double("FT");
    vector<double>* BT = data->get_double("BT");
    vector<int>* FI = data->get_int("FI");
    vector<int>* BI = data->get_int("BI");
    vector<int>* ID = data->get_int("ID");

    // define local vectors that can freely be modified
    vector<double> ft, bt;
    vector<int> fi, bi, id;
    for (int det = 0; det <= 3; det++) {
        cout << "\033[1;34m" << "\nPerforming Gauss fitting on detector " << detector_map.at(det) << "." << "\033[0m" << endl; 
        if (det == 0 || det == 1) { // S3
            tie(ft, bt, fi, bi, id) = find_coincidences(FT, BT, FI, BI, ID, det); 
        } else { // W1
            tie(ft, bt, fi, bi, id) = detector_filter(FT, BT, FI, BI, ID, det); 
        }
        vector<double> dt = calc_dt(&ft, &bt);

        // with all the data set up, we can perform the actual fits
        map<string, double> gauss = gauss_fit(&dt, plot::x_axes::gauss, plot::plot_det.at(det) && plot::gauss);
        a[det] = -n_sigma*gauss["sigma"];
        b[det] = n_sigma*gauss["sigma"];
        sigma[det] = gauss["sigma"];
        if (plot::W && plot::gauss_sigma) {
            hist_info info(plot::x_axes::centered, plot::y_axis.at(det), (format("gauss_cut_before_%1%") % detector_map.at(det)).str(), "Gauss fit to dt", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, gauss["sigma"], n_sigma, info);
        }
    }

    data->dt_filter(a, b);

    if (plot::gauss) {
        std::cout << "\nThe following plots shows the result of the Gaussian dt filter." << endl;
        for (int det = 0; det <= 3; det++) {
            vector<double> ft, bt;
            vector<int> fi, bi, id;
            if (det == 0 || det == 1) { // S3
                tie(ft, bt, fi, bi, id) = find_coincidences(FT, BT, FI, BI, ID, det); 
            } else { // W1
                tie(ft, bt, fi, bi, id) = detector_filter(FT, BT, FI, BI, ID, det); 
            }

            vector<double> dt = calc_dt(&ft, &bt);

            // std::cout << format("Plotting filtered dt histogram for %1%.") % detector_map.at(det) << endl;
            // hist_info info(plot::x_axes::gauss, plot::y_axis.at(det), "Raw data", "dt", "Count");
            // hist1D(&dt, info);
            
            std::cout << format("\nPlotting filtered 2D histogram for %1%.") % detector_map.at(det) << endl;
            hist_info info = hist_info(plot::x_axes::centered, plot::y_axis.at(det), (format("gauss_cut_after_%1%") % detector_map.at(det)).str(), "Gauss fit to dt", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, sigma[det], n_sigma, info);
        }
    }
}

/* 
    The main changes in this analysis compared to just calibrate.cpp is that all functions must respect that a whole lot of entries in FT and BT may be 0.
    To allow for this, I've had to make a few minor changes in most functions, all of which I've only made available locally in this file
*/
// int main(int argc, char *argv[]) {
//     data_container data;
//     prepare_data(argc, argv, &data, "mul==-1");

//     // define the plot colour scheme
//     gStyle->SetPalette(kBird);
//     gStyle->SetOptStat(0);
//     gStyle->SetOptTitle(0);
//     gROOT->ForceStyle();

//     // start a ROOT application window such that the plots can actually be shown
//     TApplication *app = new TApplication("ROOT window", 0, 0);
//     // save(&data, "output/raw_data.root");

//     // attempt to repair the broken peaks
//     repair_peaks(&data);

//     // attempt to align the peaks since each detector is offset slightly from the others
//     vector<vector<int>> ft_peaks = {{65400, 65500}, {65500, 65600}};
//     vector<vector<int>> bt_peaks = {{65300, 65380}, {65380, 65460}, {65460, 65520}, {65520, 65600}};
//     align_peaks(&data, ft_peaks, bt_peaks);

//     // imposes a Gaussian filter on FT and BT, to remove outliers. 
//     vector<double> ft_peak = {100, 65400, 65500};
//     vector<double> bt_peak = {200, 65250, 65450};
//     gauss_filter_custom(&data, 3, ft_peak, bt_peak);

//     auto[deltaF, deltaB, offset] = load_calibration();

//     // perform a tdc calibration on the data
//     apply_tdc_calibration(&data, deltaF, deltaB);

//     // impose a Gaussian filter on DT, to remove outliers
//     gauss_cut_custom(&data, 3);
    
//     // save the data_container
//     save(&data, "output/reconstructed_data.root");
//     return 0;
// }