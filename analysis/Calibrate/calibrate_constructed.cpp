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
        hist_info info_ft({1000.0, 65000, 66000}, "gauss_cut_ft", "title", "ft", "count");
        hist1D(&ft, ft_result.at("mean"), n_sigma*ft_result.at("sigma"), info_ft);

        hist_info info_bt({1000.0, 65000, 66000}, "gauss_cut_bt", "title", "bt", "count");
        hist1D(&bt, bt_result.at("mean"), n_sigma*bt_result.at("sigma"), info_bt);

        auto dt = calc_dt(&ft, &bt);
        hist_info info_dt({1000.0, 65000, 66000}, "gauss_cut_dt", "title", "dt", "count");
        hist1D(&dt, info_dt);
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

    std::cout << "Calibration loaded." << endl;
    return tuple(doubleF, doubleB, offset);
}

// applies the offsets from all detectors to bt
void apply_offset(data_container* data, const vector<vector<double>> offset) {
    data->add_dt();
    const vector<int> &bi = *data->get_int("BI");
    const vector<int> &id = *data->get_int("ID");
    vector<double> &bt = *data->get_double("BT");
    
    for (int i = 0; i < bt.size(); i++) {
        bt[i] -= offset[id[i]][bi[i]-1];
    }
    return;
}

// applies deltaF and deltaB to FT and BT according to the fitted equation
// this method applies the fit on data from all detectors simultaneously
void apply_fit(data_container* data, const vector<vector<double>> deltaF, const vector<vector<double>> deltaB) {
    vector<double> &FTref = *data->get_double("FT");
    vector<double> &BTref = *data->get_double("BT");
    const vector<int> &FIref = *data->get_int("FI");
    const vector<int> &BIref = *data->get_int("BI");
    const vector<int> &IDref = *data->get_int("ID");

    for (int i = 0; i < FTref.size(); i++) {
        FTref[i] += deltaF[IDref[i]][FIref[i]-1];
        BTref[i] += deltaB[IDref[i]][BIref[i]-1];
    }
    return;
}

// apply a given tdc calibration to the data
tuple<vector<bool>, vector<bool>> apply_tdc_calibration(data_container* data, vector<vector<double>> deltaF, vector<vector<double>> deltaB, vector<vector<double>> offset) {
    // dereference all the necessary data from the container
    vector<double>& ft = *data->get_double("FT");
    vector<double>& bt = *data->get_double("BT");
    vector<int>& fi = *data->get_int("FI");
    vector<int>& bi = *data->get_int("BI");
    vector<int>& id = *data->get_int("ID");
    
    // we need to keep track of which entries were zero
    vector<bool> zero_ft(ft.size());
    vector<bool> zero_bt(bt.size());
    for (int i = 0; i < ft.size(); i++) {
        zero_ft[i] = (ft[i] == 0);
        zero_bt[i] = (bt[i] == 0);
    }

    apply_offset(data, offset);
    apply_fit(data, deltaF, deltaB);
    return tuple(zero_ft, zero_bt);
}

// fits a Gaussian to each dt histogram, and cuts it at a given x*sigma. this way we remove the outliers in a statistically sound manner
void gauss_cut_custom(data_container* data, int n_sigma, vector<bool>* zero_ft, vector<bool>* zero_bt) {
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
    vector<double> ft, bt, dE;
    vector<int> fi, bi, id;
    for (int det = 0; det <= 3; det++) {
        cout << "\033[1;34m" << "\nPerforming Gauss fitting on detector " << detector_map.at(det) << "." << "\033[0m" << endl; 
        if (det == 0 || det == 1) { // S3
            tie(ft, bt, fi, bi, id) = find_coincidences(FT, BT, FI, BI, ID, det); 
            remove_0(&ft, &bt, &fi, &bi, &id, det);
        } else { // W1
            tie(ft, bt, fi, bi, id) = detector_filter(FT, BT, FI, BI, ID, det); 
            remove_0(&ft, &bt, &fi, &bi, &id, det);
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

    data->dt_filter(a, b, true);

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

int main(int argc, char *argv[]) {
    data_container data;
    prepare_data(argc, argv, &data, "mul==-1");

    // define the plot colour scheme
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gROOT->ForceStyle();

    // start a ROOT application window such that the plots can actually be shown
    TApplication *app = new TApplication("ROOT window", 0, 0);
    // save(&data, "output/raw_data.root");

    // attempt to repair the broken peaks
    repair_peaks(&data);

    // attempt to align the peaks since each detector is offset slightly from the others
    vector<vector<int>> ft_peaks = {{65400, 65500}, {65500, 65600}};
    vector<vector<int>> bt_peaks = {{65300, 65380}, {65380, 65460}, {65460, 65520}, {65520, 65600}};
    align_peaks(&data, ft_peaks, bt_peaks);

    // imposes a Gaussian filter on FT and BT, to remove outliers. 
    vector<double> ft_peak = {100, 65400, 65500};
    vector<double> bt_peak = {200, 65250, 65450};
    gauss_filter_custom(&data, 3, ft_peak, bt_peak);

    auto[deltaF, deltaB, offset] = load_calibration();

    // perform a tdc calibration on the data
    auto[zero_ft, zero_bt] = apply_tdc_calibration(&data, deltaF, deltaB, offset);

    // impose a Gaussian filter on DT, to remove outliers
    gauss_cut_custom(&data, 3, &zero_ft, &zero_bt);
    
    // save the data_container    
    // save(&data, "output/corrected_data.root"); // save the data after the above methods have been imposed upon it
    return 0;
}