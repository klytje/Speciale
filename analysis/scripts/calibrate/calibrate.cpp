// optimization stuff
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

// other root stuff
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
#include <iostream>
#include <fstream>

// my own stuff
#include "plots.cpp"
#include "utility.cpp"
#include "data_fix.cpp"

using namespace std;
using boost::format;

// if your data is very concentrated on a few strips, you may need to increase this value
static int mem_scaler = 1;

// filter the input to elements satisfying a < dt < b
template <typename T>
vector<T> dt_filter(const vector<T> *input, const vector<double> *DT, double a, double b) {
    int n = DT->size();
    const vector<double> &DTref = *DT;
    const vector<T> &inputref = *input;
    vector<T> result(n);
    int m = 0; // current index in result
    double dt;
    for (int i = 0; i < n; i++) {
        dt = DTref[i];
        if (a < dt && dt < b) {
            result[m] = inputref[i];
            m++;
        }
    }
    result.resize(m); // remove all the unset entries
    return result;
}

// a variant of the dt filter which filters FT, BT, FI, BI, ID simultaneously
tuple<vector<double>, vector<double>, vector<int>, vector<int>, vector<int>> dt_filter(const vector<double> *FT, const vector<double> *BT, const vector<double> *DT, const vector<int> *FI, const vector<int> *BI, const vector<int> *ID, double a, double b) {
    // dereference input
    int n = ID->size(); 
    const vector<double> &FTref = *FT;
    const vector<double> &BTref = *BT;
    const vector<double> &DTref = *DT;
    const vector<int> &FIref = *FI;
    const vector<int> &BIref = *BI;
    const vector<int> &IDref = *ID; 

    // initialize result vectors
    vector<double> ft(n);
    vector<double> bt(n);
    vector<int> fi(n);
    vector<int> bi(n);
    vector<int> id(n);

    int m = 0; // current index in results
    double dt;
    for (int i = 0; i < n; i++) {
        dt = DTref[i];
        if (a < dt && dt < b) {
            ft[m] = FTref[i];
            bt[m] = BTref[i];
            fi[m] = FIref[i];
            bi[m] = BIref[i];
            id[m] = IDref[i];
            m++;
        }
    }
    // remove all unset entries and return the tuple
    ft.resize(m); bt.resize(m); fi.resize(m); bi.resize(m); id.resize(m);
    return tuple(ft, bt, fi, bi, id);
}

// a class storing and handling everything related to the fit
class fit_func {
    public: 
        fit_func(vector<double> *FT, vector<double> *BT, vector<int> *FI, vector<int> *BI, int det) {
            this->FT = *FT;
            this->BT = *BT;
            this->FI = *FI;
            this->BI = *BI;
            this->det = det;
            this->size = FT->size();

            // set the number of strips
            if (det == 0 || det == 1) { // S3
                f_strips = 16; 
                b_strips = 24; 
            } else { // W1
                f_strips = 16; 
                b_strips = 16; 
            }

            // initialize the grid of strip ids
            ft = vector<vector<vector<double>>>(b_strips, vector<vector<double>>(f_strips));
            bt = vector<vector<vector<double>>>(b_strips, vector<vector<double>>(f_strips));

            // fill the grid with data
            for (int i = 1; i <= b_strips; i++) {
                separate_strip(i);
            }
            // debug_plot(6); // plot strip for debug
        }

        void set_fid(int fid) {
            this->fid = fid-1;
        }

        void set_bid(int bid) {
            this->bid = bid-1;
        }

        // fits all front strips from a single back strip
        // note that both fid and bid must be set before it can be used
        // note also that it assumes deltaB = 0, so it can only be used *once*
        double single_back_fit (const double *x) {
            double deltaF = x[0];
            double sum = 0;
            for (int i = 0; i < ft[bid][fid].size(); i++) {
                sum += abs(bt[bid][fid][i] - ft[bid][fid][i] - deltaF);
            }
            return sum;
        }

        // sets the result of the single_back_fit minimization
        // this *must* be done before a single_front_fit can be attempted
        void set_single_back_fit_result(vector<double> result) {
            this->single_back_fit_result = result;
        }

        // fits a back strip from all front strips
        // note that bid must be set before it can be used
        double single_front_fit (const double *x) {
            double deltaB = x[0];
            double sum = 0;

            for (int i = 0; i < f_strips; i++) { // iterate over front strips
                for (int j = 0; j < ft[bid][i].size(); j++) { // iterate over all entries
                    sum += abs(bt[bid][i][j] + deltaB - ft[bid][i][j] - single_back_fit_result[i]);
                }
            }
            // for (int i = 0; i < ft[bid][fid].size(); i++) {
            //     sum += abs(bt[bid][fid][i] + deltaB - ft[bid][fid][i] - single_back_fit_result[fid]);
            // }
            return sum;            
        }

        // simultaneous 32 variable fit x[0:16] = deltaB, x[16:32] = deltaF, x[16] = 0
        double fit (const double *x) {
            double sum = 0;
            for (int i = 0; i < b_strips; i++) { // loop over back strips
                for (int j = 0; j < f_strips; j++) { // loop over front strips
                    for (int k = 0; k < ft[i][j].size(); k++) { // loop over data
                        sum += abs(bt[i][j][k] + x[i] - ft[i][j][k] - x[16+j]);
                    }
                }
            }
            return sum;
        }

        // apply a fit result to the internal data
        void apply_fit(vector<double> deltaF, vector<double> deltaB) {
            for (int i = 0; i < b_strips; i++) { // loop over back strips
                for (int j = 0; j < f_strips; j++) { // loop over front strips
                    for (int k = 0; k < ft[i][j].size(); k++) { // loop over data
                        ft[i][j][k] += deltaF[j];
                        bt[i][j][k] += deltaB[i];
                    }
                }
            }
            debug_plot();
        }

    private:
        vector<double> &FT = *(new vector<double>()); // references to the actual data
        vector<double> &BT = *(new vector<double>());
        vector<int> &BI = *(new vector<int>());
        vector<int> &FI = *(new vector<int>());
        vector<vector<vector<double>>> ft; // local data split into each of the 16x16 fitting regions
        vector<vector<vector<double>>> bt; 
        vector<double> single_back_fit_result;
        int det;
        int size;
        int fid; // current front strip
        int bid; // current back strip
        int f_strips; // number of front strips (16 for W1, 16 for S3)
        int b_strips; // number of back strips (16 for W1, 24 for S3)

        // separates data from back strip "bid" into its "strips" FI components and stores the result in ft/bt
        void separate_strip(int bid) {
            vector<int> mf(f_strips); // current index of the rf vectors
            vector<int> mb(f_strips); // current index of the rb vectors
            vector<vector<double>> rf(f_strips, vector<double>(mem_scaler*size/f_strips)); // resulting ft[0]. size/strips is probably alright
            vector<vector<double>> rb(f_strips, vector<double>(mem_scaler*size/f_strips)); // resulting bt[0]

            int fid;
            for (int i = 0; i < size; i++) {
                if (BI[i] == bid) {
                    fid = FI[i]-1;
                    rf[fid][mf[fid]] = FT[i]; // mf[fid] is the current index of rf[fid]
                    rb[fid][mb[fid]] = BT[i]; // mb[fid] is the current index of rb[fid]
                    mf[fid]++;
                    mb[fid]++;
                }
            }
            
            // resize the vectors to free all unused space
            for (int i = 0; i < f_strips; i++) {
                rf[i].resize(mf[i]);            
                rb[i].resize(mb[i]);
            }

            // set the internally stored representation
            ft[bid-1] = rf;
            bt[bid-1] = rb;
        }

        // plots the data from "strip" bid (or everything if nothing is given)
        // used for debugging and nothing else
        void debug_plot(int strip = 0) {
            double x_axis[3] = {100, -50, 50}; // bins, xmin, xmax
            double y_axis[3] = {double(b_strips), 1, b_strips+1.0};    // bins, ymin, ymax
            TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);
            // TH1D *h = new TH1D("raw input histogram", "2d histo", int(x_axis[0]), x_axis[1], x_axis[2]); 
            TH2D *h = new TH2D("raw input histogram", "2d histo", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]); 
            
            int n = 0;
            for (int i = 0; i < f_strips; i++) {
                n += ft[0][i].size();
            }
            vector<double> flat_ft(n);
            vector<double> flat_bt(n);

            int start = strip-1;
            int stop = strip;
            if (stop == 0) {
                start = 0;
                stop = b_strips;
            }
            for (int i = start; i < stop; i++) {
                for (int j = 0; j < f_strips; j++) {
                    for (int k = 0; k < ft[i][j].size(); k++) {
                        h->Fill(bt[i][j][k] - ft[i][j][k], i+1);
                    }
                }
            }
            std::cout << "Debug histogram filled with " << h->GetEntries() << " data points." << endl;

            // set up the figure
            h->GetXaxis()->SetTitle("Time [ns]");
            h->GetYaxis()->SetTitle("Strip");
            h->GetXaxis()->CenterTitle();
            h->GetYaxis()->CenterTitle();
            h->Draw("colz");
            canvas->SetLogz();
            canvas->SetRightMargin(0.15);
            canvas->Modified(); canvas->Update();
            canvas->WaitPrimitive();
            canvas->Close();
            h->Delete();
        }
};

// the fitting function. it handles everything related to the fit process. 
pair<vector<double>, vector<double>> fit(vector<double> *FT, vector<double> *BT, vector<int> *FI, vector<int> *BI, vector<double> *DT, int det) {
        std::cout << "Performing a fit on the data." << endl;
        vector<double> &ft = *FT;
        vector<double> &bt = *BT;
        vector<double> &dt = *DT;
        vector<int> &fi = *FI;
        vector<int> &bi = *BI;

        // check for disabled strips. returns a boolean vector with an entry for each strip
        auto check_disabled = [] (vector<int> *BI, int nb) {
            vector<int> &bi = *BI;
            vector<int> counter(nb, 0); // contains the number of entries for each bid
            vector<bool> disabled(nb, false); // the result vector
            for (int i = 0; i < bi.size(); i++) {
                counter[bi[i]-1]++;
            }

            for (int i = 0; i < nb; i++) {
                if (counter[i] == 0) {
                    cout << "Strip " << i+1 << " appears to be disabled." << endl; 
                    disabled[i] = true;
                }
            }
            return disabled;
        };

        // determine the sizes of the fitting vectors
        int nf = 16;
        int nb, back_fit_start;
        if (det == 0 || det == 1) {
            nb = 24; // S3
            back_fit_start = 1; // the first back strip is not fixed for the S3's
        } else {
            nb = 16; // W1
            back_fit_start = 2; // we fix the first back strip to 0 for the W1's
        }

        // check for disabled strips
        vector<bool> disabled = check_disabled(BI, nb);

        // create the actual minimizer, see https://root.cern.ch/doc/v612/NumericalMinimization_8C.html for options
        ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetMaxFunctionCalls(100000); // maybe increase this
        double step = 0.01;
        double start = 0;

        fit_func f(FT, BT, FI, BI, det); // initialize the fitting class
        // pick the first back strip that's not disabled
        for (int i = 1; i < nb+1; i++) {
            if (!disabled[i-1]) {
                f.set_bid(i);
                break;
            }
        }
        vector<double> deltaB(nb);
        vector<double> deltaF(nf);

        // check if we have to do the back fit. we shouldn't do this for the S3 detectors
        if (det == 2 || det == 3) {
            // fit a single back strip containing 16 front strips to get an estimate of each deltaF[i]
            for (int i = 1; i < nf+1; i++) {
                f.set_fid(i); // set the fid we're trying to fit
                auto g = bind(&fit_func::single_back_fit, f, placeholders::_1);
                ROOT::Math::Functor func(g, 1);
                minimizer->SetFunction(func);
                minimizer->SetVariable(0, "deltaF", start, step);
                minimizer->Minimize();
                auto res = minimizer->X();
                deltaF[i-1] = res[0];
                // cout << "Front strip fit " << i << " with result " << res[0] << endl;
            }
        }

        // use the previous result to fit each individual back strip to get an estimate of deltaB[i]
        f.set_single_back_fit_result(deltaF);
        for (int i = back_fit_start; i < nb+1; i++) {
            if (disabled[i-1]) continue; // skip all disabled strips
            f.set_bid(i);
            auto g = bind(&fit_func::single_front_fit, f, placeholders::_1);
            ROOT::Math::Functor func(g, 1);
            minimizer->SetFunction(func);
            minimizer->SetVariable(0, "deltaB", start, step);
            minimizer->Minimize();
            auto res = minimizer->X();
            deltaB[i-1] = res[0];
            // cout << "Back strip fit " << i << " with result " << res[0] << endl;
        }

        //*** simultaneous 32-variable fit ***//
        // auto g = bind(&fit_func::fit, f, placeholders::_1);
        // ROOT::Math::Functor func(g, 1);
        // minimizer->SetFunction(func);

        // // set the optimization variables (two entries are pulled out of the loop so we can fix deltaB[0] = 0)
        // minimizer->SetFixedVariable(0, "deltaB_0", 0);
        // minimizer->SetVariable(16, "deltaF_0", deltaF[0], step);
        // for (int i = 1; i < 16; i++) {
        //     minimizer->SetVariable(i, (format("deltaB_%1%") % i).str(), deltaB[i], step);
        //     minimizer->SetVariable(16+i, (format("deltaF_%1%") % i).str(), deltaF[i], step);
        // }
        // minimizer->Minimize();
        // auto res = minimizer->X();
        // for (int i = 0; i < 16; i++) {
        //     deltaF[i] = res[i];
        //     deltaB[i] = res[16+i];
        //     std::cout << deltaF[i] << ", " << deltaB[i] << endl;
        // }

        // apply the fit and extract the data
        apply_fit(FT, BT, FI, BI, deltaF, deltaB);
        dt = calc_dt(FT, BT);
        return pair(deltaF, deltaB);
}

// performs a tdc-calibration of the data and returns the fit results
// note that it also centers the input around 0 and cuts everything outside the x_axis_raw bounds
tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> tdc_calibration(data_container* data) {
    print_title("*** TDC CALIBRATION INITIALIZED ***");

    // we need these four arrays to impose the fit on the actual data (input references) at the end
    vector<vector<double>> deltaF(4); // the delta_f fit parameters
    vector<vector<double>> deltaB(4); // the delta_b fit parameters
    vector<vector<double>> offset(4);

    // dereference all the necessary data from the container
    vector<double>* FT = data->get_double("FT");
    vector<double>* BT = data->get_double("BT");
    vector<int>* FI = data->get_int("FI");
    vector<int>* BI = data->get_int("BI");
    vector<int>* ID = data->get_int("ID");

    // W1 detectors
    for (int det = 2; det <= 3; det++) {
        cout << "\033[1;34m" << "\nStarting analysis of detector " << detector_map.at(det) << "." << "\033[0m" << endl; 
        auto[ft, bt, fi, bi, id] = detector_filter(FT, BT, FI, BI, ID, det);
        remove_0(&ft, &bt, &fi, &bi, &id, det);
        vector<double> dt = calc_dt(&ft, &bt);
        if (plot::W && plot::raw) {
            hist_info info(plot::x_axes::standard, plot::y_axis.at(det), (format("tdc_calibration_raw_%1%") % detector_map.at(det)).str(), "Raw data", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, info);
        }

        // center the columns on 0
        offset[det] = center_dt(&dt, &bt, &bi, 16);
        if (plot::W && plot::centered) {
            hist_info info(plot::x_axes::centered, plot::y_axis.at(det), (format("tdc_calibration_centered_%1%") % detector_map.at(det)).str(), "Raw data", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, info);
        }

        // perform the actual fit
        tie(deltaF[det], deltaB[det]) = fit(&ft, &bt, &fi, &bi, &dt, det);
        if (plot::W && plot::fit) {
            hist_info info(plot::x_axes::centered, plot::y_axis.at(det), (format("tdc_calibration_fit_%1%") % detector_map.at(det)).str(), "Raw data", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, info);
        }
    }

    // S3 detectors
    // only back side strips are connected to the TDC
    for (int det = 0; det <= 1; det++) {
        cout << "\033[1;34m" << "\nStarting analysis of detector " << detector_map.at(det) << "." << "\033[0m" << endl; 
        auto[ft, bt, fi, bi, id] = find_coincidences(FT, BT, FI, BI, ID, det);
        remove_0(&ft, &bt, &fi, &bi, &id, det);
        vector<double> dt = calc_dt(&ft, &bt);
        cout << "Found " << dt.size() << " multi-detector coincidences." << endl;
        if (plot::S && plot::raw) {
            hist_info info(plot::x_axes::standard, plot::y_axis.at(det), (format("tdc_calibration_raw_%1%") % detector_map.at(det)).str(), "Raw data", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, info);
        }

        // center the columns on 0
        offset[det] = center_dt(&dt, &bt, &bi, 24);
        if (plot::S && plot::centered) {
            hist_info info(plot::x_axes::centered, plot::y_axis.at(det), (format("tdc_calibration_centered_%1%") % detector_map.at(det)).str(), "Raw data", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, info);
        }

        // perform the previous W1 fit on the data
        // this creates a slightly worse net result, but it should be more accurate
        apply_fit(&ft, &fi, &id, deltaF);

        // perform the actual fit
        tie(deltaF[det], deltaB[det]) = fit(&ft, &bt, &fi, &bi, &dt, det);
        if (plot::S && plot::fit) {
            hist_info info(plot::x_axes::centered, plot::y_axis.at(det), (format("tdc_calibration_fit_%1%") % detector_map.at(det)).str(), "Raw data", "dt [ns]", "Back strip id");
            hist2D(&dt, &bi, info);
        }
    }
    // apply the fit on the actual data. note that these two operations are lossless, and cannot delete any entries. 
    // this means that we still have the (alpha_1 | alpha_2 | alpha_3) structure
    apply_offset(data, offset);
    apply_fit(data, deltaF, deltaB);

    return tuple(deltaF, deltaB, offset);
}

// fits a Gaussian to each dt histogram, and cuts it at a given x*sigma. this way we remove the outliers in a statistically sound manner
void gauss_cut(data_container* data, int n_sigma) {
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

    for (int det = 0; det <= 3; det++) {
        vector<double> ft, bt;
        vector<int> fi, bi, id;
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

void save_calibration(vector<vector<double>> doubleF, vector<vector<double>> doubleB, vector<vector<double>> offset) {
    print_title("*** SAVING CALIBRATION RESULTS ***");
    ofstream file;
    file.open("analysis/fit_info.txt");
    file << "doubleF" << endl;
    for (int i = 0; i < doubleF.size(); i++) {
        string line = "";
        for (int j = 0; j < doubleF[i].size(); j++) {
            line += to_string(doubleF[i][j]) + " ";
        }
        file << line << endl;
    }

    file << "\n\ndoubleB" << endl;
    for (int i = 0; i < doubleB.size(); i++) {
        string line = "";
        for (int j = 0; j < doubleB[i].size(); j++) {
            line += to_string(doubleB[i][j]) + " ";
        }
        file << line << endl;
    }

    file << "\n\noffset" << endl;
    for (int i = 0; i < offset.size(); i++) {
        string line = "";
        for (int j = 0; j < offset[i].size(); j++) {
            line += to_string(offset[i][j]) + " ";
        }
        file << line << endl;
    }
    file.close();
    std::cout << "Calibration saved." << endl;
    return;
}