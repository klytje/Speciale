#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>

#include "../../calibrate/calibrate.cpp"

using namespace std;

int main(int argc, char *argv[]) {
    data_container data;
    prepare_data(argc, argv, &data, "mul==3");

    // define the plot colour scheme
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gROOT->ForceStyle();
    gROOT->SetBatch(kTRUE); // no graphics display

    // start a ROOT application window such that the plots can actually be shown
    TApplication *app = new TApplication("ROOT window", 0, 0);

    // save the raw data before any cuts or modifications are performed
    save(&data, "output/true_raw.root");

    // set the x_axis for FT and BT plots
    plot::x_axes::FT = {400, 13900, 14300};
    plot::x_axes::BT = {400, 13900, 14300};
    plot::path += "true/";
    // ensures that the file path exists
    setup();

    // attempt to align the peaks since each detector is offset slightly from the others
    vector<vector<int>> ft_peaks = {{14000, 14150}}; // I manually read these off my data 
    vector<vector<int>> bt_peaks = {{13930, 14020}, {14020, 14070}, {14070, 14160}, {14160, 14250}};
    align_peaks(&data, ft_peaks, bt_peaks);

    save(&data, "true_aligned_peaks.root");

    // imposes a Gaussian filter on FT and BT, to remove outliers. 
    vector<double> ft_peak = {100, 14000, 14150}; // the area where we expect the peak to be. this is to help the fitting algorithm
    vector<double> bt_peak = {100, 13930, 14020};
    gauss_filter(&data, 5, ft_peak, bt_peak);

    // perform a tdc calibration on the data
    auto[deltaF, deltaB, offset] = tdc_calibration(&data);

    // save the results of the calibration
    save_calibration(deltaF, deltaB, offset);

    // impose a Gaussian filter on DT, to remove outliers
    gauss_cut(&data, 3);
    
    // save the data_container
    save(&data, "output/true_events.root");
    return 0;
}
