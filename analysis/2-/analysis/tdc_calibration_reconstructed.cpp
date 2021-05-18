#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>

#include "../../scripts/calibrate/calibrate_constructed.cpp"

using namespace std;

/* 
    The main changes in this analysis compared to just calibrate.cpp is that all functions must respect that a whole lot of entries in FT and BT may be 0.
    To allow for this, I've had to make a few minor changes in most functions, all of which I've only made available locally in this file
*/
int main(int argc, char *argv[]) {
    data_container data;
    prepare_data(argc, argv, &data, "mul==-1");

    gROOT->SetBatch(kTRUE); // no graphics display

    // start a ROOT application window such that the plots can actually be shown
    TApplication *app = new TApplication("ROOT window", 0, 0);

    // set the x_axis for FT and BT plots
    plot::x_axes::standard = {100, -150, -50};
    plot::x_axes::FT = {400, 13900, 14300};
    plot::x_axes::BT = {400, 13900, 14300};
    plot::path += "reconstructed/";
    
    // ensures that the file path exists
    setup();

    // attempt to align the peaks since each detector is offset slightly from the others
    vector<vector<int>> ft_peaks = {{14000, 14150}}; // I manually read these off my data 
    vector<vector<int>> bt_peaks = {{13930, 14020}, {14020, 14070}, {14070, 14160}, {14160, 14250}};
    align_peaks(&data, ft_peaks, bt_peaks);

    // save(&data, "reconstructed_aligned_peaks.root");

    // imposes a Gaussian filter on FT and BT, to remove outliers. 
    vector<int> ft_peak = {ft_peaks[0][1] - ft_peaks[0][0], ft_peaks[0][0], ft_peaks[0][1]}; 
    vector<int> bt_peak = {bt_peaks[0][1] - bt_peaks[0][0], bt_peaks[0][0], bt_peaks[0][1]};
    gauss_filter_custom(&data, 5, ft_peak, bt_peak);

    auto[deltaF, deltaB, offset] = load_calibration();

    // perform a tdc calibration on the data
    apply_tdc_calibration(&data, deltaF, deltaB);

    // impose a Gaussian filter on DT, to remove outliers
    gauss_cut_custom(&data, 3);
    
    // save the data_container
    save(&data, "output/reconstructed_events.root");
    return 0;
}
