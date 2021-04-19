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

    gROOT->SetBatch(kTRUE); // no graphics display. if kFALSE, you must click on each figure as it appears to continue the script
    plot::path += "reconstructed/";

    // start a ROOT application window such that the plots can actually be shown (probably not necessary in batch mode)
    TApplication *app = new TApplication("ROOT window", 0, 0);

    // ensures that the file path exists
    setup();

    // attempt to repair the broken peaks
    repair_peaks(&data);

    // attempt to align the peaks since each detector is offset slightly from the others
    vector<vector<int>> ft_peaks = {{65400, 65500}, {65500, 65600}};
    vector<vector<int>> bt_peaks = {{65300, 65380}, {65380, 65460}, {65460, 65520}, {65520, 65600}};
    align_peaks(&data, ft_peaks, bt_peaks);

    // save(&data, "reconstructed_aligned_peaks.root");

    // imposes a Gaussian filter on FT and BT, to remove outliers. 
    vector<double> ft_peak = {100, 65400, 65500};
    vector<double> bt_peak = {200, 65250, 65450};
    gauss_filter_custom(&data, 3, ft_peak, bt_peak);

    auto[deltaF, deltaB, offset] = load_calibration();

    // perform a tdc calibration on the data
    apply_tdc_calibration(&data, deltaF, deltaB);

    // impose a Gaussian filter on DT, to remove outliers
    gauss_cut_custom(&data, 3);
    
    // save the data_container
    save(&data, "output/reconstructed_events.root");
    return 0;
}
