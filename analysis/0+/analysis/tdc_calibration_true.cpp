#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>

#include "../../scripts/calibrate/calibrate.cpp"

using namespace std;

int main(int argc, char *argv[]) {
    data_container data;
    prepare_data(argc, argv, &data, "mul==3");

    gROOT->SetBatch(kTRUE); // no graphics display. if kFALSE, you must click on each figure as it appears to continue the script. 
    plot::path += "true/";
    plot::x_axes::standard = {200, -360, -200};

    // start a ROOT application window such that the plots can actually be shown (probably not necessary in batch mode)
    TApplication *app = new TApplication("ROOT window", 0, 0);

    // ensures that the file path exists
    setup();

    // attempt to repair the broken peaks
    repair_peaks(&data);

    // attempt to align the peaks since each detector is offset slightly from the others
    vector<vector<int>> ft_peaks = {{65400, 65500}, {65500, 65600}}; // I manually read these off my data (just run it without anything to see them)
    vector<vector<int>> bt_peaks = {{65100, 65180}, {65180, 65250}, {65300, 65380}, {65380, 65460}, {65460, 65520}, {65520, 65600}};
    align_peaks(&data, ft_peaks, bt_peaks);

    // save(&data, "true_aligned_peaks.root");

    // imposes a Gaussian filter on FT and BT, to remove outliers. 
    vector<double> ft_peak = {1000, 65300, 66000}; // the area where we expect the peak to be. this is to help the fitting algorithm
    vector<double> bt_peak = {1000, 65000, 65500};
    gauss_filter(&data, 3, ft_peak, bt_peak);

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
