#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>

#include "../../scripts/calibrate/calibrate.cpp"

using namespace std;

int main(int argc, char *argv[]) {
    data_container data;
    prepare_data(argc, argv, &data, "mul==3");

    gROOT->SetBatch(kTRUE); // no graphics display. if kFALSE, you must click on each figure as it appears to continue the script. 

    // start a ROOT application window such that the plots can actually be shown (probably not necessary in batch mode)
    TApplication *app = new TApplication("ROOT window", 0, 0);

    // set the x_axis for FT and BT plots
    plot::x_axes::standard = {100, -150, -50};
    plot::x_axes::FT = {400, 14060, 14120};
    plot::x_axes::BT = {400, 13950, 14150};
    plot::path += "true/";

    // ensures that the file path exists
    setup();

    // attempt to align the peaks since each detector is offset slightly from the others
    vector<vector<int>> ft_peaks = {{14060, 14120}};
    vector<vector<int>> bt_peaks = {{13950, 14050}, {14050, 14130}};
    align_peaks(&data, ft_peaks, bt_peaks);

    // imposes a Gaussian filter on FT and BT, to remove outliers. 
    vector<double> ft_peak = {100, 14060, 14120};
    vector<double> bt_peak = {200, 13950, 14050};
    gauss_filter_custom(&data, 3, ft_peak, bt_peak);

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
