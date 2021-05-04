// ROOT stuff
#include <TCanvas.h>
#include <TH1D.h>

// my own stuff
#include "../calibrate/data_fix.cpp"
#include "../calibrate/utility.cpp"
#include "../calibrate/plots.cpp"

using namespace std;

#include <TRandom.h>

int main(int argc, char *argv[])
{
    // convert the data to vector format & extract the relevant vectors
    data_container data;
    prepare_data(argc-1, argv, &data, "mul==3");
    vector<double>* FT = data.get_double("FT");
    vector<double>* BT = data.get_double("BT");
    vector<int>* ID = data.get_int("ID");

    // change some settings
    gROOT->SetBatch(kTRUE); // no graphics display
    plot::save = false;
    plot::path = string(argv[argc-1]);
    setup();

    // makes the whole plot. used twice; once before "data" is repaired, and once after.
    auto plot = [&] (TCanvas* c) {
        c->Divide(4, 2, 0, 0); // split into 2x2 plot
        for (int det = 0; det < 4; det++) {
            TH1D ft_plot("", "", int(plot::x_axes::FT[0]), plot::x_axes::FT[1], plot::x_axes::FT[2]);
            TH1D bt_plot("", "", int(plot::x_axes::BT[0]), plot::x_axes::BT[1], plot::x_axes::BT[2]);
            vector<double> ft = detector_filter(FT, ID, det);
            vector<double> bt = detector_filter(BT, ID, det);

            for (int i = 0; i < ft.size(); i++) {
                ft_plot.Fill(ft[i]);
                bt_plot.Fill(bt[i]);
            }

            // y axis labels
            if (det == 0) {
                c->cd(1);
                bt_plot.GetYaxis()->SetTitle("BT");
                bt_plot.GetYaxis()->SetTitleSize(0.1);
                bt_plot.GetYaxis()->SetTitleOffset(0.3);
                bt_plot.GetYaxis()->CenterTitle();

                c->cd(5);
                ft_plot.GetYaxis()->SetTitle("FT");
                ft_plot.GetYaxis()->SetTitleSize(0.1);
                ft_plot.GetYaxis()->SetTitleOffset(0.3);
                ft_plot.GetYaxis()->CenterTitle();
            }

            // BT plot
            c->cd(det+1);
            bt_plot.SetLabelSize(0, "X");
            bt_plot.SetTickLength(0, "X");
            bt_plot.SetLabelSize(0, "Y");
            bt_plot.SetTickLength(0, "Y");
            bt_plot.DrawClone();

            // FT plot
            c->cd(det+5); // det + 4 + 1
            ft_plot.GetXaxis()->SetTitle(detector_map.at(det).c_str());
            ft_plot.GetXaxis()->SetTitleSize(0.1);
            ft_plot.GetXaxis()->SetTitleOffset(0.3);
            ft_plot.GetXaxis()->CenterTitle();

            // x axis labels
            ft_plot.SetLabelSize(0, "X");
            ft_plot.SetTickLength(0, "X");
            ft_plot.SetLabelSize(0, "Y");
            ft_plot.SetTickLength(0, "Y");
            ft_plot.DrawClone();

            // clean up after ourselves
            ft_plot.Delete();
            bt_plot.Delete();
        }
    };

    // before
    TCanvas *cb = new TCanvas("before", "", 1200, 600);
    plot(cb);
    string path = plot::path + "tdc_plot_before" + plot::format;
    cb->SaveAs(path.c_str());

    repair_peaks(&data);

    // after
    TCanvas *ca = new TCanvas("after", "", 1200, 600);
    plot(ca);
    path = plot::path + "tdc_plot_after" + plot::format;
    ca->SaveAs(path.c_str());
}