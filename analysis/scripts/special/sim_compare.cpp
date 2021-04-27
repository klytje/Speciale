// ROOT stuff
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>

// other stuff
#include <filesystem>
#include <boost/format.hpp>
#include <iostream>

// my stuff
#include "../plot_style.cpp"

using namespace std;
using namespace ROOT;
using boost::format;

int main(int argc, char *argv[]) {
    if (!(argc == 4 || argc == 5)) {
        std::cout << "\033[1;31m" << "Received too many inputs: " << argc << ". Usage: ./sim_compare <data> <simulated data> <destination> <optional: cut>" << "\033[0m" << endl;
        exit(1);
    }

    // prepare the dataframes
    ROOT::RDF::RNode data = RDataFrame("tree", argv[1]);
    ROOT::RDF::RNode sim = RDataFrame("tree", argv[2]);
    string dest = argv[3];

    int cut = 0; // this means we do not make a cut
    if (argc == 5) {
        cut = atoi(argv[4]);
        dest += "sim_compare_cut.pdf";
    } else {
        dest += "sim_compare_raw.pdf";    
    }

    // set the axes    
    double x_axis[] = {200, -1.3, 1.3};
    double y_axis[] = {200, -1.3, 1.3};

    auto prepare_dataframe = [] (ROOT::RDF::RNode* data) {
        ROOT::RDF::RNode &df = *data;

        // define sorting methods
        auto max = [] (double e1, double e2, double e3) {return std::max({e1, e2, e3});};
        auto min = [] (double e1, double e2, double e3) {return std::min({e1, e2, e3});};
        auto mid = [] (double e1, double e2, double e3) {
            if (e1 > e2) {
                if (e2 > e3) {
                    return e2;
                } else if (e1 > e3) {
                    return e3;
                } else {
                    return e1;
                }
            } else {
                if (e1 > e3) {
                    return e1;
                } else if (e2 > e3) {
                    return e3;
                } else {
                    return e2;
                }
            }
        };

        // create the necessary columns
        df = df.Filter("abs(deltaE)<200")
               .Define("E_tot", "E_cm[0] + E_cm[1] + E_cm[2]")
               .Define("e_cm_1", "E_cm[0]/E_tot") // normalized such that e1 + e2 + e3 = 1
               .Define("e_cm_2", "E_cm[1]/E_tot")
               .Define("e_cm_3", "E_cm[2]/E_tot")
               .Define("e_1", max, {"e_cm_1", "e_cm_2", "e_cm_3"}) // we want e1 > e2 > e3
               .Define("e_2", mid, {"e_cm_1", "e_cm_2", "e_cm_3"})
               .Define("e_3", min, {"e_cm_1", "e_cm_2", "e_cm_3"})
               .Define("X","sqrt(3)*(e_2 - e_3)")
               .Define("Y","3*e_1 - 1")
               .Filter("Y < 0.93") // remove ground state decay
               .Filter("pow(X, 2) + pow(Y, 2) < 0.99"); // remove some irrelevent differences at the borders
    };

    // define the needed columns in both dataframes
    prepare_dataframe(&data);
    prepare_dataframe(&sim);

    //*** PLOT ***//
    setup_style();

    TCanvas* canvas = new TCanvas("sim_compare", "sim_compare", 600, 600);
    TH2D* hist = new TH2D("h1", "sim_compare", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);

    // generate the Dalitz plots for both the actual and simulated data
    TH2D* hdata = new TH2D("hdata", "hdata", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    TH2D* hsim = new TH2D("hsim", "hsim", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    int perms[] = {1, 2, 3};
    do {
        int i = perms[0];
        int j = perms[1];
        int k = perms[2];

        // create the current slice of the Dalitz plot for both the actual and simulated data
        TH2D hdata_slice = data.Define("x", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                          .Define("y", (format("3*e_%1% - 1") % i).str())
                          .Histo2D({"hdata_slice", "sim_compare", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x", "y").GetValue();

        TH2D hsim_slice = sim.Define("x", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                        .Define("y", (format("3*e_%1% - 1") % i).str())
                        .Histo2D({"hsim_slice", "sim_compare", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x", "y").GetValue();

        // add them to the global histogram
        hdata->Add(&hdata_slice); 
        hsim->Add(&hsim_slice);
    } while (std::next_permutation(perms, perms+3)); // repeat for all 3! = 6 permutations of i, j, k
    
    // normalize the histograms
    hdata->Scale(1/hdata->Integral()); // "1" means we add all hdata bins to hist
    hsim->Scale(1/hsim->Integral()); // "-1" means we subtract all hsim bins from hist

    // divide hdata by hsim binwise
    hist->Divide(hdata, hsim);

    // hdata->Chi2Test(hsim, "UU NORM P");

    // any bin higher than cut will be set to cut. this is used to enhance the overall figure when only a few known sources of high count are present
    // this also allows easy comparison with sim_compare figures from other models, since they will have the same z scale
    if (cut != 0) {
        TH2D &h = *hist;
        for (int i = 0; i < h.GetNbinsX(); i++) {
            for (int j = 0; j < h.GetNbinsY(); j++) {
                auto val = h.GetBinContent(i, j);
                h.SetBinContent(i, j, cut < val ? cut : val);
            }
        }
    }

    hist->GetXaxis()->SetTitle("X");
    hist->GetYaxis()->SetTitle("Y");
    hist->Draw("colz");

    //canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->SaveAs(dest.c_str());
}
