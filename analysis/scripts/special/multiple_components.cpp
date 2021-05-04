// ROOT stuff
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

// other stuff
#include <boost/format.hpp>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using namespace ROOT;
using boost::format;

int main(int argc, char *argv[]) {
    if (argc != 5) {
        cout << "Usage: ./sim_compare <output path> <0+ data> <2- data> <real data>" << endl;
        exit(1);
    }

    setup_style();
    // double delta = 2*M_PI*0.59; // 0.11
    double delta = 0; //2*M_PI*0.11
    double k = 0.630681; // ratio of 2- l = 1 & 2- l = 3 (from ang_cor_fit)
    double k_0 = 0.35; // ratio of 0+ to 2- (total)
    double k_2 = 1 - k_0;
    int bins = 200;

    ROOT::RDF::RNode dsim0 = RDataFrame("tree", argv[2]);
    ROOT::RDF::RNode dsim2 = RDataFrame("tree", argv[3]);
    ROOT::RDF::RNode ddat = RDataFrame("tree", argv[4]);
    filter(&dsim0);
    filter(&dsim2);
    filter(&ddat);

    //*** DALITZ PLOT ***//
    TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
    vector<double> x_axis = {double(bins), -1, 1};
    vector<double> y_axis = {double(bins), -1, 1};

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
    prepare_dataframe(&ddat);
    prepare_dataframe(&dsim0);
    prepare_dataframe(&dsim2);

    // check if we are dealing with sim3a_i data
    if (string(argv[3]).find("_i") != string::npos) {
        cout << "File name contains \"_i\", assuming sim3a_i data..." << endl;
        auto weights = [&k, &delta] (vector<vector<double>> f, double wU) { // weights defined by eq 42 in Morten's thesis
            return wU*(k*f[0][0]+(1-k)*f[0][1] + 2*sqrt(k*(1-k))*(f[0][2]*cos(delta) + f[0][3]*sin(delta)));
        };
        dsim2 = dsim2.Define("w", weights, {"f", "wU"});
    } else {
        dsim2 = dsim2.Define("w", "1");
    }

    TH2D* hsim2 = new TH2D("hsim2", "Dalitz plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    TH2D* hsim0 = new TH2D("hsim0", "Dalitz plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    
    // we can get the other slices simply by permutating i, j, k
    int perms[] = {1, 2, 3};
    do {
        int i = perms[0];
        int j = perms[1];
        int k = perms[2];
        TH2D htemp0 = dsim0.Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                    .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                    .Histo2D({"h1", "temp", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp").GetValue();
        TH2D htemp2 = dsim2.Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                    .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                    .Histo2D({"h1", "temp", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp", "w").GetValue();
        hsim0->Add(&htemp0);
        hsim2->Add(&htemp2);
    } while (std::next_permutation(perms, perms+3)); // repeat for each of the 3! = 6 permutations of {1, 2, 3}

    // the bins are meaningless, so we simply normalize it such that the maximum value is 1
    hsim0->Scale(k_0/hsim0->GetMaximum());
    hsim2->Scale(k_2/hsim2->GetMaximum());
    hsim2->Add(hsim0);
    hsim2->Scale(1/hsim2->GetMaximum());

    hsim2->GetXaxis()->SetTitle("x");
    hsim2->GetXaxis()->CenterTitle();
    hsim2->GetYaxis()->SetTitle("y");
    hsim2->GetYaxis()->CenterTitle();
    hsim2->Draw("colz");

    string path = string(argv[1]) + "dalitz_multiple_components.pdf";
    c1->SetLogz();
    c1->SetRightMargin(0.15);
    c1->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;


    //*** RADIAL COMPARISON ***//
    x_axis = {100, 0, 1};
    TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);

    TH1D sim0_rho = dsim0.Define("x", "sqrt(pow(X,2)+pow(Y,2))").Histo1D({"sim_rho", "sim_rho", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();
    TH1D sim2_rho = dsim2.Define("x", "sqrt(pow(X,2)+pow(Y,2))").Histo1D({"sim_rho", "sim_rho", int(x_axis[0]), x_axis[1], x_axis[2]}, "x", "w").GetValue();
    TH1D dat_rho = ddat.Define("x", "sqrt(pow(X,2)+pow(Y,2))").Histo1D({"dat_rho", "dat_rho", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();
    sim0_rho.Scale(k_0/sim0_rho.GetMaximum());
    sim2_rho.Scale(k_2/sim2_rho.GetMaximum());
    dat_rho.Scale(1/dat_rho.GetMaximum());

    sim2_rho.Add(&sim0_rho);
    sim2_rho.Scale(1/sim2_rho.GetMaximum());

    dat_rho.GetXaxis()->SetTitle("\\rho");
    dat_rho.GetXaxis()->CenterTitle();
    dat_rho.GetYaxis()->SetTitle("Arbitrary units");
    dat_rho.GetYaxis()->CenterTitle();
    
    sim2_rho.SetLineColor(kOrange+1);
    dat_rho.SetLineColor(kBlack);
    sim2_rho.SetLineWidth(2);
    dat_rho.SetLineWidth(2);
    dat_rho.Draw("HIST L");
    sim2_rho.Draw("HIST L SAME");

    path = string(argv[1]) + "rho.pdf";
    c2->SetLeftMargin(0.15);
    c2->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;


    //*** ANGULAR COMPARISON ***//
    x_axis = {100, 0, M_PI/3};

    TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);
    TH1D sim0_ang = dsim0.Define("x", "atan2(X,Y)").Histo1D({"sim0_ang", "sim_ang", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();
    TH1D sim2_ang = dsim2.Define("x", "atan2(X,Y)").Histo1D({"sim2_ang", "sim_ang", int(x_axis[0]), x_axis[1], x_axis[2]}, "x", "w").GetValue();
    TH1D dat_ang = ddat.Define("x", "atan2(X,Y)").Histo1D({"dat_ang", "dat_ang", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();
    sim0_ang.Scale(k_0/sim0_ang.GetMaximum());
    sim2_ang.Scale(k_2/sim2_ang.GetMaximum());
    dat_ang.Scale(1/dat_ang.GetMaximum());

    sim2_ang.Add(&sim0_ang);
    sim2_ang.Scale(1/sim2_ang.GetMaximum());

    dat_ang.GetXaxis()->SetTitle("\\phi");
    dat_ang.GetXaxis()->CenterTitle();
    dat_ang.GetYaxis()->SetTitle("Arbitrary units");
    dat_ang.GetYaxis()->CenterTitle();
    
    sim2_ang.SetLineColor(kOrange+1);
    dat_ang.SetLineColor(kBlack);
    sim2_ang.SetLineWidth(2);
    dat_ang.SetLineWidth(2);
    dat_ang.Draw("HIST L");
    sim2_ang.Draw("HIST L SAME");

    path = string(argv[1]) + "phi.pdf";
    c3->SetLeftMargin(0.15);
    c3->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;


    //*** ENERGY COMPARISON ***//
    x_axis = {100, 0, 7000};
    TCanvas* c4 = new TCanvas("c4", "c4", 600, 600);
    TH1D* sim0_E = new TH1D("sim0_E", "sim_E", int(x_axis[0]), x_axis[1], x_axis[2]);
    TH1D* sim2_E = new TH1D("sim2_E", "sim_E", int(x_axis[0]), x_axis[1], x_axis[2]);
    TH1D* dat_E = new TH1D("dat_E", "dat_E", int(x_axis[0]), x_axis[1], x_axis[2]);

    for (int i = 0; i < 3; i++) {
        TH1D sim0_temp = dsim0.Define("x", (format("E_cm[%1%]") % i).str()).Histo1D({"sim0_temp", "sim_temp", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();
        TH1D sim2_temp = dsim2.Define("x", (format("E_cm[%1%]") % i).str()).Histo1D({"sim2_temp", "sim_temp", int(x_axis[0]), x_axis[1], x_axis[2]}, "x", "w").GetValue();
        TH1D dat_temp = ddat.Define("x", (format("E_cm[%1%]") % i).str()).Histo1D({"dat_temp", "dat_temp", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();
        sim0_E->Add(&sim0_temp);
        sim2_E->Add(&sim2_temp);
        dat_E->Add(&dat_temp);
    }
    sim0_E->Scale(k_0/sim0_E->GetMaximum());
    sim2_E->Scale(k_2/sim2_E->GetMaximum());
    dat_E->Scale(1/dat_E->GetMaximum());

    sim2_E->Add(sim0_E);
    sim2_E->Scale(1/sim2_E->GetMaximum());

    dat_E->GetXaxis()->SetTitle("E_{cm}");
    dat_E->GetXaxis()->CenterTitle();
    dat_E->GetYaxis()->SetTitle("Arbitrary units");
    dat_E->GetYaxis()->CenterTitle();
    
    sim2_E->SetLineColor(kOrange+1);
    dat_E->SetLineColor(kBlack);
    sim2_E->SetLineWidth(2);
    dat_E->SetLineWidth(2);
    dat_E->Draw("HIST L");
    sim2_E->Draw("HIST L SAME");

    path = string(argv[1]) + "E_cm.pdf";
    c4->SetLeftMargin(0.15);
    c4->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;
}