// ROOT stuff
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
    if (argc != 4) {
        cout << "Usage: ./sim_compare <output path> <output/X.root> <comparison file>" << endl;
        exit(1);
    }

    setup_style();
    int bins = 200;

    ROOT::RDF::RNode dsim = RDataFrame("tree", argv[2]);
    ROOT::RDF::RNode ddat = RDataFrame("tree", argv[3]);
    filter(&dsim);
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
    prepare_dataframe(&dsim);

    // check if we are dealing with sim3a_i data
    if (string(argv[2]).find("_i") != string::npos) {
        double delta; // phase difference, essentially controls the amount of interference between the two terms
        double k; // ratio of l = 1 : l = 3
        if (string(argv[2]).find("0+") != string::npos) {
            delta = 2*M_PI*0.11; // 0.59 is also good
            k = 0.630681;
            cout << "File name contains \"0+\" and \"_i\", assuming 0+ sim3a_i data..." << endl;
        } else if (string(argv[2]).find("2-") != string::npos) {
            delta = 2*M_PI*0.691;
            k = 0.186;
            cout << "File name contains \"2-\" and \"_i\", assuming 2- sim3a_i data..." << endl;
        } else {
            cout << "\033[1;31m" << "File name contains \"_i\", but not a recognized state. Correct this in ../scripts/special/sim_compare.cpp" << "\033[0m" << endl;
            exit(1);
        }
        auto weights = [&k, &delta] (vector<vector<double>> f, double wU) { // weights defined by eq 42 in Morten's thesis
            return wU*(k*f[0][0]+(1-k)*f[0][1] + 2*sqrt(k*(1-k))*(f[0][2]*cos(delta) + f[0][3]*sin(delta)));
        };
        dsim = dsim.Define("w", weights, {"f", "wU"});
    } else {
        dsim = dsim.Define("w", "1");
    }

    TH2D* hsim = new TH2D("h1", "Dalitz plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    
    // we can get the other slices simply by permutating i, j, k
    int perms[] = {1, 2, 3};
    do {
        int i = perms[0];
        int j = perms[1];
        int k = perms[2];
        TH2D htemp = dsim.Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                    .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                    .Histo2D({"h1", "temp", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp", "w").GetValue();
        hsim->Add(&htemp);
    } while (std::next_permutation(perms, perms+3)); // repeat for each of the 3! = 6 permutations of {1, 2, 3}

    // the bins are meaningless, so we simply normalize it such that the maximum value is 1
    hsim->Scale(1/hsim->GetMaximum());

    hsim->GetXaxis()->SetTitle("x");
    hsim->GetYaxis()->SetTitle("y");
    hsim->GetXaxis()->CenterTitle();
    hsim->GetYaxis()->CenterTitle();
    hsim->GetXaxis()->SetNdivisions(2);
    hsim->GetYaxis()->SetNdivisions(2);
    hsim->Draw("colz");

    string path = string(argv[1]) + "dalitz.pdf";
    c1->SetLogz();
    c1->SetRightMargin(0.15);
    c1->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;


    //*** RADIAL COMPARISON ***//
    x_axis = {100, 0, 1};
    TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);

    TH1D sim_rho = dsim.Define("x", "sqrt(pow(X,2)+pow(Y,2))").Histo1D({"sim_rho", "sim_rho", int(x_axis[0]), x_axis[1], x_axis[2]}, "x", "w").GetValue();
    TH1D dat_rho = ddat.Define("x", "sqrt(pow(X,2)+pow(Y,2))").Histo1D({"dat_rho", "dat_rho", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();
    sim_rho.Scale(1/sim_rho.GetMaximum());
    dat_rho.Scale(1/dat_rho.GetMaximum());

    dat_rho.GetXaxis()->SetTitle("\\rho");
    dat_rho.GetYaxis()->SetTitle("Arbitrary units");
    dat_rho.GetXaxis()->CenterTitle();
    dat_rho.GetYaxis()->CenterTitle();
    dat_rho.GetXaxis()->SetNdivisions(205);
    dat_rho.GetYaxis()->SetNdivisions(203);
    
    sim_rho.SetLineColor(kOrange+1);
    dat_rho.SetLineColor(kBlack);
    sim_rho.SetLineWidth(2);
    dat_rho.SetLineWidth(2);
    dat_rho.Draw("HIST L");
    sim_rho.Draw("HIST L SAME");

    path = string(argv[1]) + "rho.pdf";
    c2->SetLeftMargin(0.15);
    c2->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;


    //*** ANGULAR COMPARISON ***//
    x_axis = {100, 0, M_PI/3};

    TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);
    TH1D sim_ang = dsim.Define("x", "atan2(X,Y)").Histo1D({"sim_ang", "sim_ang", int(x_axis[0]), x_axis[1], x_axis[2]}, "x", "w").GetValue();
    TH1D dat_ang = ddat.Define("x", "atan2(X,Y)").Histo1D({"dat_ang", "dat_ang", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();
    sim_ang.Scale(1/sim_ang.GetMaximum());
    dat_ang.Scale(1/dat_ang.GetMaximum());

    dat_ang.GetXaxis()->SetTitle("\\varphi");
    dat_ang.GetYaxis()->SetTitle("Arbitrary units");
    dat_ang.GetXaxis()->CenterTitle();
    dat_ang.GetYaxis()->CenterTitle();
    dat_ang.GetXaxis()->SetNdivisions(206);
    dat_ang.GetYaxis()->SetNdivisions(203);
    
    sim_ang.SetLineColor(kOrange+1);
    dat_ang.SetLineColor(kBlack);
    sim_ang.SetLineWidth(2);
    dat_ang.SetLineWidth(2);
    dat_ang.Draw("HIST L");
    sim_ang.Draw("HIST L SAME");

    path = string(argv[1]) + "phi.pdf";
    c3->SetLeftMargin(0.15);
    c3->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;


    //*** ENERGY COMPARISON ***//
    x_axis = {100, 0, 7000};
    TCanvas* c4 = new TCanvas("c4", "c4", 600, 600);
    TH1D* sim_E = new TH1D("sim_E", "sim_E", int(x_axis[0]), x_axis[1], x_axis[2]);
    TH1D* dat_E = new TH1D("dat_E", "dat_E", int(x_axis[0]), x_axis[1], x_axis[2]);

    for (int i = 0; i < 3; i++) {
        TH1D sim_temp = dsim.Define("x", (format("E_cm[%1%]") % i).str()).Histo1D({"sim_temp", "sim_temp", int(x_axis[0]), x_axis[1], x_axis[2]}, "x", "w").GetValue();
        TH1D dat_temp = ddat.Define("x", (format("E_cm[%1%]") % i).str()).Histo1D({"dat_temp", "dat_temp", int(x_axis[0]), x_axis[1], x_axis[2]}, "x").GetValue();
        sim_E->Add(&sim_temp);
        dat_E->Add(&dat_temp);
    }
    sim_E->Scale(1/sim_E->GetMaximum());
    dat_E->Scale(1/dat_E->GetMaximum());

    dat_E->GetXaxis()->SetTitle("E_{cm}");
    dat_E->GetYaxis()->SetTitle("Arbitrary units");
    dat_E->GetXaxis()->CenterTitle();
    dat_E->GetYaxis()->CenterTitle();
    dat_E->GetXaxis()->SetNdivisions(205);
    dat_E->GetYaxis()->SetNdivisions(203);
    
    sim_E->SetLineColor(kOrange+1);
    dat_E->SetLineColor(kBlack);
    sim_E->SetLineWidth(2);
    dat_E->SetLineWidth(2);
    dat_E->Draw("HIST L");
    sim_E->Draw("HIST L SAME");

    path = string(argv[1]) + "E_cm.pdf";
    c4->SetLeftMargin(0.15);
    c4->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;
}