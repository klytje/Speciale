// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>

// other stuff
#include <boost/format.hpp>
#include <iostream>
#include <math.h>

// my stuff
#include "../plot_style.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[]) {
    if (argc != 5) {
        cout << "Usage: ./ang_compare_single <output path> <nuclear state> <l> <file>" << endl;
        exit(1);
    }

    //*** PREPARE CORRELATION FUNCTION ***//
    function<double(Double_t*, Double_t*)> ang_corr; 
    string state = string(argv[2]);
    string l = string(argv[3]);

    // these two values are determined by eye for all possibilities
    vector<double> bounds; // y axis bounds
    double interference_point; // point of maximum interference

    // they are all calculated with my own angular_correlation.py script
    if (state == "0+" && l == "2") {
        bounds = {0.35, 0.45};
        interference_point = 0.7;
        ang_corr = [] (Double_t* x, Double_t* par) {
            double d = par[0]; 
            double maxval = par[1];
            double xp = M_PI/2*(1 + x[0]/d); // beta = 180 - theta', and theta' = pi/2(1 - x/d)
            double res = 2.25*pow(cos(xp), 4) - 1.5*pow(cos(xp), 2) + 0.25;
            return maxval*res;
        };
    } else if (state == "2-" && l == "1") {
        bounds = {0.28, 0.37};
        interference_point = 0.56;
        ang_corr = [] (Double_t* x, Double_t* par) {
            double d = par[0]; 
            double maxval = par[1];
            double xp = M_PI/2*(1 + x[0]/d); // beta = 180 - theta', and theta' = pi/2(1 - x/d)
            double res = 1 - pow(cos(xp), 2);
            return maxval*res;
        };
    } else if (state == "3-" && l == "1") {
        bounds = {0.35, 0.45};
        interference_point = 0.6;
        ang_corr = [] (Double_t* x, Double_t* par) {
            double d = par[0]; 
            double maxval = par[1];
            double xp = M_PI/2*(1 + x[0]/d); // beta = 180 - theta', and theta' = pi/2(1 - x/d)
            double res = 1./3*pow(cos(xp), 2) + 2./3;
            return maxval*res;
        };
    } else {
        cout << format("The given state %1% l = %2% is not implemented yet! \nYou can do it pretty easily yourself in ../scripts/special/ang_compare.cpp") % state % l << endl;
    }
    
    //*** DALITZ PLOT SETUP (copied from ../basic/dalitz.cpp) ***//
    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", argv[4]);
    double bins = 200;
    double x_axis[] = {bins, -1, 1};
    double y_axis[] = {bins, -1, 1};

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

    // define the necessary variables
    df = df.Define("E_tot","E_cm[0] + E_cm[1] + E_cm[2]")
           .Define("e_cm_1", "E_cm[0]/E_tot") // normalized such that e1 + e2 + e3 = 1
           .Define("e_cm_2", "E_cm[1]/E_tot")
           .Define("e_cm_3", "E_cm[2]/E_tot")
           .Define("e_1", max, {"e_cm_1", "e_cm_2", "e_cm_3"}) // we want e1 > e2 > e3
           .Define("e_2", mid, {"e_cm_1", "e_cm_2", "e_cm_3"})
           .Define("e_3", min, {"e_cm_1", "e_cm_2", "e_cm_3"})
           .Define("x","sqrt(3)*(e_2 - e_3)")
           .Define("y","3*e_1 - 1")
           .Filter("pow(x,2) + pow(y,2) < 1.0")
           .Filter("y < 0.93");
    
    //*** PLOT ***//
    setup_style();
    TCanvas* c = new TCanvas("c", "c", 600, 600);
    TH2D* hist = new TH2D("h1", "Dalitz plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);

    // we can get the other slices simply by permutating i, j, k
    int perms[] = {1, 2, 3};
    do {
        int i = perms[0];
        int j = perms[1];
        int k = perms[2];
        TH2D htemp = df.Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                       .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                       .Histo2D({"h1", "temp", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp").GetValue();
        hist->Add(&htemp);
    } while (std::next_permutation(perms, perms+3)); // repeat for each of the 3! = 6 permutations of {1, 2, 3}
    
    hist->GetXaxis()->SetTitle("x");
    hist->GetYaxis()->SetTitle("y");
    hist->Draw("colz");

    TLine* bot = new TLine(-1, bounds[0], 1, bounds[0]);
    TLine* top = new TLine(-1, bounds[1], 1, bounds[1]);
    top->SetLineWidth(2);
    bot->SetLineWidth(2);
    top->SetLineColor(kRed);
    bot->SetLineColor(kRed);
    top->Draw();
    bot->Draw();

    string path = string(argv[1]) + "ang_compare_raw.pdf";
    c->SetLogz();
    c->SetRightMargin(0.15);
    c->SaveAs(path.c_str());

    //*** PROJECTION ***//
    TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
    int biny1 = 0;
    int biny2 = 0;
    int b = 0;
    for (double i = -1; i < 1; i+=2./bins) {
        b++;
        if (bounds[0] < i && biny1 == 0) biny1 = b;
        if (bounds[1] < i && biny2 == 0) biny2 = b;
    }
    cout << format("bins: %1%, %2%") % biny1 % biny2 << endl;
    TH1D* h1 = hist->ProjectionX("px", biny1, biny2);
    h1->Draw();

    // angular correlation calculated with my own tool
    double d = sqrt(1 - pow((bounds[1] + bounds[0])/2, 2)); // distance across the circle at y = mean(y1, y2)

    // we need to determine the max value of the histogram, but it cannot be at the interference point (about 0.6)
    int maxval = 0;
    double location = 0;
    double width = 0.1;
    double bad_area[4] = {-(interference_point+width), -(interference_point-width), interference_point-width, interference_point+width};
    for (int i = 0; i < bins; i++) {
        int count = h1->GetBinContent(i);
        double loc = h1->GetBinCenter(i);
        if (maxval < count) {
            if (!(bad_area[0] < loc && loc < bad_area[1]) && !(bad_area[2] < loc && loc < bad_area[3])) {
                maxval = count;
                location = loc;
            }
        }
    }
    cout << format("Maximum histogram value found at bin %1%") % location << endl;
    TLine* l1 = new TLine(bad_area[0], 0, bad_area[0], maxval);
    TLine* l2 = new TLine(bad_area[1], 0, bad_area[1], maxval);
    TLine* l3 = new TLine(bad_area[2], 0, bad_area[2], maxval);
    TLine* l4 = new TLine(bad_area[3], 0, bad_area[3], maxval);
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);
    l3->SetLineWidth(2);
    l4->SetLineWidth(2);

    l1->SetLineColor(kRed);
    l2->SetLineColor(kRed);
    l3->SetLineColor(kRed);
    l4->SetLineColor(kRed);

    l1->Draw();
    l2->Draw();
    l3->Draw();
    l4->Draw();

    TF1 *corr = new TF1("corr", ang_corr, -1, 1, 2);
    corr->SetParameter(0, d);
    corr->SetParameter(1, maxval);
    corr->Draw("same");

    path = string(argv[1]) + "ang_compare_projection.pdf";
    c2->SaveAs(path.c_str());
}