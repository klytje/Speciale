// ROOT stuff
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>
#include <TVector3.h>

// other stuff
#include <boost/format.hpp>
#include <math.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using boost::format;

const double m_alpha = 3.72737*1e3;
const double mu23 = 1./2; // should be m_alpha/2, but that's just a constant away
const int n_sigma = 3; // number of sigmas to include in the peak fits

// from email correspondence with Oliver: 
// "Ved Ep = 2.00 MeV har jeg eff[Be(gs)] = 6.6(3)% og eff[Be(2+)] = 2.10(44)%"
const double gs_eff = 0.066;
const double ex_eff = 0.021;

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        cout << "Usage: ./theta_proj <output path> <input>" << endl;
        exit(1);
    }

    setup_style();
    // auto max = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
    //     TVector3 p1 = {x[0], y[0], z[0]};
    //     TVector3 p2 = {x[1], y[1], z[1]};
    //     TVector3 p3 = {x[2], y[2], z[2]};
    //     vector<double> mags = {p1.Mag(), p2.Mag(), p3.Mag()};
    //     sort(mags.begin(), mags.end());
    //     if (p1.Mag() > mags[1]) {
    //         return 0;
    //     } else if (p2.Mag() > mags[1]) {
    //         return 1;
    //     } else {
    //         return 2;
    //     }
    // };

    // auto mid = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
    //     TVector3 p1 = {x[0], y[0], z[0]};
    //     TVector3 p2 = {x[1], y[1], z[1]};
    //     TVector3 p3 = {x[2], y[2], z[2]};
    //     vector<double> mags = {p1.Mag(), p2.Mag(), p3.Mag()};
    //     sort(mags.begin(), mags.end());
    //     if (p1.Mag() > mags[1]) {
    //         if (mags[1] > p2.Mag()) {
    //             return 2;
    //         } else {
    //             return 1;
    //         }
    //     } else if (p2.Mag() > mags[1]) {
    //         if (mags[1] > p3.Mag()) {
    //             return 0;
    //         } else {
    //             return 2;
    //         }
    //     } else {
    //         if (mags[1] > p1.Mag()) {
    //             return 1;
    //         } else {
    //             return 0;
    //         }
    //     }
    // };

    // auto theta = [] (double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double r3x, double r3y, double r3z) {
    //     TVector3 p23 = {r2x-r3x, r2y-r3y, r2z-r3z};
    //     TVector3 p1 = {r1x, r1y, r1z};
    //     return cos(p1.Angle(p23));

    // };

    auto theta = [] (double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double r3x, double r3y, double r3z) {
        TVector3 p1 = {r1x, r1y, r1z};
        TVector3 p2 = {r2x, r2y, r2z};
        TVector3 p3 = {r3x, r3y, r3z};
        TVector3 p23 = p2-p3;
        auto num = p1.Dot(p23);
        auto denom = p1.Mag()*p23.Mag();
        return num/denom;
    };

    auto max = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
        TVector3 p1 = {x[0], y[0], z[0]};
        TVector3 p2 = {x[1], y[1], z[1]};
        TVector3 p3 = {x[2], y[2], z[2]};
        vector<double> m = {p1.Mag(), p2.Mag(), p3.Mag()};
        double maxv = std::max({m[0], m[1], m[2]});
        if (maxv == m[0]) {
            return 0;
        } else if (maxv == m[1]) {
            return 1;
        } else {
            return 2;
        }
    };

    auto min = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
        TVector3 p1 = {x[0], y[0], z[0]};
        TVector3 p2 = {x[1], y[1], z[1]};
        TVector3 p3 = {x[2], y[2], z[2]};
        vector<double> m = {p1.Mag(), p2.Mag(), p3.Mag()};
        double minv = std::min({m[0], m[1], m[2]});
        if (minv == m[0]) {
            return 0;
        } else if (minv == m[1]) {
            return 1;
        } else {
            return 2;
        }
    };    

    auto mid = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
        TVector3 p1 = {x[0], y[0], z[0]};
        TVector3 p2 = {x[1], y[1], z[1]};
        TVector3 p3 = {x[2], y[2], z[2]};
        vector<double> m = {p1.Mag(), p2.Mag(), p3.Mag()};
        if (m[0] > m[1]) {
            if (m[1] > m[2]) { // m0 > m1 > m2
                return 1;
            } else if (m[0] > m[2]) { // m0 > m2 > m1
                return 2;
            } else { // m2 > m0 > m1
                return 0;
            }
        } else if (m[1] > m[2]) { // m1 > m0, m2
            if (m[0] > m[2]) { // m1 > m0 > m2
                return 0;
            } else { // m1 > m2 > m0
                return 2;
            }
        } else { // m2 > m1 > m0
            return 1;
        }
    };

    auto E23 = [] (double r2x, double r2y, double r2z, double r3x, double r3y, double r3z) {
        TVector3 p23 = {r2x-r3x, r2y-r3y, r2z-r3z};
        return p23.Mag()/(2*mu23);
    };

    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", argv[2]);
    filter(&df);
    df = df.Define("i_max", max, {"px", "py", "pz"})
            .Define("i_min", min, {"px", "py", "pz"})
            .Define("i_mid", "3 - i_max - i_min") // i_min + i_mid + i_max = 3
            // .Define("px1", "px[0]")
            // .Define("py1", "py[0]")
            // .Define("pz1", "pz[0]")
            // .Define("px2", "px[1]")
            // .Define("py2", "py[1]")
            // .Define("pz2", "pz[1]")
            // .Define("px3", "px[2]")
            // .Define("py3", "py[2]")
            // .Define("pz3", "pz[2]");
            .Define("px1", "px[i_max]")
            .Define("py1", "py[i_max]")
            .Define("pz1", "pz[i_max]")
            .Define("px2", "px[i_mid]")
            .Define("py2", "py[i_mid]")
            .Define("pz2", "pz[i_mid]")
            .Define("px3", "px[i_min]")
            .Define("py3", "py[i_min]")
            .Define("pz3", "pz[i_min]");

    // TH1D h1 = df.Define("x", "cos(phi_cm/180)").Histo1D({"h1", "h", 100, 0, 1}, "x").GetValue();
    // TH1D h2 = df.Define("x", "cos(theta_cm/180)").Histo1D({"h2", "h", 100, 0, 1}, "x").GetValue();
    // h1.SetLineColor(kOrange+1);
    // h1.SetLineWidth(2);
    // h1.Scale(1./h1.GetMaximum());
    // h2.SetLineColor(kAzure+1);
    // h2.SetLineWidth(2);
    // h2.Scale(1./h2.GetMaximum());
    // h1.Draw("HIST L SAME");
    // h2.Draw("HIST L SAME");

    // auto func = [] (double* x, double* par) {
    //     double scale = par[0];
    //     double c1 = 0.368245, c3 = 0.630681;
    //     double beta = x[0];
    //     return scale*(c1*correlation_functions.at("2- 1")(beta) + c3*correlation_functions.at("2- 3")(beta));
    // };
    // TF1* prediction = new TF1("mix", func, 0, 3.14, 1);
    // prediction->SetParameter(0, h1.GetMaximum());

    // prediction->Draw("same");

//*** ANGLE PLOT FOR GS, EX ***//
    // set the axes
    double bins = 200;
    vector<double> x_axis;
    vector<double> y_axis;
    x_axis = {bins, -1, 1};
    y_axis = {bins, -1, 1};

    df = df.Define("E_tot","E_cm[0] + E_cm[1] + E_cm[2]")
            .Define("e_cm_1", "E_cm[0]/E_tot") // normalized such that e1 + e2 + e3 = 1
            .Define("e_cm_2", "E_cm[1]/E_tot")
            .Define("e_cm_3", "E_cm[2]/E_tot")
            .Define("e_1", emax, {"e_cm_1", "e_cm_2", "e_cm_3"}) // we want e1 > e2 > e3
            .Define("e_2", emid, {"e_cm_1", "e_cm_2", "e_cm_3"})
            .Define("e_3", emin, {"e_cm_1", "e_cm_2", "e_cm_3"})
            .Define("x","sqrt(3)*(e_2 - e_3)")
            .Define("y","3*e_1 - 1");

    TCanvas* c1 = new TCanvas("c1", "c", 600, 600);
    TH1D h_theta = df.Define("theta", theta, {"px1", "py1", "pz1", "px2", "py2", "pz2", "px3", "py3", "pz3"}).Histo1D({"h_theta", "h", 100, -1, 1}, "theta").GetValue();
    TH1D h_theta_gs = df.Filter("y >= 0.93").Define("theta", theta, {"px1", "py1", "pz1", "px2", "py2", "pz2", "px3", "py3", "pz3"}).Histo1D({"h_theta_gs", "h", 100, -1, 1}, "theta").GetValue();
    TH1D h_theta_ex = df.Filter("y < 0.93").Define("theta", theta, {"px1", "py1", "pz1", "px2", "py2", "pz2", "px3", "py3", "pz3"}).Histo1D({"h_theta_ex", "h", 100, -1, 1}, "theta").GetValue();
    h_theta_gs.SetLineColor(kOrange+1);
    h_theta_gs.SetLineWidth(2);
    h_theta_gs.Scale(1./h_theta.GetMaximum());
    h_theta_ex.SetLineColor(kAzure+1);
    h_theta_ex.SetLineWidth(2);
    h_theta_ex.Scale(1./h_theta.GetMaximum());
    h_theta.SetLineColor(kBlack);
    h_theta.SetLineWidth(2);
    h_theta.Scale(1./h_theta.GetMaximum());

    h_theta.Draw("HIST L");
    h_theta_gs.Draw("SAME HIST L");
    h_theta_ex.Draw("SAME HIST L");

    TLegend* legend = new TLegend(0.7, 0.65, 0.9, 0.9);
    legend->AddEntry("h_theta", "All transitions", "l");
    legend->AddEntry("h_theta_gs", "Ground state", "l");
    legend->AddEntry("h_theta_ex", "Excited state", "l");
    legend->Draw();

    string path = string(argv[1]) + "breakup_angle_data.pdf";
    c1->SaveAs(path.c_str());

//*** E_23 PLOT ***//
    df = df.Define("E_23", E23, {"px2", "py2", "pz2", "px3", "py3", "pz3"});
    TCanvas* c2 = new TCanvas("c2", "c", 600, 600);
    TH1D h2 = df.Histo1D({"h2", "h", 200, 0, 250000}, "E_23").GetValue();
    h2.SetLineColor(kBlack);
    h2.SetLineWidth(2);
    h2.Scale(1./h2.GetMaximum());

    double midpoint = 90000;
    double pars[6];
    TF1* p1 = new TF1("p1", "gaus", 0, midpoint);
    TF1* p2 = new TF1("p2", "gaus", midpoint, 250000);
    h2.Fit(p1, "LQR");
    h2.Fit(p2, "LQR+");
    p1->GetParameters(&pars[0]);
    p2->GetParameters(&pars[3]);

    TF1* pcomb = new TF1("p3", "gaus(0) + gaus(3)", 0, 250000);
    pcomb->SetParameters(pars);
    h2.Fit(pcomb, "QR");

    double mu_gs = pcomb->GetParameter(1), mu_ex = pcomb->GetParameter(4);
    double sigma_gs = pcomb->GetParameter(2), sigma_ex = pcomb->GetParameter(5);
    vector<vector<double>> peaks = {{mu_gs - n_sigma*sigma_gs, mu_gs + n_sigma*sigma_gs}, {mu_ex - n_sigma*sigma_ex, mu_ex + n_sigma*sigma_ex}};
    TLine* l1 = new TLine(peaks[0][0], 0, peaks[0][0], 1);
    TLine* l2 = new TLine(peaks[0][1], 0, peaks[0][1], 1);
    TLine* l3 = new TLine(peaks[1][0], 0, peaks[1][0], 1);
    TLine* l4 = new TLine(peaks[1][1], 0, peaks[1][1], 1);

    h2.Draw("HIST L");
    pcomb->Draw("same");
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");

    path = string(argv[1]) + "E23.pdf";
    c2->SaveAs(path.c_str());

//*** DALITZ PLOT OF EACH PEAK ***//
    // calculate and print the ratio of gs vs ex events
    double c_gs = df.Filter("y > 0.93").Count().GetValue();
    double c_ex = df.Count().GetValue() - c_gs;
    cout << "\nRATIO BASED ON DALITZ Y=0.93 CUT:" << endl;
    cout << "    Ratio of ground state vs excited detections: " << c_gs/c_ex << endl;
    cout << "    Ratio corrected for detection efficiency: " << (c_gs/gs_eff)/(c_ex/ex_eff) << endl;

    ROOT::RDF::RNode df_gs = df.Filter((format("%1% < E_23 && E_23 < %2%") % peaks[0][0] % peaks[0][1]).str());
    ROOT::RDF::RNode df_ex = df.Filter((format("!(%1% < E_23 && E_23 < %2%)") % peaks[0][0] % peaks[0][1]).str());
    TH2D* hgs = new TH2D("hgs", "Dalitz plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    TH2D* hex = new TH2D("hex", "Dalitz plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);

    // we can get the other slices simply by permutating i, j, k
    int perms[] = {1, 2, 3};
    do {
        int i = perms[0];
        int j = perms[1];
        int k = perms[2];
        TH2D htgs = df_gs.Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                    .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                    .Histo2D({"htgs", "", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp").GetValue();
        TH2D htex = df_ex.Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                    .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                    .Histo2D({"htex", "", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp").GetValue();
        hgs->Add(&htgs);
        hex->Add(&htex);
    } while (std::next_permutation(perms, perms+3)); // repeat for each of the 3! = 6 permutations of {1, 2, 3}

    c_gs = df_gs.Count().GetValue();
    c_ex = df.Count().GetValue() - c_gs;
    cout << "\nRATIO BASED ON E23 FIT (check that this looks good!):" << endl;
    cout << "    Ratio of ground state vs excited detections: " << c_gs/c_ex << endl;
    cout << "    Ratio corrected for detection efficiency: " << (c_gs/gs_eff)/(c_ex/ex_eff) << endl;

    TCanvas* c3 = new TCanvas("c3", "c", 600, 600);
    hgs->Draw("colz");
    path = string(argv[1]) + "gs.pdf";
    c3->SaveAs(path.c_str());

    TCanvas* c4 = new TCanvas("c4", "c", 600, 600);
    hex->Draw("colz");
    path = string(argv[1]) + "ex.pdf";
    c4->SaveAs(path.c_str());
}