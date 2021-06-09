// ROOT stuff
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>
#include <TVector3.h>
#include <TLatex.h>

// other stuff
#include <boost/format.hpp>
#include <math.h>
#include <filesystem>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using boost::format;

const double mu23 = m_alpha/2; 
const double mu1_23 = m_alpha*2./3;
const int n_sigma = 3; // number of sigmas to include in the peak fits

int main(int argc, char const *argv[]) {
    double gs_eff = 0;
    double ex_eff = 0;
    if (argc != 3) {
        cout << "Usage: ./theta_proj <output path> <input>" << endl;
        exit(1);
    }

    if (string(argv[1]).find("0+") != string::npos) {
        cout << "Output path contains \"0+\", assuming we are dealing with that state." << endl;
        // from email correspondence with Oliver: 
        // "Ved Ep = 2.00 MeV har jeg eff[Be(gs)] = 6.6(3)% og eff[Be(2+)] = 2.10(44)%"
        gs_eff = 0.066;
        ex_eff = 0.021;
    } else if (string(argv[1]).find("3-") != string::npos) {
        cout << "Output path contains \"3-\", assuming we are dealing with that state." << endl;
        // from email correspondence with Oliver: 
        // "Ved Ep = 2.64 MeV har jeg eff[Be(gs)] = 4.32(22)% og eff[Be(2+)] = 2.12(45)%"
        gs_eff = 0.043;
        ex_eff = 0.021;
    } else {
        cout << "\033[1;31m" << "Could not deduce the state based on output path. Branching ratio cannot be calculated." << "\033[0m" << endl;
    } 
    setup_style();
    string dir = string(argv[1]) + "momentum_analysis/";
    filesystem::create_directories(dir);

    // auto theta = [] (double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double r3x, double r3y, double r3z) {
    //     TVector3 p23 = {r2x-r3x, r2y-r3y, r2z-r3z};
    //     TVector3 p1 = {r1x, r1y, r1z};
    //     return cos(p1.Angle(p23));

    // };

    auto theta = [] (double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double r3x, double r3y, double r3z) {
        TVector3 p1 = {r1x, r1y, r1z};
        TVector3 p2 = {r2x, r2y, r2z};
        TVector3 p3 = {r3x, r3y, r3z};
        TVector3 p23 = 1./2*(p2-p3);
        TVector3 p1_23 = 2./3*p1 - 2./3*(p2+p3);
        auto num = p1_23.Dot(p23);
        auto denom = p1_23.Mag()*p23.Mag();
        return num/denom;
    };

    auto theta2 = [] (double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double r3x, double r3y, double r3z) {
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
        TVector3 p2 = {r2x, r2y, r2z};
        TVector3 p3 = {r3x, r3y, r3z};
        TVector3 p23 = p2-p3;
        return 1./(4*m_alpha)*p23.Mag2();
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

    h_theta.GetXaxis()->SetTitle("cos \\theta");
    h_theta.GetXaxis()->CenterTitle();
    h_theta.GetXaxis()->SetNdivisions(3);
    h_theta.GetYaxis()->SetTitle("Normalized count");
    h_theta.GetYaxis()->CenterTitle();
    h_theta.GetYaxis()->SetNdivisions(2);
    h_theta.Draw("HIST L");
    h_theta_gs.Draw("SAME HIST L");
    h_theta_ex.Draw("SAME HIST L");

    TLegend* legend = new TLegend(0.6, 0.65, 0.9, 0.9);
    legend->AddEntry("h_theta", "All transitions", "l");
    legend->AddEntry("h_theta_gs", "Ground state", "l");
    legend->AddEntry("h_theta_ex", "Excited state", "l");
    legend->SetTextSize(0.04);
    legend->Draw();

    string path = dir + "breakup_angle_data.pdf";
    c1->SaveAs(path.c_str());

    TCanvas* c6 = new TCanvas("c6", "c", 600, 600);
    h_theta_ex.GetXaxis()->SetTitle("cos \\theta");
    h_theta_ex.GetXaxis()->CenterTitle();
    h_theta_ex.GetXaxis()->SetNdivisions(3);
    h_theta_ex.GetYaxis()->SetTitle("Normalized count");
    h_theta_ex.GetYaxis()->CenterTitle();
    h_theta_ex.GetYaxis()->SetNdivisions(2);

    h_theta_ex.SetAxisRange(-1, 0, "X");
    h_theta_ex.Scale(1./h_theta_ex.GetMaximum());

    // perform fit with only 2- states (to get accurate error estimates)
    auto f13 = [] (double* x, double* par) {
        double theta = acos(x[0]);
        double c = par[0];
        double scale = par[1];
        double v = c*correlation_functions.at("2- 1")(theta) + (1-c)*correlation_functions.at("2- 3")(theta);
        return v*scale;
    };
    TF1* tf13 = new TF1("mix", f13, -1, 0, 3);
    tf13->SetParameter(0, 0.5); // c
    tf13->SetParameter(1, 1); // scale

    tf13->SetParLimits(0, 0, 1);
    h_theta_ex.Fit(tf13, "QRL");
    double c13 = tf13->GetParameter(0);
    double e13 = tf13->GetParError(0);
    cout << format("Fitted values: 2- 1 = %1% (%2%), 2- 3 = %3% (%4%)") % c13 % e13 % (1-c13) % e13 << endl;

    h_theta_ex.Draw("HIST L");
    tf13->Draw("SAME");
    path = dir + "breakup_angle_2-_fit.pdf";
    cout << path << endl;
    c6->SaveAs(path.c_str());

    TCanvas* c7 = new TCanvas("c7", "c", 600, 600);
    auto f = [] (double* x, double* par) {
        double theta = acos(x[0]);
        double c1 = par[0];
        double c2 = par[1];
        double scale = par[2];
        double v = (1-c2)*c1*correlation_functions.at("2- 1")(theta) + (1-c2)*(1-c1)*correlation_functions.at("2- 3")(theta) 
                    + c2*correlation_functions.at("0+ 2")(theta);
        return v*scale;
    };
    TF1* tf = new TF1("mix", f, -1, 0, 3);
    tf->SetParameter(0, 0.5); // c1
    tf->SetParameter(1, 0.4); // c2
    tf->SetParameter(2, 1); // scale

    tf->SetParLimits(0, 0, 1); // c1
    tf->SetParLimits(1, 0, 1); // c2
    h_theta_ex.Fit(tf, "QRL");
    double k1 = tf->GetParameter(0);
    double k2 = tf->GetParameter(1);
    double e1 = tf->GetParError(0);
    double e2 = tf->GetParError(1);
    cout << format("Fitted values: 0+ 2 = %1% (%2%), 2- 1 = %3% (%4%), 2- 3 = %5% (%6%)") % k2 % e2 % ((1-k2)*k1) % ((e1/k1+e2/k2)*((1-k2)*k1)) % ((1-k2)*(1-k1)) % ((e1/k1+e2/k2)*((1-k2)*k1)) << endl;

    h_theta_ex.Draw("HIST L");
    tf->Draw("SAME");
    path = dir + "breakup_angle_2-_0+_fit.pdf";
    cout << path << endl;
    c7->SaveAs(path.c_str());

//*** E_23 PLOT ***//
    double pleft[] = {0, 200};
    double pright[] = {1500, 4500};
    df = df.Define("E_23", E23, {"px2", "py2", "pz2", "px3", "py3", "pz3"});

    TCanvas *c2 = new TCanvas("c2", "c2", 900, 600);
    TH1D h2 = df.Histo1D({"h2", "h", 200, 0, 6000}, "E_23").GetValue();
    double pars[6];
    TF1* tf_gs = new TF1("tf_gs", "gaus", pleft[0], pleft[1]);
    TF1* tf_ex = new TF1("tf_ex", "gaus", pright[0], pright[1]);
    h2.Fit(tf_gs, "LQR");
    h2.Fit(tf_ex, "LQR+");
    tf_gs->GetParameters(&pars[0]);
    tf_ex->GetParameters(&pars[3]);

    TF1* tf_both = new TF1("tf_both", "gaus(0) + gaus(3)", pleft[0], pright[1]);
    tf_both->SetParameters(pars);
    h2.Fit(tf_both, "LQR");

    double mu_gs = tf_both->GetParameter(1), mu_ex = tf_both->GetParameter(4);
    double sigma_gs = tf_both->GetParameter(2), sigma_ex = tf_both->GetParameter(5);
    vector<vector<double>> peaks = {{mu_gs - n_sigma*sigma_gs, mu_gs + n_sigma*sigma_gs}, {mu_ex - n_sigma*sigma_ex, mu_ex + n_sigma*sigma_ex}};
    TLine* l1 = new TLine(peaks[0][0], 0, peaks[0][0], h2.GetMaximum());
    TLine* l2 = new TLine(peaks[0][1], 0, peaks[0][1], h2.GetMaximum());
    TLine* l3 = new TLine(peaks[1][0], 0, peaks[1][0], h2.GetMaximum());
    TLine* l4 = new TLine(peaks[1][1], 0, peaks[1][1], h2.GetMaximum());
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);
    l3->SetLineWidth(2);
    l4->SetLineWidth(2);
    l3->SetLineStyle(11);
    l4->SetLineStyle(11);
    h2.GetXaxis()->SetTitle("Energy [keV]");
    h2.GetXaxis()->CenterTitle();
    h2.GetYaxis()->SetTitle("Count");
    h2.GetYaxis()->CenterTitle();

    h2.SetNdivisions(205, "X");
    h2.SetNdivisions(205, "Y");
    h2.SetAxisRange(-1000, 5000, "X");
    h2.GetYaxis()->SetMaxDigits(3);
    h2.GetYaxis()->SetTickLength(0);
    h2.GetXaxis()->SetTitleOffset(0.75);
    h2.Draw("HIST L");
    tf_both->Draw("same");
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");

    path = dir + "E23_raw.pdf";
    c2->SaveAs(path.c_str());

//*** PRETTY E23 ***//
    TCanvas *c3 = new TCanvas("c3", "c3", 900, 600);
    TPad *p1 = new TPad("p1", "p1", 0.1, 0.1, 0.5, 0.9);
    p1->SetRightMargin(0.); 
    p1->SetBottomMargin(0.1);           
    p1->SetBorderMode(0);
    p1->Draw();

    TPad *p2 = new TPad("p2", "p2", 0.5, 0.1, 0.9, 0.9); 
    p2->SetLeftMargin(0.);   
    p2->SetBottomMargin(0.1);           
    p2->SetBorderMode(0);
    p2->Draw();

    int scale = 5; // relates the bin count in one panel to the other. needed for consistent y counts
    TH1D hleft = df.Histo1D({"hleft", "h", int(pleft[1]-pleft[0])/4, pleft[0], pleft[1]}, "E_23").GetValue();
    hleft.SetLineColor(kBlack);
    hleft.SetLineWidth(2);

    TH1D hright = df.Histo1D({"hright", "h", int(pright[1]-pright[0])/(4*scale), pright[0], pright[1]}, "E_23").GetValue();
    hright.SetLineColor(kBlack);
    hright.SetLineWidth(2);

    // scale each bin of the right panel for consistent counts
    for (int i = 0; i < hright.GetNbinsX()+1; i++) { // n+1th bin is the overflow bin
        hright.SetBinContent(i, hright.GetBinContent(i)/scale);
    }

    TLatex* x = new TLatex(0.5, 0.08, "Energy [keV]");
    x->SetNDC();
    x->SetTextSize(0.07);
    x->SetTextAlign(22);
    x->Draw();

    TLatex* yl = new TLatex(0.07, 0.5, "^{8}Be Count");
    yl->SetNDC();
    yl->SetTextAngle(90);
    yl->SetTextSize(0.07);
    yl->SetTextAlign(22);
    yl->Draw();

    TLatex* yr = new TLatex(0.93, 0.5, "^{8}Be* Count");
    yr->SetNDC();
    yr->SetTextAngle(270);
    yr->SetTextSize(0.07);
    yr->SetTextAlign(22);
    yr->Draw();

    p1->cd();
    hleft.SetNdivisions(205, "X");
    hleft.SetNdivisions(205, "Y");
    hleft.GetYaxis()->SetMaxDigits(3);
    // hleft.GetYaxis()->SetLabelSize(0.07); 
    hleft.Draw("HIST L");
    
    p2->cd();
    hright.SetNdivisions(205, "X");
    hright.SetNdivisions(205, "Y");
    hright.GetYaxis()->SetMaxDigits(3);
    // hright.GetYaxis()->SetTickLength(0.05);
    hright.Draw("HIST L Y+");

    // pad to cover the line separating the two panels
    c3->cd();
    TPad *b = new TPad("b" , "b", 0.48, 0.1, 0.52, 0.813);
    b->SetBorderMode(0);
    b->Draw();
    b->cd();

    // extension of the axes
    TLine *lleft = new TLine(0, 0.111, 0.4, 0.111);
    lleft->Draw();
    TLine *lright = new TLine(0.5860092, 0.111, 1, 0.111);
    lright->Draw();

    // vertical break lines
    TLine *lvright = new TLine(0.5143349, 0.076639, 0.6863532, 0.1524797);
    lvright->Draw();
    TLine *lvleft = new TLine(0.3423165, 0.076639, 0.5143349, 0.1524797);
    lvleft->Draw();

    path = dir + "E23.pdf";
    p2->SetRightMargin(0.12);
    c3->SetRightMargin(0.12);
    c3->SaveAs(path.c_str());

//*** DALITZ PLOT OF EACH PEAK ***//
    // calculate and print the ratio of gs vs ex events
    double c_gs = df.Filter("y > 0.93").Count().GetValue();
    double c_ex = df.Count().GetValue() - c_gs;
    cout << "\nRATIO BASED ON DALITZ Y=0.93 CUT:" << endl;
    cout << "    c_gs: " << c_gs << ", c_ex: " << c_ex << endl;
    cout << "    Ratio of ground state vs excited detections: " << c_gs/c_ex << endl;
    cout << "    Ratio corrected for detection efficiency: " << (c_gs/gs_eff)/(c_ex/ex_eff) << endl;

    // all events within 3 sigma of the ground state are assumed to belong to that branch. likewise with the excited state. 
    // when this cut is used practically, this same filter is used to extract only the excited state decays
    ROOT::RDF::RNode df_gs = df.Filter((format("%1% < E_23 && E_23 < %2%") % peaks[0][0] % peaks[0][1]).str());
    ROOT::RDF::RNode df_ex = df.Filter((format("%1% < E_23 && E_23 < %2%") % peaks[1][0] % peaks[1][1]).str());
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
    c_ex = df_ex.Count().GetValue();
    cout << "\nRATIO BASED ON E23 FIT (check that this looks good!):" << endl;
    cout << "    c_gs: " << c_gs << ", c_ex: " << c_ex << endl;
    cout << "    Ratio of ground state vs excited detections: " << c_gs/c_ex << endl;
    cout << "    Ratio corrected for detection efficiency: " << (c_gs/gs_eff)/(c_ex/ex_eff) << endl;

    TCanvas* c4 = new TCanvas("c4", "c", 600, 600);
    hgs->Draw("colz");
    path = dir + "gs.pdf";
    c4->SaveAs(path.c_str());

    TCanvas* c5 = new TCanvas("c5", "c", 600, 600);
    hex->Draw("colz");
    path = dir + "ex.pdf";
    c5->SaveAs(path.c_str());
}