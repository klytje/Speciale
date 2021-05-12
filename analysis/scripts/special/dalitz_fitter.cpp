// ROOT stuff
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>
#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <TCanvas.h>

// other stuff
#include <filesystem>
#include <boost/format.hpp>
#include <math.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using namespace ROOT::Math;

// setup for the Dalitz plots
const int bins = 100;
const double range[] = {-1, 1};

class container {
public:
    container() {};
    container(vector<double>* x, vector<double>* y, vector<vector<vector<double>>>* f, vector<double>* wU) {
        this->x = x;
        this->y = y;
        this->f = f;
        this->wU = wU;
    }
    container(vector<double>* x, vector<double>* y) {
        this->x = x;
        this->y = y;
    }
    vector<double>* x; 
    vector<double>* y; 
    vector<double>* wU;
    vector<vector<vector<double>>>* f;
};

// weights defined by eq 42 in Morten's thesis
double calc_weight(vector<vector<double>> f, double wU, double k, double delta) { 
    return wU*(k*f[0][0]+(1-k)*f[0][1] + 2*sqrt(k*(1-k))*(f[0][2]*cos(delta) + f[0][3]*sin(delta)));
}

class MinuitFitter : public ROOT::Math::IMultiGenFunction {
public:
    MinuitFitter() : fit_type(-1) {}

    MinuitFitter(TH2D* hdata, container* sim, int ndim) : Ndim(ndim), fit_type(0) {
        this->hdata = hdata;
        this->sim1 = sim;
        sim1x = *(*sim).x;
        sim1y = *(*sim).y;
        f = *(*sim).f;
        wU = *(*sim).wU;
    }

    MinuitFitter(TH2D* hdata, container* sim1, container* sim2, int ndim) : Ndim(ndim), fit_type(1) {
        this->hdata = hdata;
        this->sim1 = sim1;
        this->sim2 = sim2;
        sim1x = *(*sim1).x;
        sim1y = *(*sim1).y;
        sim2x = *(*sim2).x;
        sim2y = *(*sim2).y;
        f = *(*sim1).f;
        wU = *(*sim1).wU;
    }
    
    ROOT::Math::IBaseFunctionMultiDim *Clone() const override {
        if (fit_type == 1) {
            return new MinuitFitter(hdata, sim1, Ndim);
        } else if (fit_type == 2) {
            return new MinuitFitter(hdata, sim1, sim2, Ndim);
        }
        return nullptr;
    }

    unsigned int NDim() const override {
        return Ndim;
    }

    double DoEval(const double *params) const override {
        const double k = params[0];
        const double delta = params[1];
        const double c = fit_type == 1 ? 1 : params[2];

        vector<double> w = vector<double>(sim1x.size());
        for (int i = 0; i < w.size(); i++) {
            w[i] = calc_weight(f[i], wU[i], k, delta);
        }
        TH2D* hsim = dalitz_slice(sim1x, sim1y, bins, w);
        hsim->Scale(1./hsim->GetMaximum()); // since we scale the histograms, the chi2 value is meaningless

        double chi = 0;
        for (int i = 1; i < hsim->GetNbinsX(); i++) { // start from 1 to exclude overflow bins
            for (int j = 1; j < hsim->GetNbinsY(); j++) {
                if (hdata->GetBinContent(i, j) != 0) {
                    chi += pow(hsim->GetBinContent(i, j) + hdata->GetBinContent(i, j), 2)/hdata->GetBinContent(i, j);
                }
            }
        }
        hsim->Delete();
        return chi;
        // return hdata->Chi2Test(hsim, "UW NORM");
    }

private:
    const int fit_type;
    TH2D* hdata;
    container *sim1, *sim2; // only stored for cloning
    vector<vector<vector<double>>> f;
    vector<double> wU;
    vector<double> sim1x, sim2x;
    vector<double> sim1y, sim2y;
    int Ndim;
};

// the idea is to take two simulated data sets as input and fit them to the data through their Dalitz plots
int main(int argc, char const *argv[]) {
    if (!(argc == 4 || argc == 5)) {
        cout << "Two modes are supported: " << endl;
        cout << "./dalitz_fitter <output path> <data> <sim3a_i data> <sim data>" << endl;
        cout << "./dalitz_fitter <output path> <data> <sim3a_i data>" << endl;
        exit(1);
    }
    setup_style();

//*** PREPARE THE DATA ***//
    ROOT::RDF::RNode data = ROOT::RDataFrame("tree", argv[2]);
    ROOT::RDF::RNode sim1 = ROOT::RDataFrame("tree", argv[3]);
    filter(&data); // perform energy and momentum cut
    filter(&sim1);
    setup_dataframe(&data); // define dalitz coordinates
    setup_dataframe(&sim1);
    // cut_edges(&data); // cut edges
    // cut_edges(&sim1);
    cut_gs(&data); // cut the ground state decay events
    TH2D* hdata;
    {
        vector<double> datax = data.Take<double>("x").GetValue();
        vector<double> datay = data.Take<double>("y").GetValue();
        hdata = dalitz_slice(datax, datay, bins);
        hdata->Scale(1./hdata->GetMaximum());

        // debug plot
        // TCanvas* c = new TCanvas("c", "c", 600, 600);
        // TH2D* test = dalitz(datax, datay, 2*bins);
        // test->Draw("colz");
        // string path = string(argv[1]) + "tmp.pdf";
        // c->SetLogz();
        // c->SetRightMargin(0.15);
        // c->SaveAs(path.c_str());
    }

    container* csim1;
    container* csim2;

    vector<double> sim1x, sim1y, sim2x, sim2y, sim1wU; 
    vector<vector<vector<double>>> sim1f;
    {
        sim1x = sim1.Take<double>("x").GetValue();
        sim1y = sim1.Take<double>("y").GetValue();
        if (sim1.HasColumn("f")) {
            sim1f = sim1.Take<vector<vector<double>>>("f").GetValue();
            sim1wU = sim1.Take<double>("wU").GetValue();
        } else {
            cout << "\033[1;31m" << "ERROR: First simulation data set does not appear to be from sim3a_i" << "\033[0m" << endl;
            exit(1);
        }
        csim1 = new container(&sim1x, &sim1y, &sim1f, &sim1wU);
    }

    MinuitFitter* MinFitter;
    if (argc == 5) {
        ROOT::RDF::RNode sim2 = ROOT::RDataFrame("tree", argv[4]);
        filter(&sim2);
        setup_dataframe(&sim2);
        sim2x = sim2.Take<double>("x").GetValue();
        sim2y = sim2.Take<double>("y").GetValue();
        csim2 = new container(&sim2x, &sim2y);
        MinFitter = new MinuitFitter(hdata, csim1, csim2, 3);
    } else {
        MinFitter = new MinuitFitter(hdata, csim1, 2);
    }
    
//*** PERFORM FIT ***//
    auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan"); 
    auto functor = ROOT::Math::Functor(MinFitter, &MinuitFitter::DoEval, MinFitter->NDim());
    minimizer->SetFunction(functor);
    minimizer->SetPrintLevel(2);
    minimizer->SetLimitedVariable(0, "k", 0.186, 0.01, 0, 1);
    minimizer->SetLimitedVariable(1, "delta", 2*M_PI*0.691, 0.02, 0, 2*M_PI);
    // minimizer->SetFixedVariable(1, "delta", 0);
    minimizer->Minimize();

//*** GENERATE FIGURES ***//
    string path = string(argv[1]) + "dalitz_fit/";
    filesystem::create_directories(path);

    const double* res = minimizer->X();
    const double k = res[0], delta = res[1];
    auto weights = [&k, &delta] (vector<vector<double>> f, double wU) { 
        return wU*(k*f[0][0]+(1-k)*f[0][1] + 2*sqrt(k*(1-k))*(f[0][2]*cos(delta) + f[0][3]*sin(delta)));
    };
    sim1 = sim1.Define("w", weights, {"f", "wU"});

    // RADIAL PROJECTION
    vector<double> x_axis = {100, 0, 1};
    TCanvas* c1 = new TCanvas("c1", "c", 600, 600);
    TH1D sim_rho = sim1.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))")
                        .Histo1D({"sim_rho", "sim_rho", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp", "w").GetValue();
    TH1D dat_rho = data.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))").Histo1D({"dat_rho", "dat_rho", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp").GetValue();
    sim_rho.Scale(1/sim_rho.GetMaximum());
    dat_rho.Scale(1/dat_rho.GetMaximum());
    setup_compare_plot(&dat_rho, &sim_rho, "\\rho", "Arbitrary units");

    path = string(argv[1]) + "dalitz_fit/rho.pdf";
    c1->SetLeftMargin(0.15);
    c1->SaveAs(path.c_str());

    // ANGULAR PROJECTION
    x_axis = {100, 0, M_PI/3};
    TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
    TH1D sim_ang = sim1.Define("tmp", "atan2(x, y)").Histo1D({"sim_ang", "sim_ang", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp", "w").GetValue();
    TH1D dat_ang = data.Define("tmp", "atan2(x, y)").Histo1D({"dat_ang", "dat_ang", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp").GetValue();
    sim_ang.Scale(1/sim_ang.GetMaximum());
    dat_ang.Scale(1/dat_ang.GetMaximum());
    setup_compare_plot(&dat_ang, &sim_ang, "\\varphi", "Arbitrary units");

    path = string(argv[1]) + "dalitz_fit/phi.pdf";
    c2->SetLeftMargin(0.15);
    c2->SaveAs(path.c_str());

    // ENERGY COMPARISON
    x_axis = {100, 0, 7000};
    TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);
    TH1D* sim_E = new TH1D("sim_E", "sim_E", int(x_axis[0]), x_axis[1], x_axis[2]);
    TH1D* dat_E = new TH1D("dat_E", "dat_E", int(x_axis[0]), x_axis[1], x_axis[2]);
    for (int i = 0; i < 3; i++) {
        TH1D sim_temp = sim1.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"sim_temp", "sim_temp", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp", "w").GetValue();
        TH1D dat_temp = data.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"dat_temp", "dat_temp", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp").GetValue();
        sim_E->Add(&sim_temp);
        dat_E->Add(&dat_temp);
    }
    sim_E->Scale(1/sim_E->GetMaximum());
    dat_E->Scale(1/dat_E->GetMaximum());
    setup_compare_plot(dat_E, sim_E, "E_{cm}", "Arbitrary units");

    path = string(argv[1]) + "dalitz_fit/E_cm.pdf";
    c3->SetLeftMargin(0.15);
    c3->SaveAs(path.c_str());
}