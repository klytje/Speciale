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
#include <boost/format.hpp>
#include <math.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using namespace ROOT::Math;

// setup for the Dalitz plots
const int bins = 200;
const double range[] = {-1, 1};

class container {
public:
    container() {};
    container(vector<double>* x, vector<double>* y, vector<vector<double>>* f, vector<double>* wU) {
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
    vector<vector<double>>* f;
};

class MinuitFitter : public ROOT::Math::IMultiGenFunction {
public:
    MinuitFitter() = default;

    MinuitFitter(TH2D* hdata, container* sim, int ndim, string algo) : Ndim(ndim), AlgoType(algo) {
        this->hdata = hdata;
        this->sim1 = sim;
        fit_type = 1;
        sim1x = *(*sim).x;
        sim1y = *(*sim).y;
        f = *(*sim).f;
        wU = *(*sim).wU;
    }

    MinuitFitter(TH2D* hdata, container* sim1, container* sim2, int ndim, string algo) : Ndim(ndim), AlgoType(algo) {
        this->hdata = hdata;
        this->sim1 = sim1;
        this->sim2 = sim2;
        fit_type = 2;
        sim1x = *(*sim1).x;
        sim1y = *(*sim1).y;
        sim2x = *(*sim2).x;
        sim2y = *(*sim2).y;
        f = *(*sim1).f;
        wU = *(*sim1).wU;
    }
    
    ROOT::Math::IBaseFunctionMultiDim *Clone() const override {
        if (fit_type == 1) {
            return new MinuitFitter(hdata, sim1, Ndim, AlgoType);
        } else if (fit_type == 2) {
            return new MinuitFitter(hdata, sim1, sim2, Ndim, AlgoType);
        }
        return nullptr;
    }

    unsigned int NDim() const override {
        return Ndim;
    }

    void setFunction(std::function<double(double *, const double *)> f) {
        func = f;
    }

    double DoEval(const double *params) const override {
        const double k = params[0];
        const double delta = params[1];
        const double c = fit_type == 0 ? 1 : params[2];

        auto weights = [&k, &delta] (vector<vector<double>> f, double wU) { // weights defined by eq 42 in Morten's thesis
            return wU*(k*f[0][0]+(1-k)*f[0][1] + 2*sqrt(k*(1-k))*(f[0][2]*cos(delta) + f[0][3]*sin(delta)));
        };

        // dsim = dsim.Define("w", weights, {"f", "wU"});


        double chi = 0;

        return 0;
    }

private:
    int fit_type = 0;
    TH2D* hdata;
    container *sim1, *sim2; // only stored for cloning
    vector<vector<double>> f;
    vector<double> wU;
    vector<double> sim1x, sim2x;
    vector<double> sim1y, sim2y;
    vector<Double_t> Xdata, Ydata, error;
    int Ndim{};
    string AlgoType;
    std::function<double(double *, const double *)> func;

public:
    void setYData(vector<Double_t> ydata) {
        Ydata = ydata;
    }

    void setXData(vector<Double_t> xdata) {
        Xdata = xdata;
    }

    void setData(vector<Double_t> xdata, vector<Double_t> ydata) {
        Xdata = xdata;
        Ydata = ydata;
    }

    void setAlgoType(string algoType) {
        if (algoType == "Neyman" || algoType == "Pearson" || algoType == "ML") {
            AlgoType = algoType;
        } else {
            cout << "That Algorithm is not implemented." << endl;
            cout << "Try Neyman, Pearson, ML or UserError" << endl;
        }
    }
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
    string algorithm = "";

//*** PREPARE THE DATA ***//
    ROOT::RDF::RNode data = ROOT::RDataFrame("tree", argv[2]);
    filter(&data); // perform energy and momentum cut
    setup_dataframe(&data); // define dalitz coordinates
    // cut_edges(&data); // cut edges
    cut_gs(&data); // cut the ground state decay events
    TH2D* hdata;
    {
        vector<double> datax = data.Take<double>("x").GetValue();
        vector<double> datay = data.Take<double>("y").GetValue();
        hdata = dalitz_slice(datax, datay, bins);

        TCanvas* c = new TCanvas("c", "c", 600, 600);
        TH2D* test = dalitz(datax, datay, bins);
        test->Draw("colz");
        string path = string(argv[1]) + "tmp.pdf";
        c->SetLogz();
        c->SetRightMargin(0.15);
        c->SaveAs(path.c_str());
    }

    container* csim1;
    container* csim2;

    vector<double> sim1x, sim1y, sim2x, sim2y, sim1wU; 
    vector<vector<double>> sim1f;
    {
        ROOT::RDF::RNode sim1 = ROOT::RDataFrame("tree", argv[3]);
        filter(&sim1);
        setup_dataframe(&sim1);
        sim1x = sim1.Take<double>("x").GetValue();
        sim1y = sim1.Take<double>("y").GetValue();
        if (sim1.HasColumn("f")) {
            sim1f = sim1.Take<vector<double>>("f").GetValue();
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
        MinFitter = new MinuitFitter(hdata, csim1, csim2, 3, algorithm);
    } else {
        MinFitter = new MinuitFitter(hdata, csim1, 2, algorithm);
    }
    
//*** SETUP FIT ***//
    auto minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
    auto functor = ROOT::Math::Functor(MinFitter, &MinuitFitter::DoEval, MinFitter->NDim());
    minimum->SetFunction(functor);
    minimum->Minimize();
}
