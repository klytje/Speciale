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
#include <fstream>
#include <boost/format.hpp>
#include <math.h>
#include <time.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using namespace ROOT::Math;

// setup for the Dalitz plots
const int bins = 50; // THIS VALUE AFFECTS THE QUALITY OF THE FIT!

// this container was a good idea at the start before I refactored the code to use histograms for the data and sim2
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

class Dalitz_fitter {
public:
    Dalitz_fitter() = default;

    // fit_type 1 constructor
    Dalitz_fitter(TH2D* hdata, container* sim) {
        this->hdata = hdata;
        sim1x = *(*sim).x;
        sim1y = *(*sim).y;
        f1 = *(*sim).f;
        wU1 = *(*sim).wU;
    }

    // fit_type 2 constructor
    Dalitz_fitter(TH2D* hdata, container* sim1, TH2D* hsim2) {
        this->hdata = hdata;
        this->hsim2 = hsim2;
        sim1x = *(*sim1).x;
        sim1y = *(*sim1).y;
        f1 = *(*sim1).f;
        wU1 = *(*sim1).wU;
    }
    
    // fit_type 3 constructor
    Dalitz_fitter(TH2D* hdata, container* sim1, container* sim2) {
        this->hdata = hdata;

        sim1x = *(*sim1).x;
        sim1y = *(*sim1).y;
        f1 = *(*sim1).f;
        wU1 = *(*sim1).wU;

        sim2x = *(*sim2).x;
        sim2y = *(*sim2).y;
        f2 = *(*sim2).f;
        wU2 = *(*sim2).wU;
    }

    // type 1 evaluation method. we only fit one sim3a_i data set
    double eval_type1(const double *params) const {
        const double k = params[0];
        const double delta = params[1]*2*M_PI;

        // some of the minimizers do not respect our limits, which results in nan values
        if (k > 1) { 
            return 1e9; // I don't know how the minimizer stores its maxval, but DBL_MAX was apparently too high. I think this value is fine though
        }

        // calculate the weight
        vector<double> w = vector<double>(sim1x.size());
        for (int i = 0; i < w.size(); i++) {
            w[i] = calc_weight(f1[i], wU1[i], k, delta);
        }

        // generate the dalitz slices
        TH2D* hsim = dalitz_slice(sim1x, sim1y, bins, w);
        hsim->Scale(1./hsim->GetMaximum()); // since we scale the histograms, the chi2 value is meaningless

        // calculate chi
        double chi = maximum_likelihood(hsim);
        hsim->Delete();
        return chi;
    }

    // type 2 evaluation method. we fit both a sim3a_i and sim3a data set
    double eval_type2(const double *params) const {
        const double k = params[0];
        const double delta = params[1]*2*M_PI;
        const double c = params[2];

        // some of the minimizers do not respect our limits, which results in nan values
        if (k > 1 || c > 1) { 
            return 1e9; // I don't know how the minimizer stores its maxval, but DBL_MAX was apparently too high. I think this value is fine though
        }

        // calculate the weight
        vector<double> w = vector<double>(sim1x.size());
        for (int i = 0; i < w.size(); i++) {
            w[i] = calc_weight(f1[i], wU1[i], k, delta);
        }

        // generate the dalitz slices
        TH2D* hsim = dalitz_slice(sim1x, sim1y, bins, w);
        TH2D* hsim2c = (TH2D*) hsim2->Clone(); // we clone the locally stored sim2 histogram to avoid modifying it

        // scale & add the two simulation histograms according to their ratio c
        hsim->Scale(c/hsim->GetMaximum()); // since we scale the histograms, the chi2 value is meaningless
        hsim2c->Scale((1-c)/hsim2c->GetMaximum()); // since we scale the histograms, the chi2 value is meaningless
        hsim->Add(hsim2c);

        // calculate chi
        double chi = maximum_likelihood(hsim);
        hsim->Delete(); hsim2c->Delete(); // clean up after ourselves
        return chi;
    }

    // type 3 evaluation method. we fit two sim3a_i datasets
    double eval_type3(const double *params) const {
        const double k1 = params[0];
        const double delta1 = params[1]*2*M_PI;
        const double c = params[2];
        const double k2 = params[3];
        const double delta2 = params[4]*2*M_PI;

        // some of the minimizers do not respect our limits, which results in nan values
        if (k1 > 1 || k2 > 1 || c > 1) { 
            return 1e9; // I don't know how the minimizer stores its maxval, but DBL_MAX was apparently too high. I think this value is fine though
        }

        // calculate the weights
        vector<double> w1 = vector<double>(sim1x.size());
        vector<double> w2 = vector<double>(sim2x.size());
        for (int i = 0; i < w1.size(); i++) {
            w1[i] = calc_weight(f1[i], wU1[i], k1, delta1);
        }
        for (int i = 0; i < w2.size(); i++) {
            w2[i] = calc_weight(f2[i], wU2[i], k2, delta2);
        }

        // generate the dalitz slices
        TH2D* hsim1 = dalitz_slice(sim1x, sim1y, bins, w1);
        TH2D* hsim2 = dalitz_slice(sim2x, sim2y, bins, w2);

        // scale & add the two simulation histograms according to their ratio c
        hsim1->Scale(c/hsim1->GetMaximum()); // since we scale the histograms, the chi2 value is meaningless
        hsim2->Scale((1-c)/hsim2->GetMaximum()); // since we scale the histograms, the chi2 value is meaningless
        hsim1->Add(hsim2);

        // calculate chi
        double chi = maximum_likelihood(hsim1);
        hsim1->Delete(); hsim2->Delete(); // clean up after ourselves
        return chi;
    }

    // calculate the log likelihood for a given evaluation
    double maximum_likelihood(TH2D* h) const {
        double chi = 0, eps = 1e-9;
        for (int i = 1; i < h->GetNbinsX(); i++) { // start from 1 to exclude overflow bins
            for (int j = 1; j < h->GetNbinsY(); j++) {
                if (hdata->GetBinContent(i, j) != 0) {
                    double mu = hdata->GetBinContent(i, j); // model
                    double m = h->GetBinContent(i, j); // observed
                    chi += m*log(max(m, eps)/max(mu, eps)) + mu - m; // ML
                }
            }
        }
        return chi;
    }

    // does not work since mu can be 0
    double chi2(TH2D* h) {
        double chi = 0;
        for (int i = 1; i < h->GetNbinsX(); i++) { // start from 1 to exclude overflow bins
            for (int j = 1; j < h->GetNbinsY(); j++) {
                if (hdata->GetBinContent(i, j) != 0) {
                    double mu = hdata->GetBinContent(i, j); // model
                    double m = h->GetBinContent(i, j); // observed
                    chi += pow(m + mu, 2)/mu;
                }
            }
        }
        return chi;
    }

private:
    TH2D *hdata, *hsim2;
    vector<vector<vector<double>>> f1, f2;
    vector<double> wU1, wU2;
    vector<double> sim1x, sim1y, sim2x, sim2y;
};

// the idea is to take two simulated data sets as input and fit them to the data through their Dalitz plots
int main(int argc, char const *argv[]) {
    /* 
        Which type & algorithm do we use? Based on my own investigation, the parameter space has a multitude of local minima, so any algorithm based on 
        first derivates are useless (essentially the whole Minuit2 library). I've had some luck with the GSLMultiMin BFGS algorithm, which is reasonably
        fast and appears to find the actual minima. BFGS2 was also good, but found a wrong minima very close in value to the actual. 
        I think the best option is GSLSimAn, which takes a hell of a time to run (~15k function calls), but found the correct minimum. 
    */
    string type = "GSLSimAn", algorithm = "";
    // note that you can also change the bin count of the Dalitz plot at the top of this file, and that it may affect the quality of the fit

    if (argc < 3) {
        cout << "Two modes are supported: " << endl;
        cout << "./dalitz_fitter <data> <sim3a_i data>" << endl; // fit_type = 1
        cout << "./dalitz_fitter <data> <sim3a_i data> <sim3a data>" << endl; // fit_type = 2
        cout << "./dalitz_fitter <data> <sim3a_i data> <sim3a_i data>" << endl; // fit_type = 3
        cout << "Figures are automatically written to figures/dalitz_fit/" << endl;
        cout << "Only the name of the files should be provided, e.g. output/true_events.root --> true_events" << endl;
        cout << "Guess values can be supplied after the files. Format: .k X .delta Y .c Z" << endl;
        cout << "The algorithm can also be specified similarly to the guess values with .type A .algo B" << endl;
        exit(1);
    }
    setup_style();

    // parse arguments
    int fit_type, guess_pars = 0;
    double guess[] = {0.5, 0.5, 0.2};
    bool fix_delta = false;
    string args[argc-1];
    for (int i = 0; i < argc; i++) {
        args[i] = argv[i];
    }
    for (int i = 3; i < argc; i++) {
        if (args[i].find(".k") != string::npos) {
            guess[0] = atof(args[i+1].c_str());
            guess_pars += 2;
        }
        if (args[i].find(".delta") != string::npos) {
            if (args[i+1].substr(6) == "fixed") {
                fix_delta = true;
            } else {
                guess[1] = atof(args[i+1].c_str());
            }
            guess_pars += 2;
        }
        if (args[i].find(".c") != string::npos) {
            guess[2] = atof(args[i+1].c_str());
            guess_pars += 2;
        }
        if (args[i].find(".type") != string::npos) {
            type = args[i+1];
            if (type == "GSLSimAn") {
                cout << "Warning: Using GSL simulated annealing, expect around 15k function calls." << endl;
            }
            guess_pars += 2;
        }
        if (args[i].find(".algo") != string::npos) {
            algorithm = args[i+1];
            guess_pars += 2;
        }
    }
    if (argc == 3 + guess_pars) {
        fit_type = 1;
    } else if (argc == 4 + guess_pars) {
        fit_type = 2;
    } else {
        cout << "\033[1;31m" << "Invalid number of arguments." << "\033[0m" << endl;
        exit(1);
    }

    // create folder
    time_t curr_t = time(0);
    tm *t = localtime(&curr_t);
    string folder = "figures/dalitz_fit/" + to_string(t->tm_mon) + "." + to_string(t->tm_mday) + " " + to_string(t->tm_hour) + ":" + to_string(t->tm_min) + ":" + to_string(t->tm_sec) + "/";
    filesystem::create_directories(folder);

//*** PREPARE THE DATA ***//
    ROOT::RDF::RNode data = ROOT::RDataFrame("tree", "output/" + args[1] + ".root");
    ROOT::RDF::RNode sim1 = ROOT::RDataFrame("tree", "output/" + args[2] + ".root");
    ROOT::RDF::RNode sim2 = ROOT::RDataFrame(0); // empty dataframe
    filter(&data); // perform energy and momentum cut
    filter(&sim1);
    setup_dataframe(&data); // define dalitz coordinates
    setup_dataframe(&sim1);
    cut_gs(&data); // cut the ground state decay events
    TH2D *hdata, *hsim2;
    {
        hdata = dalitz_slice(&data, bins, true);
        hdata->Scale(1./hdata->GetMaximum());
    }

    // extract the data from the first simulation data set
    container *csim1, *csim2;
    vector<double> sim1x, sim1y, sim2x, sim2y, sim1wU, sim2wU; 
    vector<vector<vector<double>>> sim1f, sim2f;
    {
        if (!sim1.HasColumn("f")) {
            cout << "\033[1;31m" << "ERROR: First simulation data set does not appear to be from sim3a_i" << "\033[0m" << endl;
            exit(1);        
        }
        sim1x = sim1.Take<double>("x").GetValue();
        sim1y = sim1.Take<double>("y").GetValue();
        sim1f = sim1.Take<vector<vector<double>>>("f").GetValue();
        sim1wU = sim1.Take<double>("wU").GetValue();
        csim1 = new container(&sim1x, &sim1y, &sim1f, &sim1wU);
    }

    // extract the data from the second simulation data set
    if (fit_type == 2) {
        sim2 = ROOT::RDataFrame("tree", "output/" + args[3] + ".root");
        filter(&sim2);
        setup_dataframe(&sim2);
        if (sim2.HasColumn("f")) { // check if we are dealing with sim3a_i data
            cout << "Second simulation data set appears to be from sim3a_i." << endl;
            fit_type == 3;

            sim2x = sim2.Take<double>("x").GetValue();
            sim2y = sim2.Take<double>("y").GetValue();
            sim2f = sim2.Take<vector<vector<double>>>("f").GetValue();
            sim2wU = sim2.Take<double>("wU").GetValue();
            csim2 = new container(&sim1x, &sim1y, &sim1f, &sim1wU);
        } else {
            hsim2 = dalitz_slice(&sim2, bins, true);
            hsim2->Scale(1./hdata->GetMaximum());
        }
    }

    // prepare the fitting algorithm
    Dalitz_fitter* fitter;
    ROOT::Math::Functor functor;
    int pars = 0;
    if (fit_type == 1) {
        pars = 2; // k & delta
        fitter = new Dalitz_fitter(hdata, csim1);
        functor = ROOT::Math::Functor(fitter, &Dalitz_fitter::eval_type1, pars);
    } else if (fit_type == 2) {
        pars = 3; // k, delta, c
        fitter = new Dalitz_fitter(hdata, csim1, hsim2);
        functor = ROOT::Math::Functor(fitter, &Dalitz_fitter::eval_type2, pars);
    } else if (fit_type == 3) {
        pars = 5; // k1, delta1, c, k2, delta2
        fitter = new Dalitz_fitter(hdata, csim1, csim2);
        functor = ROOT::Math::Functor(fitter, &Dalitz_fitter::eval_type3, pars);
    }

    // a small function to set the fitting parameters. I know it's longer than it strictly needs to be, but I think this makes it easier to read
    auto set_params = [&fit_type, &fix_delta, &guess] (ROOT::Math::Minimizer* m) {
        if (fit_type == 1) {
            m->SetLimitedVariable(0, "k", guess[0], 0.01, 0, 1); // Morten found k = 0.186, delta = 0.691 for 2-
            if (fix_delta)
                m->SetFixedVariable(1, "delta", 0);
            else 
                m->SetLimitedVariable(1, "delta", guess[1], 0.01, 0, 1);
        }
        if (fit_type == 2) {
            m->SetLimitedVariable(0, "k", guess[0], 0.01, 0, 1);
            if (fix_delta)
                m->SetFixedVariable(1, "delta", 0);
            else 
                m->SetLimitedVariable(1, "delta", guess[1], 0.01, 0, 1);
            m->SetLimitedVariable(2, "c", guess[2], 0.01, 0, 1); // c is the ratio of sim1 to sim2, i.e. c*sim1 + (1-c)*sim2
        }
        if (fit_type == 3) { // order: k1, d1, c, k2, c2
            m->SetLimitedVariable(0, "k1", guess[0], 0.01, 0, 1);
            m->SetLimitedVariable(3, "k2", guess[3], 0.01, 0, 1);
            m->SetLimitedVariable(2, "c", guess[2], 0.01, 0, 1); // c is the ratio of sim1 to sim2, i.e. c*sim1 + (1-c)*sim2
            if (fix_delta) {
                m->SetFixedVariable(1, "delta1", 0);
                m->SetFixedVariable(4, "delta2", 0);
            } else {
                m->SetLimitedVariable(1, "delta1", guess[1], 0.01, 0, 1);
                m->SetLimitedVariable(4, "delta2", guess[4], 0.01, 0, 1);
            }
        }
    };

    // first minimization
    auto minimizer = ROOT::Math::Factory::CreateMinimizer(type.c_str(), algorithm.c_str()); 
    minimizer->SetFunction(functor);
    minimizer->SetPrintLevel(2); // level 2 shows some additional live fitting information, which is nice so we know it is progressing
    set_params(minimizer);
    minimizer->Minimize();
    const double* res = minimizer->X();

    // prints MINOS errors
    std::ofstream file(folder + "info.txt"); // create a file to store all of the fit information
    auto minos_errs = [&fit_type, &res, &file] (ROOT::Math::Minimizer* m) {
        // const double* res2 = m->X();
        double err_up[5], err_low[5];
        m->GetMinosError(0, err_low[0], err_up[0]);
        m->GetMinosError(1, err_low[1], err_up[1]);
        if (fit_type == 2 || fit_type == 3) {
            m->GetMinosError(2, err_low[2], err_up[2]);
        }
        if (fit_type == 3) {
            m->GetMinosError(3, err_low[3], err_up[3]);
            m->GetMinosError(4, err_low[4], err_up[4]);
        }
        file << "\nEstimating errors with MINOS: " << endl;
        file << format("k: %1% | +%2% | %3%") % res[0] % err_up[0] % err_low[0] << endl;
        file << format("delta: %1% | +%2% | %3%") % res[1] % err_up[1] % err_low[1] << endl;
        if (fit_type == 2 || fit_type == 3) {
            file << format("c: %1% | +%2% | %3%") % res[2] % err_up[2] % err_low[2] << endl;
        }
        if (fit_type == 3) {
            file << format("k2: %1% | +%2% | %3%") % res[3] % err_up[3] % err_low[3] << endl;
            file << format("delta2: %1% | +%2% | %3%") % res[4] % err_up[4] % err_low[4] << endl;
        }
    };

    // print a bunch of information to the file
    file << "*** FIT REPORT ***" << endl;
    file << format("Type: %1%, algorithm: %2%") % type % algorithm << endl;
    file << format("Fitting %1% with %2%") % args[1] % args[2];
    if (fit_type == 2 || fit_type == 3) {
        file << " and " << args[3];
    }
    file << format("\nUsing bin width: %1%") % bins << endl;
    file << format("Function value at minimum: FVAL = %1%") % minimizer->MinValue() << endl;
    file << "\nROOT Minimizer report: " << endl;
    auto *coutbuf = std::cout.rdbuf();

    // unless the first fit was a Migrad, we perform a quick second fit with it to estimate the errors.
    if (algorithm != "Migrad") {
        file << "    Results from first fit: " << endl;
        file << format("        k = %1% (FREE)") % res[0] << endl;
        file << format("        delta = %1% (%2%)") % res[1] % (fix_delta ? "FIXED" : "FREE") << endl;
        if (fit_type == 2) {
            file << format("        c = %1% (FREE)") % res[2] << endl;
        }
        cout << "    \nEstimating errors with Migrad" << endl;

        // update our guess values
        for (int i = 0; i < pars; i++) guess[i] = res[i];

        auto minimize2 = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimize2->SetFunction(functor);
        set_params(minimize2);
        minimize2->Minimize(); 
        res = minimize2->X();
        std::cout.rdbuf(file.rdbuf()); // redirect std::cout to fit.txt
        minimize2->PrintResults();
        minos_errs(minimize2);
    } else {
        std::cout.rdbuf(file.rdbuf());
        minimizer->PrintResults();
        minos_errs(minimizer);
    }
    std::cout.rdbuf(coutbuf); // remove the redirection
    file.close();

//*** GENERATE FIGURES ***//
    cout << "\nGenerating figures (this may take a while)" << endl;
    const double k = res[0]; 
    const double delta = res[1]*2*M_PI;
    const double c = fit_type == 2 ? res[2] : 1;
    const double k2 = fit_type == 3 ? res[3] : 1;
    const double delta2 = fit_type == 3 ? res[4] : 1;

    auto weights = [] (vector<vector<double>> f, double wU, double k, double delta) { 
        return wU*(k*f[0][0]+(1-k)*f[0][1] + 2*sqrt(k*(1-k))*(f[0][2]*cos(delta) + f[0][3]*sin(delta)));
    };
    auto w1 = [&k, &delta, &weights] (vector<vector<double>> f, double wU) {return weights(f, wU, k, delta);};
    auto w2 = [&k2, &delta2, &weights] (vector<vector<double>> f, double wU) {return weights(f, wU, k2, delta2);};
    sim1 = sim1.Define("w", weights, {"f", "wU"});

    bool sim2_weighted = false;
    if (fit_type == 3) {
        sim2 = sim2.Define("w", weights, {"f", "wU"});
        sim2_weighted = true;
    }

    // DALITZ PLOTS
    TCanvas* c1 = new TCanvas("c1", "c", 900, 600);
    c1->Divide(2, 1, 0, 0);
    
    c1->cd(1);
    TH2D* dalitz_data = dalitz(&data, 2*bins, true);
    double data_scale = dalitz_data->GetMaximum();
    dalitz_data->Scale(1./data_scale);
    setup_dalitz_plot(dalitz_data, "col");

    c1->cd(2);
    TH2D* dalitz_sim = dalitz(&sim1, 2*bins, true, true);
    dalitz_sim->Scale(c/dalitz_sim->GetMaximum());
    if (fit_type == 2) { // we need to add the two simulations with their respective ratios
        TH2D* dalitz_sim2 = dalitz(&sim2, 2*bins, true, sim2_weighted);
        dalitz_sim2->Scale((1-c)/dalitz_sim2->GetMaximum());
        dalitz_sim->Add(dalitz_sim2);
    }
    setup_dalitz_plot(dalitz_sim, "col");
    
    string path = folder + "dalitz.pdf";
    c1->SetLogz();
    c1->SetRightMargin(0.15);
    c1->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

    // DALITZ DIFFERENCE PLOT
    TCanvas* c2 = new TCanvas("c2", "c", 600, 600);
    dalitz_data->Scale(data_scale); // scale both the data and simulation so the z-axis makes some kind of sense
    dalitz_sim->Scale(data_scale);

    dalitz_data->Add(dalitz_sim, -1); // subtract the two plots
    dalitz_data->Draw("colz");

    path = folder + "dalitz_diff.pdf";
    c2->SetRightMargin(0.15);
    c2->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

    // RADIAL PROJECTION
    vector<double> x_axis = {100, 0, 1};
    TCanvas* c3 = new TCanvas("c3", "c", 600, 600);
    TH1D dat_rho = data.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))").Histo1D({"dat_rho", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp").GetValue();
    TH1D sim_rho = sim1.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))").Histo1D({"sim_rho", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp", "w").GetValue();
    dat_rho.Scale(1./dat_rho.GetMaximum());
    sim_rho.Scale(c/sim_rho.GetMaximum());
    if (fit_type == 2) {
        TH1D sim2_rho;
        if (fit_type == 3) {
            sim2_rho = sim2.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))")
                                .Histo1D({"sim2_rho", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp", "w").GetValue();
        } else {
            sim2_rho = sim2.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))")
                        .Histo1D({"sim2_rho", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp").GetValue();
        }
        sim2_rho.Scale((1-c)/sim2_rho.GetMaximum());
        sim_rho.Add(&sim2_rho);
    }
    setup_compare_plot(&dat_rho, &sim_rho, "\\rho", "Arbitrary units");

    path = folder + "rho.pdf";
    c3->SetLeftMargin(0.15);
    c3->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

    // ANGULAR PROJECTION
    x_axis = {100, 0, M_PI/3};
    TCanvas* c4 = new TCanvas("c4", "c", 600, 600);
    TH1D dat_ang = data.Define("tmp", "atan2(x, y)").Histo1D({"dat_ang", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp").GetValue();
    TH1D sim_ang = sim1.Define("tmp", "atan2(x, y)").Histo1D({"sim_ang", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp", "w").GetValue();
    dat_ang.Scale(1./dat_ang.GetMaximum());
    sim_ang.Scale(c/sim_ang.GetMaximum());
    if (fit_type == 2) {
        TH1D sim2_ang;
        if (fit_type == 3) {
            sim2_ang = sim2.Define("tmp", "atan2(x, y)").Histo1D({"sim2_ang", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp", "w").GetValue();
        } else {
            sim2_ang = sim2.Define("tmp", "atan2(x, y)").Histo1D({"sim2_ang", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp").GetValue();
        }
        sim2_ang.Scale((1-c)/sim2_ang.GetMaximum());
        sim_ang.Add(&sim2_ang);
    }
    setup_compare_plot(&dat_ang, &sim_ang, "\\varphi", "Arbitrary units");

    path = folder + "phi.pdf";
    c4->SetLeftMargin(0.15);
    c4->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

    // ENERGY COMPARISON
    x_axis = {100, 0, 7000};
    TCanvas* c5 = new TCanvas("c5", "c", 600, 600);
    TH1D* dat_E = new TH1D("dat_E", "h", int(x_axis[0]), x_axis[1], x_axis[2]);
    TH1D* sim_E = new TH1D("sim_E", "h", int(x_axis[0]), x_axis[1], x_axis[2]);
    for (int i = 0; i < 3; i++) {
        TH1D dat_temp = data.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"dat_temp", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp").GetValue();
        TH1D sim_temp = sim1.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"sim_temp", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp", "w").GetValue();
        dat_E->Add(&dat_temp);
        sim_E->Add(&sim_temp);
    }
    dat_E->Scale(1/dat_E->GetMaximum());
    sim_E->Scale(c/sim_E->GetMaximum());
    if (fit_type == 2) {
        TH1D* sim2_E = new TH1D("sim2_E", "h", int(x_axis[0]), x_axis[1], x_axis[2]);
        for (int i = 0; i < 3; i++) {
            TH1D sim2_temp;
            if (fit_type == 3) {
                sim2_temp= sim2.Define("tmp", (format("E_cm[%1%]") % i).str())
                                .Histo1D({"sim2_temp", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp", "w").GetValue();
            } else {
                sim2_temp= sim2.Define("tmp", (format("E_cm[%1%]") % i).str())
                                .Histo1D({"sim2_temp", "h", int(x_axis[0]), x_axis[1], x_axis[2]}, "tmp").GetValue();
            }
            sim2_E->Add(&sim2_temp);
        }
        sim2_E->Scale((1-c)/sim2_E->GetMaximum());
        sim_E->Add(sim2_E);
    }
    setup_compare_plot(dat_E, sim_E, "E_{cm}", "Arbitrary units");

    path = folder + "E_cm.pdf";
    c5->SetLeftMargin(0.15);
    c5->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;
}