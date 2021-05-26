// ROOT stuff
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>
#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <TCanvas.h>
#include <TGraph2D.h>

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

int bins;
double y_cut;
vector<vector<double>> p_evals;

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

    // general constructor supporting any number of sim3a_i simulations (untested beyond 2)
    Dalitz_fitter(container* cdata, vector<container*> sim) {
        std::tie(d, std::ignore) = hist(cdata);
        f = vector<vector<vector<vector<double>>>>(sim.size());
        wU = vector<vector<double>>(sim.size());
        for (int j = 0; j < sim.size(); j++) {
            f[j] = *(*sim[j]).f;
            wU[j] = *(*sim[j]).wU;
        }
        if (sim.size() == 1) {
            type = 1;
        } else if (sim.size() == 2) {
            type = 3;
        } else {
            cout << "Critical failure: more than two simulation files passed to constructor. Implement your own evaluation function to support this." << endl;
            exit(1);
        }
        setup_vectors(sim);
    }

    // specific constructor for 1 sim3a_i simulation, and 1 sim3a simulation
    Dalitz_fitter(container* cdata, container* sim3ai, container* sim3a) {
        std::tie(d, std::ignore) = hist(cdata);
        f = vector<vector<vector<vector<double>>>>(1);
        wU = vector<vector<double>>(1);
        f[0] = *(*sim3ai).f;
        wU[0] = *(*sim3ai).wU;
        type = 2;
        setup_vectors({sim3ai, sim3a});
    }

    // type 1 evaluation method for the case of a single sim3a_i simulation
    double eval_type1(const double *params) {
        vector<double> k = {params[0]};
        vector<double> delta = {params[1]*2*M_PI};
        vector<double> P = {1};
        return eval(k, delta, P);
    }

    // type 3 evaluation method for the case of two sim3a_i simulations
    double eval_type3(const double *params) {
        vector<double> k = {params[0], params[3]};
        vector<double> delta = {params[1]*2*M_PI, params[4]*2*M_PI};
        vector<double> P = {params[2], 1-params[2]};
        return eval(k, delta, P);
    }

    // general evaluation function for the modified maximum likelihood
    double eval(const vector<double> k, const vector<double> delta, const vector<double> P) {
        // some of the minimizers do not respect our limits, which results in nan values
        for (int j = 0; j < m; j++) {
            if (k[j] > 1 || P[j] > 1) { 
                return 1e9; // I don't know how the minimizer stores its maxval, but DBL_MAX was apparently too high. I think this value is fine though
            }
        }

        // calculate the weight
        vector<vector<double>> w(m);
        for (int j = 0; j < m; j++) {
            w[j] = vector<double>(entry_map[j].size());
            for (int i = 0; i < w[j].size(); i++) {
                w[j][i] = calc_weight(f[j][i], wU[j][i], k[j], delta[j]);
            }
        }
        vector<vector<double>> avg_w = calc_avg_w(w);
        return modified_maximum_likelihood(&avg_w, &P);
    }

    double eval_type2(const double *params) {
        double k = params[0];
        double delta = params[1]*2*M_PI;
        // vector<double> P = {params[2], 1-params[2]};
        vector<double> P = {1, 0};
        
        // calculate the weight
        vector<vector<double>> w(m);
        w[0] = vector<double>(entry_map[0].size());
        w[1] = vector<double>(entry_map[1].size(), 1);
        for (int i = 0; i < w[0].size(); i++) {
            w[0][i] = calc_weight(f[0][i], wU[0][i], k, delta);
        }
        vector<vector<double>> avg_w = calc_avg_w(w);
        return modified_maximum_likelihood(&avg_w, &P);
    }

    // calculate chi^2 for a given fit result
    double chisquare(const double* params) {
        double chi = 0;
        if (type == 1) {
            chi = eval_type1(params);
        } else if (type == 2) {
            chi = eval_type2(params);
        } else if (type == 3) {
            chi = eval_type3(params);
        } else {
            cout << "Critical failure in evaluating chisquare: unknown evaluation type " << type << endl;
            exit(1);
        }
        for (int i = 0; i < bins2; i++) {
            if (skip_bin[i]) {
                continue;
            }
            if (d[i] != 0) {
                chi += 2*(d[i]*log(d[i]) - d[i]);
            }
            for (int j = 0; j < m; j++) {
                if (a[j][i] != 0) {
                    chi += 2*(a[j][i]*log(a[j][i]) - a[j][i]);
                }
            }
        }
        return chi;
    }

    // get the vector designating the skipped bins
    vector<bool> get_skipped_bins() {
        return skip_bin;
    }

    // get the contribution of each bin to the chi2
    vector<double> binwise_chisquare(const double* params) const {
        vector<double> chi = binwise_chi;

        // subtract the likelihood of the data to get the likelihood ratio
        for (int i = 0; i < bins2; i++) {
            if (d[i] != 0 && !skip_bin[i]) {
                chi[i] += 2*(d[i]*log(d[i]) - d[i]);
            }
            for (int j = 0; j < m; j++) {
                if (a[j][i] != 0 && !skip_bin[i]) {
                    chi[i] += 2*(a[j][i]*log(a[j][i]) - a[j][i]);
                }
            }
        }
        return chi;
    }

    // return the actual number of bins used, i.e. bins^2 - (all bins outside the allowed regions in the Dalitz plot)
    int get_bins() {
        int skipped = 0;
        for (int i = 0; i < bins2; i++) {
            if (skip_bin[i]) {
                skipped++;
            }
        }
        return bins2 - skipped;
    }

    double modified_maximum_likelihood(const vector<vector<double>>* weights, const vector<double>* P_ratio) {
        const vector<vector<double>> &w = *weights;
        const vector<double> &P = *P_ratio;
        vector<double> p(m); // normalized fractions
        vector<double> t(bins2);
        vector<vector<double>> A(m, vector<double>(bins2));
        double chi = 0; // likelihood

        vector<double>Nc(m);
        for (int i = 0; i < bins2; i++) {
            for (int j = 0; j < m; j++) {
                Nc[j] += w[j][i]*a[j][i];
            }
        }

        // calculate the normalized fractions
        for (int j = 0; j < m; j++) {
            p[j] = P[j]*Nd/Nc[j];
        }

        // definition of f for the bisection method
        auto f = [&] (double t, int i) {
            double val = 0;
            for (int j = 0; j < 1; j++) {
                val += p[j]*w[j][i]*a[j][i]/(1 + p[j]*w[j][i]*t);
            }
            val -= d[i]/(1 - t);
            return val;
        };

        // evaluate chi of the ith bin and add it to the running total
        auto eval_chi = [&] (double i) {
            double prev_chi = chi;

            // determine fi
            double fi = 0;
            if (t[i] == 1) {
                for (int j = 0; j < m; j++) {
                    fi += p[j]*w[j][i]*A[j][i];
                }
                // if fi = 0, it would imply wji or aji are also 0, which is only possible outside the boundaries of our Dalitz plot
                if (fi == 0) { 
                    skip_bin[i] = true;
                    return;
                }
            } else {
                fi = d[i]/(1 - t[i]);
            }

            chi += d[i]*log(fi) - fi;
            for (int j = 0; j < m; j++) {
                chi += a[j][i]*log(A[j][i]) - A[j][i];
            }

            // check for nan values
            if (isnan(chi-prev_chi)) {
                cout << "Critical error in bin " << i << ", change in chi is nan (old chi = -inf?)" << endl;
                cout << format("f = %1%, d = %2%, t = %3%") % fi % d[i] % t[i] << endl;
                for (int j = 0; j < m; j++) {
                    cout << format("\ta%1% = %2%, A%1% = %3%, w%1% = %4%") % j % a[j][i] % A[j][i] % w[j][i] << endl;
                }
                exit(1);
            }

            binwise_chi[i] = -2*(chi-prev_chi);
        };

        double tol = 1e-15; // tolerance on ti for bisection method
        int max_evals = 100; // max evaluations before stopping the bisection method
        for (int i = 0; i < bins2; i++) {
            if (skip_bin[i]) { // skips any bin outside x^2 + y^2 = 1 and with y > y_cut
                continue;
            }

        //*** CASE di = 0 ***//
            if (d[i] == 0) {
                // cout << "case di == 0" << endl;
                t[i] = 1;
                for (int j = 0; j < m; j++) {
                    A[j][i] = a[j][i]/(1 + p[j]*w[j][i]*t[i]);
                }
                eval_chi(i);
                continue;
            }

        //*** CASE a_ji = 0 ***//
            int k = -1;
            double pw = -1;
            for (int j = 0; j < m; j++) {
                if (a[j][i] == 0) {
                    A[j][i] = 0;
                    if (p[j]*w[j][i] > pw) {
                        pw = p[j]*w[j][i];
                        k = j;
                    }
                }
            }
            if (pw == 0) { // only possible if all our weights are 0, i.e. no counts have been registered
                skip_bin[i] = true;
                continue;
            }

            // calculate A_ki
            if (k != -1) {
                // cout << "case aji == 0" << endl;
                t[i] = -1./pw;
                A[k][i] = d[i]/(1 + pw);
                for (int j = 0; j < m; j++) {
                    if (j != k) {
                        A[k][i] -= p[j]*w[j][i]*a[j][i]/(p[k] - p[j]);
                    }
                }
                eval_chi(i);
                continue;
            }

        //*** NORMAL CASE ***//
            if (d[i] != 0) {
                // cout << "normal case" << endl;
                for (int j = 0; j < m; j++) {
                    if (p[j]*w[j][i] > pw) {
                        pw = p[j]*w[j][i];
                    }
                }
                double upper_bound = 1;
                double lower_bound = -1./pw;

                // the value of ti where fi changes sign must be between upper_bound and lower_bound
                // we use a binary search (or bisection method, if you prefer) to find the value
                t[i] = 0;
                for (int n = 0; n < max_evals; n++) {
                    double val = f(t[i], i);
                    if (abs(val) < tol) { // break when |f| < tol
                        break;
                    }
                    else if (val < 0) { // sign is unchanged, so we move the upper bound
                        upper_bound = t[i];
                    } else { // sign is changed, so we move the lower bound
                        lower_bound = t[i];
                    }
                    t[i] = (upper_bound + lower_bound)/2;
                    if (isnan(val)) {
                        cout << "Critical error in bin " << i << endl;
                        cout << "f = " << val << ", di = " << d[i] << ", a0i = " << a[0][i] << ", ti = " << t[i] << ", w0i = " << w[0][i] << endl;
                        exit(1);
                    }
                }
                for (int j = 0; j < m; j++) {
                    A[j][i] = a[j][i]/(1 + p[j]*w[j][i]*t[i]);
                }

                eval_chi(i);
                continue;
            }
        }
        return -2*chi;
    }

private:
    // for calculating the weights
    vector<vector<vector<vector<double>>>> f;
    vector<vector<double>> wU;
    int type = -1;

    // everything related to the modified_maximum_likelihood evaluation
    int bins2 = pow(bins, 2); // bins^2
    int m; // number of simulated data sets
    int Nd; // number of data events
    vector<int> N; // number of simulated events
    vector<vector<int>> a; // number of simulated counts for source j in bin i
    vector<int> d; // expected number of events
    vector<double> p; // normalized fractions
    vector<vector<int>> entry_map; // a map from entry --> bin for each source
    vector<bool> skip_bin; // skipped entries (outside the Dalitz plot)
    vector<double> binwise_chi; // contains the chi2 value for each bin

    // calculate the average weight for each bin
    vector<vector<double>> calc_avg_w(vector<vector<double>> weights) const {
        vector<vector<double>> w(m, vector<double>(bins2, 0));
        for (int j = 0; j < m; j++) {
            vector<double> counts(bins2, 0);
            for (int n = 0; n < N[j]; n++) {
                int i = entry_map[j][n];
                if (i != -1) { // -1 means the entry is invalid (either outside the Dalitz plot or outside the y_cut)
                    counts[i] += 1;
                    w[j][i] = (w[j][i]*(counts[i]-1) + weights[j][n])/counts[i];
                }
            }
        }
        return w;
    }

    void setup_vectors(vector<container*> sim) {
        m = sim.size();

        // define the other vectors needed
        Nd = 0; // number of data events
        N = vector<int>(m, 0); // number of simulated events
        a = vector<vector<int>>(m); // number of simulated events for source j in bin i
        entry_map = vector<vector<int>>(m);
        skip_bin = vector<bool>(bins2, false);
        binwise_chi = vector<double>(bins2, 0);

        // define the raw count histograms
        for (int j = 0; j < m; j++) {
            std::tie(a[j], entry_map[j]) = hist(sim[j]);
        }

        // calculate counts
        for (int i = 0; i < bins2; i++) {
            Nd += d[i];
            for (int j = 0; j < m; j++) {
                N[j] += a[j][i];
            }
        }

        // determine which bins we should skip
        for (int i = 0; i < bins2; i++) {
            int zeros = 0;
            for (int j = 0; j < m; j++) {
                if (a[j][i] == 0) { // check if all simulations are zero
                    zeros++;
                }
            }
            if (zeros == m && d[i] == 0) { // if all counts are zero, skip this bin. this should remove any bin outside the dalitz slice
                skip_bin[i] = true;
            }
            else if (i > y_cut*bins2) { // if the bin is above the y_cut, skip it
                skip_bin[i] = true;
            }
        }
    }

    // a simple binning method, which also creates an index mapping each event to its bin
    tuple<vector<int>, vector<int>> hist(container* data) {
        vector<double> &x = *data->x, &y = *data->y;
        vector<int> counts(bins2);
        vector<int> entry_map(x.size());
        double step = 1./bins; // interval between each bin
        for (int i = 0; i < x.size(); i++) {
            if (pow(x[i], 2) + pow(y[i], 2) <= 1 && y[i] < y_cut) {
                int binx = floor(x[i]/step);
                int biny = floor(y[i]/step);
                int bin = binx + bins*biny;
                entry_map[i] = bin;
                counts[bin]++;
            } else {
                entry_map[i] = -1;
            }
        }
        return make_tuple(counts, entry_map);
    }
};

// the idea is to take two simulated data sets as input and fit them to the data through their Dalitz plots
int main(int argc, char const *argv[]) {
    /* 
        Which type & algorithm do we use? Based on my own investigation, the parameter space has a multitude of local minima, so any algorithm based on 
        first derivates are useless (essentially the whole Minuit2 library). I've had some luck with the GSLMultiMin BFGS algorithm, which is reasonably
        fast and appears to find the actual minima. BFGS2 was also good, but found a wrong minima very close in value to the actual. 
        I think the best option is GSLSimAn, which takes a hell of a time to run (~15k function calls), but found the correct minimum. 
    */
    if (argc < 3) {
        cout << "Two modes are supported: " << endl;
        cout << "\t./dalitz_fitter <data> <sim3a_i data>" << endl; // fit_type = 1
        cout << "\t./dalitz_fitter <data> <sim3a_i data> <sim3a data>" << endl; // fit_type = 2
        cout << "\t./dalitz_fitter <data> <sim3a_i data> <sim3a_i data>" << endl; // fit_type = 3
        cout << "Figures are automatically written to figures/dalitz_fit/" << endl;
        cout << "Only the name of the files should be provided, e.g. output/true_events.root --> true_events" << endl;
        cout << "The following list of parameters can be supplied like \".X Y\"" << endl;
        cout << "\tk:      guess for the l : l' ratio of the first simulated data set" << endl;
        cout << "\tdelta:  guess for the delta value of the first simulated data set" << endl;
        cout << "\tk2:     guess for the l : l' ratio of the second simulated data set" << endl;
        cout << "\tdelta2: guess for the delta value of the second simulated data set" << endl;
        cout << "\ttype:   fitting class used (e.g. GSLMultiMin, GSLSimAn, Genetic)" << endl;
        cout << "\talgo:   fitting algorithm used (e.g. BFGS2, Migrad)" << endl;
        cout << "\tbins:   number of bins. This may affect the quality of the fit" << endl;
        cout << "\tycut:   imposes a cut on the Dalitz y-coordinate at this value" << endl;
        cout << "Either delta can also be set to \"fixed\", in which case they will be fixed to 0" << endl;
        exit(1);
    }
    setup_style();

    // parse arguments
    string type = "GSLSimAn", algorithm = "";
    bins = 100; // THIS VALUE AFFECTS THE QUALITY OF THE FIT! this is the number of bins for each axis on a Dalitz *slice*, so the total plot has twice this number
    y_cut = 1;
    int fit_type, guess_pars = 0;
    double guess[] = {0.5, 0.5, 0.2, 0.5, 0.5};
    bool fix_delta = false;
    string args[argc-1];
    for (int i = 0; i < argc; i++) {
        args[i] = argv[i];
    }
    cout << "\nOptional parameters parsed: " << endl;
    for (int i = 3; i < argc; i++) {
        if (args[i].find(".k") != string::npos || args[i].find(".k1") != string::npos) {
            guess[0] = atof(args[i+1].c_str());
            cout << "\tk guess: " << guess[0] << endl;
            guess_pars += 2;
        }
        else if (args[i].find(".delta") != string::npos || args[i].find(".delta") != string::npos) {
            if (args[i+1] == "fixed") {
                fix_delta = true;
            } else {
                guess[1] = atof(args[i+1].c_str());
                cout << "\tdelta guess: " << guess[1] << endl;
            }
            guess_pars += 2;
        }
        else if (args[i].find(".c") != string::npos) {
            guess[2] = atof(args[i+1].c_str());
            guess_pars += 2;
            cout << "\tc guess: " << guess[2] << endl;
        }
        else if (args[i].find(".k2") != string::npos) {
            guess[3] = atof(args[i+1].c_str());
            guess_pars += 2;
            cout << "\tk2 guess: " << guess[3] << endl;
        }
        else if (args[i].find(".delta2") != string::npos) {
            if (args[i+1] == "fixed") {
                fix_delta = true;
            } else {
                guess[4] = atof(args[i+1].c_str());
                cout << "\tdelta2 guess: " << guess[4] << endl;
            }
            guess_pars += 2;
        }
        else if (args[i].find(".type") != string::npos) {
            type = args[i+1];
            guess_pars += 2;
        }
        else if (args[i].find(".algo") != string::npos) {
            algorithm = args[i+1];
            guess_pars += 2;
        }
        else if (args[i].find(".bins") != string::npos) {
            bins = atof(args[i+1].c_str());
            guess_pars += 2;
            cout << "\tnumber of bins set to " << bins << endl;
        }
        else if (args[i].find(".yc") != string::npos) {
            y_cut = atof(args[i+1].c_str());
            guess_pars += 2;
            cout << "\tperforming cut at y = " << y_cut << endl;
        }
    }
    if (fix_delta) {
        cout << "\tAll \u03B4 are fixed at 0, and will not be varied." << endl;
    }
    if (argc == 3 + guess_pars) {
        fit_type = 1;
        cout << format("\nFitting %1% with %2%") % args[1] % args[2] << endl;
    } else if (argc == 4 + guess_pars) {
        fit_type = 2;
        cout << format("\nFitting %1% with %2% and %3%") % args[1] % args[2] % args[3] << endl;
    } else {
        cout << "\033[1;31m" << "Invalid number of arguments." << "\033[0m" << endl;
        exit(1);
    }
    if (type == "GSLSimAn") {
        cout << "Warning: Using GSL simulated annealing, expect around 15k function calls." << endl;
    } else if (type == "Genetic") {
        cout << "Warning: Using Genetic algorithm, expect around an hour of computation time." << endl;
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
    cut_circle(&data); // cut everything outside the unit circle
    cut_circle(&sim1);
    if (y_cut != 1) {
        cut_y(&data, y_cut); // perform a cut on the y axis
        cut_y(&sim1, y_cut);
    }
    cut_gs(&data); // cut the ground state decay events

    // extract the data and put it into a container
    container *cdata, *csim1, *csim2;
    vector<double> datax, datay, sim1x, sim1y, sim2x, sim2y; 
    datax = data.Take<double>("x").GetValue();
    datay = data.Take<double>("y").GetValue();
    cdata = new container(&datax, &datay);

    // extract the data from the first simulation data set
    vector<double> sim1wU, sim2wU;
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
        cut_circle(&sim2);
        if (y_cut != 1) {
            cut_y(&sim2, y_cut);
        }
        sim2x = sim2.Take<double>("x").GetValue();
        sim2y = sim2.Take<double>("y").GetValue();
        if (sim2.HasColumn("f")) { // check if we are dealing with sim3a_i data
            cout << "Second simulation data set appears to be from sim3a_i." << endl;
            fit_type = 3;
            sim2f = sim2.Take<vector<vector<double>>>("f").GetValue();
            sim2wU = sim2.Take<double>("wU").GetValue();
            csim2 = new container(&sim1x, &sim1y, &sim1f, &sim1wU);
        } else {
            csim2 = new container(&sim1x, &sim1y);
        }
    }

    // prepare the fitting algorithm
    Dalitz_fitter* fitter;
    ROOT::Math::Functor functor;
    int pars = 0;
    if (fit_type == 1) {
        pars = 2; // k & delta
        fitter = new Dalitz_fitter(cdata, {csim1});
        functor = ROOT::Math::Functor(fitter, &Dalitz_fitter::eval_type1, pars);
    } else if (fit_type == 2) {
        pars = 3; // k, delta, c
        fitter = new Dalitz_fitter(cdata, csim1, csim2);
        functor = ROOT::Math::Functor(fitter, &Dalitz_fitter::eval_type2, pars);
    } else if (fit_type == 3) {
        pars = 5; // k1, delta1, c, k2, delta2
        fitter = new Dalitz_fitter(cdata, {csim1, csim2});
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
            if (fix_delta) // variables *must* be supplied in order, so we need 2 if-else statements
                m->SetFixedVariable(1, "delta1", 0);
            else
                m->SetLimitedVariable(1, "delta1", guess[1], 0.01, 0, 1);
            m->SetLimitedVariable(2, "c", guess[2], 0.01, 0, 1); // c is the ratio of sim1 to sim2, i.e. c*sim1 + (1-c)*sim2
            m->SetLimitedVariable(3, "k2", guess[3], 0.01, 0, 1);
            if (fix_delta)
                m->SetFixedVariable(4, "delta2", 0);
            else 
                m->SetLimitedVariable(4, "delta2", guess[4], 0.01, 0, 1);
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
        file << "    Estimating errors with MINOS: " << endl;
        file << format("        k: %1% | +%2% | %3%") % res[0] % err_up[0] % err_low[0] << endl;
        file << format("        delta: %1% | +%2% | %3%") % res[1] % err_up[1] % err_low[1] << endl;
        if (fit_type == 2 || fit_type == 3) {
            file << format("        c: %1% | +%2% | %3%") % res[2] % err_up[2] % err_low[2] << endl;
        }
        if (fit_type == 3) {
            file << format("        k2: %1% | +%2% | %3%") % res[3] % err_up[3] % err_low[3] << endl;
            file << format("        delta2: %1% | +%2% | %3%") % res[4] % err_up[4] % err_low[4] << endl;
        }
        file << endl;
    };

    // print the result of a fit
    auto fit_report = [&fit_type, &res, &fix_delta, &file, &fitter] (ROOT::Math::Minimizer* m) {
        file << format("        FVAL = %1%") % m->MinValue() << endl;
        file << format("        chi2 = %1%") % fitter->chisquare(res) << endl;
        file << format("        k = %1% (FREE)") % res[0] << endl;
        file << format("        delta = %1% (%2%)") % res[1] % (fix_delta ? "FIXED" : "FREE") << endl;
        if (fit_type == 2 || fit_type == 3) {
            file << format("        c = %1% (FREE)") % res[2] << endl;
        } 
        if (fit_type == 3) {
            file << format("        k2 = %1% (FREE)") % res[3] << endl;
            file << format("        delta2 = %1% (%2%)") % res[4] % (fix_delta ? "FIXED" : "FREE") << endl;
        } 
        file << endl;
    };

    // print a bunch of information to the file
    file << "*** FIT REPORT ***" << endl;
    file << format("Type: %1%, algorithm: %2%") % type % algorithm << endl;
    file << format("Fitting %1% with %2%") % args[1] % args[2];
    if (fit_type == 2 || fit_type == 3) {
        file << " and " << args[3];
    }
    file << format("\nNumber of events in files: %1% | %2%") % data.Count().GetValue() % sim1.Count().GetValue();
    if (fit_type == 2 || fit_type == 3) {
        file << " | " << to_string(sim2.Count().GetValue());
    }
    file << format("\nTotal number of bins: %1%, where only %2% are used.") % pow(bins, 2) % fitter->get_bins() << endl;
    if (y_cut != 1) {
        file << format("Performed a cut on y = %1%") % y_cut << endl;
    }
    file << "\nROOT Minimizer report: " << endl;
    file << "    Initial guess values:" << endl;
    file << "        k = " << guess[0] << endl;
    file << "        delta = " << guess[1] << endl;
    if (fit_type == 2 || fit_type == 3) {
        file << "        c = " << guess[2] << endl;
    }
    if (fit_type == 3) {
        file << "        k2 = " << guess[3] << endl;
        file << "        delta2 = " << guess[4] << endl;
    }
    file << "\n    Results from first fit: " << endl;
    fit_report(minimizer);

    // unless the first fit was a Migrad, we perform a quick second fit with it to estimate the errors.
    auto *coutbuf = std::cout.rdbuf();
    if (algorithm != "Migrad") {
        // update our guess values
        for (int i = 0; i < pars; i++) guess[i] = res[i];

        cout << "    \nEstimating errors with Migrad" << endl;
        auto minimize2 = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimize2->SetFunction(functor);
        set_params(minimize2);
        minimize2->Minimize(); 
        res = minimize2->X();
        std::cout.rdbuf(file.rdbuf()); // redirect std::cout to fit.txt
        file << "    Results from second fit: " << endl;
        fit_report(minimize2);
        minos_errs(minimize2);
    } else {
        std::cout.rdbuf(file.rdbuf());
        fit_report(minimizer);
        minos_errs(minimizer);
    }
    std::cout.rdbuf(coutbuf); // remove the redirection

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
    sim1 = sim1.Define("w", w1, {"f", "wU"});

    bool sim2_weighted = false;
    if (fit_type == 3) {
        sim2 = sim2.Define("w", w2, {"f", "wU"});
        sim2_weighted = true;
    }

//*** DALITZ PLOTS ***//
    TCanvas* c1 = new TCanvas("c1", "c", 1200, 600);
    c1->Divide(2, 1, 0); // 0 x margin
    
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
    c1->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

//*** DALITZ DIFFERENCE PLOT ***//
    TCanvas* c2 = new TCanvas("c2", "c", 600, 600);
    dalitz_data->Scale(data_scale); // scale both the data and simulation so the z-axis makes some kind of sense
    dalitz_sim->Scale(data_scale);

    dalitz_data->Add(dalitz_sim, -1); // subtract the two plots
    dalitz_data->Draw("colz");

    path = folder + "dalitz_diff.pdf";
    c2->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

//*** RADIAL PROJECTION ***//
    vector<double> x_axis = {0, 1};
    TCanvas* c3 = new TCanvas("c3", "c", 600, 600);
    TH1D dat_rho = data.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))").Histo1D({"dat_rho", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
    TH1D sim_rho = sim1.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))").Histo1D({"sim_rho", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
    dat_rho.Scale(1./dat_rho.GetMaximum());
    sim_rho.Scale(c/sim_rho.GetMaximum());
    if (fit_type == 2 || fit_type == 3) {
        TH1D sim2_rho;
        if (fit_type == 3) {
            sim2_rho = sim2.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))")
                                .Histo1D({"sim2_rho", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
        } else {
            sim2_rho = sim2.Define("tmp", "sqrt(pow(x, 2) + pow(y, 2))")
                        .Histo1D({"sim2_rho", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
        }
        sim2_rho.Scale((1-c)/sim2_rho.GetMaximum());
        sim_rho.Add(&sim2_rho);
    }
    setup_compare_plot(&dat_rho, &sim_rho, "\\rho", "Arbitrary units");

    path = folder + "rho.pdf";
    c3->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

//*** ANGULAR PROJECTION ***//
    x_axis = {0, M_PI/3};
    TCanvas* c4 = new TCanvas("c4", "c", 600, 600);
    TH1D dat_ang = data.Define("tmp", "atan2(x, y)").Histo1D({"dat_ang", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
    TH1D sim_ang = sim1.Define("tmp", "atan2(x, y)").Histo1D({"sim_ang", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
    dat_ang.Scale(1./dat_ang.GetMaximum());
    sim_ang.Scale(c/sim_ang.GetMaximum());
    if (fit_type == 2 || fit_type == 3) {
        TH1D sim2_ang;
        if (fit_type == 3) {
            sim2_ang = sim2.Define("tmp", "atan2(x, y)").Histo1D({"sim2_ang", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
        } else {
            sim2_ang = sim2.Define("tmp", "atan2(x, y)").Histo1D({"sim2_ang", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
        }
        sim2_ang.Scale((1-c)/sim2_ang.GetMaximum());
        sim_ang.Add(&sim2_ang);
    }
    setup_compare_plot(&dat_ang, &sim_ang, "\\varphi", "Arbitrary units");

    path = folder + "phi.pdf";
    c4->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

//*** ENERGY COMPARISON ***//
    x_axis = {0, 7000};
    TCanvas* c5 = new TCanvas("c5", "c", 600, 600);
    TH1D* dat_E = new TH1D("dat_E", "h", bins, x_axis[0], x_axis[1]);
    TH1D* sim_E = new TH1D("sim_E", "h", bins, x_axis[0], x_axis[1]);
    for (int i = 0; i < 3; i++) {
        TH1D dat_temp = data.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"dat_temp", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
        TH1D sim_temp = sim1.Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"sim_temp", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
        dat_E->Add(&dat_temp);
        sim_E->Add(&sim_temp);
    }
    dat_E->Scale(1/dat_E->GetMaximum());
    sim_E->Scale(c/sim_E->GetMaximum());
    if (fit_type == 2 || fit_type == 3) {
        TH1D* sim2_E = new TH1D("sim2_E", "h", bins, x_axis[0], x_axis[1]);
        for (int i = 0; i < 3; i++) {
            TH1D sim2_temp;
            if (fit_type == 3) {
                sim2_temp= sim2.Define("tmp", (format("E_cm[%1%]") % i).str())
                                .Histo1D({"sim2_temp", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
            } else {
                sim2_temp= sim2.Define("tmp", (format("E_cm[%1%]") % i).str())
                                .Histo1D({"sim2_temp", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
            }
            sim2_E->Add(&sim2_temp);
        }
        sim2_E->Scale((1-c)/sim2_E->GetMaximum());
        sim_E->Add(sim2_E);
    }
    setup_compare_plot(dat_E, sim_E, "E_{cm}", "Arbitrary units");

    path = folder + "E_cm.pdf";
    c5->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

//*** BINWISE CHI2 DALITZ PLOT ***//
    gStyle->SetPalette(kThermometer); // change to blue/red color scheme to indicate differences
    gROOT->ForceStyle();
    TCanvas* c6 = new TCanvas("c2", "c", 600, 600);
    vector<double> binwise_chi = fitter->binwise_chisquare(res);
    vector<bool> skip_bin = fitter->get_skipped_bins();
    TH2D* chi2 = new TH2D("chi", "h", bins, 0, 1, bins, 0, 1);
    for (int y = 0; y < bins; y++) {
        for (int x = 0; x < bins; x++) {
            int my_bin = x + y*bins; // the logical choice of bins
            int root_bin = (x + 1) + (y + 1)*(bins + 2); // index 0 is underflow bin and index bins+1 is overflow bin

            if (skip_bin[my_bin]) {
                chi2->SetBinContent(root_bin, -1);
            } else {
                chi2->SetBinContent(root_bin, binwise_chi[my_bin]);
            }
            chi2->SetMinimum(0);
        }
    }
    setup_dalitz_plot(chi2); // setup axes, labels etc
    chi2->Draw("colz");

    path = folder + "binwise_chi2.pdf";
    c6->SetRightMargin(0.15);
    c6->SaveAs(path.c_str());
    cout << "Created " << path << "." << endl;

//*** CHI2 VALS FOR PROJECTION PLOTS ***//
    double chi_1 = 0, chi_2 = 0, chi_3 = 0;
    int bins1 = bins, bins2 = bins, bins3 = bins;
    for (int i = 1; i < bins; i++) {
        double mu_1 = dat_rho.GetBinContent(i); // model
        double m_1 = sim_rho.GetBinContent(i); // observed
        if (mu_1 != 0) {
            // chi_1 += 2*(mu_1 - m_1 + m_1*log(m_1/mu_1));
            chi_1 += pow(m_1 + mu_1, 2)/mu_1;
        } else {
            bins1--;
        }

        double mu_2 = dat_ang.GetBinContent(i); // model
        double m_2 = sim_ang.GetBinContent(i); // observed
        if (mu_2 != 0) {
            // chi_2 += 2*(mu_2 - m_2 + m_2*log(m_2/mu_2));
            chi_2 += pow(m_2 + mu_2, 2)/mu_2;
        } else {
            bins2--;
        }

        double mu_3 = dat_E->GetBinContent(i); // model
        double m_3 = sim_E->GetBinContent(i); // observed
        if (mu_3 != 0) {
            // chi_3 += 2*(mu_3 - m_3 + m_3*log(m_3/mu_3));
            chi_3 += pow(m_3 + mu_3, 2)/mu_3;
        } else {
            bins3--;
        }
    }
    file << format("Radial projection:  chi2 = %1%, bins = %2%") % chi_1 % bins1 << endl;
    file << format("Angular projection: chi2 = %1%, bins = %2%") % chi_2 % bins2 << endl;
    file << format("Energy projection:  chi2 = %1%, bins = %2%") % chi_3 % bins3 << endl;
    file.close();
}