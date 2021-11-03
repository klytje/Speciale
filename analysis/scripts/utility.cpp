#pragma once

// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TROOT.h>
#include <TVector3.h>
#include <TF1.h>

// other stuff
#include <boost/format.hpp>

using namespace std;
using boost::format;

const static map<string, function<double(double)>> correlation_functions = {
    {"0+ 2", [] (double theta) {return 2.25*pow(cos(theta), 4) - 1.5*pow(cos(theta), 2) + 0.25;}},
    {"1- 1", [] (double theta) {return 0.75*pow(cos(theta), 2) + 0.25;}},
    {"1+ 2", [] (double theta) {return -4*pow(cos(theta), 4) + 4*pow(cos(theta), 2);}},
    {"1- 3", [] (double theta) {return 1.25*pow(cos(theta), 4) - 0.5*pow(cos(theta), 2) + 0.25;}},
    {"2- 1", [] (double theta) {return 1 - pow(cos(theta), 2);}},
    {"2+ 2", [] (double theta) {return 2.25*pow(cos(theta), 4) - 2.25*pow(cos(theta), 2) + 1;}},
    {"2- 3", [] (double theta) {return -3.529*pow(cos(theta), 4) + 3.294*pow(cos(theta), 2) + 0.235;}},
    {"3- 1", [] (double theta) {return 1./3*pow(cos(theta), 2) + 2./3;}},
    {"3+ 2", [] (double theta) {return -0.6*pow(cos(theta), 4) - 0.4*pow(cos(theta), 2) + 1.0;}},
    {"3- 3", [] (double theta) {return 2.368*pow(cos(theta), 4) - 2.526*pow(cos(theta), 2) + 1;}}
};

// returns the correlation function as a lambda for a given state and l
tuple<function<double(Double_t*, Double_t*)>, vector<double>, double> get_angular_correlation_function(string state, string l) { 
    // these two values are determined by eye for all possibilities
    vector<double> bounds; // y axis bounds
    double interference_point; // point of maximum interference
    function<double(Double_t*, Double_t*)> ang_corr;
    double E_tot;

    // these are all calculated with my own angular_correlation.py script
    if (state == "0+" && l == "2") {
        bounds = {0.35, 0.45}; // the bounds on the beam of interest
        interference_point = 0.7; // rough location of the points of interference. These will not be used for scaling the angular correlation function
        ang_corr = [] (Double_t* x, Double_t* par) {
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("0+ 2")(beta); // get the correlation function
            return maxval*corr;
        };
    } else if (state == "1-" && l == "1") {
        std::cout << "\033[1;31m" << "WARNING: Specified state is not configured. Set bounds and interference_point in ../scripts/utility.cpp" << "\033[0m" << endl;
        bounds = {0, 0};
        interference_point = 0;
        ang_corr = [] (Double_t* x, Double_t* par) {            
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("1- 1")(beta); // get the correlation function
            return maxval*corr;
        };
    } else if (state == "1+" && l == "2") {
        std::cout << "\033[1;31m" << "WARNING: Specified state is not configured. Set bounds and interference_point in ../scripts/utility.cpp" << "\033[0m" << endl;
        bounds = {0, 0};
        interference_point = 0;
        ang_corr = [] (Double_t* x, Double_t* par) {            
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("1+ 2")(beta); // get the correlation function
            return maxval*corr;
        };
    } else if (state == "1-" && l == "3") {
        std::cout << "\033[1;31m" << "WARNING: Specified state is not configured. Set bounds and interference_point in ../scripts/utility.cpp" << "\033[0m" << endl;
        bounds = {0, 0};
        interference_point = 0;
        ang_corr = [] (Double_t* x, Double_t* par) {            
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("1- 3")(beta); // get the correlation function
            return maxval*corr;
        };
    } else if (state == "2-" && l == "1") {
        bounds = {0.28, 0.37};
        interference_point = 0.56;
        ang_corr = [] (Double_t* x, Double_t* par) {            
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("2- 1")(beta); // get the correlation function
            return maxval*corr;
        };
    } else if (state == "2+" && l == "2") {
        std::cout << "\033[1;31m" << "WARNING: Specified state is not configured. Set bounds and interference_point in ../scripts/utility.cpp" << "\033[0m" << endl;
        bounds = {0, 0};
        interference_point = 0;
        ang_corr = [] (Double_t* x, Double_t* par) {            
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("2+ 2")(beta); // get the correlation function
            return maxval*corr;
        };
    } else if (state == "2-" && l == "3") {
        std::cout << "\033[1;31m" << "WARNING: Specified state is not configured. Set bounds and interference_point in ../scripts/utility.cpp" << "\033[0m" << endl;
        bounds = {0, 0};
        interference_point = 0;
        ang_corr = [] (Double_t* x, Double_t* par) {            
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("2- 3")(beta); // get the correlation function
            return maxval*corr;
        };
    } else if (state == "3-" && l == "1") {
        bounds = {0.4, 0.47};
        interference_point = 0.78;
        ang_corr = [] (Double_t* x, Double_t* par) {
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("3- 1")(beta); // get the correlation function
            return maxval*corr;
        };
    } else if (state == "3+" && l == "2") {
        std::cout << "\033[1;31m" << "WARNING: Specified state is not configured. Set bounds and interference_point in ../scripts/utility.cpp" << "\033[0m" << endl;
        bounds = {0, 0};
        interference_point = 0;
        ang_corr = [] (Double_t* x, Double_t* par) {
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("3+ 2")(beta); // get the correlation function
            return maxval*corr;
        };
    } else if (state == "3-" && l == "3") {
        std::cout << "\033[1;31m" << "WARNING: Specified state is not configured. Set bounds and interference_point in ../scripts/utility.cpp" << "\033[0m" << endl;
        bounds = {0, 0};
        interference_point = 0;
        ang_corr = [] (Double_t* x, Double_t* par) {            
            double maxval = par[0];
            double y = par[1];
            
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta; // beta = 180 - theta
            double corr = correlation_functions.at("3- 3")(beta); // get the correlation function
            return maxval*corr;
        };
    } else {
        cout << format("The given state %1% l = %2% is not specified yet! \nYou can do it pretty easily yourself in ../scripts/utility.cpp") % state % l << endl;
        exit(1);
    }
    return make_tuple(ang_corr, bounds, interference_point);
}

// define sorting methods
const auto emax = [] (double e1, double e2, double e3) {return std::max({e1, e2, e3});};
const auto emin = [] (double e1, double e2, double e3) {return std::min({e1, e2, e3});};
const auto emid = [] (double e1, double e2, double e3) {
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

// weights defined by eq 42 in Morten's thesis
double calc_weight(vector<vector<double>> f, double wU, double k, double delta) { 
    return wU*(k*f[0][0]+(1-k)*f[0][1] + 2*sqrt(k*(1-k))*(f[0][2]*cos(delta) + f[0][3]*sin(delta)));
}

// this filter is supposed to be applied before any figure is made. I've placed it here so it's easy to change later (or even remove completely)
void filter(ROOT::RDF::RNode* df) {
    *df = df->Filter("abs(deltaE) < 300")
            .Filter("p_tot < 40e3");
    return;
}

// Defines the Dalitz coordinates x and y
void setup_dataframe(ROOT::RDF::RNode* df) {
    *df = df->Define("E_tot","E_cm[0] + E_cm[1] + E_cm[2]")
            .Define("e_cm_1", "E_cm[0]/E_tot") // normalized such that e1 + e2 + e3 = 1
            .Define("e_cm_2", "E_cm[1]/E_tot")
            .Define("e_cm_3", "E_cm[2]/E_tot")
            .Define("e_1", emax, {"e_cm_1", "e_cm_2", "e_cm_3"}) // we want e1 > e2 > e3
            .Define("e_2", emid, {"e_cm_1", "e_cm_2", "e_cm_3"})
            .Define("e_3", emin, {"e_cm_1", "e_cm_2", "e_cm_3"})
            .Define("x","sqrt(3)*(e_2 - e_3)")
            .Define("y","3*e_1 - 1");
    return;
}

// cuts a Dalitz plot at y = 0.93 or whatever is input
void cut_y(ROOT::RDF::RNode* df, double y = 0.93) {
    *df = df->Filter("y < " + to_string(y));
}

// removes any events outside a circle of radius 1
void cut_circle(ROOT::RDF::RNode* df) {
    *df = df->Filter("pow(x, 2) + pow(y, 2) < 1.0");
}

const double m_alpha = 3.72737*1e6; // in keV
const auto E23 = [] (double r2x, double r2y, double r2z, double r3x, double r3y, double r3z) {
    TVector3 p2 = {r2x, r2y, r2z};
    TVector3 p3 = {r3x, r3y, r3z};
    TVector3 p23 = p2-p3;
    return 1./(4*m_alpha)*p23.Mag2();
};

const auto Emax = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
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

const auto Emin = [] (ROOT::VecOps::RVec<Double_t> x, ROOT::VecOps::RVec<Double_t> y, ROOT::VecOps::RVec<Double_t> z) {
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

// fit and remove the ground state decay events
void cut_gs_alt(ROOT::RDF::RNode* df) {
    double E_exc[] = {0, 200};
        *df = df->Define("i_max", Emax, {"px", "py", "pz"})
            .Define("i_min", Emin, {"px", "py", "pz"})
            .Define("i_mid", "3 - i_max - i_min") // i_min + i_mid + i_max = 3
            .Define("px2", "px[i_mid]")
            .Define("py2", "py[i_mid]")
            .Define("pz2", "pz[i_mid]")
            .Define("px3", "px[i_min]")
            .Define("py3", "py[i_min]")
            .Define("pz3", "pz[i_min]")
            .Define("E_23", E23, {"px2", "py2", "pz2", "px3", "py3", "pz3"});
    *df = df->Filter("200 < E_23");
}

// fit and extract only the excited decay events
void cut_gs(ROOT::RDF::RNode* df) {
    double E_exc[] = {1000, 4500};
        *df = df->Define("i_max", Emax, {"px", "py", "pz"})
            .Define("i_min", Emin, {"px", "py", "pz"})
            .Define("i_mid", "3 - i_max - i_min") // i_min + i_mid + i_max = 3
            .Define("px2", "px[i_mid]")
            .Define("py2", "py[i_mid]")
            .Define("pz2", "pz[i_mid]")
            .Define("px3", "px[i_min]")
            .Define("py3", "py[i_min]")
            .Define("pz3", "pz[i_min]")
            .Define("E_23", E23, {"px2", "py2", "pz2", "px3", "py3", "pz3"});
    TH1D hist = (*df).Histo1D({"cut_gs", "h", 200, E_exc[0], E_exc[1]}, "E_23").GetValue();
    TF1* tf_ex = new TF1("tf_ex", "gaus", E_exc[0], E_exc[1]);
    hist.Fit(tf_ex, "LQR+");

    double mu = tf_ex->GetParameter(1);
    double sigma = tf_ex->GetParameter(2);
    *df = df->Filter((format("%1% < E_23 && E_23 < %2%") % (mu-3*sigma) % (mu+3*sigma)).str());

    // debug plot
    // TCanvas* c1 = new TCanvas("c1", "c", 600, 600);
    // TLine* l1 = new TLine(mu-3*sigma, 0, mu-3*sigma, hist.GetMaximum());
    // TLine* l2 = new TLine(mu+3*sigma, 0, mu+3*sigma, hist.GetMaximum());
    // hist.Draw();
    // tf_ex->Draw("same");
    // l1->Draw();
    // l2->Draw();
    // c1->SaveAs("tmp.pdf");
}

// makes a Dalitz plot from a dataframe. it assumes the Dalitz coordinates x and y are already defined. 
// if bounded, it cuts away all events outside the unit circle. if weighted, it uses the w column as weights
TH2D* dalitz(ROOT::RDF::RNode* df, const int bins = 200, const bool bounded = false, const bool weighted = false) {
    if (!(df->HasColumn("x") && df->HasColumn("y"))) {
        std::cout << "\033[1;31m" << "ERROR: The Dalitz coordinates x and y are not defined in the input dataframe!" << "\033[0m" << endl;
        exit(1);
    } 
    if (weighted && !df->HasColumn("w")) {
        std::cout << "\033[1;31m" << "ERROR: weighted was specified, but no w is defined in the input dataframe!" << "\033[0m" << endl;
        exit(1);
    } 

    double lim = bounded ? 1 : 1.3;
    vector<double> x_axis = {double(bins), -lim, lim};
    vector<double> y_axis = {double(bins), -lim, lim};
    TH2D* hist = new TH2D("hist", "h", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    int perms[] = {1, 2, 3}; // we can get the other slices simply by permutating i, j, k
    do {
        int i = perms[0];
        int j = perms[1];
        int k = perms[2];
        TH2D htemp;
        if (weighted) {
            htemp = df->Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                        .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                        .Histo2D({"h1", "temp", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp", "w").GetValue();
        } else {
            htemp = df->Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                        .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                        .Histo2D({"h1", "temp", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp").GetValue();
        }
        hist->Add(&htemp);
    } while (std::next_permutation(perms, perms+3)); // repeat for each of the 3! = 6 permutations of {1, 2, 3}
    return hist;
}

// makes a complete Dalitz plot based on the input file, including various cuts
TH2D* dalitz(const char* file, int bins = 200, bool bounded = false, bool remove_gs = true) {
    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", file);
    filter(&df);
    setup_dataframe(&df);
    if (bounded) {
        cut_circle(&df);
    }
    if (remove_gs) {
        cut_gs(&df);
    }

    return dalitz(&df, bins, bounded);
}

// makes a Dalitz plot from vectors. this is slightly more complicated since we do not have access to the energies
TH2D* dalitz(const vector<double> x, const vector<double> y, const int bins = 200, vector<double> w = {}) {
    // if no weights are supplied, create a new array filled with 1s
    if (w.size() == 0) {
        w = vector<double>(x.size(), 1);
    }
    TH2D* hist = new TH2D("hist", "h", bins, -1, 1, bins, -1, 1);
    const vector<double> cost = {cos(0), cos(2./3*M_PI), cos(4./3*M_PI)};
    const vector<double> sint = {sin(0), sin(2./3*M_PI), sin(4./3*M_PI)};

    // in this first loop, each event is mirrored across the radial line at 60 degrees
    for (int i = 0; i < x.size(); i++) {
        double phi = atan2(y[i], x[i]);
        double X = x[i]*cos(2*phi-M_PI/3) + y[i]*sin(2*phi-M_PI/3);
        double Y = -x[i]*sin(2*phi-M_PI/3) + y[i]*cos(2*phi-M_PI/3);

        // in this second loop, both the original and mirrored event are duplicated and rotated by 120 degrees
        // the rotation matrix is precalculated for efficiency
        for (int r = 0; r < 3; r++) {
            double x1 = x[i]*cost[r] + y[i]*sint[r];
            double x2 = X*cost[r] + Y*sint[r];
            double y1 = -x[i]*sint[r] + y[i]*cost[r];
            double y2 = -X*sint[r] + Y*cost[r];
            hist->Fill(x1, y1, w[i]);
            hist->Fill(x2, y2, w[i]);
        }
    }
    return hist;
}

// makes a Dalitz slice plot from vectors
TH2D* dalitz_slice(const vector<double> x, const vector<double> y, const int bins = 100, vector<double> w = {}) {
    // if no weights are supplied, create a new array filled with 1s
    if (w.size() == 0) {
        w = vector<double>(x.size(), 1);
    }
    TH2D* hist = new TH2D("hist", "h", bins, 0, 1, bins, 0, 1);
    for (int i = 0; i < x.size(); i++) {
        hist->Fill(x[i], y[i], w[i]);
    }
    return hist;
}


// makes a Dalitz slice from a given file input
TH2D* dalitz_slice(ROOT::RDF::RNode* df, const int bins = 200, const bool bounded = false, const bool weighted = false) {
    if (!(df->HasColumn("x") && df->HasColumn("y"))) {
        std::cout << "\033[1;31m" << "ERROR: The Dalitz coordinates x and y are not defined in the input dataframe!" << "\033[0m" << endl;
        exit(1);
    } 
    if (weighted && !df->HasColumn("w")) {
        std::cout << "\033[1;31m" << "ERROR: weighted was specified, but no w is defined in the input dataframe!" << "\033[0m" << endl;
        exit(1);
    } 

    double lim = bounded ? 1 : 1.3;
    vector<double> x_axis = {double(bins), 0, lim};
    vector<double> y_axis = {double(bins), 0, lim};
    TH2D tmp; // .GetPtr() and direct pointer assignment (*hist = ...) crashes ROOT, so we have to do it in this roundabout fashion
    if (weighted) {
        tmp = df->Histo2D({"tmp", "h", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x", "y", "w").GetValue();
    } else {
        tmp = df->Histo2D({"tmp", "h", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x", "y").GetValue();
    }
    TH2D* hist = new TH2D("hist", "h", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    hist->Add(&tmp);
    return hist;
}

/**
 *  Create a Dalitz slice for the input file.
 * @param file The path to the target file.
 * @param bins The number of bins used for both the x and y axes.
 * @param cut Whether to impose cuts on the plot or not.
 * @return A TH2D* containing the Dalitz slice
 */
TH2D* dalitz_slice(const char* file, const int bins = 100, bool cut = true) {
    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", file);
    filter(&df);
    setup_dataframe(&df);
    if (cut) {
        cut_circle(&df);
        cut_gs(&df);
    }

    return dalitz_slice(&df, bins, cut);
}

/**
    Energy projection for a single simulation or data file. 
    @param df The input dataframe
    @param weighted Determines whether each bin will be weighted by the "w" column, which must thus already be defined. 
    @param name The name used internally in ROOT. Change this when creating multiple projections. 
    @param scale The scaling imposed on the histogram. Options: "max", "integral", "none"
    @return A TH1D* ROOT histogram containing the energy projection. 
*/
TH1D* energy_projection(ROOT::RDF::RNode* df, bool weighted = false, string name = "proj", string scale = "none") {
    if (weighted && !df->HasColumn("w")) {
        std::cout << "\033[1;31m" << "ERROR: weighted was specified, but no w is defined in the input dataframe!" << "\033[0m" << endl;
        exit(1);
    } 
    vector<double> x_axis = {0, 7000};
    int bins = 100;

    TH1D* proj = new TH1D(name.c_str(), "h", bins, x_axis[0], x_axis[1]);
    for (int i = 0; i < 3; i++) {
        TH1D temp;
        if (weighted) {
            temp = df->Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"temp", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
        } else {
            temp = df->Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"temp", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
        }
        proj->Add(&temp);
    }
    if (scale == "max") {
        proj->Scale(1/proj->GetMaximum());
    } else if (scale == "integral") {
        proj->Scale(1/proj->Integral());
    }
    return proj;
}

/**
    Energy projection for two simulation files. 
    @param df The two dataframes containing the simulated events.
    @param ratio The fitted ratio between the two components. 
    @param weighted Determines whether each bin will be weighted by the "w" column, which must thus already be defined. 
    @param name The name used internally in ROOT. Change this when creating multiple projections. 
    @param scale The scaling imposed on the histogram. Options: "max", "integral", "none"
    @return A TH1D* ROOT histogram containing the energy projection. 
*/
TH1D* energy_projection(ROOT::RDF::RNode* df1, ROOT::RDF::RNode* df2, double ratio, bool weighted1 = false, bool weighted2 = false, string name = "proj", string scale = "none") {
    if (weighted1 && !df1->HasColumn("w") || weighted2 && !df2->HasColumn("w")) {
        std::cout << "\033[1;31m" << "ERROR: weighted was specified, but no w is defined in the input dataframe!" << "\033[0m" << endl;
        exit(1);
    } 
    vector<double> x_axis = {0, 7000};
    int bins = 100;

    TH1D* proj = new TH1D(name.c_str(), "h", bins, x_axis[0], x_axis[1]);
    TH1D* proj2 = new TH1D("proj2", "h", bins, x_axis[0], x_axis[1]);
    for (int i = 0; i < 3; i++) {
        TH1D temp1, temp2;
        if (weighted1) {
            temp1 = df1->Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"temp", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
        } else {
            temp1 = df1->Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"temp", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
        }
        if (weighted2) {
            temp2 = df2->Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"temp", "h", bins, x_axis[0], x_axis[1]}, "tmp", "w").GetValue();
        } else {
            temp2 = df2->Define("tmp", (format("E_cm[%1%]") % i).str()).Histo1D({"temp", "h", bins, x_axis[0], x_axis[1]}, "tmp").GetValue();
        }
        proj->Add(&temp1);
        proj2->Add(&temp2);
    }

    proj->Scale(ratio/proj->Integral());
    proj2->Scale((1-ratio)/proj2->Integral());
    double data_scale = proj->GetMaximum();
    proj->Scale(1/data_scale);
    proj2->Scale(1/data_scale);
    proj->Add(proj2);

    if (scale == "max") {
        proj->Scale(1/proj->GetMaximum());
    } else if (scale == "integral") {
        proj->Scale(1/proj->Integral());
    }
    return proj;
}

// I ended up needing these specific plot options quite a lot, so I relocated it here as a common method. note that it also draws the histograms
void setup_compare_plot(TH1D* hdat, TH1D* hsim, string xlabel, string ylabel) {
    hdat->GetXaxis()->SetTitle(xlabel.c_str());
    hdat->GetYaxis()->SetTitle(ylabel.c_str());
    hdat->GetXaxis()->CenterTitle();
    hdat->GetYaxis()->CenterTitle();
    hdat->GetXaxis()->SetNdivisions(205);
    hdat->GetYaxis()->SetNdivisions(203);
    
    hsim->SetLineColor(kOrange+1);
    hdat->SetLineColor(kBlack);
    hsim->SetLineWidth(2);
    hdat->SetLineWidth(2);
    hdat->Draw("HIST L");
    hsim->Draw("HIST L SAME");
}

// I also ended up needing this quite a bit
void setup_dalitz_plot(TH2D* h, string draw_options = "") {
    h->GetXaxis()->SetTitle("x");
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetNdivisions(2);
    h->GetYaxis()->SetTitle("y");
    h->GetYaxis()->CenterTitle();
    h->GetYaxis()->SetNdivisions(2);

    if (draw_options == "") {
        draw_options = "colz";
    }
    h->Draw(draw_options.c_str());
}