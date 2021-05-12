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

void cut_edges(ROOT::RDF::RNode* df) {
    *df = df->Filter("pow(x,2) + pow(y,2) < 1.0")
            .Filter("y < 0.93");
}

const double m_alpha = 3.72737*1e6;
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

// fit and extract only the excited decay events
void cut_gs(ROOT::RDF::RNode* df) {
    double E_exc[] = {1500, 4500};
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
    TH1D hist = (*df).Histo1D({"h2", "h", 200, E_exc[0], E_exc[1]}, "E_23").GetValue();
    TF1* tf_ex = new TF1("tf_ex", "gaus", E_exc[0], E_exc[1]);
    hist.Fit(tf_ex, "LQR+");

    // TCanvas* c = new TCanvas("c", "c", 600, 600);
    // hist.Draw();
    // string path = "tmp.pdf";
    // c->SaveAs(path.c_str());

    double mu = tf_ex->GetParameter(1);
    double sigma = tf_ex->GetParameter(2);
    *df = df->Filter((format("%1% < E_23 && E_23 < %2%") % (mu-3*sigma) % (mu+3*sigma)).str());
}

// I ended up needing this in a whole lot of different files, so I moved it here to a common #include file instead
TH2D* dalitz(const char* file, int bins = 200, bool bounded = false) {
    // set the axes
    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", file);
    filter(&df);
    setup_dataframe(&df);

    vector<double> x_axis;
    vector<double> y_axis;
    if (bounded) {
        x_axis = {double(bins), -1, 1};
        y_axis = {double(bins), -1, 1};
        cut_edges(&df);
    } else {
        x_axis = {double(bins), -1.3, 1.3};
        y_axis = {double(bins), -1.3, 1.3};
    }

    TH2D* hist = new TH2D("h1", "Dalitz plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    int perms[] = {1, 2, 3}; // we can get the other slices simply by permutating i, j, k
    do {
        int i = perms[0];
        int j = perms[1];
        int k = perms[2];
        TH2D htemp = df.Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                    .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                    .Histo2D({"h1", "temp", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp").GetValue();
        hist->Add(&htemp);
    } while (std::next_permutation(perms, perms+3)); // repeat for each of the 3! = 6 permutations of {1, 2, 3}
    return hist;
}

TH2D* dalitz(const vector<double> x, const vector<double> y, const int bins = 200, vector<double> w = {}) {
    // if no weights are supplied, create a new array filled with 1s
    if (w.size() == 0) {
        w = vector<double>(x.size(), 1);
    }
    TH2D* hist = new TH2D("h1", "Dalitz plot", bins, -1, 1, bins, -1, 1);
    const vector<double> cost = {cos(0), cos(2./3*M_PI), cos(4./3*M_PI)};
    const vector<double> sint = {sin(0), sin(2./3*M_PI), sin(4./3*M_PI)};
    for (int i = 0; i < x.size(); i++) {
        double phi = atan2(y[i], x[i]);
        double X = x[i]*cos(2*phi-M_PI/3) + y[i]*sin(2*phi-M_PI/3);
        double Y = -x[i]*sin(2*phi-M_PI/3) + y[i]*cos(2*phi-M_PI/3);
        for (int r = 0; r < 3; r++) {
            double x1 = x[i]*cost[r] + y[i]*sint[r];
            double x2 = X*cost[r] + Y*sint[r];
            double y1 = -x[i]*sint[r] + y[i]*cost[r];
            double y2 = -X*sint[r] + Y*cost[r];
            hist->Fill(x1, y1, w[i]);
            hist->Fill(x2, y2, w[i]);
        }
    }
    for (int r = 0; r < 6; r++) {
        // precalculate cos(theta) and sin(theta) to reduce computation time
        for (int i = 0; i < x.size(); i++) {
        }
    }
    return hist;
}

TH2D* dalitz_slice(const vector<double> x, const vector<double> y, const int bins = 100, vector<double> w = {}) {
    // if no weights are supplied, create a new array filled with 1s
    if (w.size() == 0) {
        w = vector<double>(x.size(), 1);
    }
    TH2D* hist = new TH2D("h1", "Dalitz plot", bins, -1, 1, bins, -1, 1);
    for (int i = 0; i < x.size(); i++) {
        hist->Fill(x[i], y[i], w[i]);
    }
    return hist;
}

TH2D dalitz_slice(const char* file, const int bins = 100, bool bounded = true) {
    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", file);
    filter(&df);
    setup_dataframe(&df);

    vector<double> x_axis;
    vector<double> y_axis;
    if (bounded) {
        x_axis = {double(bins), 0, 1};
        y_axis = {double(bins), 0, 1};
        cut_edges(&df);
    } else {
        x_axis = {double(bins), 0, 1.3};
        y_axis = {double(bins), 0, 1.3};
    }

    TH2D hist = df.Histo2D({"h2", "Dalitz slice", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x", "y").GetValue();
    return hist;
}

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