// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TROOT.h>

// other stuff
#include <boost/format.hpp>

using namespace std;
using boost::format;

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

            double corr = 2.25*pow(cos(theta), 4) - 1.5*pow(cos(theta), 2) + 0.25; // the correlation function
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
            double corr = 0.75*pow(cos(theta), 2) + 0.25;
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
            double corr = -4*pow(cos(theta), 4) + 4*pow(cos(theta), 2);
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
            double corr = 1.25*pow(cos(theta), 4) - 0.5*pow(cos(theta), 2) + 0.25;
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
            double corr = 1 - pow(cos(theta), 2);
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
            double corr = 2.25*pow(cos(theta), 4) - 2.25*pow(cos(theta), 2) + 1;
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
            double corr = -3.529*pow(cos(theta), 4) + 3.294*pow(cos(theta), 2) + 0.235;
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

            double corr = 1./3*pow(cos(theta), 2) + 2./3;
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
            double corr = 2.368*pow(cos(theta), 4) - 2.526*pow(cos(theta), 2) + 1;
            return maxval*corr;
        };
    } else {
        cout << format("The given state %1% l = %2% is not specified yet! \nYou can do it pretty easily yourself in ../scripts/utility.cpp") % state % l << endl;
        exit(1);
    }
    return make_tuple(ang_corr, bounds, interference_point);
}

// I ended up needing this in a whole lot of different files, so I moved it here to a common #include file instead
TH2D* dalitz(const char* file, int bins = 200, bool show_borders = false) {
    // set the axes
    vector<double> x_axis;
    vector<double> y_axis;
    if (show_borders) {
        x_axis = {double(bins), -1.3, 1.3};
        y_axis = {double(bins), -1.3, 1.3};
    } else {
        x_axis = {double(bins), -1, 1};
        y_axis = {double(bins), -1, 1};
    }

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

    ROOT::RDF::RNode df = ROOT::RDataFrame("tree", file);
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
    return hist;
}