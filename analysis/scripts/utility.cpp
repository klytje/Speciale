// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TROOT.h>

// other stuff
#include <boost/format.hpp>

using namespace std;
using boost::format;

function<double(double, double)> calc_theta(string state) {
    double E_tot;
    if (state == "0+") {
        E_tot = (17.76 - 10.39) + (10.39 - 7.27);
    } else if (state == "2-") {
        E_tot = (16.67 - 10.39) + (10.39 - 7.27);
    } else if (state == "3-") {
        E_tot = (18.350 - 10.39) + (10.39 - 7.27);
    } else {
        cout << format("The given state %1% is not specified yet! \nYou can do it pretty easily yourself in ../scripts/utility.cpp") % state << endl;
        exit(1);
    }
    return [&E_tot] (double x, double y) {
        double Q1 = (y+1)*E_tot/2;
        double Q2 = E_tot - Q1;
        return acos(x*E_tot/(2*sqrt(Q1*Q2)));
    };
}

// returns the correlation function for a given state and l
tuple<function<double(Double_t*, Double_t*)>, vector<double>, double> get_angular_correlation_function(string state, string l) { 
    // these two values are determined by eye for all possibilities
    vector<double> bounds; // y axis bounds
    double interference_point; // point of maximum interference
    function<double(Double_t*, Double_t*)> ang_corr;

    // these are all calculated with my own angular_correlation.py script
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
        // ang_corr = [] (Double_t* x, Double_t* par) {
        //     double d = par[0]; 
        //     double maxval = par[1];
        //     double xp = M_PI/2*(1 + x[0]/d); // beta = 180 - theta', and theta' = pi/2(1 - x/d)
        //     double res = 2.25*pow(cos(xp), 4) - 1.5*pow(cos(xp), 2) + 0.25;
        //     return maxval*res;
        // };
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
        bounds = {0.4, 0.47};
        interference_point = 0.78;
        ang_corr = [] (Double_t* x, Double_t* par) {
            double d = par[0]; 
            double maxval = par[1];
            double xp = M_PI/2*(1 + x[0]/d); // beta = 180 - theta', and theta' = pi/2(1 - x/d)
            double res = 1./3*pow(cos(xp), 2) + 2./3;
            return maxval*res;
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