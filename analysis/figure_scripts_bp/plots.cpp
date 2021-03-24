// ROOT stuff
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <TROOT.h>
#include <TStyle.h>

// other stuff
#include <filesystem>
#include <boost/format.hpp>
#include <iostream>

using namespace std;
using namespace ROOT;
using boost::format;

const static map<string, int> plot_map = {{"dalitz_fig", 1}, 
                                          {"dalitz_forfit", 2}, 
                                          {"E_cm", 3}, 
                                          {"Excitation_C12", 4},
                                          {"fynbo_plot", 5},
                                          {"kinematic", 6},
                                          {"kinematic_lab", 7},
                                          {"phi", 8},
                                          {"rho", 9}};

void dalitz_fig(RDataFrame* df) {}
void dalitz_forfit(RDataFrame* df) {}
void E_cm(RDataFrame* df) {}
void Excitation_C12(RDataFrame* df) {}
void fynbo_plot(RDataFrame* df) {}
void kinematic(RDataFrame* df) {}
void kinematic_lab(RDataFrame* df) {}
void phi(RDataFrame* df) {}
void rho(RDataFrame* df) {}

int main(int argc, char *argv[]) {
    // prepare the data and plot type
    string plot = argv[1]; // first argument is the plot type
    RDataFrame df("tree", argv[2]);

    // define the plot colour scheme
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gROOT->ForceStyle();

    switch (plot_map.at(plot)) {
        case 1:

    }
}