#include <TStyle.h>
#include <TROOT.h>

void setup_style() {
    gStyle->SetPalette(kViridis);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gROOT->ForceStyle();
}