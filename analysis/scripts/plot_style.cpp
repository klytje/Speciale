#include <TStyle.h>
#include <TROOT.h>

void setup_style() {
    gStyle->SetPalette(kViridis); // set the global color scheme of figures
    gStyle->SetOptStat(0); // hide legends
    gStyle->SetOptTitle(0); // hide titles
    gROOT->ForceStyle();
}
