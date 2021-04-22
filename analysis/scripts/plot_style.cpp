#include <TStyle.h>
#include <TROOT.h>

void setup_style() {
    gStyle->SetLineStyleString(11, "20 10"); // smaller dashes than the standard
    gStyle->SetPalette(kViridis); // set the global color scheme of figures
    gStyle->SetOptStat(0); // hide legends
    gStyle->SetOptTitle(0); // hide titles
    gROOT->ForceStyle();
}
