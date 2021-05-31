#include <TStyle.h>
#include <TROOT.h>

void setup_style() {
    gStyle->SetLabelSize(0.05, "X");
    gStyle->SetLabelSize(0.05, "Y");
    gStyle->SetLabelSize(0.05, "Z");
    gStyle->SetTitleSize(0.07, "X");
    gStyle->SetTitleSize(0.07, "Y");
    gStyle->SetTitleOffset(0.65, "X");
    gStyle->SetTitleOffset(0.65, "Y");
    // gStyle->SetTitleXSize(0.04);
    // gStyle->SetTitleYSize(0.04);
    // gStyle->SetTickLength(0.04);
    gStyle->SetLineStyleString(11, "20 10"); // smaller dashes than the standard
    gStyle->SetLineStyleString(12, "20 20"); // smaller, more spread out dashes than the standard
    gStyle->SetPalette(kViridis); // set the global color scheme of figures
    gStyle->SetOptStat(0); // hide legends
    gStyle->SetOptTitle(0); // hide titles
    gROOT->ForceStyle();
}
