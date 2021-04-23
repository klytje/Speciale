// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>

// other stuff
#include <boost/format.hpp>
#include <iostream>
#include <math.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[])
{
    TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
    double y = 0;
    double maxval = 1;
    vector<double> x_bounds = {-sqrt(1-pow(y, 2)), sqrt(1-pow(y, 2))};

    TF1 *dummy = new TF1("dummy", "1", -1, 1);
    dummy->SetTitle("Angular correlation functions");
    dummy->SetMaximum(2);
    dummy->SetMinimum(0);
    dummy->SetLineColor(kBlack);
    dummy->Draw();

    // plot multiple correlation functions
    vector<string> states = {"0+", "1-", "1-", "2-", "2-", "3-", "3-"}; // parity doesn't actually matter
    vector<string> ls =     {"2" , "1" , "3" , "1" , "3" , "1" , "3"};
    vector<int> color = {kRed, kYellow, kGreen, kViolet, kRose, kSienna, kRust};
    function<double(Double_t*, Double_t*)> ang_corr;

    auto legend = new TLegend(0.1, 0.6, 0.4, 0.9);
    for (int i = 0; i < states.size(); i++) {
        tie(ang_corr, std::ignore, std::ignore) = get_angular_correlation_function(states[i], ls[i]);

        string label = (format("%1% %2%") % states[i] % ls[i]).str();
        TF1 *corr = new TF1(label.c_str(), ang_corr, x_bounds[0], x_bounds[1], 2);
        corr->SetParameter(0, maxval);
        corr->SetParameter(1, y);
        corr->SetLineColor(color[i]);
        corr->DrawClone("same");

        legend->AddEntry(label.c_str(), label.c_str(), "l");
    }
    legend->Draw();
    
    string path = string(argv[1]) + "corr_funcs.pdf";
    c2->SaveAs(path.c_str());
}