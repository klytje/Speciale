// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>

// other stuff
#include <boost/format.hpp>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[]) {
    setup_style();
    TCanvas* c = new TCanvas("c", "c", 900, 600);
    double y = 0;
    double maxval = 1;
    vector<double> x_bounds = {-sqrt(1-pow(y, 2)), sqrt(1-pow(y, 2))};

    TF1 *dummy = new TF1("dummy", "1", -1, 1);
    dummy->SetTitle("Angular correlation functions");
    dummy->GetXaxis()->SetNdivisions(3);
    dummy->GetXaxis()->SetTitle("x");
    dummy->GetXaxis()->CenterTitle();

    dummy->GetYaxis()->SetNdivisions(2);
    dummy->GetYaxis()->SetTitle("Angular correlation");
    dummy->GetYaxis()->CenterTitle();

    dummy->SetMaximum(1.2);
    dummy->SetMinimum(0);
    dummy->SetLineColor(kBlack);
    dummy->Draw();

    // plot multiple correlation functions
    vector<string> states = {"0+", "1-", "1+", "1-", "2-", "2+", "2-", "3-", "3+", "3-"}; // parity doesn't actually matter
    vector<string> ls =     {"2" , "1" , "2" , "3" , "1" , "2" , "3" , "1" , "2" , "3"};
    vector<int> color = {kBlack, kPink-8, kOrange+1, kYellow-7, kSpring-1, kTeal+3, kCyan+1, kAzure+1, kBlue-9, kViolet+1};
    function<double(Double_t*, Double_t*)> ang_corr;

    auto legend = new TLegend(0.1, 0.77, 0.9, 0.9);
    legend->SetNColumns(7);
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
    c->SaveAs(path.c_str());
}