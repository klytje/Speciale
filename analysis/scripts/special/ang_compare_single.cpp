// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>

// other stuff
#include <boost/format.hpp>
#include <iostream>
#include <math.h>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[]) {
    if (argc != 5) {
        cout << "Usage: ./ang_compare_single <output path> <nuclear state> <l> <file>" << endl;
        exit(1);
    }

    function<double(Double_t*, Double_t*)> ang_corr; // angular correlation function
    vector<double> bounds; // y axis bounds
    double interference_point; // point of maximum interference
    tie(ang_corr, bounds, interference_point) = get_angular_correlation_function(argv[2], argv[3]);
    
    setup_style();
    int bins = 200;
    TCanvas* c = new TCanvas("c", "c", 600, 600);
    TH2D* hist = dalitz(argv[4], bins);

    hist->GetXaxis()->SetTitle("x");
    hist->GetYaxis()->SetTitle("y");
    hist->Draw("colz");

    TLine* bot = new TLine(-1, bounds[0], 1, bounds[0]);
    TLine* top = new TLine(-1, bounds[1], 1, bounds[1]);
    top->SetLineWidth(2);
    bot->SetLineWidth(2);
    top->SetLineColor(kRed);
    bot->SetLineColor(kRed);
    top->Draw();
    bot->Draw();

    string path = string(argv[1]) + "ang_compare_raw.pdf";
    c->SetLogz();
    c->SetRightMargin(0.15);
    c->SaveAs(path.c_str());

    //*** PROJECTION ***//
    TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
    int biny1 = 0;
    int biny2 = 0;
    int b = 0;
    for (double i = -1; i < 1; i+=2./bins) {
        b++;
        if (bounds[0] < i && biny1 == 0) biny1 = b;
        if (bounds[1] < i && biny2 == 0) biny2 = b;
    }
    cout << format("bins: %1%, %2%") % biny1 % biny2 << endl;
    TH1D* h1 = hist->ProjectionX("px", biny1, biny2);
    h1->Draw();

    // angular correlation calculated with my own tool
    double d = sqrt(1 - pow((bounds[1] + bounds[0])/2, 2)); // distance across the circle at y = mean(y1, y2)

    // we need to determine the max value of the histogram, but it cannot be at the interference point (about 0.6)
    int maxval = 0;
    double location = 0;
    double width = 0.1;
    double bad_area[4] = {-(interference_point+width), -(interference_point-width), interference_point-width, interference_point+width};
    for (int i = 0; i < bins; i++) {
        int count = h1->GetBinContent(i);
        double loc = h1->GetBinCenter(i);
        if (maxval < count) {
            if (!(bad_area[0] < loc && loc < bad_area[1]) && !(bad_area[2] < loc && loc < bad_area[3])) {
                maxval = count;
                location = loc;
            }
        }
    }
    cout << format("Maximum histogram value found at bin %1%") % location << endl;
    TLine* l1 = new TLine(bad_area[0], 0, bad_area[0], maxval);
    TLine* l2 = new TLine(bad_area[1], 0, bad_area[1], maxval);
    TLine* l3 = new TLine(bad_area[2], 0, bad_area[2], maxval);
    TLine* l4 = new TLine(bad_area[3], 0, bad_area[3], maxval);
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);
    l3->SetLineWidth(2);
    l4->SetLineWidth(2);

    l1->SetLineColor(kRed);
    l2->SetLineColor(kRed);
    l3->SetLineColor(kRed);
    l4->SetLineColor(kRed);

    l1->Draw();
    l2->Draw();
    l3->Draw();
    l4->Draw();

    TF1 *corr = new TF1("corr", ang_corr, -1, 1, 2);
    corr->SetParameter(0, d);
    corr->SetParameter(1, maxval);
    corr->Draw("same");

    path = string(argv[1]) + "ang_compare_projection.pdf";
    c2->SaveAs(path.c_str());
}