// ROOT stuff
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>

// other stuff
#include <boost/format.hpp>
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
    string state = argv[2];
    string l = argv[3];
    bool draw_excluded_area = true; // draw the areas which are excluded from the maximum value search (used to scale correlation function)

    function<double(Double_t*, Double_t*)> ang_corr; // angular correlation function
    vector<double> y_bounds; // y axis bounds
    double interference_point; // point of maximum interference
    tie(ang_corr, y_bounds, interference_point) = get_angular_correlation_function(state, l);

    setup_style();
    int bins = 200;
    TCanvas* c = new TCanvas("c", "c", 600, 600);
    TH2D* hist = dalitz(argv[4], bins);

    hist->GetXaxis()->SetTitle("x");
    hist->GetYaxis()->SetTitle("y");
    hist->Draw("colz");

    TLine* bot = new TLine(-1, y_bounds[0], 1, y_bounds[0]);
    TLine* top = new TLine(-1, y_bounds[1], 1, y_bounds[1]);
    top->SetLineWidth(2);
    bot->SetLineWidth(2);
    top->SetLineColor(kRed);
    bot->SetLineColor(kRed);
    top->Draw();
    bot->Draw();

    string path = string(argv[1]) + "ang_compare_cut.pdf";
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
        if (y_bounds[0] < i && biny1 == 0) biny1 = b;
        if (y_bounds[1] < i && biny2 == 0) biny2 = b;
    }
    cout << format("bins: %1%, %2%") % biny1 % biny2 << endl;
    TH1D* h1 = hist->ProjectionX("px", biny1, biny2);
    h1->Draw();

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
    if (draw_excluded_area) {
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

        l1->DrawClone();
        l2->DrawClone();
        l3->DrawClone();
        l4->DrawClone();
    }

    double y = (y_bounds[1] + y_bounds[0])/2; // we use the middle of the bounded area as the y coordinate
    vector<double> x_bounds = {-sqrt(1-pow(y, 2)), sqrt(1-pow(y, 2))};
    TF1 *corr = new TF1("corr", ang_corr, x_bounds[0], x_bounds[1], 2);
    corr->SetParameter(0, maxval);
    corr->SetParameter(1, y);
    corr->Draw("same");

    path = string(argv[1]) + "ang_compare_projection.pdf";
    c2->SaveAs(path.c_str());
}