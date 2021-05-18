// ROOT stuff
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>

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

    vector<double> y_bounds; // y axis bounds
    double interference_point; // point of maximum interference
    tie(std::ignore, y_bounds, interference_point) = get_angular_correlation_function(state, l);

    // since the real data is slightly different from the simulations, I specify other values of y_bounds for them
    bool real_data = false;
    if (string(argv[argc-1]).find("events") != string::npos) { // if we are dealing with real data
        cout << "File name contains \"events\", assuming real data..." << endl;
        real_data = true;
        if (state == "0+") {
            y_bounds = {0.38, 0.54};
        } else if (state == "2-") {
            y_bounds = {0.30, 0.52};
        } else if (state == "3-") {
            y_bounds = {0.3, 0.6};
        }
    }

    setup_style();
    int bins = 200;
    TCanvas* c = new TCanvas("c", "c", 600, 600);
    TH2D* hist = dalitz(argv[4], bins, true);

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

    string path = string(argv[1]) + "ang_cor_fit_cut.pdf";
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
    TH1D* h1 = hist->ProjectionX("px", biny1, biny2);
    h1->GetXaxis()->SetNdivisions(3);
    h1->GetXaxis()->SetTitle("x");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->SetNdivisions(2);
    h1->GetYaxis()->SetTitle("Angular correlation");
    h1->GetYaxis()->CenterTitle();
    h1->SetLineColor(kBlack);
    h1->SetLineWidth(3);
    h1->Draw("HIST L");

    // we need to determine the max value of the histogram, but it cannot be at the interference point (about 0.6)
    double maxval = 0;
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
    double y = (y_bounds[1] + y_bounds[0])/2; // we use the middle of the bounded area as the y coordinate
    // the x bounds needs to be sliiightly smaller than one would predict due to floating point errors (theta is undefined for x-values outside this bound)
    vector<double> x_bounds = {-sqrt(1-pow(y, 2))+0.0001, sqrt(1-pow(y, 2))-0.0001};

    // plot multiple correlation functions
    vector<string> states = {"0+", "1-", "1+", "1-", "2-", "2+", "2-", "3-", "3+", "3-"}; // parity doesn't actually matter
    vector<string> ls =     {"2" , "1" , "2" , "3" , "1" , "2" , "3" , "1" , "2" , "3"};
    vector<int> color = {kBlack, kPink-8, kOrange+1, kYellow-7, kSpring-1, kTeal+3, kCyan+1, kAzure+1, kBlue-9, kViolet+1};
    function<double(Double_t*, Double_t*)> ang_corr;

    maxval = 1./h1->GetMaximum()*maxval;
    h1->Scale(1./h1->GetMaximum());
    h1->SetMaximum(1.2); // make space for the legend
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
    
    path = string(argv[1]) + "ang_cor_fit.pdf";
    c2->SaveAs(path.c_str());

    //*** PLOT ONLY RELEVANT FUNCTIONS ***//
    TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);
    h1->Draw("HIST L");
    auto legend3 = new TLegend(0.1, 0.75, 0.2, 0.9);
    int i = 0;
    for (auto const& [key, val] : correlation_functions) {
        string Jp = state.substr(0, 1) + "+";
        string Jm = state.substr(0, 1) + "-";
        string s = "";

        bool plot = false;
        if (key.find(Jp) != string::npos) {
            plot = true;
            s = Jp;
        } else if (key.find(Jm) != string::npos) {
            plot = true;
            s = Jm;
        }
        if (plot) {
            i++;
            string ll = key.substr(3);
            tie(ang_corr, std::ignore, std::ignore) = get_angular_correlation_function(s, ll);

            string label = (format("%1% %2%") % s % ll).str();
            TF1 *corr = new TF1(label.c_str(), ang_corr, x_bounds[0], x_bounds[1], 2);
            corr->SetParameter(0, maxval);
            corr->SetParameter(1, y);
            corr->SetLineColor(color[i]);
            corr->DrawClone("same");

            legend3->AddEntry(label.c_str(), label.c_str(), "l");
        }
    }
    legend3->Draw();
    path = string(argv[1]) + "ang_cor_relevant_funcs.pdf";
    cout << path << endl;
    c3->SaveAs(path.c_str());


    if (!real_data) {
        exit(0);
    }

    //*** FIT ***//
    // attempt to fit the data (very specific to my thesis)
    if (state == "0+" || state == "2-") {
        double guess = 0; // guess value for c
        if (state == "0+") {
            guess = 1./3;
        } else if (state == "2-") {
            guess = 0;
        }

        TCanvas* c4 = new TCanvas("c4", "c4", 600, 600);
        h1->SetMaximum(1);
        auto f13 = [] (double* x, double* par) {
            double y = par[0];
            double c = par[1];
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta;
            return c*correlation_functions.at("2- 1")(beta) + (1-c)*correlation_functions.at("2- 3")(beta);
        };
        TF1* tf13 = new TF1("mix", f13, x_bounds[0], x_bounds[1], 2);
        tf13->FixParameter(0, y);
        tf13->SetParameter(1, guess);
        tf13->SetParLimits(1, 0, 1);
        h1->Fit(tf13, "QRL");
        double c13 = tf13->GetParameter(1);
        cout << format("Fitted values: c1 = %1%, c3 = %2%") % c13 % (1-c13) << endl;

        auto f0 = [] (double* x, double* par) {
            double y = par[0];
            double c1 = par[1];
            double c2 = par[2];
            double theta = acos(x[0]/sqrt(1-pow(y, 2)));
            double beta = M_PI - theta;
            return (1-c2)*c1*correlation_functions.at("2- 1")(beta) + (1-c2)*(1-c1)*correlation_functions.at("2- 3")(beta) 
                    + c2*correlation_functions.at("0+ 2")(beta);
        };
        TF1* tf0 = new TF1("mix", f0, x_bounds[0], x_bounds[1], 3);
        tf0->FixParameter(0, y);
        tf0->SetParameter(1, c13); // c1
        tf0->SetParameter(2, 0.1); // c2
        tf0->SetParLimits(1, 0, 1);
        tf0->SetParLimits(2, 0, 1);
        h1->Fit(tf0, "QRLM");
        double pars[3];
        tf0->GetParameters(pars);
        double c1 = pars[1], c2 = pars[2];
        double c0 = tf0->GetParameter(2);
        cout << format("Fitted values: 0+ 2 = %1%, 2- 1 = %2%, 2- 3 = %3%") % c2 % ((1-c2)*c1) % ((1-c2)*(1-c1)) << endl;

        h1->Draw("HIST L");
        tf0->Draw("SAME");
        path = string(argv[1]) + "ang_cor_fitted.pdf";
        cout << path << endl;
        c4->SaveAs(path.c_str());
    } 
}