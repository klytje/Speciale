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

// my stuff
#include "../plot_style.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[]) {
    //*** PREPARE CORRELATION FUNCTION ***//
    // they are all calculated with my own angular_correlation.py script
    function<double(Double_t*, Double_t*)> ang_corr; 
    string state = string(argv[2]);
    string l = string(argv[3]);
    if (state == "0+" && l == "2") {
        ang_corr = [] (Double_t* x, Double_t* par) {
            double d = par[0]; 
            double maxval = par[1];
            double xp = M_PI/2*(1 + x[0]/d); // beta = 180 - theta', and theta' = pi/2(1 - x/d)
            double res = 2.25*pow(cos(xp), 4) - 1.5*pow(cos(xp), 2) + 0.25;
            return maxval*res;
        };
    } else if (state == "2-" && l == "1") {
        ang_corr = [] (Double_t* x, Double_t* par) {
            double d = par[0]; 
            double maxval = par[1];
            double xp = M_PI/2*(1 + x[0]/d); // beta = 180 - theta', and theta' = pi/2(1 - x/d)
            double res = 1 - pow(cos(xp), 2);
            return maxval*res;
        };
    } else if (state == "3-" && l == "1") {
        ang_corr = [] (Double_t* x, Double_t* par) {
            double d = par[0]; 
            double maxval = par[1];
            double xp = M_PI/2*(1 + x[0]/d); // beta = 180 - theta', and theta' = pi/2(1 - x/d)
            double res = 1./3*pow(cos(xp), 2) + 2./3;
            return maxval*res;
        };
    } else {
        cout << format("The given state %1% l = %2% is not implemented yet! \nYou can do it pretty easily yourself in ../scripts/special/ang_compare.cpp") % state % l << endl;
    }

    //*** DALITZ PLOT SETUP ***//
    double bins = 200;
    // set the axes    
    double x_axis[] = {bins, -1, 1};
    double y_axis[] = {bins, -1, 1};

    // define sorting methods
    auto max = [] (double e1, double e2, double e3) {return std::max({e1, e2, e3});};
    auto min = [] (double e1, double e2, double e3) {return std::min({e1, e2, e3});};
    auto mid = [] (double e1, double e2, double e3) {
        if (e1 > e2) {
            if (e2 > e3) {
                return e2;
            } else if (e1 > e3) {
                return e3;
            } else {
                return e1;
            }
        } else {
            if (e1 > e3) {
                return e1;
            } else if (e2 > e3) {
                return e3;
            } else {
                return e2;
            }
        }
    };

    //*** DALITZ PLOT (copied from ../basic/dalitz.cpp) ***//
    auto dalitz = [&] (const char* file) {
        ROOT::RDF::RNode df = ROOT::RDataFrame("tree", file);
        df = df.Define("E_tot","E_cm[0] + E_cm[1] + E_cm[2]")
                .Define("e_cm_1", "E_cm[0]/E_tot") // normalized such that e1 + e2 + e3 = 1
                .Define("e_cm_2", "E_cm[1]/E_tot")
                .Define("e_cm_3", "E_cm[2]/E_tot")
                .Define("e_1", max, {"e_cm_1", "e_cm_2", "e_cm_3"}) // we want e1 > e2 > e3
                .Define("e_2", mid, {"e_cm_1", "e_cm_2", "e_cm_3"})
                .Define("e_3", min, {"e_cm_1", "e_cm_2", "e_cm_3"})
                .Define("x","sqrt(3)*(e_2 - e_3)")
                .Define("y","3*e_1 - 1")
                .Filter("pow(x,2) + pow(y,2) < 1.0")
                .Filter("y < 0.93");

        TH2D* hist = new TH2D("h1", "Dalitz plot", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);
    
        // we can get the other slices simply by permutating i, j, k
        int perms[] = {1, 2, 3};
        do {
            int i = perms[0];
            int j = perms[1];
            int k = perms[2];
            TH2D htemp = df.Define("x_temp", (format("sqrt(3)*(e_%1% - e_%2%)") % j % k).str())
                        .Define("y_temp", (format("3*e_%1% - 1") % i).str())
                        .Histo2D({"h1", "temp", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x_temp", "y_temp").GetValue();
            hist->Add(&htemp);
        } while (std::next_permutation(perms, perms+3)); // repeat for each of the 3! = 6 permutations of {1, 2, 3}
        return hist;
    };

    //*** ACTUAL CODE ***//
    setup_style();
    // call convention is ./ang_compare <output path> <nuclear state> <l1> <sim3a files> <simX files>
    int sim3a_start = 4; // number of arguments preceding the root files
    int fileno = 0; // used for a simple progress bar
    string x_labels[8]; // contains the name of each column
    double bounds[2] = {0.35, 0.45}; // y axis bounds where the cut is imposed
    
    auto plot = [&] (const char* files[], int no) {
        // prepare the plots
        TCanvas *simX = new TCanvas("simX", "c", 1800, 900);
        TCanvas *sim3a = new TCanvas("sim3a", "c", 1800, 900);
        TCanvas *c = new TCanvas("c", "c", 1800, 900);
        simX->Divide(4, 2, 0, 0);
        sim3a->Divide(4, 2, 0, 0);
        c->Divide(4, 2, 0, 0);
        
        // loop over each file
        for (int file = 0; file < 8; file++) {
            fileno++;
            cout << format("Analyzing file %1%/%2%.") % (fileno) % (argc-sim3a_start) << "\r";
            TH2D* h = dalitz(files[file]);

            // this small section extracts gamma or Gamma from the file names
            string name = string(files[file]);
            int start = 10; // we want to start at X in output/0+_X
            int len = 0; // length of the number
            for (; start+len < name.size(); len++) {
                if (!isdigit(name.at(start+len))) {
                    if (name.at(start+len) == '.') { // check if it is a dot
                        if (isdigit(name.at(start+len+1))) { // check if it is a decimal separator (a number)
                            continue;
                        }
                    }
                    break;
                }
            }
            x_labels[file] = name.substr(start, len).c_str();

            // calculate the projection plot
            int biny1 = 0;
            int biny2 = 0;
            int b = 0;
            for (double i = -1; i < 1; i+=2./bins) {
                b++;
                if (bounds[0] < i && biny1 == 0) biny1 = b;
                if (bounds[1] < i && biny2 == 0) biny2 = b;
            }
                
            TH1D* hp = h->ProjectionX("px", biny1, biny2);
            TF1 *corr; // defined outside the lambda expression because we need it afterwards

            //*** FIGURE SETUP ***//
            // y labels for the leftmost plots
            auto set_ylabel = [] (const char title[], TH1* hist) {
                hist->GetYaxis()->SetTitle(title);
                hist->GetYaxis()->SetTitleSize(0.1);
                hist->GetYaxis()->SetTitleOffset(0.3);
                hist->GetYaxis()->CenterTitle();
            };

            if (file == 0) {
                set_ylabel("Dalitz", h);
                set_ylabel("Projection", hp);
            } else if (file == 4) {
                set_ylabel("Dalitz", h);
                set_ylabel("Projection", hp);
            }

            // hide ticks for all plots
            h->SetLabelSize(0, "X");
            h->SetTickLength(0, "X");
            h->SetLabelSize(0, "Y");
            h->SetTickLength(0, "Y");

            hp->SetLabelSize(0, "X");
            hp->SetTickLength(0, "X");
            hp->SetLabelSize(0, "Y");
            hp->SetTickLength(0, "Y");

            // x labels for all bottom row plots
            auto set_xlabel = [] (const char title[], TH1* hist) {
                hist->GetXaxis()->SetTitle(title);
                hist->GetXaxis()->SetTitleSize(0.1);
                hist->GetXaxis()->SetTitleOffset(0.3);
                hist->GetXaxis()->CenterTitle();
            };
            set_xlabel((x_labels[file]).c_str(), hp);

            // makes a Dalitz plot with two lines showing where the data is cut
            auto plot_top = [&] () {
                TLine* bot = new TLine(-1, bounds[0], 1, bounds[0]);
                TLine* top = new TLine(-1, bounds[1], 1, bounds[1]);
                top->SetLineWidth(2);
                bot->SetLineWidth(2);
                top->SetLineColor(kRed);
                bot->SetLineColor(kRed);
                h->DrawClone("colz");
                top->DrawClone("same");
                bot->DrawClone("same");
            };

            // makes a projection of the histogram down on the x-axis for the cut region. Also plots the correlation function
            auto plot_bot = [&] () {
                double maxval = hp->GetMaximum(); // maximum value in the histogram, used to scale the correlation function
                double d = sqrt(1 - pow((bounds[1] + bounds[0])/2, 2)); // distance from 0 to the border of the circle at y = mean(y1, y2)

                corr = new TF1("corr", ang_corr, -1, 1, 2);
                corr->SetParameter(0, d);
                corr->SetParameter(1, maxval);

                hp->DrawClone();
                corr->DrawCopy("same");
            };

            //*** MAKE THE sim3a PLOT ***//
            if (file < 4) {
                // top row
                sim3a->cd(file+1);
                plot_top();

                // bottom row
                sim3a->cd(file+5);
                plot_bot();
            } 

            //*** MAKE THE simX PLOT ***//
            else { // 4 < file
                // top row
                simX->cd(file-3);
                plot_top();

                // bottom row
                simX->cd(file+1);
                plot_bot();
            }

            //*** MIXED PLOT ***//
            c->cd(file+1);
            if (file < 4) {
                set_ylabel("sim3a", hp);
            }
            else {
                set_xlabel((x_labels[file-4]).c_str(), hp); // we want the xlabels to be small gamma
                set_ylabel("simX", hp);
            }

            hp->DrawClone();
            corr->DrawCopy("same");
        }
        // save the figures
        string path = string(argv[1]) + (format("ang_compare_%1%_sim3a.pdf") % no).str();
        sim3a->SetLogz();
        sim3a->SaveAs(path.c_str());

        path = string(argv[1]) + (format("ang_compare_%1%_simX.pdf") % no).str();
        simX->SetLogz();
        simX->SaveAs(path.c_str());

        path = string(argv[1]) + (format("ang_compare_%1%.pdf") % no).str();
        c->SaveAs(path.c_str());
    };
    
    if ((argc-sim3a_start) % 8 != 0) {
        cout << "Number of sim3a and simX files must be an integer multiple of 8! Received " << (argc-sim3a_start) << " files" << endl;
        cout << "Calling convention is: ./ang_compare <output path> <nuclear state> <l1> <sim3a files> <simX files>" << endl;
        exit(1);
    }

    // we need to convert the input to a more workable format. The idea is to create sets of 4 sim3a and their corresponding simX files, stored in
    // the first 4 entries and last 4 entries of files[8], respectively. 
    int simX_start = sim3a_start + (argc-sim3a_start)/2; // starting location of the simX files. sim3a files start at "sim3a_start"
    int plots = argc/8;
    int offset = 0;
    for (int i = 0; i < plots; i++) {
        const char* files[8];
        for (int j = 0; j < 4; j++) {
            files[j] = argv[sim3a_start+offset];
            files[4+j] = argv[simX_start+offset];
            offset += 1;
        }
        plot(files, i);
    }
}
