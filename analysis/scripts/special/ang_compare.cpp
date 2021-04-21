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
//W(θ) = 0.25*(3*sin(θ)**2 - 2)**2

double ang_corr(Double_t* x, Double_t* par) {
        double d = par[0]; 
        double maxval = par[1];
        double xp = M_PI/2*(1 + x[0]/d); // beta = 180 - theta', and theta' = pi/2(1 - x/d)
        double res = 2.25*pow(cos(xp), 4) - 1.5*pow(cos(xp), 2) + 0.25;
        return 150*res;
}

int main(int argc, char const *argv[]) {
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

            // setup labels and general plot style
            if (file == 0) {
                h->GetYaxis()->SetTitle("sim3a");
                h->GetYaxis()->SetTitleSize(0.1);
                h->GetYaxis()->SetTitleOffset(0.3);
                h->GetYaxis()->CenterTitle();
            } else if (file == 4) {
                h->GetYaxis()->SetTitle("simX");
                h->GetYaxis()->SetTitleSize(0.1);
                h->GetYaxis()->SetTitleOffset(0.3);
                h->GetYaxis()->CenterTitle();
            }

            if (file >= 4) {
                string title = to_string(5*(file+1-4) + 20*(no-1));
                h->GetXaxis()->SetTitle(title.c_str());
                h->GetXaxis()->SetTitleSize(0.1);
                h->GetXaxis()->SetTitleOffset(0.3);
                h->GetXaxis()->CenterTitle();    
            }
            
            h->SetLabelSize(0, "X");
            h->SetTickLength(0, "X");
            h->SetLabelSize(0, "Y");
            h->SetTickLength(0, "Y");
            h->SetLabelSize(0, "Z");
            h->SetTickLength(0, "Z");

            if (file < 4) {
                sim3a->cd(file+1);
            } else {
                simX->cd(file+1);
            }
            
            // top row of Dalitz plot + Y cut
            double bounds[2] = {0.35, 0.45}; // y axis bounds
            TLine* bot = new TLine(-1, bounds[0], 1, bounds[0]);
            TLine* top = new TLine(-1, bounds[1], 1, bounds[1]);
            top->SetLineWidth(2);
            bot->SetLineWidth(2);
            top->SetLineColor(kRed);
            bot->SetLineColor(kRed);
            h->DrawClone("colz");
            top->DrawClone("same");
            bot->DrawClone("same");

            // bottom row of projection plot + angular correlation function
            if (file < 4) {
                sim3a->cd(file+5);
            } else {
                simX->cd(file+5);
            }

            int biny1 = 0;
            int biny2 = 0;
            int b = 0;
            for (double i = -1; i < 1; i+=2./bins) {
                b++;
                if (bounds[0] < i && biny1 == 0) biny1 = b;
                if (bounds[1] < i && biny2 == 0) biny2 = b;
            }
            
            // distance from 0 to the border of the circle at y = mean(y1, y2)
            double d = sqrt(1 - pow((bounds[1] + bounds[0])/2, 2)); 

            TF1 *corr = new TF1("corr", ang_corr, -1, 1, 1);
            TH1D* hp = h->ProjectionX("px", biny1, biny2);
            corr->SetParameter(0, d);

            hp->DrawClone();
            corr->DrawCopy("same");

            // mixed plot
            c->cd(file+1);
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
