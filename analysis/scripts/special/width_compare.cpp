// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>

// other stuff
#include <boost/format.hpp>
#include <iostream>

// my stuff
#include "../plot_style.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[])
{
    //*** DALITZ PLOT SETUP ***//
    // set the axes    
    double x_axis[] = {200, -1.3, 1.3};
    double y_axis[] = {200, -1.3, 1.3};

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

    //*** ACTUAL CODE ***//
    setup_style();

    int fileno = 0; // used for a simple progress bar
    auto plot = [&] (const char* files[], int no) {
        TCanvas *c = new TCanvas("c", "c", 1800, 900);
        c->Divide(4, 2, 0, 0);
        for (int file = 0; file < 8; file++) {
            fileno++;
            cout << format("Analyzing file %1%/%2%.") % (fileno) % (argc-2) << "\r";
            c->cd(file+1);
            ROOT::RDF::RNode df = ROOT::RDataFrame("tree", files[file]);

            // Dalitz plot copied from ../basic/dalitz.cpp
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

            if (file == 0) {
                hist->GetYaxis()->SetTitle("sim3a");
                hist->GetYaxis()->SetTitleSize(0.1);
                hist->GetYaxis()->SetTitleOffset(0.3);
                hist->GetYaxis()->CenterTitle();
            } else if (file == 4) {
                hist->GetYaxis()->SetTitle("simX");
                hist->GetYaxis()->SetTitleSize(0.1);
                hist->GetYaxis()->SetTitleOffset(0.3);
                hist->GetYaxis()->CenterTitle();
            }

            if (file >= 4) {
                string title = to_string(5*(file+1-4) + 20*(no-1));
                hist->GetXaxis()->SetTitle(title.c_str());
                hist->GetXaxis()->SetTitleSize(0.1);
                hist->GetXaxis()->SetTitleOffset(0.3);
                hist->GetXaxis()->CenterTitle();    
            }
            
            hist->SetLabelSize(0, "X");
            hist->SetTickLength(0, "X");
            hist->SetLabelSize(0, "Y");
            hist->SetTickLength(0, "Y");
            hist->SetLabelSize(0, "Z");
            hist->SetTickLength(0, "Z");
            hist->DrawClone("col");
        }

        string path = string(argv[1]) + (format("width_compare_%1%.pdf") % no).str();
        c->SetLogz();
        c->SaveAs(path.c_str());
    };
    const char* first[8]; 
    const char* second[8]; 
    for (int i = 0; i < 4; i++) {
        first[i] = argv[2+i];
        second[i] = argv[2+4+i];
        first[4+i] = argv[2+8+i];
        second[4+i] = argv[2+12+i];
    }
    plot(first, 1);
    plot(second, 2);
}
