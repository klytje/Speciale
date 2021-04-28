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

int main(int argc, char const *argv[]) {
    //*** DALITZ PLOT SETUP ***//
    // set the axes    
    double x_axis[] = {200, -1, 1};
    double y_axis[] = {200, -1, 1};

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
    int sim3a_start = 2; // number of arguments preceding the root files
    int fileno = 0; // used for a simple progress bar

    string x_labels[4]; // contains the name of each column
    auto plot = [&] (const char* files[], int no) {
        TCanvas *c = new TCanvas("c_raw", "c", 1800, 900);
        c->Divide(4, 2, 0, 0);

        for (int file = 0; file < 8; file++) {
            fileno++;
            cout << format("Analyzing file %1%/%2%.") % (fileno) % (argc-sim3a_start) << "\r";
            TH2D* hsim = dalitz(files[file]);

            if (file == 0) {
                hsim->GetYaxis()->SetTitle("sim3a");
                hsim->GetYaxis()->SetTitleSize(0.1);
                hsim->GetYaxis()->SetTitleOffset(0.3);
                hsim->GetYaxis()->CenterTitle();
            } else if (file == 4) {
                hsim->GetYaxis()->SetTitle("simX");
                hsim->GetYaxis()->SetTitleSize(0.1);
                hsim->GetYaxis()->SetTitleOffset(0.3);
                hsim->GetYaxis()->CenterTitle();
            }

            // this small section extracts gamma from the sim3a file names, and uses them as xlabels
            if (file < 4) {
                string name = string(files[file]);
                int start = 10; // we want to start at X in output/0+_X
                int len = 0; // length of the number
                for (; start+len < name.size(); len++) {
                    if (!isdigit(name.at(start+len))) {
                        break;
                    }
                }
                x_labels[file] = name.substr(start, len).c_str();
            } else {
                hsim->GetXaxis()->SetTitle(x_labels[file-4].c_str());
                hsim->GetXaxis()->SetTitleSize(0.1);
                hsim->GetXaxis()->SetTitleOffset(0.3);
                hsim->GetXaxis()->CenterTitle();    
            }
            
            hsim->SetLabelSize(0, "X");
            hsim->SetTickLength(0, "X");
            hsim->SetLabelSize(0, "Y");
            hsim->SetTickLength(0, "Y");
            hsim->SetLabelSize(0, "Z");
            hsim->SetTickLength(0, "Z");

            c->cd(file+1);
            hsim->DrawClone("col");
        }

        string path = string(argv[1]) + (format("depend_width_%1%.pdf") % no).str();
        c->SetLogz();
        c->SaveAs(path.c_str());
    };

    if ((argc-sim3a_start) % 8 != 0) {
        cout << "Number of sim3a and simX files must be an integer multiple of 8! Received " << (argc-sim3a_start) << " files" << endl;
        cout << "Calling convention is: ./ang_compare <output path> <sim3a files> <simX files>" << endl;
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
