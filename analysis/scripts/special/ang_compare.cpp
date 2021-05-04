// ROOT stuff
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>

// other stuff
#include <boost/format.hpp>

// my stuff
#include "../plot_style.cpp"
#include "../utility.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[]) {
    vector<double> y_bounds; // y axis bounds
    double interference_point; // point of maximum interference
    function<double(Double_t*, Double_t*)> ang_corr;
    tie(ang_corr, y_bounds, interference_point) = get_angular_correlation_function(argv[2], argv[3]);
    int bins = 200;
    
    setup_style();
    // call convention is ./ang_compare <output path> <nuclear state> <l1> <sim3a files> <simX files>
    int sim3a_start = 4; // number of arguments preceding the root files
    int fileno = 0; // used for a simple progress bar
    string x_labels[8]; // contains the name of each column
    
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
            TH2D* h = dalitz(files[file], bins);

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
                if (y_bounds[0] < i && biny1 == 0) biny1 = b;
                if (y_bounds[1] < i && biny2 == 0) biny2 = b;
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
                TLine* bot = new TLine(-1, y_bounds[0], 1, y_bounds[0]);
                TLine* top = new TLine(-1, y_bounds[1], 1, y_bounds[1]);
                top->SetLineWidth(2);
                bot->SetLineWidth(2);
                top->SetLineColor(kRed);
                bot->SetLineColor(kRed);
                
                top->SetLineStyle(11); //dashed linestyle
                bot->SetLineStyle(11); 

                h->DrawClone("colz");
                top->DrawClone("same");
                bot->DrawClone("same");
            };

            // makes a projection of the histogram down on the x-axis for the cut region. Also plots the correlation function
            auto plot_bot = [&] () {
                // we need to determine the max value of the histogram, but it cannot be at the interference point (about 0.6)
                int maxval = 0;
                double width = 0.1;
                double bad_area[4] = {-(interference_point+width), -(interference_point-width), interference_point-width, interference_point+width};
                for (int i = 0; i < bins; i++) {
                    int count = hp->GetBinContent(i);
                    double loc = hp->GetBinCenter(i);
                    if (maxval < count) {
                        if (!(bad_area[0] < loc && loc < bad_area[1]) && !(bad_area[2] < loc && loc < bad_area[3])) {
                            maxval = count;
                        }
                    }
                }

                double y = (y_bounds[1] + y_bounds[0])/2; // we use the middle of the bounded area as the y coordinate
                vector<double> x_bounds = {-sqrt(1-pow(y, 2)), sqrt(1-pow(y, 2))};
                corr = new TF1("corr", ang_corr, x_bounds[0], x_bounds[1], 2);
                corr->SetParameter(0, maxval);
                corr->SetParameter(1, y);

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
