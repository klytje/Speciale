// optimization stuff
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

// other root stuff
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>
#include <ROOT/RDFHelpers.hxx>

// other stuff
#include <filesystem>
#include <boost/format.hpp>
#include <iostream>
#include <cmath>

using namespace std;
using namespace ROOT;
using boost::format;
using ROOT::RDF::RNode;

// merge the analyzed files from output/X.root with the matched files from matched/Xm.root
void merge(int num, char *path[]) {
    for (int i = 1; i < num; i++) { // skip first entry (that's the path to this script)
        // path to analyzed root file
        const char *apath = path[i];

        // path to matched root file
        filesystem::path p(apath);
        string tmp = "match/" + p.stem().string() + "m.root"; //stem is filename without extension; stem(x/y.z) = y
        const char *mpath = tmp.c_str();

        // open the root files
        TFile *fa = TFile::Open(apath);
        TFile *fm = TFile::Open(mpath);
        TTree *ta = (TTree *)fa->Get("a");
        TTree *tm = (TTree *)fm->Get("a101");

        // define variables from analyzed tree
        int mul, N;
        double pt, deltaE, E_cm[3], E_lab[3], exC12, theta_lab[3], theta_cm[3], phi_lab[3], phi_cm[3];
        ta->SetBranchAddress("N", &N);
        ta->SetBranchAddress("mul", &mul);
        ta->SetBranchAddress("pt", &pt);
        ta->SetBranchAddress("deltaE", &deltaE);
        ta->SetBranchAddress("E_cm", &E_cm);
        ta->SetBranchAddress("E_lab", &E_lab);
        ta->SetBranchAddress("exC12", &exC12);
        ta->SetBranchAddress("theta_cm", &theta_cm);
        ta->SetBranchAddress("phi_cm", &phi_cm);
        ta->SetBranchAddress("theta_lab", &theta_lab);
        ta->SetBranchAddress("phi_lab", &phi_lab);

        // define variables from matched tree
        double FT[3], BT[3], FI[3], BI[3], id[3];
        tm->SetBranchAddress("FI", &FI);
        tm->SetBranchAddress("FT", &FT);
        tm->SetBranchAddress("BI", &BI);
        tm->SetBranchAddress("BT", &BT);
        tm->SetBranchAddress("id", &id);

        // define destination tree and set its branches
        string dest = "merged/" + p.filename().string(); // file destination
        TFile f(dest.c_str(), "recreate");
        TTree t("tree", "merged tree for TDC calibration");
        //TFile *f = new TFile(dest.c_str(), "recreate");
        //TTree *t = new TTree("tree", "merged tree for TDC calibration");
        t.Branch("FI", &FI, "FI[3]/I");
        t.Branch("BI", &BI, "BI[3]/I");
        t.Branch("FT", &FT, "FT[3]/I");
        t.Branch("BT", &BT, "BT[3]/I");
        t.Branch("id", &id, "id[3]/I");
        t.Branch("E_cm", &E_cm, "E_cm[3]/D");
        t.Branch("E_lab", &E_lab, "E_lab[3]/D");
        t.Branch("theta_cm", &theta_cm, "theta_cm[3]/D");
        t.Branch("phi_cm", &phi_cm, "phi_cm[3]/D");
        t.Branch("theta_lab", &theta_lab, "theta_lab[3]/D");
        t.Branch("phi_lab", &phi_lab, "phi_lab[3]/D");
        t.Branch("pt", &pt, "pt/D");
        t.Branch("deltaE", &deltaE, "deltaE/D");
        t.Branch("mul", &mul, "mul/I");
        t.Branch("exC12", &exC12, "exC12/D");

        // loop over every event in the analyzed tree ta
        for (double i = 0; i < ta->GetEntries(); i++) {
            try {
                int status = ta->GetEntry(i); // try to get entry i
                if (status == 0) // if unsuccesful, print error
                    throw 1;
            }
            catch (int e) {
                cout << "Error in accessing entry; program will exit prematurely." << endl;
                exit(1);
            }
            tm->GetEntry(N); // N is the corresponding index in tm for each event in ta
            t.Fill(); // fill the entries into the merged tree t
        }
        t.Write(); // write to disk
    }
}

void print_title(string title) {
    cout << "\n\033[1;32m" + title + "\033[0m" << endl;
    return;
}

void tdcraw(RDataFrame *df) {
    print_title("*** TDC RAW FIGURE ***");
    // axis properties
    double x_axis[3] = {4097, 0, 4096}; // bins, xmin, xmax
    double y_axis[3] = {4097, 0, 4096}; // bins, ymin, ymax

    // W detectors (id == 2 or 3) //
    TCanvas *canvas = new TCanvas("tdc channels", "tdc channels", 600, 600);
    TH2D *h = new TH2D("bidi_h", "2d histo", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);    

    double n = 0; // total number of entries
    cout << "Binning data..." << endl;    
    for (int i = 0; i < 3; i++) {
        TH2D htemp = df->Filter((format("id[%1%] == 2 || id[%1%] == 3") % i).str())
                .Filter("mul == 3 && pt < 50e3 && abs(deltaE) < 200") // filter from his tdcmerger.C
                .Filter((format("FT[%1%] > 1000") % i).str())
                .Filter((format("BT[%1%] > 1000") % i).str())
                .Define("x", (format("(FT[%1%]%%651264)/100") % i).str()) // %% is formatted as a single %
                .Define("y", (format("(BT[%1%]%%651264)/100") % i).str())
                .Histo2D({"h2", "title", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "x", "y")
                .GetValue();
        n += h->GetEntries();
        h->Add(&htemp);
    }
    cout << "Finished binning data. Found a total of " << n << " events for W1 & W2." << endl;
    h->Draw("colz");
    canvas->SetLogz();
    canvas->Modified(); canvas->Update();
    canvas->WaitPrimitive();
    canvas->Close();
    return;
}

// external class for individual strip minimization. 
// this is necessary because we need to pass more information than just the pure fitting variables, which the ROOT minimizer doesn't support
class fit_func {
    public: 
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> *data;
        vector<double> FT[3]; // the FT branch in vector format
        vector<double> BT[3]; // the BT branch in vector format
        map<string, int> detector_map = {{"S3_1", 0}, {"S3_2", 1}, {"W1_1", 2}, {"W1_2", 3}};
        int det;
        int size; // length of the FT/BT vectors

        fit_func(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> *df, string detector, double offset[], int FID, int BID) {
            data = df;
            det = detector_map[detector];
            setup(offset[FID], FID, BID);
        }

        virtual double operator()(const double *x) {return 0;}

    private: 

        // extract the required data from the dataframe and put it into a vector
        // this is infinitely faster than operating directly on the dataframe (seconds versus milliseconds)
        void setup(double offset, int FID, int BID) {
            for (int i = 0; i < 3; i++) {
                // auto filter = [&i, &FID, &BID] (vector<int> FI, vector<int> BI) {return FI[i] == FID && BI[i] == BID;};
                string filter = (format("ID[%1%] == %2% && FI[%1%] == %3% && BI[%1%] == %4%") % i % det % FID % BID).str();
                FT[i] = data->Filter(filter).Define("x", (format("FT[%1%]*1e-3 - %2%") % i % offset).str()).Take<double>("x").GetValue();
                BT[i] = data->Filter(filter).Define("x", (format("BT[%1%]*1e-3") % i).str()).Take<double>("x").GetValue();
            }
            size = FT[0].size();
            cout << "Preparing fit to " << size*3 << " data points." << endl;
        }
};

// this fit function assumes deltaF = 0, such that there is only one free variable
class back_fit : fit_func {
    public:
        back_fit(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> *df, string detector, double offset[], int FID, int BID) : fit_func(df, detector, offset, FID, BID) {}

        double operator() (const double *x) override {
            double deltaB = x[0];
            double sum = 0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < size; j++) {
                    sum += abs(BT[i][j] + deltaB - FT[i][j]);
                }
            }
            return sum;
        }
};

void diff(RDataFrame *df, bool plotOffset = false, bool plotFit = true) {
    print_title("*** FRONT/BACK DIFFERENCE FIGURE ***");
    auto data = df->Filter("mul == 3 && pt < 50e3 && abs(deltaE) < 200"); // filter from his tdcmerger.C

    // *** HISTOGRAMS *** //
    // axis properties
    double x_axis[3] = {1000, -500, 500}; // bins, xmin, xmax
    double y_axis[3] = {16, 1, 17};    // bins, ymin, ymax

    // 1st W detector (id == 2) //
    TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);
    TH2D *h1 = new TH2D("bidi_h", "2d histo", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);    

    // iterate over each column, corresponding to each alpha particle
    cout << "Binning data for W1..." << endl;
    double n = 0; // total number of entries
    for (int i = 0; i < 3; i++) {
        TH2D htemp = data.Filter((format("id[%1%] == 2") % i).str()) // filter to events detected by W1
                         .Define("dt", (format("(BT[%1%] - FT[%1%])*1e-3") % i).str()) // define the difference between FT[i] and BT[i]
                         //.Filter("0 < dt && dt < 200")
                         .Define("x", (format("BI[%1%]") % i).str()) // define B[i] so we can use it in our histogram
                         .Histo2D({"h2", "title", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "dt", "x")
                         .GetValue(); // convert the ROOT smart pointer to the actual object        
        n += htemp.GetEntries();
        h1->Add(&htemp); // add the histogram to our initial one
    }
    cout << "Finished binning data. Found a total of " << n << " events for W1." << endl;

    // setting up the figure
    h1->GetXaxis()->SetTitle("Time [ns]");
    h1->GetYaxis()->SetTitle("Strip");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->Draw("colz");
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->Modified(); canvas->Update();
    canvas->WaitPrimitive();
    canvas->Close();


    // *** FITTING *** //    
    cout << "Starting fitting process for W1..." << endl;
    
    // we want to center all strips around 0, so we calculate the current mean
    double offset[16];
    for (int i = 0; i < 16; i++) {
        offset[i] = h1->ProjectionX("project x", i+1, i+2)->GetMean(); // remember that strip indices start from 1
    }

    // axis properties, _c for centered
    double x_axis_c[3] = {200, -100, 100}; // bins, xmin, xmax
    double y_axis_c[3] = {16, 1, 17};    // bins, ymin, ymax

    if (plotOffset) {
        // with the offsets found, we want to subtract them from each strip and redo the figure. Sadly there's no easy way of doing this
        TCanvas *canvas_offset = new TCanvas("diff", "diff", 600, 600);
        TH2D *h2 = new TH2D("bidi_h", "2d histo", int(x_axis_c[0]), x_axis_c[1], x_axis_c[2], int(y_axis_c[0]), y_axis_c[1], y_axis_c[2]); 
        n = 0; // reset the counter
        for (int i = 0; i < 3; i++) {
            // defining filters as lambda expressions makes a noticeable difference in computation time
            // we almost always have to use a .Define to extract the columns for the filters (BI[i] doesn't work), so unless it's in a loop it's probably not worth it
            auto dt_filter = [] (double dt) {return -100 < dt && dt < 100;}; 
            auto data_tmp = data.Filter((format("id[%1%] == 2") % i).str()) // filter to events detected by W1
                                .Define("BI_i", (format("BI[%1%]") % i).str());
            for (int j = 1; j < 17; j++) {
                auto BI_filter = [&j] (int BI) {return BI == j;};
                auto htemp = data_tmp.Filter(BI_filter, {"BI_i"}) // filter to events detected by the jth strip 
                                    .Define("dt", (format("(FT[%1%] - BT[%1%])*1e-3 - %2%") % i % offset[j-1]).str())
                                    .Filter(dt_filter, {"dt"})
                                    .Define("x", (format("BI[%1%]") % i).str())
                                    .Histo2D({"h2", "title", int(x_axis_c[0]), x_axis_c[1], x_axis_c[2], int(y_axis_c[0]), y_axis_c[1], y_axis_c[2]}, "dt", "x")
                                    .GetValue();
                n += htemp.GetEntries();
                h2->Add(&htemp); // add the histogram to our initial one
                cout << "Generated histogram " << i*16 + j << " of 48." << endl;
            }
        }
        cout << "Finished binning centered data. Found a total of " << n << " events for W1." << endl;
        
        // setting up the figure
        h2->GetXaxis()->SetTitle("Time [ns]");
        h2->GetYaxis()->SetTitle("Strip");
        h2->GetXaxis()->CenterTitle();
        h2->GetYaxis()->CenterTitle();
        h2->Draw("colz");
        canvas_offset->SetLogz();
        canvas_offset->SetRightMargin(0.15);
        canvas_offset->Modified(); canvas->Update();
        canvas_offset->WaitPrimitive();
        canvas_offset->Close();
    }

    // create the actual minimizer, see https://root.cern.ch/doc/v612/NumericalMinimization_8C.html for options
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(10000); // maybe increase this

    // INDIVIDUAL MINIMIZATION //
    // we want to minimize each strip by itself to obtain some good starting points for the global minimization

    // 1st fit //
    // we have 32 independent variables. by fixing one of them, the rest should be uniquely determined. I pick deltaF1 = 0
    double result[16];
    for (int i = 1; i < 17; i++) {
        double step = 0.01;
        double start = 0;
        back_fit f(&data, "W1_1", offset, 1, i); // FID = 1, BID = i
        ROOT:Math::Functor func(f, 1);
        minimizer->SetFunction(func);
        minimizer->SetVariable(0, "deltaB", start, step);
        minimizer->SetPrintLevel(1);
        minimizer->Minimize();
        auto res = minimizer->X();
        result[i-1] = res[0];
        cout << "Fit complete with result: " << result[0] << endl;
    }

    if (plotFit) {
        TCanvas *canvas_offset = new TCanvas("diff", "diff", 600, 600);
        TH2D *h2 = new TH2D("bidi_h", "2d histo", int(x_axis_c[0]), x_axis_c[1], x_axis_c[2], int(y_axis_c[0]), y_axis_c[1], y_axis_c[2]); 
        n = 0; // reset the counter
        for (int i = 0; i < 3; i++) {
            auto dt_filter = [] (double dt) {return -100 < dt && dt < 100;}; 
            auto data_tmp = data.Filter((format("id[%1%] == 2") % i).str()) // filter to events detected by W1
                                .Define("BI_i", (format("BI[%1%]") % i).str());
            for (int j = 1; j < 17; j++) {
                auto BI_filter = [&j] (int BI) {return BI == j;};
                auto htemp = data_tmp.Filter(BI_filter, {"BI_i"}) // filter to events detected by the jth strip 
                                    .Define("dt", (format("(FT[%1%] - BT[%1%])*1e-3 - %2% - %3%") % i % offset[j-1]).str())
                                    .Filter(dt_filter, {"dt"})
                                    .Define("x", (format("BI[%1%]") % i).str())
                                    .Histo2D({"h2", "title", int(x_axis_c[0]), x_axis_c[1], x_axis_c[2], int(y_axis_c[0]), y_axis_c[1], y_axis_c[2]}, "dt", "x")
                                    .GetValue();
                n += htemp.GetEntries();
                h2->Add(&htemp); // add the histogram to our initial one
                cout << "Generated histogram " << i*16 + j << " of 48." << endl;
            }
        }
        cout << "Finished binning centered data. Found a total of " << n << " events for W1." << endl;
        
        // setting up the figure
        h2->GetXaxis()->SetTitle("Time [ns]");
        h2->GetYaxis()->SetTitle("Strip");
        h2->GetXaxis()->CenterTitle();
        h2->GetYaxis()->CenterTitle();
        h2->Draw("colz");
        canvas_offset->SetLogz();
        canvas_offset->SetRightMargin(0.15);
        canvas_offset->Modified(); canvas->Update();
        canvas_offset->WaitPrimitive();
        canvas_offset->Close();
    }

//     double step[2] = {0.01, 0.01}; 
//     double start[2] = {0, 0}; 
//     for (int i = 1; i < 17; i++) {
//         strip_fitter f(data, 1, i);
//         ROOT:Math::Functor func(f, 2);
//         minimizer->SetFunction(func);
//         minimizer->SetVariable(0, "deltaF", start, step);
//         minimizer->SetVariable(1, "deltaB", start, step);
//         minimizer->Minimize();
//         auto result = minimizer->X();
//         cout << "Fit complete with result: " << result[0] << endl;
//  }
   
    // // GLOBAL MINIMIZATION //
    // // total minimization function
    // auto dt_sum = [] (RDataFrame* df, double deltaF[], double deltaB[]) {
    //     double sum;
    //     for (int i = 0; i < 3; i++) {
    //         for (int j = 1; j < 17; j++) {
    //             sum = df->Define("f", (format("dt_%1%_%2% + %3% - %4%") % i % j % deltaF[j] % deltaB[j]).str()).Sum("f").GetValue();
    //         }            
    //     }
    //     return sum;
    // };

//    ROOT::Math::Functor total_func(&dt_sum, 32);

    //minimizer->SetFunction(func);
    return;
}

void debug(RDataFrame *df) {
    print_title("*** FRONT/BACK DIFFERENCE FIGURE ***");
    auto data = df->Filter("mul == 3 && pt < 50e3 && abs(deltaE) < 200"); // filter from his tdcmerger.C

    auto reduced = data.Range(100);
    
    vector<double> *FT = new vector<double>(), *BT = new vector<double>();
    vector<int> *FI = new vector<int>(), *BI = new vector<int>(), *ID = new vector<int>(); 
    // append y to x
    auto append = [](auto x, auto y) {
        x->insert(x->end(), y.begin(), y.end());
    };

    for (int i = 0; i < 3; i++) {
        append(FT, reduced.Define("x", (format("FT[%1%]*1e-3") % i).str()).Take<double>("x").GetValue());
        append(BT, reduced.Define("x", (format("BT[%1%]*1e-3") % i).str()).Take<double>("x").GetValue());
        append(FI, reduced.Define("x", (format("FI[%1%]") % i).str()).Take<int>("x").GetValue());
        append(BI, reduced.Define("x", (format("BI[%1%]") % i).str()).Take<int>("x").GetValue());
        append(ID, reduced.Define("x", (format("id[%1%]") % i).str()).Take<int>("x").GetValue());
    }

    for (int i = 0; i < 1000; i++) {
        cout << format("%1% : %2% : %3% : %4% : %5%") % FT->at(i) % FI->at(i) % BT->at(i) % BI->at(i) % ID->at(i) << endl;
    }
    
    
    // *** HISTOGRAMS *** //
    // axis properties
    double x_axis[3] = {200, -100, 100}; // bins, xmin, xmax
    double y_axis[3] = {16, 1, 17};    // bins, ymin, ymax

    // 1st W detector (id == 2) //
    TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);
    // TH2D *h1 = new TH2D("bidi_h", "2d histo", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]);    
    TH1D *h1 = new TH1D("bidi_h", "2d histo", int(x_axis[0]), x_axis[1], x_axis[2]);

    // iterate over each column, corresponding to each alpha particle
    cout << "Binning data for W1..." << endl;
    double n = 0; // total number of entries
    for (int i = 0; i < 3; i++) {
        TH1D htemp = data.Filter((format("id[%1%] == 2") % i).str()) // filter to events detected by W1
                         .Define("dt", (format("(BT[%1%] - FT[%1%])*1e-3") % i).str()) // define the difference between FT[i] and BT[i]
                         //.Filter("0 < dt && dt < 200")
                         .Define("x", (format("BT[%1%]") % i).str()) // define B[i] so we can use it in our histogram
                         .Histo1D({"h", "title", int(x_axis[0]), x_axis[1], x_axis[2]}, "x")
                        //  .Histo2D({"h2", "title", int(x_axis[0]), x_axis[1], x_axis[2], int(y_axis[0]), y_axis[1], y_axis[2]}, "dt", "x")
                         .GetValue(); // convert the ROOT smart pointer to the actual object        
        n += htemp.GetEntries();
        h1->Add(&htemp); // add the histogram to our initial one
    }
    cout << "Finished binning data. Found a total of " << n << " events for W1." << endl;
    // set up the figure
    h1->GetXaxis()->SetTitle("Time [ns]");
    h1->GetYaxis()->SetTitle("Strip");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->Draw("colz");
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    canvas->Modified(); canvas->Update();
    canvas->WaitPrimitive();
    canvas->Close();
}

void debug2(RDataFrame *df) {
    // TCanvas *canvas = new TCanvas("debug", "debug", 600, 600);
    double x_axis[3] = {61, -30, 30}; // bins, xmin, xmax
    double y_axis[3] = {32, -33, 33};    // bins, ymin, ymax

    auto data = df->Filter("mul == 3 && pt < 50e3 && abs(deltaE) < 200");  // filter from his tdcmerger.C

    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(10000); // maybe increase this
    back_fit f(&data, "W1_1", {0}, 1, 1);
    ROOT:Math::Functor func(f, 2);

    double start[2] = {-1, 1.2};
    double step[2] = {0.01, 0.01};
    minimizer->SetPrintLevel(1);
    minimizer->SetFunction(func);
    minimizer->SetVariable(0, "x", start[0], step[0]);
    minimizer->SetVariable(1, "y", start[1], step[1]);
    cout << "Starting minimization" << endl;
    minimizer->Minimize();
    auto result = minimizer->X();
    cout << "Fit complete with result: " << result[0] << endl;

    // canvas->Modified(); canvas->Update();
    // canvas->WaitPrimitive();
    // canvas->Close();
}

int main(int argc, char *argv[]) {
    // ROOT::EnableImplicitMT(8); // enable multithreading with 8 threads. NOTE: scrambles event ordering, irrelevant for us

    // check if all merged files exists, and create them if not
    string paths[argc-1]; // merged file paths
    bool mergeflag = false;
    cout << "\nChecking if all prerequisite files exists.." << endl;
    for (int i = 1; i < argc; i++) {
        filesystem::path p(argv[i]);
        string file = "merged/" + p.filename().string();
        if (!filesystem::exists(file)) {
            mergeflag = true;
            cout << "\033[1;31m" << boost::format("Missing %1%") % file << "\033[0m" << endl; // red colour
        }
        paths[i-1] = file; // loop starts at 1 while paths start at 0
    }
    if (mergeflag) {
        cout << "Files missing, recreating all merged files.." << endl;
        merge(argc, argv);
    } else {
        cout << "All merged files found." << endl;
    }

    TChain chain("tree");
    for (int i = 0; i < argc-1; i++) {
        chain.Add(paths[i].c_str());
    }
    RDataFrame df(chain);
    cout << format("Dataframe initialized with %1% entries.") % df.Count().GetValue() << endl;

    // define the plot colour scheme
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gROOT->ForceStyle();

    // start a ROOT application window such that the plots can actually be shown
    TApplication *app = new TApplication("ROOT window", 0, 0);
    // debug(&df); // various stuff for debugging
    //tdcraw(&df); // generate the "raw tdc" figure
    diff(&df); // generate the front/back detector difference figure
    app->Run(); // show all canvas
    return 0;
}