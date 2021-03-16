#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <TChain.h>

#include <iostream>

#include "plots.cpp"
#include "utility.cpp"

using namespace std;

int main(int argc, char *argv[]) {
    TChain chain("a101");
    for (int i = 1; i < argc; i++) {
        chain.Add(argv[i]);
    }
    ROOT::RDataFrame df(chain);

    TApplication *app = new TApplication("ROOT window", 0, 0);
    TCanvas *canvas = new TCanvas("diff", "diff", 600, 600);

    for (int det = 0; det < 4; det++) {
        cout << "Calculating, please wait.." << endl;
        auto h = df.Filter((format("mul==1 && id[0]==%1%") % det).str()).Define("x", "BT[0]*1e-3").Histo1D({"legend", (format("%1%, mul == 1") % detector_map.at(det)).str().c_str(), 600, 65000, 65600}, "x");
        h->Draw();
        canvas->Modified(); canvas->Update();
        canvas->WaitPrimitive();
        h->Reset("ICESM");

        cout << "Calculating, please wait.." << endl;
        h = df.Filter((format("mul==2 && id[0]==%1%") % det).str()).Define("x", "BT[0]*1e-3").Histo1D({"legend", (format("%1%, mul == 2") % detector_map.at(det)).str().c_str(), 600, 65000, 65600}, "x");
        auto h1 = df.Filter((format("mul==2 && id[1]==%1%") % det).str()).Define("x", "BT[1]*1e-3").Histo1D({"legend", (format("%1%, mul == 2") % detector_map.at(det)).str().c_str(), 600, 65000, 65600}, "x").GetValue();
        h->Add(&h1);
        h->Draw();
        canvas->Modified(); canvas->Update();
        canvas->WaitPrimitive();
        h->Reset("ICESM");
        h1.Reset("ICESM");

        cout << "Calculating, please wait.." << endl;
        h = df.Filter((format("mul==3 && id[0]==%1%") % det).str()).Define("x", "BT[0]*1e-3").Histo1D({"legend", (format("%1%, mul == 3") % detector_map.at(det)).str().c_str(), 600, 65000, 65600}, "x");
        h1 = df.Filter((format("mul==3 && id[1]==%1%") % det).str()).Define("x", "BT[1]*1e-3").Histo1D({"legend", (format("%1%, mul == 3") % detector_map.at(det)).str().c_str(), 600, 65000, 65600}, "x").GetValue();
        auto h2 = df.Filter((format("mul==3 && id[2]==%1%") % det).str()).Define("x", "BT[2]*1e-3").Histo1D({"legend", (format("%1%, mul == 3") % detector_map.at(det)).str().c_str(), 600, 65000, 65600}, "x").GetValue();
        h->Add(&h1);
        h->Add(&h2);
        h->Draw();
        canvas->Modified(); canvas->Update();
        canvas->WaitPrimitive();
        h->Reset("ICESM");
        h1.Reset("ICESM");
        h2.Reset("ICESM");
    }
    canvas->Close();
}
