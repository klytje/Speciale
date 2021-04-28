// ROOT stuff
#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <TChain.h>
#include <TText.h>
#include <TPad.h>

// my own stuff
#include "../calibrate/utility.cpp"

using namespace std;

int main(int argc, char *argv[]) {
    TChain chain("a101");
    for (int i = 1; i < argc-1; i++) {
        chain.Add(argv[i]);
    }
    ROOT::RDataFrame df(chain);

    // change some settings
    gROOT->SetBatch(kTRUE); // no graphics display
    plot::save = false;
    plot::path = string(argv[argc-1]);
    setup();

    TCanvas *c = new TCanvas("mul_compare", "", 1200, 900);
    c->Divide(4, 3, 0, 0); // split into 2x2 plot

    for (int det = 0; det < 4; det++) {
        print_title("Now plotting detector " + detector_map.at(det) + ".");
        std::cout << "Calculating mul == 1..." << endl;
        auto top = df.Filter((format("mul==1 && id[0]==%1%") % det).str()).Define("x", "BT[0]*1e-3").Histo1D({"legend", (format("%1%, mul == 1") % detector_map.at(det)).str().c_str(), int(plot::x_axes::BT[0]), plot::x_axes::BT[1], plot::x_axes::BT[2]}, "x");

        std::cout << "Calculating mul == 2..." << endl;
        auto mid = df.Filter((format("mul==2 && id[0]==%1%") % det).str()).Define("x", "BT[0]*1e-3").Histo1D({"legend", (format("%1%, mul == 2") % detector_map.at(det)).str().c_str(), int(plot::x_axes::BT[0]), plot::x_axes::BT[1], plot::x_axes::BT[2]}, "x");
        auto h1 = df.Filter((format("mul==2 && id[1]==%1%") % det).str()).Define("x", "BT[1]*1e-3").Histo1D({"legend", (format("%1%, mul == 2") % detector_map.at(det)).str().c_str(), int(plot::x_axes::BT[0]), plot::x_axes::BT[1], plot::x_axes::BT[2]}, "x").GetValue();
        mid->Add(&h1);
        h1.Reset("ICESM");

        std::cout << "Calculating mul == 3..." << endl;
        auto bot = df.Filter((format("mul==3 && id[0]==%1%") % det).str()).Define("x", "BT[0]*1e-3").Histo1D({"legend", (format("%1%, mul == 3") % detector_map.at(det)).str().c_str(), int(plot::x_axes::BT[0]), plot::x_axes::BT[1], plot::x_axes::BT[2]}, "x");
        h1 = df.Filter((format("mul==3 && id[1]==%1%") % det).str()).Define("x", "BT[1]*1e-3").Histo1D({"legend", (format("%1%, mul == 3") % detector_map.at(det)).str().c_str(), int(plot::x_axes::BT[0]), plot::x_axes::BT[1], plot::x_axes::BT[2]}, "x").GetValue();
        auto h2 = df.Filter((format("mul==3 && id[2]==%1%") % det).str()).Define("x", "BT[2]*1e-3").Histo1D({"legend", (format("%1%, mul == 3") % detector_map.at(det)).str().c_str(), int(plot::x_axes::BT[0]), plot::x_axes::BT[1], plot::x_axes::BT[2]}, "x").GetValue();
        bot->Add(&h1);
        bot->Add(&h2);
        h1.Reset("ICESM");
        h2.Reset("ICESM");

        // y axis labels
        if (det == 0) {
            c->cd(1);
            top->GetYaxis()->SetTitle("1");
            top->GetYaxis()->SetTitleSize(0.1);
            top->GetYaxis()->SetTitleOffset(0.3);
            top->GetYaxis()->CenterTitle();

            c->cd(5);
            mid->GetYaxis()->SetTitle("2");
            mid->GetYaxis()->SetTitleSize(0.1);
            mid->GetYaxis()->SetTitleOffset(0.3);
            mid->GetYaxis()->CenterTitle();

            c->cd(9);
            bot->GetYaxis()->SetTitle("3");
            bot->GetYaxis()->SetTitleSize(0.1);
            bot->GetYaxis()->SetTitleOffset(0.3);
            bot->GetYaxis()->CenterTitle();
        }

        c->cd(det+1);
        top->SetLabelSize(0, "X");
        top->SetTickLength(0, "X");
        top->SetLabelSize(0, "Y");
        top->SetTickLength(0, "Y");
        top->DrawClone();

        c->cd(det+5);
        mid->SetLabelSize(0, "X");
        mid->SetTickLength(0, "X");
        mid->SetLabelSize(0, "Y");
        mid->SetTickLength(0, "Y");
        mid->DrawClone();

        // bottom plot & x axis labels
        c->cd(det+9);
        bot->SetLabelSize(0, "X");
        bot->SetTickLength(0, "X");
        bot->SetLabelSize(0, "Y");
        bot->SetTickLength(0, "Y");

        bot->GetXaxis()->SetTitle(detector_map.at(det).c_str());
        bot->GetXaxis()->SetTitleSize(0.1);
        bot->GetXaxis()->SetTitleOffset(0.3);
        bot->GetXaxis()->CenterTitle();
        bot->DrawClone();
    }

    c->cd(5);
    TText* t = new TText(0, 0.5, "Coincidence");
    t->SetNDC();
    t->SetTextSize(0.1);
    t->SetTextAngle(90);
    t->SetTextAlign(22);
    t->DrawClone();

    string path = plot::path + "depend_mul.pdf";
    c->SaveAs(path.c_str());
    c->Close();
}
