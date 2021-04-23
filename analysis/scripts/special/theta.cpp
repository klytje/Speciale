// ROOT stuff
#include <TAxis.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH2D.h>

// other stuff
#include <boost/format.hpp>

// my stuff
#include "../plot_style.cpp"

using namespace std;
using boost::format;

int main(int argc, char const *argv[]) {
    setup_style();
    int bins = 200;
    double E_tot = (17.76 - 10.39) + (10.39 - 7.27); // Q1 + Q2
    TCanvas* c = new TCanvas("c", "c", 600, 600);
    TH2D* hist = new TH2D("h1", "Dalitz plot", bins, -1, 1, bins, -1, 1);
    
    for (double x = -1; x < 1; x+=2./bins) {
        for (double y = -1; y < 1; y+=2./bins) {
            if (pow(x, 2) + pow(y, 2) <= 1) {
                /*  
                    x = sqrt(3)(e2 - e3) = sqrt(3)*(E2-E3)/(E1+E2+E3) = sqrt(3)*2costheta*sqrt(Q1 Q2/3)/sqrt(3)/(E1+E2+E3)
                    E1+E2+E3 = Q1+Q2
                    so costheta = x*(E1+E2+E3)/(2*sqrt(Q1 Q2/3)) = x*(Q1+Q2)/(2*sqrt(Q1 Q2)) 
                */
                /*
                    y = 3e1 - 1 = 3E1/(E1+E2+E3) - 1 = 3*2/3*Q1/E_tot - 1 = 2Q1/E_tot - 1
                    so Q1 = (y+1)*E_tot/2
                */
                double Q1 = (y+1)*E_tot/2;
                double Q2 = E_tot - Q1;
                double theta = acos(x*E_tot/(2*sqrt(Q1*Q2)));
                int binx = hist->GetXaxis()->FindBin(x);
                int biny = hist->GetYaxis()->FindBin(y);
                hist->SetBinContent(hist->GetBin(binx, biny), theta);
            }
        }
    }

    auto circle = [] (double* x, double* par) {
        return par[0]*sqrt(1-pow(x[0], 2));
    };

    hist->SetContour(100);
    hist->GetXaxis()->SetNdivisions(-2);
    hist->GetYaxis()->SetNdivisions(-2);
    hist->GetZaxis()->SetNdivisions(-2);
    hist->GetZaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "0");
    hist->GetZaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "#pi/2");
    hist->GetZaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "#pi");
    hist->Draw("colz");

    string path = string(argv[1]) + "theta.pdf";
    c->SaveAs(path.c_str());
}