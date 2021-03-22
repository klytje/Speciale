void test() {
    char filter[] = "mul==3";
    
    gStyle->SetOptTitle(0);
    
    TFile *file = new TFile("aligned.root", "READ");
    TTree *t = (TTree *)file->Get("a");

    TCanvas *c = new TCanvas("c", "c");
    TH1D *h1 = new TH1D("h1", "h1", 80, 65390, 65470);
    TH1D *h2 = new TH1D("h2", "h2", 80, 65390, 65470);

    h1->SetLineColor(kAzure+1);
    h2->SetLineColor(kOrange+1);

    t->Draw("FT-20 >> h2", filter);
    
    if (!strncmp(filter, "mul==3", 6)) {
        t->Draw("BT+269 >> h1", filter);
    } else {
        t->Draw("BT+88 >> h1", filter);
    }
    
    h1->Scale( 1./h1->Integral());
    h2->Scale( 1./h2->Integral());

    h2->Draw("HIST L");
    h1->Draw("HIST L same");

    c->cd();  // c is the TCanvas
    TPad *pad = new TPad("all","all",0,0,1,1);
    pad->SetFillStyle(4000);  // transparent
    pad->Draw();
    pad->cd();
    TLatex *lat = new TLatex();
    lat->DrawLatexNDC(.4,.95,"BT & FT (aligned)");

    TCanvas *c1 = new TCanvas("c1", "c1");
    if (!strncmp(filter, "mul==3", 6)) {
        t->Draw("dt >> (70, -320, -250)", filter);
    } else {
        t->Draw("dt >> (70, -140, -70)", filter);
    }

    c1->cd();  // c is the TCanvas
    TPad *pad1 = new TPad("all","all",0,0,1,1);
    pad1->SetFillStyle(4000);  // transparent
    pad1->Draw();
    pad1->cd();
    TLatex *lat1 = new TLatex();
    lat1->DrawLatexNDC(.4,.95,"dt (aligned)");

    TCanvas *c2 = new TCanvas("c2", "c2");
    TFile *file2 = new TFile("merged/all.root", "READ");
    TTree *t2 = (TTree *)file2->Get("tree");
    t2->Draw("(BT-FT)*1e-3 >> (70, -140, -70)", filter);

    c2->cd();  // c is the TCanvas
    TPad *pad2 = new TPad("all","all",0,0,1,1);
    pad2->SetFillStyle(4000);  // transparent
    pad2->Draw();
    pad2->cd();
    TLatex *lat2 = new TLatex();
    lat2->DrawLatexNDC(.4,.95,"dt (raw)");}
