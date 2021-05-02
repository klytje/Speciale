void CopyDir(TDirectory *source) {
   //copy all objects and subdirs of directory source as a subdir of the current directory   
   source->ls();
   TDirectory *savdir = gDirectory;
   TDirectory *adir = savdir->mkdir(source->GetName());
   adir->cd();
   //loop on all entries of this directory
   TKey *key;
   TIter nextkey(source->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      const char *classname = key->GetClassName();
      TClass *cl = gROOT->GetClass(classname);
      if (!cl) continue;
      if (cl->InheritsFrom("TDirectory")) {
         source->cd(key->GetName());
         TDirectory *subdir = gDirectory;
         adir->cd();
         CopyDir(subdir);
         adir->cd();
      } else if (cl->InheritsFrom("TTree")) {
         TTree *T = (TTree*)source->Get(key->GetName());
         adir->cd();
         TTree *newT = T->CloneTree(-1,"fast");
         newT->Write();
      } else {
         source->cd();
         TObject *obj = key->ReadObj();
         adir->cd();
         obj->Write();
         delete obj;
     }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void fix_files() {
	vector<string> files = {"1828_000m.root", "1828_001m.root", "1828_002m.root", "1828_003m.root", "1828_004m.root",
			"1828_005m.root", "1828_006m.root", "1828_007m.root", "1828_008m.root", "1828_009m.root",
			"1828_010m.root", "1828_011m.root", "1828_012m.root", "1828_013m.root", "1828_014m.root",
			"1828_015m.root", "1828_016m.root", "1828_017m.root", "1828_018m.root", "1828_019m.root"};
	for (int i = 0; i < files.size(); i++) { 
		cout << "Reading file " << files[i] << endl;
		string input = "download/" + files[i];	
		TFile f_in(input.c_str(), "READ");
		TTree* t_in = (TTree*)f_in.Get("a101");

		string output = "match/" + files[i];
		TFile f_out(output.c_str(), "RECREATE");
		TTree* t_out = new TTree("a101", "Repaired match file.");

		UInt_t mul, N, FI[100], FT[100], BI[100], BT[100], id[100], clock, charge;
		Double_t FE[100], BE[100], theta, phi; 
		
		t_in->SetBranchAddress("mul", &mul);
		t_in->SetBranchAddress("FI", &FI);
		t_in->SetBranchAddress("FT", &FT);
		t_in->SetBranchAddress("FE", &FE);
		t_in->SetBranchAddress("BI", &BI);
		t_in->SetBranchAddress("BT", &BT);
		t_in->SetBranchAddress("BE", &BE);
		t_in->SetBranchAddress("theta", &theta);
		t_in->SetBranchAddress("phi", &phi);
		t_in->SetBranchAddress("id", &id);
		t_in->SetBranchAddress("CLOCK", &clock);
		t_in->SetBranchAddress("CHARGE", &charge);

		t_out->Branch("mul", &mul, "mul/i");
		t_out->Branch("___N___", &N, "___N___/i");
		t_out->Branch("FI", &FI, "FI[mul]/i");
		t_out->Branch("FT", &FT, "FT[mul]/i");
		t_out->Branch("FE", &FE, "FE[mul]/D");
		t_out->Branch("BI", &BI, "BI[mul]/i");
		t_out->Branch("BT", &BT, "BT[mul]/i");
		t_out->Branch("BE", &BE, "BE[mul]/D");
		t_out->Branch("theta", &theta, "theta[mul]/D");
		t_out->Branch("phi", &phi, "phi[mul]/D");
		t_out->Branch("id", &id, "id[mul]/i");
		t_out->Branch("CLOCK", &clock, "CLOCK[mul]/i");
		t_out->Branch("VCHARGE", &charge, "VCHARGE[mul]/i"); // rename "CHARGE" to "VCHARGE", as it should be
		
		int n = t_in->GetEntries();
		for (int j = 0; j < t_in->GetEntries(); j++) {
			N = j; // the original N are broken, but can easily be fixed
			t_in->GetEntry(j);
			t_out->Fill();
		}
		auto hash = f_in.Get("AUSALIB_HASH");
		auto branch = f_in.Get("AUSALIB_BRANCH");
		auto guid = f_in.Get("GUID");
		TObjArray* detectors = (TObjArray*) f_in.GetObjectChecked("detectors", "TObjArray");
		auto config = f_in.Get("MATCHER_CONFIG");

		f_in.cd("sortQC");
		TDirectory* sortqc = gDirectory;
		f_out.cd();
		CopyDir(sortqc);
		f_out.cd();

		hash->Write("AUSALIB_HASH");
		branch->Write("AUSALIB_BRANCH");
		guid->Write("GUID");
		f_out.WriteObject(detectors, "detectors");
		sortqc->Write("sortQC");
		config->Write("MATCHER_CONFIG");

		f_out.Write();
	}
}
