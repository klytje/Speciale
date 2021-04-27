{
	string original,output;
	cin >> original >> output;

	TFile f1(original.c_str(),"READ");
	t1 = (TTree*) f1.Get("sample");

	TFile f2(output.c_str(),"UPDATE");
	t2 = (TTree*) f2.Get("a");

	for (auto name : {"J"}) {
		std::vector<int>* vi;
		f1.GetObject(name,vi);
		f2.WriteObject(vi,name);
	}
	
	Int_t N;
	Double_t wU;
	std::vector<std::vector<double>>* factors;

	t2->SetBranchAddress("N",&N);

	t1->SetBranchAddress("f",&factors);
	t1->SetBranchAddress("wU",&wU);
	TBranch *b_f = t2->Branch("f",&factors);
	TBranch *b_wU = t2->Branch("wU",&wU);

	Long64_t nentries = t2->GetEntries();
	for (Long64_t i=0;i<nentries;i++) {
		t2->GetEntry(i);
		t1->GetEntry(N);
		b_f->Fill();
		b_wU->Fill();
	}

	f2.WriteTObject(t2);
	f2.Flush();
	f2.Close();
}
