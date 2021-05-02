void rename_CHARGE() {
	string files[] = {"match/1828_000m.root", "match/1828_001m.root", "match/1828_002m.root", "match/1828_003m.root", "match/1828_004m.root",
			"match/1828_005m.root", "match/1828_006m.root", "match/1828_007m.root", "match/1828_008m.root", "match/1828_009m.root",
			"match/1828_010m.root", "match/1828_011m.root", "match/1828_012m.root", "match/1828_013m.root", "match/1828_014m.root",
			"match/1828_015m.root", "match/1828_016m.root", "match/1828_017m.root", "match/1828_018m.root", "match/1828_019m.root"};
	for (int i = 0; i < 20; i++) {
		cout << "Reading file " << files[i] << endl;
		TFile file(files[i].c_str(), "UPDATE");
		TTree* tree = (TTree*)file.Get("a101");
		UInt_t charge, mul, N; 
		tree->SetBranchAddress("CHARGE", &charge);
		auto b_vcharge = tree->Branch("VCHARGE", &charge);
		int n = tree->GetEntries();
		for (int j = 0; j < tree->GetEntries(); j++) {
			if (j % 1000 == 0) {
				cout << j << "/" << tree->GetEntries() << endl;	
			}
			tree->GetEntry(j);
			b_vcharge->Fill();
		}
		file.Write();
	}
}
