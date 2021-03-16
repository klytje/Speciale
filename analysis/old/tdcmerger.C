{
	string fname;
	cin >> fname;

	TFile output(("output/"+fname+".root").c_str(),"READ");
	t_output = (TTree*) output.Get("a");
	
	TFile match(("match/"+fname+"m.root").c_str(),"READ");
	t_match = (TTree*) match.Get("a101");
	
	TFile merged(("merged/"+fname+".root").c_str(),"RECREATE");
	TTree t_merged("m","merged tree");


	UInt_t mul;
	Int_t N;
	double pT, deltaE, eCM[3], eLab[3], exC12, thetaLab[3], thetaCM[3], phiLab[3], phiCM[3];
	t_output->SetBranchAddress("N",&N);
	t_output->SetBranchAddress("mul",&mul);
	t_output->SetBranchAddress("pT",&pT);
	t_output->SetBranchAddress("deltaE",&deltaE);
	t_output->SetBranchAddress("eCM",&eCM);
	t_output->SetBranchAddress("eLab",&eLab);
	t_output->SetBranchAddress("exC12",&exC12);
	t_output->SetBranchAddress("thetaCM",&thetaCM);
	t_output->SetBranchAddress("phiCM",&phiCM);
	t_output->SetBranchAddress("thetaLab",&thetaLab);
	t_output->SetBranchAddress("phiLab",&phiLab);
	//t_output->SetBranchAddress("d",&id);
	
	int FI[3], BI[3];
	double t[3];
	//t_output->SetBranchAddress("FI",&FI);
	//t//_output->SetBranchAddress("BI",&BI);
	//t_//output->SetBranchAddress("time",&t);
	
	UInt_t id[3], N2, id2[3];
	int FT[3], BT[3];
	t_match->SetBranchAddress("FI",&FI);
	t_match->SetBranchAddress("FT",&FT);
	t_match->SetBranchAddress("BI",&BI);
	t_match->SetBranchAddress("BT",&BT);
	t_match->SetBranchAddress("id",&id);
	//t_match->SetBranchAddress("___N___",&N2);
	//t_match->GetEntry(123454,true);

	TBranch *b_FI = t_merged.Branch("FI",&FI,"FI[3]/I");
	TBranch *b_BI = t_merged.Branch("BI",&BI,"BI[3]/I");
	TBranch *b_FT = t_merged.Branch("FT",&FT,"FT[3]/I");
	TBranch *b_BT = t_merged.Branch("BT",&BT,"BT[3]/I");
	//TBranch *b_time = t_merged.Branch("time",&t,"time[3]/D");

	TBranch *b_id = t_merged.Branch("id",&id,"id[3]/I");
	TBranch *b_eCM = t_merged.Branch("eCM",&eCM,"eCM[3]/D");
	TBranch *b_eLab = t_merged.Branch("eLab",&eLab,"eLab[3]/D");
	TBranch *b_thetaCM = t_merged.Branch("thetaCM",&thetaCM,"thetaCM[3]/D");
	TBranch *b_phiCM = t_merged.Branch("phiCM",&phiCM,"phiCM[3]/D");
	TBranch *b_thetaLab = t_merged.Branch("thetaLab",&thetaLab,"thetaLab[3]/D");
	TBranch *b_phiLab = t_merged.Branch("phiLab",&phiLab,"phiLab[3]/D");
	TBranch *b_pT = t_merged.Branch("pT",&pT,"pT/D");
	TBranch *b_deltaE = t_merged.Branch("deltaE",&deltaE,"deltaE/D");
	TBranch *b_mul = t_merged.Branch("mul",&mul,"mul/I");
	TBranch *b_exC12 = t_merged.Branch("exC12",&exC12,"exC12/D");
	//TBranch *b_wU = t2->Branch("wU",&wU);
	

	Long64_t nentries = t_output->GetEntries();
	const long consta = nentries;
	for (Long64_t i=0;i<nentries;i++) {
		int status = t_output->GetEntry(i);
		if (status == 0) break;
		t_match->GetEntry(N);
		if (mul==3 && pT<50e3 && abs(deltaE)<200) {
			//b_FI->Fill();
			//b_BI->Fill();
			t_merged.Fill();
		}
		//b_f->Fill();
		//b_wU->Fill();
	}
	
	t_merged.Write();
	//merged.WriteTObject(t_merged);
	//merged.Flush();
	//merged.Close();
}
