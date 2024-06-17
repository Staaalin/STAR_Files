/*
  root.exe -q -b -x 'muMc.C(1e6,"../*MuDst.root")'
*/
void analysis(const Char_t *inputFile="./datalist/19.6GeV_Run19/test.list", int jobindex=999999, bool isPico = false)
{
#if !defined(__CINT__)
  std::cout << "This code cannot be compiled" << std::endl;
#else
  //  gSystem->SetFPEMask(kInvalid | kDivByZero | kOverflow );
  gROOT->LoadMacro("lMuDst.C");
  
  int year = 2014;
  Int_t N = 1000000000;

  int iJob = jobindex;
  TString outputKFParticleQA = Form("KFParticleQA_%.6d.root", iJob);
  TString output = Form("output_%.6d.root", iJob);
              
  if(isPico)
    lMuDst(-1,inputFile,"ry2016,RpicoDst,mysql,kfpAna,quiet,nodefault",outputKFParticleQA.Data());
  else
    lMuDst(-1,inputFile,"ry2016,picoEvt,RMuDst,mysql,kfpAna,quiet,nodefault",outputKFParticleQA.Data());
    
  StKFParticleAnalysisMaker* kfpAnalysis = (StKFParticleAnalysisMaker*) StMaker::GetTopChain()->Maker("KFParticleAnalysis");
  if(!isPico) 
  {
    kfpAnalysis->AnalyseMuDst();
    kfpAnalysis->ProcessSignal();
    kfpAnalysis->SetCentrality(0); // 0 for no centrality selection
    kfpAnalysis->SetMCPID(26); // Lambda geant ID, 18 for Lambda, 26 for AntiLambda, also set daughter IDs
    kfpAnalysis->SetPDGPID(-3122); // Lambda PDG ID, 3122 for Lambda, -3122 for AntiLambda, also set daughter IDs
  }
  
//   if(year == 2016)
//   {
//     kfpAnalysis->UseTMVA();
//     // D0->Kpi
//     kfpAnalysis->SetTMVABinsD0("0:2:3:4:5:6:7:8:9","-1:1000");
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_0_1_pt0_80_BDT.weights.xml", 0.075, 0);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_2_2_pt0_80_BDT.weights.xml", 0.05,  1);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_3_3_pt0_80_BDT.weights.xml", 0.05,  2);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_4_4_pt0_80_BDT.weights.xml", 0.1,   3);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_5_5_pt0_80_BDT.weights.xml", 0.1,   4);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_6_6_pt0_80_BDT.weights.xml", 0.125, 5);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_7_7_pt0_80_BDT.weights.xml", 0.125, 6);
//     kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_8_8_pt0_80_BDT.weights.xml", 0.125, 7);
    
//     kfpAnalysis->RunCentralityAnalysis();
//     kfpAnalysis->SetCentralityFile("/gpfs01/star/pwg/mzyzak/Femto/Template/Centrality/centrality_2016.txt");
//   }
//   if(year == 2014)
//   {
//     kfpAnalysis->UseTMVA();
//     kfpAnalysis->SetTMVAcutsD0(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D0.xml",    0.1);
//     kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/DPlus.xml", 0.05);
//     kfpAnalysis->SetTMVAcutsDs(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/Ds.xml",    0.075);
//     kfpAnalysis->SetTMVAcutsLc(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/Lc.xml",    0.1);
//     kfpAnalysis->SetTMVAcutsD0KK( "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D0KK.xml",  0.125);
//     kfpAnalysis->SetTMVAcutsD04(  "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D04.xml",  -0.05);
//     kfpAnalysis->SetTMVAcutsBPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/BPlus.xml",-0.1);
//     kfpAnalysis->SetTMVAcutsB0(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/B0.xml",   -0.1);
    
// //     kfpAnalysis->RunCentralityAnalysis();
// //     kfpAnalysis->SetCentralityFile("/gpfs01/star/pwg/mzyzak/Femto/Template/Centrality/centrality_2014.txt");
//   }
  
  chain->Init();

  StKFParticleInterface::instance()->CleanLowPVTrackEvents();
	//StKFParticleInterface::instance()->UseHFTTracksOnly();
	StKFParticleInterface::instance()->SetSoftKaonPIDMode();
	StKFParticleInterface::instance()->SetSoftTofPidMode();
	StKFParticleInterface::instance()->SetChiPrimaryCut(10);
	StKFParticleInterface::instance()->SetChiPrimaryCut2D(3); // default >3
	StKFParticleInterface::instance()->SetChi2Cut2D(10);      // default <10
	StKFParticleInterface::instance()->SetLCut(1.0); // default >5.0
  //StKFParticleInterface::instance()->SetLdLCut2D(0.0); // default >5.0, trying 1.0 to study efficiency trend

  //Add decays to the reconstruction list
  //StKFParticleInterface::instance()->AddDecayToReconstructionList(  310);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3122);
  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3122);
  // StKFParticleInterface::instance()->AddDecayToReconstructionList( 3312);
  // StKFParticleInterface::instance()->AddDecayToReconstructionList(-3312);
  
  Long64_t nevent = N;
  if(isPico)
  {
    StPicoDstMaker* maker = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
    if (! maker) return;
    maker->SetStatus("*",1);
    TChain *tree = maker->chain();
    Long64_t nentries = tree->GetEntries();
    if (nentries <= 0) return;
    nevent = TMath::Min(nevent,nentries);
    cout << nentries << " events in chain " << nevent << " will be read." << endl;
  }else{
  StMuDstMaker* maker = (StMuDstMaker *) StMaker::GetTopChain()->Maker("MuDst");
  if (! maker) return;
    maker->SetStatus("*",1);
    TChain *tree = maker->chain();
    Long64_t nentries = tree->GetEntries();
    if (nentries <= 0) return;
    nevent = TMath::Min(nevent,nentries);
    nevent = nentries;
    cout << nentries << " events in chain " << nevent << " will be read." << endl;
  }
  
  // chain->EventLoop(nevent);
  time_t time_start;
	time_t time_now;
	time(&time_start);

	Int_t istat = 0, i = 1;
	while(i<=nevent && istat!=2) {
		//if(i%5000==0){cout << endl; cout << "== Event " << i << " start ==" << endl;}
		chain->Clear();
		istat = chain->Make(i);
		if(i%100==0) {
			cout << endl; 
			cout << "== Event " << i << " finish == " << flush;
			time(&time_now);
			int time_diff = (int)difftime(time_now, time_start);
			cout << time_diff/60 << "min " << time_diff%60 << "s: " << 1.0*time_diff/i << "s/event" << endl;
		}
		if (istat == 2)
			cout << "Last  event processed. Status = " << istat << endl;
		if (istat == 3)
			cout << "Error event processed. Status = " << istat << endl;
		i++;
	}


	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nevent << endl;
	cout << "****************************************** " << endl;

	delete chain;
#endif
  
}
