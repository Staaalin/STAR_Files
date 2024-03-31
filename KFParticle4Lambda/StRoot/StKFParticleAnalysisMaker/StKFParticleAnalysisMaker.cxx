#include <fstream>
#include <iostream>

#include "StKFParticleAnalysisMaker.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "StMessMgr.h"
#include "TSystem.h"
#include <algorithm>
#include <TMath.h>
#include <map>

#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"
#include "KaonPID.h"
#include "ProtonPID.h"
#include "PionPID.h"

#include "StTrackHelix.h"
#include "StLambdaDecayPair.h"
#include "MyToolkit.h"

#define pi                 TMath::Pi()
#define OmegaPdgMass	   1.67245
//#define OmegaMassSigma     0.0021
#define LambdaPdgMass      1.11568
#define ProtonPdgMass      0.938272
#define PionPdgMass        0.139570
#define KaonPdgMass		   0.493677
#define LambdaPdg          3122
#define XiPdg              3312
#define phiPdg			   333
#define OmegaPdg           3334
#define KaonPdg			   321
#define ProtonPdg          2212
#define KsPdg			   310
#define Xi1530Pdg		   3324
#define PionPdg           -211

// #define DEBUGGING


//-----------------------------------------------------------------------------
ClassImp(StKFParticleAnalysisMaker)

//-----------------------------------------------------------------------------
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char* name, const char* outName)
: StMaker(name)
{
	mOutName = outName;

	mRun    =0;
	mEnergy =0;
	mListDir="./";
}

//----------------------------------------------------------------------------- 
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker() {
	SafeDelete(KFParticleInterface);
	SafeDelete(KFParticlePerformanceInterface);
}

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::openFile()
{    
	// ======= StRefMultCorr ======= //
	// set up StRefMultCorr for centrality
	mRefMultCorr = CentralityMaker::instance()->getRefMultCorr() ;
	// ======= StRefMultCorr end ======= //

	cout << "-----------------------------------" << endl;
	cout << "------- user's files loaded -------" << endl;
	cout << "-----------------------------------" << endl;

	return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::Init() {

        cout << mOutName << endl;
        const char *sDir = mOutName.Data();
        int lDir = mOutName.Length();
        int iDir = lDir-6;
        while(!std::isdigit(sDir[iDir])) iDir --;
        while(std::isdigit(sDir[iDir]))  iDir --;
        mJob = std::atoi(&sDir[iDir+1]);
        cout << "current job id: " << mJob << endl;

	PI = M_PI;
	twoPI = 2*M_PI;

	badList.clear();
	runList.clear();

	if(!readRunList())return kStFatal;
	if(!readBadList())return kStFatal;

	DeclareHistograms();
	Int_t openFileStatus = openFile();
	if(openFileStatus == kStFatal) return kStFatal;

	dcatoPV_hi = 3.0;
	// pid boundary
	pT_lo = 0.2;
	pT_hi = 2.0;
	
	// for EP
	pT_asso_lo = 0.15;
	pT_asso_hi = 2.0;
	pT_trig_lo = 0.2;
	pT_trig_hi = 2.0;
	eta_trig_cut = 1.0;
	
	// pion cut
	pion_pT_lo = 0.2;
	pion_pT_hi = 1.2;
	pion_pT_TOFth = 0.4;
	pion_m2_lo = -0.1;
	pion_m2_hi = 0.1;

	// proton cut
	proton_pT_lo = 0.3;
	proton_pT_hi = 1.8;
	proton_pT_TOFth = 0.6;
	proton_m2_lo = 0.75;
	proton_m2_hi = 1.1;

	TFile *f = GetTFile(); // These two lines need to be HERE (though I don't know /why/)- don't throw in another function
	if(f){f->cd(); BookVertexPlots();}

	return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::Finish() {
	if(mOutName!="") {
		TFile *fout = new TFile(mOutName.Data(),"RECREATE");
		fout->cd();
		WriteHistograms();
		fout->Close();
	}
	return kStOK;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::DeclareHistograms() {

	const int nRuns=runList.size();
	int mRunL=runList.at(0);
	int mRunH=runList.at(nRuns-1);
	int mDayL=(int) ((mRunL)/1000%1000);
	int mDayH=(int) ((mRunH)/1000%1000)+1;
	const int nDays=mDayH-mDayL;
	int nTmps=nDays;//>40?40:nDays;
	mStps=ceil((nDays*1.0)/(nTmps*1.0));
	const int nTims=ceil((nDays*1.0)/(mStps*1.0));

	cout<<nDays<<" "<<mDayL<<" "<<mDayH<<" "<<nTims<<" "<<mStps<<endl;

	hNRefMult = new TH1F("RefMult" , "Reference Multiplicity" , 1000, 0.0, 1000.0 ) ;
	hNRefMultA= new TH1F("RefMultA", "Reference MultiplicityA", 1000, 0.0, 1000.0 ) ;
	hNRefMultB= new TH1F("RefMultB", "Reference MultiplicityB", 1000, 0.0, 1000.0 ) ;
	hVertexXY = new TH2F("VertexXY", "Vertex XY Position", 200, -10.0, 10.0 , 200, -10., 10 ) ;
	hVertexZ  = new TH1F("VertexZ" , "Event Vertex Z Position", 200, -100.0, 100.0 ) ;
	hVertex2D = new TH2F("Vertex2D", "VertexZ vs VPD Vz", 200, -100.0, 100.0 , 200, -100., 100 ) ;
	hDiffVz   = new TH1F("VertexZdiff" , "VertexZ-VPDVz diff", 100, -10.0, 10.0 ) ;
	hcent     = new TH1F("centrality","centrality"  ,nCent,0.,nCent);
	hcentw    = new TH1F("centralityw","centralityw",nCent,0.,nCent);

	hcentRefM = new TProfile("hcentRefM","hcentRefM",nCent,0.,nCent,0,1000);
	hcentRefW = new TProfile("hcentRefW","hcentRefW",nCent,0.,nCent,0,1000);

    buffer_size = 5000000;
    hadronTree = new TTree("hadronTree", "Tree_STAR");
    // hadronTree->Branch("buffer_size"       ,&buffer_size         ,"buffer_size/I"                       );
    hadronTree->Branch("PDGMult"            ,&PDGMult             ,"PDGMult/I"                           );
    hadronTree->Branch("refMult"            ,&CrefMult            ,"refMult/I"                           );
    hadronTree->Branch("grefMult"           ,&CgrefMult           ,"grefMult/I"                          );
    hadronTree->Branch("EventID"            ,&evtID               ,"EventID/I"                           );
    hadronTree->Branch("RunID"              ,&runID               ,"RunID/I"                             );
    hadronTree->Branch("PDG"                ,&PDG                 );
    hadronTree->Branch("mix_px"             ,&px                  );
    hadronTree->Branch("mix_py"             ,&py                  );
    hadronTree->Branch("mix_pz"             ,&pz                  );
    hadronTree->Branch("InvariantMass"      ,&InvariantMass       );

	// Used for QA
    hadronTree->Branch("dEdx"               ,&QA_dEdx              );
    hadronTree->Branch("m2"                 ,&QA_m2                );
    hadronTree->Branch("dcatopv"            ,&QA_DCA_V0_PV         );
    hadronTree->Branch("hasTOF"             ,&QA_hasTOF            );
    hadronTree->Branch("nSigmaProton"       ,&QA_nSigmaProton      );
    hadronTree->Branch("nSigmaPion"         ,&QA_nSigmaPion        );
    hadronTree->Branch("nSigmaKaon"         ,&QA_nSigmaKaon        );
    hadronTree->Branch("zTOF_proton"        ,&QA_zTOF_proton       );
    hadronTree->Branch("zTOF_pion"          ,&QA_zTOF_pion         );
    hadronTree->Branch("zTOF_kaon"          ,&QA_zTOF_kaon         );
	hadronTree->Branch("IfConfuse"          ,&QA_IfConfuse         );
	hadronTree->Branch("Decay_Length"       ,&QA_Decay_Length      );
	hadronTree->Branch("Chi2"               ,&QA_Chi2              );
	hadronTree->Branch("IfBadReconstructed" ,&QA_IfBadReconstructed);
	hadronTree->Branch("DCA_Daughters"      ,&QA_DCA_Daughters     );
	
	hdEdx_pQ = new TH2D("hdEdx_p_NO_CUT","dE/dx vs. p*Q without cut",2000,10,10,2000,0,20);
	hdEdx_pQ->GetXaxis()->SetTitle("p*Q [GeV]");
	hdEdx_pQ->GetYaxis()->SetTitle("dE/dx [keV/cm]");
	
	hdEdx_pQ_1cut = new TH2D("hdEdx_p_1_CUT","dE/dx vs. p*Q HITS cut",2000,10,10,2000,0,20);
	hdEdx_pQ_1cut->GetXaxis()->SetTitle("p*Q [GeV]");
	hdEdx_pQ_1cut->GetYaxis()->SetTitle("dE/dx [keV/cm]");
	
	hdEdx_pQ_2cut = new TH2D("hdEdx_p_2_CUT","dE/dx vs. p*Q HITS & PID cut",2000,10,10,2000,0,20);
	hdEdx_pQ_2cut->GetXaxis()->SetTitle("p*Q [GeV]");
	hdEdx_pQ_2cut->GetYaxis()->SetTitle("dE/dx [keV/cm]");
	
	hLN_M = new TH2D("hLN_M","",400,0,4,600,0,60);
	hLN_M->GetXaxis()->SetTitle("Mass [GeV]");
	hLN_M->GetYaxis()->SetTitle("HM");
	
	hXY = new TH2D("h_XY","h_ X&Y of PicoDST Tracks",200,0,10,200,0,10);
	hXY->GetXaxis()->SetTitle("X [cm]");
	hXY->GetYaxis()->SetTitle("Y [cm]");

	hHXY = new TH2D("h_HXY","h_ X&Y of StHelix",     200,0,10,200,0,10);
	hHXY->GetXaxis()->SetTitle("X [cm]");
	hHXY->GetYaxis()->SetTitle("Y [cm]");

	hHM_Chi2 = new TH2D("h_HM_Chi2","Mass vs. Chi2",     1000,0,10,500,0,10);
	hHM_Chi2->GetXaxis()->SetTitle("Mass [GeV]");
	hHM_Chi2->GetYaxis()->SetTitle("Chi2");

	hHM_ParentDCA = new TH2D("hH_M_ParentDCA","The DCA between parent particles vs. Mass",     2000,0,5,500,0,5);
	hHM_ParentDCA->GetXaxis()->SetTitle("Mass [GeV]");
	hHM_ParentDCA->GetYaxis()->SetTitle("DCA [cm]");

	const int tPDGList[]      = {   3122  ,   -3122   ,   3334  ,  -3334};
	const TString tNameList[] = {"Lambda" , "Lambdab" , "Omega" , "Omegab"};
	for (int Itr = 0;Itr < PDG2NameSize;Itr++){
		PDGList[Itr] = tPDGList[Itr];NameList[Itr] = tNameList[Itr];

		TString HistName1 = "HM_";
		TString HistName2 = "The Mass of ";
		HistName1 += NameList[Itr];HistName2 += NameList[Itr];
		H_ALL_NO_CUT[Itr] = new TH1F(HistName1,HistName2,6000,0,3);
		H_ALL_NO_CUT[Itr]->GetXaxis()->SetTitle("Mass [GeV]");

		HistName1 = "HM_DaughtersDCA_";
		HistName2 = "The Mass of ";
		HistName1 += NameList[Itr];HistName2 += NameList[Itr];
		HistName2 += " those DayghtersDCA < 0.6 [cm]";
		H_DaughterDCA[Itr] = new TH1F(HistName1,HistName2,6000,0,3);
		H_DaughterDCA[Itr]->GetXaxis()->SetTitle("Mass [GeV]");
	}

	Recorded_events = 0;

	cout << "-----------------------------------------" << endl;
	cout << "------- histograms & tree claimed -------" << endl;
	cout << "-----------------------------------------" << endl;

	return;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::WriteHistograms() {

	///////////////////
	hadronTree ->Write();

	//-- Used for  test --
	// hNRefMult ->Write();  
	// hNRefMultA->Write();  
	// hNRefMultB->Write();  
	// hVertexXY ->Write();  
	// hVertexZ  ->Write();  
	// hVertex2D ->Write();
	// hDiffVz   ->Write();
	// hcent     ->Write();  
	// hcentw    ->Write();  

	// hcentRefM ->Write();
	// hcentRefW ->Write();

	// hdEdx_pQ->Write();
	// hdEdx_pQ_1cut->Write();
	// hdEdx_pQ_2cut->Write();

	// hLN_M->Write();
	// hXY->Write();
	// hHXY->Write();
	// hHM_Chi2->Write();
	// hHM_ParentDCA->Write();
	
	for (int i=0;i<PDG2NameSize;i++){
		H_ALL_NO_CUT[i]->Write();

		H_DaughterDCA[i]->Write();

	}

	return;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::setRunEnergyAndListDir(int run, double energy,char ListDir[256])
{
	mRun    =run;
	mEnergy =energy;
	mListDir=ListDir;
}

//-----------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::readBadList()
{ 
	if(0) return kTRUE;
	TString inf=mListDir + "/badList/";
	inf += Form("badrun%dList%.1f.list",mRun,mEnergy);
	ifstream inrun; 
	inrun.open(inf);
	if ( inrun.fail() ) {
		cout<< "cannot open " << inf.Data() << endl;
		return kFALSE;
	}
	Int_t runid;
	while ( inrun >> runid ) { badList.push_back(runid); }
	inrun.close();
	sort(badList.begin(),badList.end());

	vector<int>::iterator it;
	it = std::unique (badList.begin(), badList.end());
	badList.resize( std::distance(badList.begin(),it) );

	cout <<"badrun list :" <<inf.Data() << " loaded." << endl;
	cout <<"Total       :" <<badList.size()<< " bad runs. "<< endl;
	return kTRUE;
}

//------------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::readRunList()
{ 
	if(0) return kTRUE;
	TString inf=mListDir + "/runList/";
	inf += Form("run%dList%.1f.list",mRun,mEnergy);
	ifstream inrun; 
	inrun.open(inf);
	if ( inrun.fail() ) {
		cout<< "cannot open " << inf.Data() << endl;
		return kFALSE;
	}
	Int_t runid;
	while ( inrun >> runid ) { runList.push_back(runid); }
	inrun.close();
	sort(runList.begin(),runList.end());

	vector<int>::iterator it;
	it = std::unique (runList.begin(), runList.end());
	runList.resize( std::distance(runList.begin(),it) );

	cout <<"Run list :" <<inf.Data() << " loaded. "<< endl;
	cout <<"Total    :" <<runList.size()<< " runs. "<< endl;

	if(runList.size()<1){cout<<"no run number found!!!"<<endl; return kFALSE;}

	return kTRUE;
}

//-----------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::removeBadID(int runnumber)
{
	for (std::vector<int>::iterator it=badList.begin(); it!=badList.end(); ++it) { 
		if(runnumber==*it) {
			return kTRUE;
		}
	}
	return kFALSE;
}  

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::CheckrunNumber(int runnumber)
{    
	int pointer=-999; 
	int id=0; 
	for (std::vector<int>::iterator it=runList.begin(); it!=runList.end(); ++it) { 
		if(runnumber==*it) pointer=id;
		id++;
	}

	if(pointer==-999) cout << "Run number are not found! "<< runnumber << endl;
	return pointer;
} 

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::findCentrality(int mult)
{
	// NOT in use
	int centrality=0;
	float centFull[9] = {30,46,65,89,118,154,200,262,305};
	if      (mult>centFull[8]) centrality=9;
	else if (mult>centFull[7]) centrality=8;
	else if (mult>centFull[6]) centrality=7;
	else if (mult>centFull[5]) centrality=6;
	else if (mult>centFull[4]) centrality=5;
	else if (mult>centFull[3]) centrality=4;
	else if (mult>centFull[2]) centrality=3;
	else if (mult>centFull[1]) centrality=2;
	else if (mult>centFull[0]) centrality=1;

	return centrality;
}

//----------------------------------------------------------------------------- 
void StKFParticleAnalysisMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::Make() 
{
	// cout<<"Start Make"<<endl;
	PicoDst = StPicoDst::instance(); 		
	StPicoDst* mPicoDst = PicoDst;
	if(!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStOK;
	}

	//     pass event  
	/////////////////////////////////////////////////////////
	StPicoEvent* mEvent= (StPicoEvent*) mPicoDst->event(); 
	if(!mEvent)return kStOK;

	runID    = mEvent->runId();
	evtID    = mEvent->eventId();
	#ifdef DEBUGGING
	std::cout << "Reading evtID : " << mEvent->eventId() <<std::endl;
	std::cout << "Read evtID : " << evtID <<std::endl;
	#endif
	const int  refMult  = mEvent->refMult();
	#ifdef DEBUGGING
	std::cout << "Reading refMult : " << mEvent->refMult() <<std::endl;
	std::cout << "Read refMult : " << refMult <<std::endl;
	#endif
	const int grefMult  = mEvent->grefMult();
	const int  ranking  = mEvent->ranking();
	const int tofMatch  = mEvent->nBTOFMatch();

	const double magnet = mEvent->bField();

	// int SizeOf_Recorded_runID = Recorded_runID.size();
	// if (SizeOf_Recorded_runID == 0){
	// 	Recorded_runID.push_back(runID);
	// 	TString LogName = "/star/data01/pwg/svianping/JobID/id";
	// 	LogName += mJob;
	// 	LogName += ".log";
	// 	ofstream outfile(LogName, ios::out | ios::trunc);
	// 	outfile << runID << endl;
	// 	outfile.close();
	// }else{
	// 	for(int i=0;i<SizeOf_Recorded_runID;i++){
	// 		if (runID == Recorded_runID[i]) {break;}
	// 		if (i == SizeOf_Recorded_runID - 1) {
	// 			Recorded_runID.push_back(runID);
	// 			TString LogName = "/star/data01/pwg/svianping/JobID/id";
	// 			LogName += mJob;
	// 			LogName += ".log";
	// 			ofstream outfile(LogName, ios::out | ios::app);
	// 			outfile.seekp(0, ios::end);
	// 			outfile << runID << endl;
	// 			outfile.close();
	// 		}
	// 	}
	// }

	// cout<<"Here OK"<<endl;

	// if(       (!mEvent->isTrigger(610001))
	// 		&&(!mEvent->isTrigger(610011))
	// 		&&(!mEvent->isTrigger(610021))
	// 		&&(!mEvent->isTrigger(610031))
	// 		&&(!mEvent->isTrigger(610041))
	// 		&&(!mEvent->isTrigger(610051))
	// )return kStOK; 
	if ((!mEvent->isTrigger(530003))
	  &&(!mEvent->isTrigger(530002))
	  &&(!mEvent->isTrigger(530806))
	  &&(!mEvent->isTrigger(530101))
	  &&(!mEvent->isTrigger(530102))
	  &&(!mEvent->isTrigger(530201))
	  &&(!mEvent->isTrigger(530202))
	  &&(!mEvent->isTrigger(530213))
	  &&(!mEvent->isTrigger(530851))// bbc
	  &&(!mEvent->isTrigger(530852))// zdce
	  &&(!mEvent->isTrigger(530853))// vpd-10
	  &&(!mEvent->isTrigger(530854))// vpd
	  &&(!mEvent->isTrigger(530855))// zdc
	  &&(!mEvent->isTrigger(530856))// bbcnotac
	  &&(!mEvent->isTrigger(530857))// zdc-notac
	  &&(!mEvent->isTrigger(530861))// bbc
	  &&(!mEvent->isTrigger(530866))// bbcnotac
	  &&(!mEvent->isTrigger(2)) //17132063
	  &&(!mEvent->isTrigger(3))
	  &&(!mEvent->isTrigger(4)) //VPD-5
	  &&(!mEvent->isTrigger(6)) //highMult-VPD-5
	  &&(!mEvent->isTrigger(15))
	  &&(!mEvent->isTrigger(16)) //BHT2-VPD-30
	  &&(!mEvent->isTrigger(17))
	  &&(!mEvent->isTrigger(54))
	  &&(!mEvent->isTrigger(55))
	  &&(!mEvent->isTrigger(56))
	  &&(!mEvent->isTrigger(57))
	  &&(!mEvent->isTrigger(58))
	  &&(!mEvent->isTrigger(59))
	  &&(!mEvent->isTrigger(61)) //vpd-30
	  )return kStOK;

	// cout<<"Trigger OK"<<endl;

	const TVector3 Vertex3D=mEvent->primaryVertex();
	const double VertexX = Vertex3D.x(); 
	const double VertexY = Vertex3D.y(); 
	const double VertexZ = Vertex3D.z(); 
	const double VertexR = sqrt(VertexX*VertexX + VertexY*VertexY);
	const double vpdVz   = mEvent->vzVpd();

	//event cut
	//if(refMult <=2 || refMult > 1000) return kStOK;
	if(removeBadID(runID)) return kStOK;            
	if(mRefMultCorr->isBadRun(runID)) return kStOK; // reject bad run of StRefMultCorr
	if(!mRefMultCorr->passnTofMatchRefmultCut(1.*refMult, 1.*tofMatch)) return kStOK; // reject pileup of StRefMultCorr

	hVertex2D ->Fill(VertexZ,vpdVz);
	hDiffVz   ->Fill(VertexZ-vpdVz); 

	const double DVz = VertexZ-vpdVz;

	if(fabs(VertexZ) > 80) return kStOK; 
	if(sqrt(pow(VertexX,2.)+pow(VertexY,2.))>2.0) return kStOK; 
	if(fabs(VertexZ-vpdVz)>3.) return kStOK;       // no vpd cut in low energy?

	//check run number
	int runnumberPointer = -999;
	runnumberPointer = CheckrunNumber(runID);
	if(runnumberPointer==-999) return kStOK;

	int dayPointer = (int)((runID)/1000%1000);
	int mRunL=runList.at(0);
	int mDayL=(int) ((mRunL)/1000%1000);
	dayPointer -= mDayL;
	int timePointer = dayPointer/mStps;

	//StRefMultCorr
	mRefMultCorr->init(runID);
	mRefMultCorr->initEvent(refMult,VertexZ,mEvent->ZDCx());
	double refmultWght = mRefMultCorr->getWeight();
	double refmultCorr = mRefMultCorr->getRefMultCorr() ;
	int centrality     = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
	//double mWght    = 1;
	//double mult_corr= refMult;
	//int centrality  = findCentrality(mult_corr)-1;  // 0-8 
	if( centrality<0||centrality>=(nCent-1)) return kStOK;
	int cent = centrality+1;  

	double mWght = refmultWght;
	double mult_corr = refmultCorr;

	///////////////////////////
	hNRefMult ->Fill(grefMult);
	hNRefMultA->Fill(mult_corr);
	hNRefMultB->Fill(mult_corr,mWght);
	hVertexXY ->Fill(VertexX,VertexY);
	hVertexZ  ->Fill(VertexZ); 
	hcent     ->Fill(cent);         // 1-9
	hcentw    ->Fill(cent,mWght);   // 1-9

	///////////////
	hcentRefM    ->Fill(cent,mult_corr);     
	hcentRefW    ->Fill(cent,mult_corr,mWght);  
	hcentRefM    ->Fill(0.,mult_corr);     
	hcentRefW    ->Fill(0.,mult_corr,mWght);  
	///////////////

// ======= KFParticle ======= //
	std::vector<StLambdaDecayPair> KFParticleLambdaDecayPair;
	std::vector<KFParticle> ParticleVec , OmegaVec , LambdaVec;
	ParticleVec.resize(0);OmegaVec.resize(0);LambdaVec.resize(0);

	float dcatoPV_hi = 3.0; // Upper limit of DCA to PVs
	float pT_trig_lo = 0.2;
	float pT_trig_hi = 2.0;
	float eta_trig_cut = 1.0;

	// cout<<"In InterfaceCantProcessEvent"<<endl;
	SetupKFParticle();
	// if (InterfaceCantProcessEvent) return;
	// cout<<"InterfaceCantProcessEvent OK"<<endl;
	CrefMult = refMult;CgrefMult = grefMult;
	PDG.resize(0);px.resize(0);py.resize(0);pz.resize(0);InvariantMass.resize(0);
	// QA
	QA_dEdx.resize(0);QA_DCA_V0_PV.resize(0);QA_m2.resize(0);QA_nSigmaProton.resize(0);
	QA_nSigmaPion.resize(0);QA_nSigmaKaon.resize(0);QA_zTOF_proton.resize(0);QA_zTOF_pion.resize(0);QA_zTOF_kaon.resize(0);
	QA_hasTOF.resize(0);QA_IfConfuse.resize(0);QA_Decay_Length.resize(0);QA_Chi2.resize(0);QA_IfBadReconstructed.resize(0);
	QA_DCA_Daughters.resize(0);

	// std::vector<int> Constructed_KFParticle_Vec_index; Constructed_KFParticle_Vec_index.resize(0);
	// for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++){
	// 	KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
	// 	if ((particle.GetPDG() == OmegaPdg) || 
	// 		(particle.GetPDG() == LambdaPdg))
	// 	{
	// 		for(int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++){ 
	// 			const int daughterId = particle.DaughterIds()[iDaughter]; 
	// 			Constructed_KFParticle_Vec_index.push_back(daughterId);
	// 		}  // iDaughter
	// 	}
	// }

	Omega_Omegab_Num = 0;
	for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++){ 
		KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];

		bool IfWellConstrcuted = true;
		bool IfHelix = false;

		#ifdef DEBUGGING
		std::cout << "Parsing refMult : " << refMult <<std::endl;
		std::cout << "Parsed CrefMult : " << CrefMult <<std::endl;
		#endif
		
		if ((fabs(particle.GetPDG()) != OmegaPdg) && (fabs(particle.GetPDG()) != LambdaPdg)) {continue;}

		//SCHEME 1: reconstruction of V0, the parent particle
		if (particle.NDaughters() != 2){cout<<"FUCK! particle.NDaughters() = "<<particle.NDaughters()<<endl;}
		int iTrack,kTrack;
		if (particle.GetPDG() == 3334 || particle.GetPDG() == -3334) {
			cout<<"Particle ID = "<<particle.GetPDG()<<endl;
		}
		for (int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++){
			const int daughterId = particle.DaughterIds()[iDaughter]; 
			const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
			if (particle.GetPDG() == 3334 || particle.GetPDG() == -3334) {
				cout<<"daughter ID = "<<daughter.GetPDG()<<endl;
				if (daughter.GetPDG() == -1) {
					for (int jDaughter=0; jDaughter < daughter.NDaughters(); jDaughter++){
						int DdaughterId = daughter.DaughterIds()[jDaughter];
						KFParticle Ddaughter = KFParticleInterface->GetParticles()[DdaughterId];
						cout<<"Grand daughter ID = "<<Ddaughter.GetPDG()<<endl;
					}
				}
			}
			const int globalTrackId = daughter.DaughterIds()[0];
			Int_t nTracks = mPicoDst->numberOfTracks();
			Int_t iTrackStart = globalTrackId - 1;
			if (globalTrackId >= nTracks) {iTrackStart = nTracks - 1;}
			for (Int_t jTrack = iTrackStart;jTrack >= 0;jTrack--){
				StPicoTrack *track = mPicoDst->track(jTrack);
				if (track->id() == globalTrackId){
					// int TrackPDG = TrackID(track , Vertex3D , magnet , false);
					if (iDaughter == 0){iTrack = jTrack;}
					if (iDaughter == 1){kTrack = jTrack;}
					break;
				}
			}
		}

		StPicoTrack* mTrackI = (StPicoTrack*)mPicoDst->track(iTrack);
		StPicoTrack* mTrackK = (StPicoTrack*)mPicoDst->track(kTrack);
		// StPicoPhysicalHelix cTrackI = mTrackI->helix(magnet);
		// StPicoPhysicalHelix cTrackK = mTrackK->helix(magnet);
		
		// for (int Itr = 0;Itr < PDG2NameSize;Itr++){
		// 	if (particle.GetPDG() == PDGList[Itr]){
		// 		pair<Double_t , Double_t>RV = cTrackI.pathLengths(cTrackK , 0.1 , 0.1);
		// 		// TVector3 LTrackI = LocAfterTransfer(cTrackI , RV.first);
		// 		// TVector3 LTrackK = LocAfterTransfer(cTrackK , RV.second);
		// 		TVector3 LTrackI = cTrackI.at(RV.first);
		// 		TVector3 LTrackK = cTrackK.at(RV.second);
		// 		double TrackDCA = DistanceBetween(LTrackI , LTrackK);
		// 		if (TrackDCA < 0.6) {
		// 			H_ALL_NO_CUT[Itr]->Fill(particle.GetMass());
		// 			H_DaughterDCA[Itr]->Fill(particle.GetMass());
		// 			QA_IfBadReconstructed.emplace_back(0);
		// 		}
		// 		else
		// 		{
		// 			H_ALL_NO_CUT[Itr]->Fill(particle.GetMass());
		// 			IfWellConstrcuted = false;
		// 			QA_IfBadReconstructed.emplace_back(1);
		// 		}
				
		// 		hHM_ParentDCA->Fill(particle.GetMass(),TrackDCA);
		// 		QA_DCA_Daughters.emplace_back(TrackDCA);
		// 		// cout<<"DCA to Parent = "<<DistanceBetween(LTrackI , LTrackK)<<endl;
		// 		break;
		// 	}
		// }

		for (int Itr = 0;Itr < PDG2NameSize;Itr++){
			if (particle.GetPDG() == PDGList[Itr]){
				if (StKFParticleAnalysisMaker::IfGoodDaughterDCA(mPicoDst , iKFParticle , magnet , 0.6 , 0.6)){
					H_ALL_NO_CUT[Itr]->Fill(particle.GetMass());
					H_DaughterDCA[Itr]->Fill(particle.GetMass());
					QA_IfBadReconstructed.emplace_back(0);
				}
				else
				{
					H_ALL_NO_CUT[Itr]->Fill(particle.GetMass());
					IfWellConstrcuted = false;
					QA_IfBadReconstructed.emplace_back(1);
				}
				
				// hHM_ParentDCA->Fill(particle.GetMass(),TrackDCA);
				QA_DCA_Daughters.emplace_back(-1.0);
				// cout<<"DCA to Parent = "<<DistanceBetween(LTrackI , LTrackK)<<endl;
				break;
			}
		}
		
		// StPicoTrack* mTrackI = (StPicoTrack*)mPicoDst->track(iTrack);
		// StPicoTrack* mTrackK = (StPicoTrack*)mPicoDst->track(kTrack);
		TVector3 xv0, op1, op2;
		double dca1to2 = closestDistance(mTrackI, mTrackK, magnet, Vertex3D, xv0, op1, op2);
		TVector3 pv0 = op1 + op2;
		TVector3 xv0toPV = xv0 - Vertex3D;
		double rdotp = xv0toPV.Dot(pv0);
		double dcav0toPV = rdotp*rdotp/pv0.Mag2();
		dcav0toPV = sqrt(xv0toPV.Mag2() - dcav0toPV);
		double v0decaylength = xv0toPV.Mag();
		double v0cosrdotp = rdotp/v0decaylength/pv0.Mag();// cout<<"SCHEME 1: DecayLength = "<<v0decaylength<<";  ";
		//SCHEME 2:
		KFParticle tempParticle(particle);
		float l,dl;
		KFParticle pv(KFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
		// pv += particle;
		tempParticle.SetProductionVertex(pv);
		tempParticle.GetDecayLength(l, dl);// cout<<"SCHEME 2: DecayLength = "<<l<<";  ";if (fabs(v0decaylength/l)>1.15 || fabs(v0decaylength/l)<0.95){cout<<particle.GetPDG()<<"  "<<particle.GetMass()<<endl;}else{cout<<" "<<endl;}
		if (particle.GetPDG() == OmegaPdg ) { OmegaVec.push_back(particle);Omega_Omegab_Num ++;}
		if (particle.GetPDG() == -1*OmegaPdg ) {Omega_Omegab_Num ++;}
		if (particle.GetPDG() == LambdaPdg) {LambdaVec.push_back(particle);}
		ParticleVec.push_back(particle);

		if (IfWellConstrcuted) {
			PDG.emplace_back(particle.GetPDG());QA_Decay_Length.emplace_back(v0decaylength);QA_DCA_V0_PV.emplace_back(dcav0toPV);QA_IfBadReconstructed.emplace_back(0);
		}
		else{
			// continue;
			PDG.emplace_back(particle.GetPDG());QA_Decay_Length.emplace_back(v0decaylength);QA_DCA_V0_PV.emplace_back(dcav0toPV);
			QA_IfBadReconstructed.emplace_back(1);
		}
		if (IfHelix && (fabs(particle.GetPDG()) == OmegaPdg)) {

			// helix
			TVector3 MomentumOfParticle(particle.GetPx(), particle.GetPy(), particle.GetPz());
			TVector3 PositionOfParticle(particle.GetX(), particle.GetY(), particle.GetZ());
			TLorentzVector OmegaLorentz(MomentumOfParticle, particle.GetE());
			StPicoPhysicalHelix heliPositionOfParticle(MomentumOfParticle, PositionOfParticle, magnet*kilogauss, particle.GetQ());
			double pathlength = heliPositionOfParticle.pathLength(Vertex3D, false);
			TVector3 MomentumOfParticle_tb = heliPositionOfParticle.momentumAt(pathlength, magnet*kilogauss); 
			px.emplace_back(MomentumOfParticle_tb.X());
			py.emplace_back(MomentumOfParticle_tb.Y());
			pz.emplace_back(MomentumOfParticle_tb.Z());
			InvariantMass.emplace_back(OmegaLorentz.M());//cout<<"particle.GetAtProductionVertex() = "<<particle.GetAtProductionVertex()<<endl;
			float DL = 0. , eDL = 0.;particle.GetDecayLength(DL,eDL);
			QA_Decay_Length.emplace_back(DL);
			OmegaVec.push_back(particle);ParticleVec.push_back(particle);
			QA_Chi2.emplace_back(particle.GetChi2());
			hHM_Chi2->Fill(particle.GetMass(),particle.GetChi2());
			// cout<<"particle.GetPz()="<<particle.GetPz()<<", "<<"MomentumOfParticle_tb.Z()="<<MomentumOfParticle_tb.Z()<<endl; 
			// cout<<"kilogauss = "<<kilogauss<<endl;
			// cout<<"MomentumOfParticle_tb.Mag() = "<<MomentumOfParticle_tb.Mag()<<endl;
			// cout<<"MomentumOfParticle.Mag() = "<<MomentumOfParticle.Mag()<<endl;
		}
		else
		{
			px.emplace_back(particle.GetPx());
			py.emplace_back(particle.GetPy());
			pz.emplace_back(particle.GetPz());
			InvariantMass.emplace_back(particle.GetMass());
			QA_Chi2.emplace_back(particle.GetChi2());
			hHM_Chi2->Fill(particle.GetMass(),particle.GetChi2());

		}
		// cout<<"CrefMult:"<<CrefMult<<endl;
		// cout<<"PDG:"<<particle.GetPDG()<<endl; 

		hLN_M->Fill(particle.GetMass(),H_ProcessEventNum);

		// Filling QA
		QA_hasTOF.emplace_back(0);
		QA_dEdx.emplace_back(0.);QA_m2.emplace_back(0.);QA_nSigmaProton.emplace_back(0.);
		QA_nSigmaPion.emplace_back(0.);QA_nSigmaKaon.emplace_back(0.);QA_zTOF_proton.emplace_back(0.);QA_zTOF_pion.emplace_back(0.);QA_zTOF_kaon.emplace_back(0.);
		QA_IfConfuse.emplace_back(0);


		int upQ; if (particle.GetPDG() == LambdaPdg) upQ = 1; else if (particle.GetPDG() == -1*LambdaPdg) upQ = -1; else continue;
		int eLambda = -(upQ-1)/2; // 0 if Lambda, 1 if AntiLambda

		SetDaughterTrackPointers(iKFParticle);
		if (ProtonTrackIndex == -99999 || PionTrackIndex == -99999) continue; if(!ProtonTrack) continue; if(!PionTrack) continue;

		double dmass = -999; // just a placeholder
		TLorentzVector p4Pair, p4Proton; // just a placeholder
		StLambdaDecayPair TmpLambdaDecayPair(p4Pair, p4Proton, ProtonTrackIndex, PionTrackIndex, (eLambda==0), dmass);
		KFParticleLambdaDecayPair.push_back(TmpLambdaDecayPair);
	} // End loop over KFParticles

	// // HighLight Reconstructed Track
	ReCons_TrackID.resize(0);
	for (int iParticle = 0;iParticle < ParticleVec.size();iParticle++){
		KFParticle particle = ParticleVec[iParticle];
		for (int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++)
		{ 
			const int daughterId = particle.DaughterIds()[iDaughter]; 
			const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
			const int globalTrackId = daughter.DaughterIds()[0];
			ReCons_TrackID.push_back(globalTrackId);

			Int_t nTracks = mPicoDst->numberOfTracks();
			// for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
			// 	StPicoTrack *track = mPicoDst->track(iTrack);
			// 	if (track->id() == globalTrackId){
			// 		if (iTrack > track->id()){
			// 			cout<<"track location = "<<iTrack<<" , TrackId = "<<track->id()<<endl;
			// 		}
			// 		break;
			// 	}
			// }
			Int_t iTrackStart = globalTrackId - 1;
			if (globalTrackId >= nTracks) {iTrackStart = nTracks - 1;}
			for (Int_t iTrack = iTrackStart;iTrack >= 0;iTrack--){
				StPicoTrack *track = mPicoDst->track(iTrack);
				if (track->id() == globalTrackId){
					(mPicoDst->track(iTrack))->setNHitsFit(0);
					break;
				}
			}
		}
	}

	// Filling Track
	Int_t nTracks = mPicoDst->numberOfTracks();
	std::vector<int> Particle_tracks; Particle_tracks.resize(0);
	std::vector<int> track_index;
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
		StPicoTrack *track = mPicoDst->track(iTrack);
		hdEdx_pQ->Fill(1.0*track->charge()*track->gMom().Mag(),track->dEdx());
    	if (! track)            continue;
    	if (! track->charge())  continue;
    	if (  track->nHitsFit() < 15) continue;
		if (  track->nHitsDedx() < 15) continue;
		if (  track->nHitsFit()*1.0 / track->nHitsMax() < 0.52 || track->nHitsFit()*1.0 / track->nHitsMax() > 1.05) continue;
		if (  track->dEdxError() < 0.04 || track->dEdxError() > 0.12) continue; // same as kfp
		if (! track->isPrimary()) continue;
		track_index.push_back(iTrack);

		hdEdx_pQ_1cut->Fill(1.0*track->charge()*track->gMom().Mag(),track->dEdx());

		// track info
		float p = track->gMom().Mag();
		float pt = track->gMom().Perp();
		float phi = track->gMom().Phi();
		float eta = track->gMom().Eta();
		float dcatopv = track->gDCA(Vertex3D).Mag();
		float nSigmaKaon = track->nSigmaKaon();
		float nSigmaPion = track->nSigmaPion();
		float nSigmaProton = track->nSigmaProton();
		float pt_prim = track->pMom().Perp();
		float phi_prim = track->pMom().Phi();
		float eta_prim = track->pMom().Eta();

		// TOF Info
		bool hasTOF = false;
		int tofindex = track->bTofPidTraitsIndex();
		float m2 = -999.;
		float beta = -999.;
		if (tofindex >= 0) 
		{
			int tofflag = (mPicoDst->btofPidTraits(tofindex))->btofMatchFlag();
			float tof = (mPicoDst->btofPidTraits(tofindex))->btof();
			float BtofYLocal = (mPicoDst->btofPidTraits(tofindex))->btofYLocal();
			// hgbtofYlocal->Fill(BtofYLocal);
			if((tofflag >= 1) && (tof > 0) && (BtofYLocal > -1.8) && (BtofYLocal < 1.8)) hasTOF = true;
		}
		StPicoPhysicalHelix helix = track->helix(magnet);
		TVector3 pkaon = helix.momentum(magnet*kilogauss);
		if (hasTOF)
		{
			beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
			m2 = pkaon.Mag2()*(1.0 / beta / beta - 1.0);

			// some kaon QA
			//if (track->nSigmaKaon() >  6) hgptm2_largenSigmaKaon->Fill(track->gMom().Perp(), m2);
			//if (track->nSigmaKaon() < -6) hgptm2_smallnSigmaKaon->Fill(track->gMom().Perp(), m2);
			zTOF_proton = 1/beta - sqrt(ProtonPdgMass*ProtonPdgMass/pkaon.Mag2()+1);
			zTOF_pion   = 1/beta - sqrt(PionPdgMass*PionPdgMass/pkaon.Mag2()+1);
			zTOF_kaon   = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
		}

		// Fill tracks
		bool IfRecordThisTrack = false;
		// for (int i = 0; i < OmegaVec.size(); i++){ if (IsKaonOmegaDaughter(OmegaVec[i], track->id())) {
		// 	cout<<"This is event: "<<evtID<<endl;
		// 	cout<<"FUCK !!!"<<endl;
		// 	if(! track){cout<<"! track = YES"<<endl;}
		// 	else       {cout<<"! track = NO"<<endl;}
		// 	if(track->nHitsFit()<15){cout<<"nHitsFit()<15 = YES"<<endl;}
		// 	else             {cout<<"nHitsFit()<15 = NO"<<endl;}
		// 	KFParticle particle = OmegaVec[i];
		// 	for (int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++)
		// 	{ 
		// 		const int daughterId = particle.DaughterIds()[iDaughter]; 
		// 		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
		// 		const int globalTrackId = daughter.DaughterIds()[0];
		// 		cout<<"globalTrackId = "<<globalTrackId<<endl;
		// 	}
		// 	cout<<"Real TrackId = "<<track->id()<<endl;
		// 	cout<<"track->nHitsFit() = "<<track->nHitsFit()<<endl;
		// 	cout<<"End Fuck"<<endl;
		// }};

		// Test if Proton
		bool proton_cut = true;
		if (fabs(nSigmaProton) > 2) proton_cut = false;
		if (pt < pT_trig_lo || pt > pT_trig_hi) proton_cut = false; 
		if (fabs(eta_prim) > eta_trig_cut) proton_cut = false;
		if (!hasTOF && pt > proton_pT_TOFth) proton_cut = false;
		if (pt > proton_pT_TOFth && (m2 > proton_m2_hi || m2 < proton_m2_lo)) proton_cut = false;
		ProtonPID proton_pid(0., nSigmaProton, pt); // not using zTOF
		if (!proton_pid.IsProtonSimple(2., track->charge())) proton_cut = false; // only 0.2 < pt < 2.0!!!
		if (dcatopv > dcatoPV_hi) proton_cut = false;
		
		// Test if Pion
		bool pion_cut = true;
		if (fabs(nSigmaPion) > 2) pion_cut = false;
		if (pt < pT_trig_lo || pt > pT_trig_hi) pion_cut = false; // use p < 2
		if (fabs(eta_prim) > eta_trig_cut) pion_cut = false;
		PionPID pion_pid(0., nSigmaPion, pt); // not using zTOF
		if (!pion_pid.IsPionSimple(2., track->charge())) pion_cut = false; // only 0.2 < pt < 2.0!!!
		if (dcatopv > dcatoPV_hi) pion_cut = false;
		if (!hasTOF && pt > pion_pT_TOFth) pion_cut = false;
		if (pt > pion_pT_TOFth && (m2 > pion_m2_hi || m2 < pion_m2_lo)) pion_cut = false;

		// Test if Kaon
		bool kaon_cut = true;
		if (fabs(nSigmaKaon) > 2) kaon_cut = false;
		if (pt < pT_trig_lo || pt > 1.6) kaon_cut = false; // use p < 1.6
		if (fabs(eta_prim) > eta_trig_cut) kaon_cut = false;
		if (!hasTOF && pt > 0.4) kaon_cut = false;
		if (pt > 0.4 && (m2 > 0.34 || m2 < 0.15)) kaon_cut = false;
		double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
		KaonPID kaon_pid(zTOF, nSigmaKaon, pt); // not using zTOF
		if (!kaon_pid.IsKaonSimple(2., track->charge())) kaon_cut = false; // only 0.2 < pt < 2.0!!!
		if (dcatopv > 2.0) kaon_cut = false;
		// for (int i = 0; i < OmegaVec.size(); i++) if (IsKaonOmegaDaughter(OmegaVec[i], track->id())) kaon_cut = false;
		// for (int i = 0; i < OmegaVec.size(); i++) if (IsTrackParticleDaughter(OmegaVec[i], track->id())) kaon_cut = false;

		if (proton_cut + pion_cut + kaon_cut == 1) {IfRecordThisTrack = true;QA_IfConfuse.emplace_back(0);}
		if (proton_cut + pion_cut + kaon_cut > 1){IfRecordThisTrack = true;QA_IfConfuse.emplace_back(1);}

		if (IfRecordThisTrack == true) {
			hdEdx_pQ_2cut->Fill(1.0*track->charge()*track->gMom().Mag(),track->dEdx());

			px.emplace_back(track->gMom().X());
			py.emplace_back(track->gMom().Y());
			pz.emplace_back(track->gMom().Z());
			if      (proton_cut) {
				IfRecordThisTrack = true;
				if (track->charge() > 0) {PDG.emplace_back( 2212);}
				else                     {PDG.emplace_back(-2212);}
				InvariantMass.emplace_back(ProtonPdgMass);
			}
			else if (pion_cut) {
				IfRecordThisTrack = true;
				if (track->charge() > 0) {PDG.emplace_back( 211);}
				else                     {PDG.emplace_back(-211);}
				InvariantMass.emplace_back(PionPdgMass);
			}
			else if (kaon_cut) {
				IfRecordThisTrack = true;
				if (track->charge() > 0) {PDG.emplace_back( 321);}
				else                     {PDG.emplace_back(-321);}
				InvariantMass.emplace_back(KaonPdgMass);
			}

			// Filling QA
			if (hasTOF) {
				QA_hasTOF.emplace_back(1);
				QA_zTOF_proton.emplace_back(zTOF_proton);QA_zTOF_pion.emplace_back(zTOF_pion);QA_zTOF_kaon.emplace_back(zTOF_kaon);
				QA_m2.emplace_back(m2);
			}
			else{
				QA_hasTOF.emplace_back(0);
				QA_zTOF_proton.emplace_back(0.);QA_zTOF_pion.emplace_back(0.);QA_zTOF_kaon.emplace_back(0.);
				QA_m2.emplace_back(0.);
			}
			QA_nSigmaProton.emplace_back(nSigmaProton);
			QA_nSigmaPion.emplace_back(nSigmaPion);
			QA_nSigmaKaon.emplace_back(nSigmaKaon);
			QA_dEdx.emplace_back(track->dEdx());QA_DCA_V0_PV.emplace_back(dcatopv);
			QA_Decay_Length.emplace_back(-1.0);
			QA_Chi2.emplace_back(-1.0);
			QA_IfBadReconstructed.emplace_back(-1);
			QA_DCA_Daughters.emplace_back(-1.0);
		}

	}

// ======= KFParticle end ======= //

// ======= Lambda loop ======= //
	for(int j=0; j<KFParticleLambdaDecayPair.size(); j++) {
		int i = KFParticleLambdaDecayPair[j].get_idxProton();
		int k = KFParticleLambdaDecayPair[j].get_idxPion();
		if(k == i) continue;

		StPicoTrack* mTrackI = (StPicoTrack*)mPicoDst->track(i);

		int    mchgI = mTrackI->charge();
		int    mhitI = mTrackI->nHitsFit();
		double mMomI = mTrackI->gMom().Mag();
		double mp0xI = mTrackI->gMom().X();
		double mp0yI = mTrackI->gMom().Y();
		double mp0zI = mTrackI->gMom().Z();
		double mpt0I = mTrackI->gMom().Perp();
		double mphiI = mTrackI->gMom().Phi();
		double metaI = mTrackI->gMom().PseudoRapidity();
		double mdcaI = mTrackI->gDCA(Vertex3D).Mag();

		if(mphiI<0) mphiI += 2*M_PI;
		if(mphiI>=2*M_PI) mphiI -= 2*M_PI;

		StPicoTrack* mTrackK = (StPicoTrack*)mPicoDst->track(k);

		int    mchgK = mTrackK->charge();
		int    mhitK = mTrackK->nHitsFit();
		double mMomK = mTrackK->gMom().Mag();
		double mp0xK = mTrackK->gMom().X();
		double mp0yK = mTrackK->gMom().Y();
		double mp0zK = mTrackK->gMom().Z();
		double mpt0K = mTrackK->gMom().Perp();
		double mphiK = mTrackK->gMom().Phi();
		double metaK = mTrackK->gMom().PseudoRapidity();
		double mdcaK = mTrackK->gDCA(Vertex3D).Mag();

		if(mphiK<0) mphiK += 2*M_PI;
		if(mphiK>=2*M_PI) mphiK -= 2*M_PI;

		bool isPP = mchgI>0 && mchgK>0;
		bool isNN = mchgI<0 && mchgK<0;
		bool isPN = mchgI>0 && mchgK<0;
		bool isNP = mchgI<0 && mchgK>0;
		bool isSelfLambda = isPN;
		bool isAntiLambda = isNP;

		// remove the SS cases
		if(isPP || isNN) continue;

		//reconstruction of V0, the parent particle
		TVector3 xv0, op1, op2;
		double dca1to2 = closestDistance(mTrackI, mTrackK, magnet, Vertex3D, xv0, op1, op2);
		TVector3 pv0 = op1 + op2;
		TVector3 xv0toPV = xv0 - Vertex3D;
		double rdotp = xv0toPV.Dot(pv0);
		double dcav0toPV = rdotp*rdotp/pv0.Mag2();
		dcav0toPV = sqrt(xv0toPV.Mag2() - dcav0toPV);
		double v0decaylength = xv0toPV.Mag();
		double v0cosrdotp = rdotp/v0decaylength/pv0.Mag();

		TLorentzVector p4ProtonI, p4PionK;
		p4ProtonI.SetPxPyPzE(op1.X(), op1.Y(), op1.Z(), sqrt(op1.Mag2() + pmass*pmass));
		p4PionK.SetPxPyPzE(  op2.X(), op2.Y(), op2.Z(), sqrt(op2.Mag2() + pimass*pimass));

		// ProtonI & PionK
		TLorentzVector p4Pair = p4ProtonI + p4PionK;
		double massPair = p4Pair.M();
		double ptPair   = p4Pair.Pt();
		double etaPair  = p4Pair.Eta();
		double phiPair  = p4Pair.Phi();
	}
// ======= Lambda loop ends ======= //

	if (PDG.size()>0){
		PDGMult = PDG.size(); // This is multiplicity of Recorded Particles
		if (Omega_Omegab_Num != 0){
			cout<<"Found Omega"<<endl;
			// hadronTree->Fill();
		}
		hadronTree->Fill();
	}
	/////////////////////////////////////////////////////////
	return kStOK;

}

//------------------------------------------------------------
bool StKFParticleAnalysisMaker::isGoodObs(double Obs) {

	bool mGoodObs = true;
	mGoodObs = mGoodObs && !std::isnan(Obs);
	mGoodObs = mGoodObs && !std::isinf(Obs);
	//mGoodObs = mGoodObs && fabs(Obs)<2.1;

	return mGoodObs;
}

//------------------------------------------------------------
void StKFParticleAnalysisMaker::BookVertexPlots()
{
  KFParticleInterface = new StKFParticleInterface;
  bool storeMCHistograms = false;
  KFParticlePerformanceInterface = new StKFParticlePerformanceInterface(KFParticleInterface->GetTopoReconstructor(), storeMCHistograms);
}

//------------------------------------------------------------
void StKFParticleAnalysisMaker::SetupKFParticle(){
	int maxGBTrackIndex = -1; //find max global track index
	for (unsigned int iTrack = 0; iTrack < PicoDst->numberOfTracks(); iTrack++){
		StPicoTrack *track = PicoDst->track(iTrack);
		if ( !track ) continue;
		if ( track->id() > maxGBTrackIndex ) maxGBTrackIndex = track->id();
	}
	vector<KFMCTrack> mcTracks(0);
	vector<int> triggeredTracks;
	vector<int> mcIndices(maxGBTrackIndex+1);
	for (unsigned int iIndex = 0; iIndex < mcIndices.size(); iIndex++) mcIndices[iIndex] = -1;
	if (maxGBTrackIndex > 0)  KFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex+1);

	if ( !KFParticleInterface->ProcessEvent(PicoDst, triggeredTracks,H_ProcessEventNum) ) InterfaceCantProcessEvent = true; else InterfaceCantProcessEvent = false;

	trackMap.resize(maxGBTrackIndex+1, -1); //make a map from trackID to track index in global track array
	for(unsigned int iTrack = 0; iTrack < PicoDst->numberOfTracks(); iTrack++)
	{
		StPicoTrack *track = PicoDst->track(iTrack); if (!track) continue;
		int index = track->id();
		trackMap[index] = iTrack;
	}

	KFParticlePerformanceInterface->SetMCTracks(mcTracks);
	KFParticlePerformanceInterface->SetMCIndexes(mcIndices);    
	KFParticlePerformanceInterface->SetCentralityBin(-1);
	KFParticlePerformanceInterface->SetCentralityWeight(1.);
	KFParticlePerformanceInterface->SetPrintEffFrequency(100000);
	KFParticlePerformanceInterface->PerformanceAnalysis();
} // void SetupKFParticle

//------------------------------------------------------------
void StKFParticleAnalysisMaker::SetDaughterTrackPointers(int iKFParticle){ // Get Daughter Tracks to calculate decay variables manually
	const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
	int upQ; if (particle.GetPDG() == LambdaPdg) upQ = 1; else if (particle.GetPDG() == -1*LambdaPdg) upQ = -1; else upQ = 0;
	ProtonTrackIndex = -99999, PionTrackIndex = -99999, KaonTrackIndex = -99999;
	for(int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++){ 
		const int daughterId = particle.DaughterIds()[iDaughter]; 
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
		const int globalTrackId = daughter.DaughterIds()[0];
		int trackIndex = trackMap[globalTrackId];    
		if (daughter.GetPDG() == upQ*ProtonPdg) ProtonTrackIndex = trackIndex;
		else if (daughter.GetPDG() == upQ*PionPdg) PionTrackIndex = trackIndex;
		else if (daughter.GetPDG() == upQ*KaonPdg) KaonTrackIndex = trackIndex;
	}  // iDaughter
	ProtonTrack = PicoDst->track(ProtonTrackIndex); PionTrack = PicoDst->track(PionTrackIndex);	KaonTrack = PicoDst->track(KaonTrackIndex);	
} // void SetDaughterTrackPointers

bool StKFParticleAnalysisMaker::IsKaonOmegaDaughter(KFParticle particle, int kaonTrackId)
{
	if (fabs(particle.GetPDG()) != OmegaPdg) return false;
	for(int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++)
	{ 
		const int daughterId = particle.DaughterIds()[iDaughter]; 
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
		if (fabs(daughter.GetPDG()) != KaonPdg) continue;
		const int globalTrackId = daughter.DaughterIds()[0];
		
		if (globalTrackId == kaonTrackId) return true;
	}  // iDaughter
	return false;
}

bool StKFParticleAnalysisMaker::IsTrackParticleDaughter(KFParticle particle, int TrackId)
{
	for(int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++)
	{ 
		const int daughterId = particle.DaughterIds()[iDaughter]; 
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
		const int globalTrackId = daughter.DaughterIds()[0];
		
		if (globalTrackId == TrackId) return true;
	}  // iDaughter
	return false;
}

void StKFParticleAnalysisMaker::SetDaughterTrackHits(KFParticle particle)
{
	ReCons_TrackID.resize(0);
	for(int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++)
	{ 
		const int daughterId = particle.DaughterIds()[iDaughter]; 
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
		const int globalTrackId = daughter.DaughterIds()[0];
		ReCons_TrackID.push_back(globalTrackId);
	}
}

int StKFParticleAnalysisMaker::TrackID(StPicoTrack *track , TVector3 Vertex3D , double magnet , bool Track_has_tof , float m2 = -999. , float beta = -999.){

	float TrackID_pt = track->gMom().Perp();
	float TrackID_dcatopv = track->gDCA(Vertex3D).Mag();
	float TrackID_nSigmaKaon = track->nSigmaKaon();
	float TrackID_nSigmaPion = track->nSigmaPion();
	float TrackID_nSigmaProton = track->nSigmaProton();
	float TrackID_eta_prim = track->pMom().Eta();

	float dcatoPV_hi = 3.0; // Upper limit of DCA to PVs
	float pT_trig_lo = 0.2;
	float pT_trig_hi = 2.0;
	float eta_trig_cut = 1.0;
	bool hasTOF = false;

	// Test if Proton
	bool proton_cut = true;
	if (fabs(TrackID_nSigmaProton) > 3) proton_cut = false;
	if (TrackID_pt < pT_trig_lo || TrackID_pt > pT_trig_hi) proton_cut = false; 
	if (fabs(TrackID_eta_prim) > eta_trig_cut) proton_cut = false;
	if (!hasTOF && TrackID_pt > proton_pT_TOFth) proton_cut = false;
	if (TrackID_pt > proton_pT_TOFth && (m2 > proton_m2_hi || m2 < proton_m2_lo)) proton_cut = false;
	ProtonPID proton_pid(0., TrackID_nSigmaProton, TrackID_pt); // not using zTOF
	if (!proton_pid.IsProtonSimple(2., track->charge())) proton_cut = false; // only 0.2 < TrackID_pt < 2.0!!!
	if (TrackID_dcatopv > dcatoPV_hi) proton_cut = false;
	
	// Test if Pion
	bool pion_cut = true;
	if (fabs(TrackID_nSigmaPion) > 3) pion_cut = false;
	if (TrackID_pt < pT_trig_lo || TrackID_pt > pT_trig_hi) pion_cut = false; // use p < 2
	if (fabs(TrackID_eta_prim) > eta_trig_cut) pion_cut = false;
	PionPID pion_pid(0., TrackID_nSigmaPion, TrackID_pt); // not using zTOF
	if (!pion_pid.IsPionSimple(2., track->charge())) pion_cut = false; // only 0.2 < TrackID_pt < 2.0!!!
	if (TrackID_dcatopv > dcatoPV_hi) pion_cut = false;
	if (!hasTOF && TrackID_pt > pion_pT_TOFth) pion_cut = false;
	if (TrackID_pt > pion_pT_TOFth && (m2 > pion_m2_hi || m2 < pion_m2_lo)) pion_cut = false;

	// Test if Kaon
	bool kaon_cut = true;
	if (fabs(TrackID_nSigmaKaon) > 2.5) kaon_cut = false;
	if (TrackID_pt < pT_trig_lo || TrackID_pt > 1.6) kaon_cut = false; // use p < 1.6
	if (fabs(TrackID_eta_prim) > eta_trig_cut) kaon_cut = false;
	if (!hasTOF && TrackID_pt > 0.4) kaon_cut = false;
	if (TrackID_pt > 0.4 && (m2 > 0.34 || m2 < 0.15)) kaon_cut = false;
	StPicoPhysicalHelix TrackID_helix = track->helix(magnet);
	TVector3 pkaon = TrackID_helix.momentum(magnet*kilogauss);
	double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
	KaonPID kaon_pid(zTOF, TrackID_nSigmaKaon, TrackID_pt); // not using zTOF
	if (!kaon_pid.IsKaonSimple(2., track->charge())) kaon_cut = false; // only 0.2 < TrackID_pt < 2.0!!!
	if (TrackID_dcatopv > 2.0) kaon_cut = false;

	if (proton_cut + pion_cut + kaon_cut == 1){
		if (proton_cut){
			if (track->charge() > 0) {return  2212;}
			else                     {return -2212;}
		}
		if (pion_cut){
			if (track->charge() > 0) {return  211;}
			else                     {return -211;}
		}
		if (kaon_cut){
			if (track->charge() > 0) {return  321;}
			else                     {return -321;}
		}
	}
	else{
		return 0; // Failed to identify
	}
}

TVector3 StKFParticleAnalysisMaker::LocAfterTransfer(StPicoPhysicalHelix Track , Double_t Length){
	double dPhase = (Length*cos(Track.dipAngle()))/(1/Track.curvature());
	TVector3 Position(cos(Track.phase() + dPhase)*Track.curvature() + Track.xcenter(),
					  sin(Track.phase() + dPhase)*Track.curvature() + Track.ycenter(),
					  Length*sin(Track.dipAngle()) + (Track.origin()).z());
	return Position;
}

double StKFParticleAnalysisMaker::DistanceBetween(TVector3 LA , TVector3 LB){
	double rX = (LA.x() - LB.x());
	double rY = (LA.y() - LB.y());
	double rZ = (LA.z() - LB.z());
	double Dis = pow(rX*rX + rY*rY + rZ*rZ,0.5);
	return Dis;
}

bool StKFParticleAnalysisMaker::IfGoodDaughterDCA(StPicoDst* mPicoDst , int iKFParticle , double magnet , double Gen1_DCALim , double Gen2_DCALim){
	bool result = true;
	KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
	int iTrack,kTrack;
	std::vector<StPicoPhysicalHelix> cTrack;
	for (int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++){
		const int daughterId = particle.DaughterIds()[iDaughter]; 
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
		if (fabs(daughter.GetPDG()) == 211 || fabs(daughter.GetPDG()) == 321 || fabs(daughter.GetPDG()) == 2212) // meaning this daughter particle was NOT reconstructed
		{
			cout<<"daughter ID = "<<daughter.GetPDG()<<endl;
			const int globalTrackId = daughter.DaughterIds()[0];
			Int_t nTracks = mPicoDst->numberOfTracks();
			Int_t iTrackStart = globalTrackId - 1;
			if (globalTrackId >= nTracks) {iTrackStart = nTracks - 1;}
			for (Int_t jTrack = iTrackStart;jTrack >= 0;jTrack--){
				StPicoTrack *track = mPicoDst->track(jTrack);
				if (track->id() == globalTrackId){
					// int TrackPDG = TrackID(track , Vertex3D , magnet , false);
					// if (iDaughter == 0){
					// 	StPicoTrack* mTrackI = (StPicoTrack*)mPicoDst->track(jTrack);
					// 	cTrackI = mTrackI->helix(magnet);
					// }
					// if (iDaughter == 1){
					// 	StPicoTrack* mTrackK = (StPicoTrack*)mPicoDst->track(jTrack);
					// 	cTrackI = mTrackI->helix(magnet);
					// }
					cTrack.emplace_back(track->helix(magnet));
					break;
				}
			}
		}
		else
		{
			result = StKFParticleAnalysisMaker::IfGoodDaughterDCA(mPicoDst,daughterId,magnet,Gen2_DCALim,-1.0);
			if (result = false) {return result;}

			TVector3 MomentumOfParticle(daughter.GetPx(), daughter.GetPy(), daughter.GetPz());
			TVector3 PositionOfParticle(daughter.GetX(),  daughter.GetY(),  daughter.GetZ());
			StPicoPhysicalHelix heliPositionOfParticle(MomentumOfParticle, PositionOfParticle, magnet*kilogauss, particle.GetQ());
			cTrack.emplace_back(heliPositionOfParticle);
		}
	}

	pair<Double_t , Double_t>RV = cTrack[0].pathLengths(cTrack[1] , 0.05 , 0.05);
	TVector3 LTrackI = cTrack[0].at(RV.first);
	TVector3 LTrackK = cTrack[1].at(RV.second);
	double TrackDCA = StKFParticleAnalysisMaker::DistanceBetween(LTrackI , LTrackK);
	cout<<"TrackDCA = "<<TrackDCA<<endl;

	if (TrackDCA < Gen1_DCALim) {return result;}
	else {return false;}
}

// StPicoHelix StKFParticleAnalysisMaker::StPicoTrack2StPicoHelix(StPicoTrack* Track){
// 	StPicoHelix Result;
// 	Result.setParameters(Double_t c, Double_t dip, Double_t phase,const TVector3& o, Int_t h)
// }