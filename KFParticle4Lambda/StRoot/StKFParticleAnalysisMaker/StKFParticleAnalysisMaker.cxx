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
#include "StPicoEvent/StPicoETofPidTraits.h"

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
#include "TPCandTOF.h"

#include "StTrackHelix.h"
#include "StLambdaDecayPair.h"
#include "TriggerList.h"
#include "MyToolkit.h"

// #define DataName           "pAu_200_15"
// #define DataName           "AuAu_27_18"
// #define DataName           "dAu_200_16"
#define DataName           "dAu_200_21"
// #define DataName           "dAu_62_16"
// #define DataName           "dAu_39_16"
// #define DataName           "dAu_20_16"
// #define DataName           "pp_200_15"
// #define DataName           "OO_200_21"
#define pi                 TMath::Pi()
#define OmegaPdgMass	   1.67245
#define XiPdgMass	       1.3223
#define LambdaPdgMass      1.11568
#define ProtonPdgMass      0.938272
#define PionPdgMass        0.139570
#define KaonPdgMass		   0.493677
#define K0SPdgMass		   0.49794
#define PhiPdgMass		   1.01926
#define LambdaPdg          3122
#define XiPdg              3312
#define OmegaPdg           3334
#define XiPdg              3312
#define KaonPdg			   321
#define ProtonPdg          2212
#define K0SPdg			   310
#define PhiPdg			   333
#define Xi1530Pdg		   3324
#define PionPdg            211
#define LambdaPdgMassSigma 0.0014
#define XiPdgMassSigma     0.0018
#define OmegaPdgMassSigma  0.0027
#define K0SPdgMassSigma    0.0043
#define PhiPdgMassSigma    0.0031

#define IfQAMode           false // If Writing Hist of QA;
#define IfTree             true  // If Writing Tree;

// #define DEBUGGING


TPCandTOF TPCandTOF_Gen(DataName);

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
		fout = new TFile(mOutName.Data(),"RECREATE");
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

	hEventNum = new TH1D("Events_Total","Events_Total",1,0,2);

	const int APDGList[]         = {     3122     ,   -3122   ,   3334    ,  -3334    , 3312        ,  -3312      ,   310   ,   333   };
	const TString ANameList[]    = {  "Lambda"    , "Lambdab" ,   "Omega" , "Omegab"  , "Xi"        ,  "Xib"      ,  "K0S"  ,  "Phi"  };
	const int BPDGList[]         = {    321       ,   -321    ,    211    , -211      ,    2212     ,   -2212     };
	const TString BNameList[]    = {  "Kaon+"     , "Kaon-"   ,   "Pi+"   , "Pi-"     , "Proton"    , "Protonb"   };
	const float TBPDGListMass[]  = { KaonPdgMass  ,KaonPdgMass,PionPdgMass,PionPdgMass,ProtonPdgMass,ProtonPdgMass};
	const int TCPDGList[]      = {    321       ,    211    ,    2212   };
	const TString TCNameList[] = {  "Kaon"      ,  "Pion"   ,  "Proton" };
	for (int Itr = 0;Itr < PDG2NameSize3;Itr++){
		CPDGList[Itr] = TCPDGList[Itr];
		CNameList[Itr] = TCNameList[Itr];
	}
	for (int Itr = 0;Itr < PDG2NameSize3;Itr++){
		BPDGListMass[Itr]  = TBPDGListMass[Itr];
	}
	for (int Itr = 0;Itr < PDG2NameSize;Itr++){
		PDGList[Itr] = APDGList[Itr];NameList[Itr] = ANameList[Itr];
	}
	for (int Itr = PDG2NameSize;Itr < PDG2NameSize + PDG2NameSize2;Itr++){
		int Jtr = Itr - PDG2NameSize;
		PDGList[Itr] = BPDGList[Jtr];NameList[Itr] = BNameList[Jtr];
	}

	if (IfQAMode) {

		hNRefMult = new TH1F("RefMult" , "Reference Multiplicity" , 1000, 0.0, 1000.0 ) ;
		hNRefMultA= new TH1F("RefMultA", "Reference MultiplicityA", 1000, 0.0, 1000.0 ) ;
		hNRefMultB= new TH1F("RefMultB", "Reference MultiplicityB", 1000, 0.0, 1000.0 ) ;
		hVertexXY = new TH2F("VertexXY", "Vertex XY Position", 200, -10.0, 10.0 , 200, -10., 10 ) ;
		hVertexZ  = new TH1F("VertexZ" , "Event Vertex Z Position", 200, -100.0, 100.0 ) ;
		hVertex2D = new TH2F("Vertex2D", "VertexZ vs VPD Vz", 200, -100.0, 100.0 , 200, -100., 100 ) ;
		hDiffVz   = new TH1F("VertexZdiff" , "VertexZ-VPDVz diff", 100, -10.0, 10.0 ) ;
		hcent     = new TH1F("centrality","centrality"  ,nCent,0.,nCent);
		hcentw    = new TH1F("centralityw","centralityw",nCent,0.,nCent);
		hNch_per_VertexZ = new TH2F("hNch_per_VertexZ" , "Event Nch vs. Vertex Z Position", 100, -100.0, 100.0 , 200 , 0,200) ;
		hNch_per_VertexZ->GetXaxis()->SetTitle("PVz [cm]");
		hNch_per_VertexZ->GetYaxis()->SetTitle("Nch");

		hcentRefM = new TProfile("hcentRefM","hcentRefM",nCent,0.,nCent,0,1000);
		hcentRefW = new TProfile("hcentRefW","hcentRefW",nCent,0.,nCent,0,1000);
		
		hdEdx_pQ = new TH2F("hdEdx_p_NO_CUT","dE/dx vs. p*Q without cut",500,-5,5,500,0,50);
		hdEdx_pQ->GetXaxis()->SetTitle("p*Q [GeV]");
		hdEdx_pQ->GetYaxis()->SetTitle("dE/dx [keV/cm]");
		
		hdEdx_pQ_1cut = new TH2F("hdEdx_p_1_CUT","dE/dx vs. p*Q HITS cut",500,-5,5,500,0,50);
		hdEdx_pQ_1cut->GetXaxis()->SetTitle("p*Q [GeV]");
		hdEdx_pQ_1cut->GetYaxis()->SetTitle("dE/dx [keV/cm]");
		
		hdEdx_pQ_2cut = new TH2F("hdEdx_p_2_CUT","dE/dx vs. p*Q HITS & PID cut",500,-5,5,500,0,50);
		hdEdx_pQ_2cut->GetXaxis()->SetTitle("p*Q [GeV]");
		hdEdx_pQ_2cut->GetYaxis()->SetTitle("dE/dx [keV/cm]");
		
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

		H_Total_Pz = new TH1F("H_Total_Pz","Total P_z for all tracks",     400,-100,100);
		H_Total_Pz->GetXaxis()->SetTitle("P_z [GeV]");

		H_Total_Pxy = new TH2F("H_Total_Pxy","Total P_xy for all tracks",     100,-10,10,100,-10,10);
		H_Total_Pxy->GetXaxis()->SetTitle("P_x [GeV]");
		H_Total_Pxy->GetYaxis()->SetTitle("P_y [GeV]");

		H_Pt_m2 = new TH2F("H_Pt_m2","m2 vs. p_t",     400,0,10,500,-0.5,2);
		H_Pt_m2->GetXaxis()->SetTitle("p_t [GeV]");
		H_Pt_m2->GetYaxis()->SetTitle("m2 [Gev^2]");

		H_Pt_nSigmaKaon = new TH2F("H_Pt_nSigmaKaon","nSigmaKaon vs. p_t",     400,0,10,800,-10,10);
		H_Pt_nSigmaKaon->GetXaxis()->SetTitle("p_t [GeV]");
		H_Pt_nSigmaKaon->GetYaxis()->SetTitle("nSigmaKaon");

		H_Pt_nSigmaKaonTOF = new TH2F("H_Pt_nSigmaKaonTOF","nSigmaKaonTOF vs. p_t",     400,0,10,800,-10,10);
		H_Pt_nSigmaKaonTOF->GetXaxis()->SetTitle("p_t [GeV]");
		H_Pt_nSigmaKaonTOF->GetYaxis()->SetTitle("nSigmaKaonTOF");

		TriggerList Trigger_List_Data(DataName);
		std::vector<int> Trigger_List = Trigger_List_Data.GetTriggerList();
		int TriggerListLength = Trigger_List.size();
		for (int TriggerItr = 0;TriggerItr < TriggerListLength;TriggerItr++){
			TString ChargeName[3] = {"Positive","Negative","All"};
			for (int Itr = 0;Itr < 3;Itr++){
			
				TString HistName1;
				TString HistName2;
				HistName1 = "H_eta_nSigmaKaon_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "nSigmaKaon vs. eta, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				H_eta_nSigmaKaon[TriggerItr][Itr] = new TH2F(HistName1,HistName2,     200,-2,2,200,-8,8);
				H_eta_nSigmaKaon[TriggerItr][Itr]->GetXaxis()->SetTitle("eta");
				H_eta_nSigmaKaon[TriggerItr][Itr]->GetYaxis()->SetTitle("nSigmaKaon");

				HistName1 = "H_eta_nSigmaPion_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "nSigmaPion vs. eta, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				H_eta_nSigmaPion[TriggerItr][Itr] = new TH2F(HistName1,HistName2,     200,-2,2,200,-8,8);
				H_eta_nSigmaPion[TriggerItr][Itr]->GetXaxis()->SetTitle("eta");
				H_eta_nSigmaPion[TriggerItr][Itr]->GetYaxis()->SetTitle("nSigmaPion");
				
				HistName1 = "H_eta_nSigmaProton_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "nSigmaProton vs. eta, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				H_eta_nSigmaProton[TriggerItr][Itr] = new TH2F(HistName1,HistName2,     200,-2,2,200,-8,8);
				H_eta_nSigmaProton[TriggerItr][Itr]->GetXaxis()->SetTitle("eta");
				H_eta_nSigmaProton[TriggerItr][Itr]->GetYaxis()->SetTitle("nSigmaProton");
				
				HistName1 = "H_eta_m2_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "m2 vs. eta, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				H_eta_m2[TriggerItr][Itr] = new TH2F(HistName1,HistName2,     200,-2,2,200,-0.5,2);
				H_eta_m2[TriggerItr][Itr]->GetXaxis()->SetTitle("eta");
				H_eta_m2[TriggerItr][Itr]->GetYaxis()->SetTitle("m2 [GeV^2]");
				
				HistName1 = "H_eta_PVz_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "Primary Vertex Z vs. eta, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				H_eta_PVz[TriggerItr][Itr] = new TH2F(HistName1,HistName2,     200,-2,2,100,-100,100);
				H_eta_PVz[TriggerItr][Itr]->GetXaxis()->SetTitle("eta");
				H_eta_PVz[TriggerItr][Itr]->GetYaxis()->SetTitle("PVz [cm]");
				
				HistName1 = "H_eta_PVr_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "Primary Vertex R vs. eta, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				H_eta_PVr[TriggerItr][Itr] = new TH2F(HistName1,HistName2,     200,-2,2,100,0,5);
				H_eta_PVr[TriggerItr][Itr]->GetXaxis()->SetTitle("eta");
				H_eta_PVr[TriggerItr][Itr]->GetYaxis()->SetTitle("PVr [cm]");
				
				HistName1 = "H_eta_DVz_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "VertexZ-vpdVz vs. eta, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				// H_eta_DVz[TriggerItr][Itr] = new TH2F(HistName1,HistName2,     200,-2,2,100,-5,5);
				H_eta_DVz[TriggerItr][Itr] = new TH2F(HistName1,HistName2,     200,-2,2,1000,-50,50); // band
				H_eta_DVz[TriggerItr][Itr]->GetXaxis()->SetTitle("eta");
				H_eta_DVz[TriggerItr][Itr]->GetYaxis()->SetTitle("DVz [cm]");
				
				HistName1 = "H_eta_triggerBIN_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "eta distribution, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				H_eta_triggerBIN[TriggerItr][Itr] = new TH1F(HistName1,HistName2,     200,-2,2);
				H_eta_triggerBIN[TriggerItr][Itr]->GetXaxis()->SetTitle("eta");
				
				HistName1 = "H_eta_triggerBIN_hasTOF_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "hasTOF eta distribution, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				H_eta_triggerBIN_hasTOF[TriggerItr][Itr] = new TH1F(HistName1,HistName2,     200,-2,2);
				H_eta_triggerBIN_hasTOF[TriggerItr][Itr]->GetXaxis()->SetTitle("eta");
				
				HistName1 = "H_Nch_triggerBIN_";
				HistName2 = ChargeName[Itr];
				HistName2 += "Tracks ";
				HistName2 += "N_Charge distribution, Trigger: ";
				HistName1 += ChargeName[Itr];
				HistName1 += "_Trigger_";
				HistName1 += TriggerItr;
				HistName2 += Trigger_List[TriggerItr];
				H_Nch_triggerBIN[TriggerItr][Itr] = new TH1F(HistName1,HistName2,     400,0,400);
				H_Nch_triggerBIN[TriggerItr][Itr]->GetXaxis()->SetTitle("Nch");
			}
		}
		H_eta_trigger = new TH2F("H_eta_trigger","trigger vs. eta",     200,-2,2,TriggerListLength,-0.5,TriggerListLength-0.5);
		H_eta_trigger->GetXaxis()->SetTitle("eta");
		H_eta_trigger->GetYaxis()->SetTitle("trigger");

		H_m2_nSigmaKaon_Pt = new TH3F("H_m2_nSigmaKaon","m2 vs. nSigmaKaon vs. Pt",     400,-0.5,2 , 400,-10,10 , 18,0.2,2.0);
		H_m2_nSigmaKaon_Pt->GetXaxis()->SetTitle("m2 [GeV^2]");
		H_m2_nSigmaKaon_Pt->GetYaxis()->SetTitle("nSigmaKaon");
		H_m2_nSigmaKaon_Pt->GetZaxis()->SetTitle("Pt [GeV]");

		H_All_nSigmaKaon_y   = new TH2F("H_All_nSigmaKaon_y"  ,"nSigmaKaon vs. y for all tracks"  ,200,-2,2,200,0,6);
		H_All_nSigmaKaon_eta = new TH2F("H_All_nSigmaKaon_eta","nSigmaKaon vs. eta for all tracks",200,-2,2,200,0,6);
		H_All_nSigmaKaon_y->GetXaxis()->SetTitle("rapidity");
		H_All_nSigmaKaon_y->GetYaxis()->SetTitle("nSigmaKaon");
		H_All_nSigmaKaon_eta->GetXaxis()->SetTitle("eta");
		H_All_nSigmaKaon_eta->GetYaxis()->SetTitle("nSigmaKaon");

		for (int Itr = 0;Itr < PDG2NameSize;Itr++){

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

			HistName1 = "HM_WDaughter_";
			HistName2 = "The Mass of (Wrong Daughter) ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_WrongDaughter[Itr] = new TH1F(HistName1,HistName2,6000,0,3);
			H_WrongDaughter[Itr]->GetXaxis()->SetTitle("Mass [GeV]");

			HistName1 = "HM_CDaughter_";
			HistName2 = "The Mass of (Corect Daughter) ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_CrectDaughter[Itr] = new TH1F(HistName1,HistName2,6000,0,3);
			H_CrectDaughter[Itr]->GetXaxis()->SetTitle("Mass [GeV]");
			
			HistName1 = "H_Hyperon_Rap_";
			HistName2 = "The Rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_Hyperon_Rap[Itr] = new TH1F(HistName1,HistName2,120,-1.5,1.5);
			H_Hyperon_Rap[Itr]->GetXaxis()->SetTitle("y");
		}
		for (int Itr = PDG2NameSize;Itr < PDG2NameSize + PDG2NameSize2;Itr++){
			int Jtr = Itr - PDG2NameSize;

			TString HistName1 = "HY_";
			TString HistName2 = "The rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			// H_rapidity[Jtr] = new TH1F(HistName1,HistName2,30,-1.5,1.5); // KFP Setting
			H_rapidity[Jtr] = new TH1F(HistName1,HistName2,120,-1.5,1.5);
			H_rapidity[Jtr]->GetXaxis()->SetTitle("y");
			
			HistName1 = "HY_eTOF_";
			HistName2 = "The eTOF rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_rapidity_eTOF[Jtr] = new TH1F(HistName1,HistName2,120,-1.5,1.5);
			H_rapidity_eTOF[Jtr]->GetXaxis()->SetTitle("y");
			
			HistName1 = "HY_eTOF_Only";
			HistName2 = "The rapidity in eTOF contained event of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_rapidity_Only_eTOF[Jtr] = new TH1F(HistName1,HistName2,120,-1.5,1.5);
			H_rapidity_Only_eTOF[Jtr]->GetXaxis()->SetTitle("y");

			HistName1 = "HP_";
			HistName2 = "The momentum of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_P[Jtr] = new TH1F(HistName1,HistName2,100,0,10);
			H_P[Jtr]->GetXaxis()->SetTitle("p [GeV]");

			HistName1 = "HPt_";
			HistName2 = "The momentum_t of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_Pt[Jtr] = new TH1F(HistName1,HistName2,100,0,10);
			H_Pt[Jtr]->GetXaxis()->SetTitle("Pt [GeV]");

			HistName1 = "H_y_Pt_";
			HistName2 = "The y vs. momentum_t of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_Pt[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,0,5);
			H_y_Pt[Jtr]->GetXaxis()->SetTitle("y");
			H_y_Pt[Jtr]->GetYaxis()->SetTitle("Pt [GeV]");

			HistName1 = "H_y_P_";
			HistName2 = "The y vs. momentum of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_P[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,0,5);
			H_y_P[Jtr]->GetXaxis()->SetTitle("y");
			H_y_P[Jtr]->GetYaxis()->SetTitle("P [GeV]");

			HistName1 = "H_y_m2_";
			HistName2 = "The y vs. m2 of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_m2[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,0,1);
			H_y_m2[Jtr]->GetXaxis()->SetTitle("y");
			H_y_m2[Jtr]->GetYaxis()->SetTitle("m^2 [GeV^2]");

			HistName1 = "H_y_nSigmaKaon_";
			HistName2 = "The y vs. nSigmaKaon of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_nSigmaKaon[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-5,5);
			H_y_nSigmaKaon[Jtr]->GetXaxis()->SetTitle("y");
			H_y_nSigmaKaon[Jtr]->GetYaxis()->SetTitle("nSigmaKaon");

			HistName1 = "H_y_nSigmaPion_";
			HistName2 = "The y vs. nSigmaPion of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_nSigmaPion[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-5,5);
			H_y_nSigmaPion[Jtr]->GetXaxis()->SetTitle("y");
			H_y_nSigmaPion[Jtr]->GetYaxis()->SetTitle("nSigmaPion");

			HistName1 = "H_y_nSigmaElectron_";
			HistName2 = "The y vs. nSigmaElectron of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_nSigmaElectron[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-5,5);
			H_y_nSigmaElectron[Jtr]->GetXaxis()->SetTitle("y");
			H_y_nSigmaElectron[Jtr]->GetYaxis()->SetTitle("nSigmaElectron");

			HistName1 = "H_y_nHitsFit_";
			HistName2 = "The y vs. nHitsFit of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_nHitsFit[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,70,0,70);
			H_y_nHitsFit[Jtr]->GetXaxis()->SetTitle("y");
			H_y_nHitsFit[Jtr]->GetYaxis()->SetTitle("nHitsFit");

			HistName1 = "H_y_nHitsDedx_";
			HistName2 = "The y vs. nHitsDedx of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_nHitsDedx[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,70,0,70);
			H_y_nHitsDedx[Jtr]->GetXaxis()->SetTitle("y");
			H_y_nHitsDedx[Jtr]->GetYaxis()->SetTitle("nHitsDedx");

			HistName1 = "H_y_nHitsFit2nHitsMax_";
			HistName2 = "The y vs. nHitsFit/nHitsMax of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_nHitsFit2nHitsMax[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,40,0.4,1.4);
			H_y_nHitsFit2nHitsMax[Jtr]->GetXaxis()->SetTitle("y");
			H_y_nHitsFit2nHitsMax[Jtr]->GetYaxis()->SetTitle("nHitsFit/nHitsMax");

			HistName1 = "H_y_eta_";
			HistName2 = "The y vs. eta of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_eta[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-2,2);
			H_y_eta[Jtr]->GetXaxis()->SetTitle("y");
			H_y_eta[Jtr]->GetYaxis()->SetTitle("eta");
			
			HistName1 = "H_dEdx_p";
			HistName2 = "The dEdx vs. momentum of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_dEdx_p[Jtr] = new TH2F(HistName1,HistName2,2000,-10,10,1000,0,20);
			H_dEdx_p[Jtr]->GetXaxis()->SetTitle("P [GeV]");
			H_dEdx_p[Jtr]->GetYaxis()->SetTitle("dE/dx [keV/cm]");

			HistName1 = "hgbtofYlocal";
			HistName2 = "The btofYlocal vs. rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			hgbtofYlocal[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-5,5);
			hgbtofYlocal[Jtr]->GetXaxis()->SetTitle("y");
			hgbtofYlocal[Jtr]->GetYaxis()->SetTitle("btofYlocal");

			HistName1 = "H_y_Vz";
			HistName2 = "The Event Primary Vertex Z vs. rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_Vz[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,160,-80,80);
			H_y_Vz[Jtr]->GetXaxis()->SetTitle("y");
			H_y_Vz[Jtr]->GetYaxis()->SetTitle("PVz [cm]");
			
			HistName1 = "H_y_Pz";
			HistName2 = "The P_z vs. rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_Pz[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,200,-10,10);
			H_y_Pz[Jtr]->GetXaxis()->SetTitle("y");
			H_y_Pz[Jtr]->GetYaxis()->SetTitle("p_z [GeV]");
			
			HistName1 = "H_y_Nch";
			HistName2 = "The Nch vs. rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_Nch[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,200,0,200);
			H_y_Nch[Jtr]->GetXaxis()->SetTitle("y");
			H_y_Nch[Jtr]->GetYaxis()->SetTitle("Nch");
			
			HistName1 = "H_Pz_Nch";
			HistName2 = "The Nch vs. Pz of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_Pz_Nch[Jtr] = new TH2F(HistName1,HistName2,200,-10,10,200,0,200);
			H_Pz_Nch[Jtr]->GetXaxis()->SetTitle("Pz [GeV]");
			H_Pz_Nch[Jtr]->GetYaxis()->SetTitle("Nch");

			for (int Ktr=0;Ktr < PDG2NameSize3;Ktr++) {
				HistName1 = "H_Pt_nSigma";
				HistName2 = "The P_t vs. nSigma";
				HistName1 += CNameList[Ktr];HistName2 += CNameList[Ktr];
				HistName1 += "_";HistName1 += NameList[Itr];
				HistName2 += " of ";HistName2 += NameList[Itr];
				H_Pt_nSigma[Jtr][Ktr] = new TH2F(HistName1,HistName2,1000,-10,10,1000,0,8);
				HistName1 = "nSigma";HistName1 += CNameList[Ktr];
				H_Pt_nSigma[Jtr][Ktr]->GetXaxis()->SetTitle(HistName1);
				H_Pt_nSigma[Jtr][Ktr]->GetYaxis()->SetTitle("pT [GeV]");
			}

			HistName1 = "H_DCAtoPV";
			HistName2 = "The DCA(to PV) of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_DCAtoPV[Jtr] = new TH1F(HistName1,HistName2,500,0,5);
			H_DCAtoPV[Jtr]->GetXaxis()->SetTitle("DCA [cm]");

			HistName1 = "H_eta";
			HistName2 = "The eta of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_eta[Jtr] = new TH1F(HistName1,HistName2,200,-2,2);
			H_eta[Jtr]->GetXaxis()->SetTitle("eta");

			HistName1 = "H_nHitsFit";
			HistName2 = "The nHitsFit of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_nHitsFit_p[Jtr] = new TH2F(HistName1,HistName2,60,0,60,500,0,10);
			H_nHitsFit_p[Jtr]->GetXaxis()->SetTitle("nHitsFit");
			H_nHitsFit_p[Jtr]->GetYaxis()->SetTitle("p [GeV]");

			HistName1 = "H_nHitsFit_nHitsMax";
			HistName2 = "The nHitsFit/nHitsMax of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_nHitsFit_nHitsMax[Jtr] = new TH1F(HistName1,HistName2,100,0,1);
			H_nHitsFit_nHitsMax[Jtr]->GetXaxis()->SetTitle("nHitsFit/nHitsMax");

			HistName1 = "H_ndEdx";
			HistName2 = "The ndEdx of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_ndEdx[Jtr] = new TH1F(HistName1,HistName2,50,0,50);
			H_ndEdx[Jtr]->GetXaxis()->SetTitle("ndEdx");

			HistName1 = "H_nSigmaTOF_p_";
			HistName2 = "The nSigmaTOFKaon vs p of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_nSigmaTOF_p[Jtr] = new TH2F(HistName1,HistName2,250,-10,10,250,0,10);
			H_nSigmaTOF_p[Jtr]->GetXaxis()->SetTitle("nSigmaTOFKaon");
			H_nSigmaTOF_p[Jtr]->GetYaxis()->SetTitle("p [GeV]");
			
			HistName1 = "H_y_nSigmaTOFKaon_";
			HistName2 = "The nSigmaTOFKaon vs rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_y_nSigmaTOFKaon[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,250,-10,10);
			H_y_nSigmaTOFKaon[Jtr]->GetXaxis()->SetTitle("y");
			H_y_nSigmaTOFKaon[Jtr]->GetYaxis()->SetTitle("nSigmaTOFKaon");
			
			HistName1 = "H_m2_nSigmaTOFKaon_";
			HistName2 = "The nSigmaTOFKaon vs m^2 of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_m2_nSigmaTOFKaon[Jtr] = new TH2F(HistName1,HistName2,100,-0.5,2,100,-10,10);
			H_m2_nSigmaTOFKaon[Jtr]->GetXaxis()->SetTitle("m2");
			H_m2_nSigmaTOFKaon[Jtr]->GetYaxis()->SetTitle("nSigmaTOFKaon");
			
			HistName1 = "H_Pxy_";
			HistName2 = "The Px vs Py of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_Pxy[Jtr] = new TH2F(HistName1,HistName2,100,-10,10,100,-10,10);
			H_Pxy[Jtr]->GetXaxis()->SetTitle("Px [GeV]");
			H_Pxy[Jtr]->GetYaxis()->SetTitle("Py [GeV]");
			
			HistName1 = "H_Pz_";
			HistName2 = "The Pz of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_Pz[Jtr] = new TH1F(HistName1,HistName2,400,-10,10);
			H_Pz[Jtr]->GetXaxis()->SetTitle("Pz [GeV]");

			// KFP PID QA
			HistName1 = "H_KFP_Y_";
			HistName2 = "KFP: The rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_rapidity[Jtr] = new TH1F(HistName1,HistName2,120,-1.5,1.5);
			H_KFP_rapidity[Jtr]->GetXaxis()->SetTitle("y");

			HistName1 = "H_KFP_P_";
			HistName2 = "KFP: The momentum of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_P[Jtr] = new TH1F(HistName1,HistName2,100,0,10);
			H_KFP_P[Jtr]->GetXaxis()->SetTitle("p [GeV]");

			HistName1 = "H_KFP_Pt_";
			HistName2 = "KFP: The momentum_t of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_Pt[Jtr] = new TH1F(HistName1,HistName2,100,0,10);
			H_KFP_Pt[Jtr]->GetXaxis()->SetTitle("Pt [GeV]");

			HistName1 = "H_KFP_y_Pt_";
			HistName2 = "KFP: The y vs. momentum_t of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_Pt[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,0,5);
			H_KFP_y_Pt[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_Pt[Jtr]->GetYaxis()->SetTitle("Pt [GeV]");

			HistName1 = "H_KFP_y_m2_";
			HistName2 = "KFP: The y vs. m2 of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_m2[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,0,1);
			H_KFP_y_m2[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_m2[Jtr]->GetYaxis()->SetTitle("m^2 [GeV^2]");

			HistName1 = "H_KFP_y_nSigmaKaon_";
			HistName2 = "KFP: The y vs. nSigmaKaon of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_nSigmaKaon[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-5,5);
			H_KFP_y_nSigmaKaon[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_nSigmaKaon[Jtr]->GetYaxis()->SetTitle("nSigmaKaon");

			HistName1 = "H_KFP_y_nSigmaPion_";
			HistName2 = "KFP: The y vs. nSigmaPion of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_nSigmaPion[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-5,5);
			H_KFP_y_nSigmaPion[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_nSigmaPion[Jtr]->GetYaxis()->SetTitle("nSigmaPion");

			HistName1 = "H_KFP_y_nSigmaElectron_";
			HistName2 = "KFP: The y vs. nSigmaElectron of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_nSigmaElectron[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-5,5);
			H_KFP_y_nSigmaElectron[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_nSigmaElectron[Jtr]->GetYaxis()->SetTitle("nSigmaElectron");

			HistName1 = "H_KFP_y_nHitsFit_";
			HistName2 = "KFP: The y vs. nHitsFit of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_nHitsFit[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,70,0,70);
			H_KFP_y_nHitsFit[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_nHitsFit[Jtr]->GetYaxis()->SetTitle("nHitsFit");

			HistName1 = "H_KFP_y_nHitsDedx_";
			HistName2 = "KFP: The y vs. nHitsDedx of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_nHitsDedx[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,70,0,70);
			H_KFP_y_nHitsDedx[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_nHitsDedx[Jtr]->GetYaxis()->SetTitle("nHitsDedx");

			HistName1 = "H_KFP_y_nHitsFit2nHitsMax_";
			HistName2 = "KFP: The y vs. nHitsFit/nHitsMax of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_nHitsFit2nHitsMax[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,40,0.4,1.4);
			H_KFP_y_nHitsFit2nHitsMax[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_nHitsFit2nHitsMax[Jtr]->GetYaxis()->SetTitle("nHitsFit/nHitsMax");

			HistName1 = "H_KFP_y_eta_";
			HistName2 = "KFP: The y vs. eta of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_eta[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-2,2);
			H_KFP_y_eta[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_eta[Jtr]->GetYaxis()->SetTitle("eta");
			
			HistName1 = "H_KFP_dEdx_p";
			HistName2 = "KFP: The dEdx vs. momentum of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_dEdx_p[Jtr] = new TH2F(HistName1,HistName2,2000,-10,10,1000,0,20);
			H_KFP_dEdx_p[Jtr]->GetXaxis()->SetTitle("P [GeV]");
			H_KFP_dEdx_p[Jtr]->GetYaxis()->SetTitle("dE/dx [keV/cm]");

			HistName1 = "h_KFP_tofYlocal";
			HistName2 = "KFP: The btofYlocal vs. rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			h_KFP_tofYlocal[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,100,-5,5);
			h_KFP_tofYlocal[Jtr]->GetXaxis()->SetTitle("y");
			h_KFP_tofYlocal[Jtr]->GetYaxis()->SetTitle("btofYlocal");

			HistName1 = "H_KFP_y_Vz";
			HistName2 = "KFP: The Event Primary Vertex Z vs. rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_Vz[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,160,-80,80);
			H_KFP_y_Vz[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_Vz[Jtr]->GetYaxis()->SetTitle("PVz [cm]");
			
			HistName1 = "H_KFP_y_Pz";
			HistName2 = "KFP: The P_z vs. rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_Pz[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,200,-10,10);
			H_KFP_y_Pz[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_Pz[Jtr]->GetYaxis()->SetTitle("p_z [GeV]");
			
			HistName1 = "H_KFP_y_Nch";
			HistName2 = "KFP: The Nch vs. rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_Nch[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,200,0,200);
			H_KFP_y_Nch[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_Nch[Jtr]->GetYaxis()->SetTitle("Nch");
			
			HistName1 = "H_KFP_Pz_Nch";
			HistName2 = "KFP: The Nch vs. Pz of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_Pz_Nch[Jtr] = new TH2F(HistName1,HistName2,200,-10,10,200,0,200);
			H_KFP_Pz_Nch[Jtr]->GetXaxis()->SetTitle("Pz [GeV]");
			H_KFP_Pz_Nch[Jtr]->GetYaxis()->SetTitle("Nch");

			for (int Ktr=0;Ktr < PDG2NameSize3;Ktr++) {
				HistName1 = "H_KFP_Pt_nSigma";
				HistName2 = "KFP: The P_t vs. nSigma";
				HistName1 += CNameList[Ktr];HistName2 += CNameList[Ktr];
				HistName1 += "_";HistName1 += NameList[Itr];
				HistName2 += " of ";HistName2 += NameList[Itr];
				H_KFP_Pt_nSigma[Jtr][Ktr] = new TH2F(HistName1,HistName2,1000,-10,10,1000,0,8);
				HistName1 = "nSigma";HistName1 += CNameList[Ktr];
				H_KFP_Pt_nSigma[Jtr][Ktr]->GetXaxis()->SetTitle(HistName1);
				H_KFP_Pt_nSigma[Jtr][Ktr]->GetYaxis()->SetTitle("pT [GeV]");
			}

			HistName1 = "H_KFP_DCAtoPV";
			HistName2 = "KFP: The DCA(to PV) of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_DCAtoPV[Jtr] = new TH1F(HistName1,HistName2,500,0,5);
			H_KFP_DCAtoPV[Jtr]->GetXaxis()->SetTitle("DCA [cm]");

			HistName1 = "H_KFP_eta";
			HistName2 = "KFP: The eta of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_eta[Jtr] = new TH1F(HistName1,HistName2,200,-2,2);
			H_KFP_eta[Jtr]->GetXaxis()->SetTitle("eta");

			HistName1 = "H_KFP_nHitsFit";
			HistName2 = "KFP: The nHitsFit of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_nHitsFit_p[Jtr] = new TH2F(HistName1,HistName2,60,0,60,500,0,10);
			H_KFP_nHitsFit_p[Jtr]->GetXaxis()->SetTitle("nHitsFit");
			H_KFP_nHitsFit_p[Jtr]->GetYaxis()->SetTitle("p [GeV]");

			HistName1 = "H_KFP_nHitsFit_nHitsMax";
			HistName2 = "KFP: The nHitsFit/nHitsMax of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_nHitsFit_nHitsMax[Jtr] = new TH1F(HistName1,HistName2,100,0,1);
			H_KFP_nHitsFit_nHitsMax[Jtr]->GetXaxis()->SetTitle("nHitsFit/nHitsMax");

			HistName1 = "H_KFP_ndEdx";
			HistName2 = "KFP: The ndEdx of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_ndEdx[Jtr] = new TH1F(HistName1,HistName2,50,0,50);
			H_KFP_ndEdx[Jtr]->GetXaxis()->SetTitle("ndEdx");

			HistName1 = "H_KFP_nSigmaTOF_p_";
			HistName2 = "KFP: The nSigmaTOFKaon vs p of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_nSigmaTOF_p[Jtr] = new TH2F(HistName1,HistName2,250,-10,10,250,0,10);
			H_KFP_nSigmaTOF_p[Jtr]->GetXaxis()->SetTitle("nSigmaTOFKaon");
			H_KFP_nSigmaTOF_p[Jtr]->GetYaxis()->SetTitle("p [GeV]");
			
			HistName1 = "H_KFP_y_nSigmaTOFKaon_";
			HistName2 = "KFP: The nSigmaTOFKaon vs rapidity of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_y_nSigmaTOFKaon[Jtr] = new TH2F(HistName1,HistName2,100,-2,2,250,-10,10);
			H_KFP_y_nSigmaTOFKaon[Jtr]->GetXaxis()->SetTitle("y");
			H_KFP_y_nSigmaTOFKaon[Jtr]->GetYaxis()->SetTitle("nSigmaTOFKaon");
			
			HistName1 = "H_KFP_m2_nSigmaTOFKaon_";
			HistName2 = "KFP: The nSigmaTOFKaon vs m^2 of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_m2_nSigmaTOFKaon[Jtr] = new TH2F(HistName1,HistName2,100,-0.5,2,100,-10,10);
			H_KFP_m2_nSigmaTOFKaon[Jtr]->GetXaxis()->SetTitle("m2");
			H_KFP_m2_nSigmaTOFKaon[Jtr]->GetYaxis()->SetTitle("nSigmaTOFKaon");
			
			HistName1 = "H_KFP_Pxy_";
			HistName2 = "KFP: The Px vs Py of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_Pxy[Jtr] = new TH2F(HistName1,HistName2,100,-10,10,100,-10,10);
			H_KFP_Pxy[Jtr]->GetXaxis()->SetTitle("Px [GeV]");
			H_KFP_Pxy[Jtr]->GetYaxis()->SetTitle("Py [GeV]");
			
			HistName1 = "H_KFP_Pz_";
			HistName2 = "KFP: The Pz of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_Pz[Jtr] = new TH1F(HistName1,HistName2,400,-10,10);
			H_KFP_Pz[Jtr]->GetXaxis()->SetTitle("Pz [GeV]");
			
			HistName1 = "H_KFP_Pt_m2_";
			HistName2 = "KFP: The m2 vs pt of ";
			HistName1 += NameList[Itr];HistName2 += NameList[Itr];
			H_KFP_Pt_m2[Jtr] = new TH2F(HistName1,HistName2,100,0,5,100,-0.5,2);
			H_KFP_Pt_m2[Jtr]->GetXaxis()->SetTitle("Pt [GeV]");
			H_KFP_Pt_m2[Jtr]->GetYaxis()->SetTitle("m2 [GeV^2]");
			// KFP PID QA End
		}
	}

	if (IfTree){
		buffer_size = 5000000;
		hadronTree = new TTree("hadronTree", "Tree_STAR");
		// hadronTree->Branch("buffer_size"       ,&buffer_size         ,"buffer_size/I"                       );
		hadronTree->Branch("PDGMult"            ,&PDGMult             ,"PDGMult/I"                           );
		hadronTree->Branch("refMult"            ,&CrefMult            ,"refMult/I"                           );
		hadronTree->Branch("grefMult"           ,&CgrefMult           ,"grefMult/I"                          );
		hadronTree->Branch("EventID"            ,&evtID               ,"EventID/I"                           );
		hadronTree->Branch("RunID"              ,&runID               ,"RunID/I"                             );
		hadronTree->Branch("TriggerID"          ,&TriggerID           ,"TriggerID/I"                         );
		hadronTree->Branch("Nch"                ,&Nch                 ,"Nch/I"                               );
		hadronTree->Branch("PDG"                ,&PDG                 );
		hadronTree->Branch("mix_px"             ,&px                  );
		hadronTree->Branch("mix_py"             ,&py                  );
		hadronTree->Branch("mix_pz"             ,&pz                  );
		hadronTree->Branch("QA_eta"             ,&QA_eta                 );

		// Used for PID QA
		hadronTree->Branch("dEdx"               ,&QA_dEdx              );
		hadronTree->Branch("m2"                 ,&QA_m2                );
		hadronTree->Branch("dcatopv"            ,&QA_DCA_V0_PV         );
		hadronTree->Branch("nSigmaProton"       ,&QA_nSigmaProton      );
		hadronTree->Branch("nSigmaPion"         ,&QA_nSigmaPion        );
		hadronTree->Branch("nSigmaKaon"         ,&QA_nSigmaKaon        );
		
		// Used for Reconstruction QA
		hadronTree->Branch("InvariantMass"      ,&InvariantMass       );
		hadronTree->Branch("Decay_Length"       ,&QA_Decay_Length      );
		hadronTree->Branch("Chi2"               ,&QA_Chi2              );

		// Used for restore corralated information
		hadronTree->Branch("ParentA"       ,&ParentA      );
		hadronTree->Branch("ParentB"       ,&ParentB      );
		hadronTree->Branch("ParentC"       ,&ParentC      );
		hadronTree->Branch("ParentD"       ,&ParentD      );
		hadronTree->Branch("ParentE"       ,&ParentE      );

	}


	cout << "-----------------------------------------" << endl;
	cout << "------- histograms & tree claimed -------" << endl;
	cout << "-----------------------------------------" << endl;

	return;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::WriteHistograms() {

	hEventNum->Write();

	if (IfTree){
		hadronTree ->Write();
	}

	if (IfQAMode){

		folder_EventQA  = fout->mkdir("Event_QA");
		folder_PIDQA    = fout->mkdir("PID_QA");
		folder_ReconsQA = fout->mkdir("Reconstruction_QA");
		KFPPIDQA        = fout->mkdir("KFP_PID_QA");

		////////////////////////////////////// Event-QA //////////////////////////////////////


		folder_EventQA->cd();
		//-- Used for  test --
		// hNRefMult ->Write();  
		// hNRefMultA->Write();  
		// hNRefMultB->Write();  
		hVertexXY ->Write();  
		hVertexZ  ->Write();  
		hVertex2D ->Write();
		hDiffVz   ->Write();
		// hcent     ->Write();  
		// hcentw    ->Write();  

		// hcentRefM ->Write();
		// hcentRefW ->Write();
		H_Total_Pz->Write();
		H_Total_Pxy->Write();
		hNch_per_VertexZ->Write();

		fout->cd();


		////////////////////////////////////// PID-QA //////////////////////////////////////

		folder_PIDQA->cd();
		// hdEdx_pQ->Write();
		// hdEdx_pQ_1cut->Write();
		// hdEdx_pQ_2cut->Write();

		// hXY->Write();
		// hHXY->Write();
		// hHM_Chi2->Write();
		// hHM_ParentDCA->Write();

		// H_All_nSigmaKaon_y  ->Write();
		// H_All_nSigmaKaon_eta->Write();
		// H_Pt_m2->Write();
		// H_Pt_nSigmaKaon->Write();
		// H_Pt_nSigmaKaonTOF->Write();
		H_m2_nSigmaKaon_Pt->Write();
		TriggerList Trigger_List_Data(DataName);
		std::vector<int> Trigger_List = Trigger_List_Data.GetTriggerList();
		int TriggerListLength = Trigger_List.size();
		for (int TriggerItr = 0;TriggerItr < TriggerListLength;TriggerItr++){
			for (int Itr = 0;Itr < 3;Itr++){
				if (Itr != 2) continue; // Do not write specific charge bin
				if (H_eta_nSigmaPion[TriggerItr][Itr]->Integral() == 0){continue;}
				H_eta_nSigmaKaon         [TriggerItr][Itr]->Write();
				H_eta_nSigmaPion         [TriggerItr][Itr]->Write();
				H_eta_nSigmaProton       [TriggerItr][Itr]->Write();
				H_eta_m2                 [TriggerItr][Itr]->Write();
				H_eta_PVz                [TriggerItr][Itr]->Write();
				// H_eta_PVr                [TriggerItr][Itr]->Write();
				H_eta_DVz                [TriggerItr][Itr]->Write();
				H_eta_triggerBIN         [TriggerItr][Itr]->Write();
				H_eta_triggerBIN_hasTOF  [TriggerItr][Itr]->Write();
				H_Nch_triggerBIN         [TriggerItr][Itr]->Write();
			}
		}
		H_eta_trigger     ->Write();

		//////////////////////////////////// Used for test //////////////////////////////////////////////////////////////////////////////////////
		for (int Itr = PDG2NameSize;Itr < PDG2NameSize + PDG2NameSize2;Itr++){
			int Jtr = Itr - PDG2NameSize;

			// if (abs(PDGList[Itr]) != KaonPdg) {continue;}

			PID_Tracks[Jtr]  = folder_PIDQA->mkdir(NameList[Itr]);
			PID_Tracks[Jtr] -> cd();

			H_rapidity[Jtr] -> Write();
			H_rapidity_eTOF[Jtr] -> Write();

			H_P[Jtr] -> Write();

			H_Pt[Jtr] -> Write();
			H_dEdx_p[Jtr]->Write();
			for (int Ktr=0;Ktr < PDG2NameSize3;Ktr++) {
				H_Pt_nSigma[Jtr][Ktr]->Write();
			}
			H_DCAtoPV[Jtr]->Write();
			H_eta[Jtr]->Write();
			H_nHitsFit_p[Jtr]->Write();
			H_nHitsFit_nHitsMax[Jtr]->Write();
			H_ndEdx[Jtr]->Write();
			H_nSigmaTOF_p[Jtr]->Write();
			H_y_Pt[Jtr]->Write();
			H_y_m2[Jtr]->Write();
			H_y_nSigmaKaon[Jtr]->Write();
			H_y_nSigmaPion[Jtr]->Write();
			H_y_nSigmaElectron[Jtr]->Write();
			H_y_nHitsFit[Jtr]->Write();
			H_y_nHitsDedx[Jtr]->Write();
			H_y_nHitsFit2nHitsMax[Jtr]->Write();
			H_y_eta[Jtr]->Write();
			hgbtofYlocal[Jtr]->Write();
			H_y_Pz[Jtr]->Write();
			H_y_P[Jtr]->Write();
			H_Pxy[Jtr]->Write();
			H_Pz[Jtr]->Write();
			H_y_Nch[Jtr]->Write();
			H_Pz_Nch[Jtr]->Write();
			H_y_nSigmaTOFKaon[Jtr]->Write();
			H_y_Vz[Jtr]->Write();
			H_m2_nSigmaTOFKaon[Jtr]->Write();
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		fout->cd();

		// KFP PID
		KFPPIDQA->cd();
		for (int Itr = PDG2NameSize;Itr < PDG2NameSize + PDG2NameSize2;Itr++){
			int Jtr = Itr - PDG2NameSize;
			if (abs(PDGList[Itr]) != KaonPdg) {continue;}
			KFPPID[Jtr]  = KFPPIDQA->mkdir(NameList[Itr]);
			KFPPID[Jtr] -> cd();

			H_KFP_rapidity[Jtr] -> Write();
			H_KFP_P[Jtr] -> Write();
			H_KFP_Pt[Jtr] -> Write();
			H_KFP_dEdx_p[Jtr]->Write();
			for (int Ktr=0;Ktr < PDG2NameSize3;Ktr++) {
				H_KFP_Pt_nSigma[Jtr][Ktr]->Write();
			}
			H_KFP_DCAtoPV[Jtr]->Write();
			H_KFP_eta[Jtr]->Write();
			H_KFP_nHitsFit_p[Jtr]->Write();
			H_KFP_nHitsFit_nHitsMax[Jtr]->Write();
			H_KFP_ndEdx[Jtr]->Write();
			H_KFP_nSigmaTOF_p[Jtr]->Write();
			H_KFP_y_Pt[Jtr]->Write();
			H_KFP_y_m2[Jtr]->Write();
			H_KFP_y_nSigmaKaon[Jtr]->Write();
			H_KFP_y_nSigmaPion[Jtr]->Write();
			H_KFP_y_nSigmaElectron[Jtr]->Write();
			H_KFP_y_nHitsFit[Jtr]->Write();
			H_KFP_y_nHitsDedx[Jtr]->Write();
			H_KFP_y_nHitsFit2nHitsMax[Jtr]->Write();
			H_KFP_y_eta[Jtr]->Write();
			h_KFP_tofYlocal[Jtr]->Write();
			H_KFP_y_Pz[Jtr]->Write();
			H_KFP_Pxy[Jtr]->Write();
			H_KFP_Pz[Jtr]->Write();
			H_KFP_y_Nch[Jtr]->Write();
			H_KFP_Pz_Nch[Jtr]->Write();
			H_KFP_y_nSigmaTOFKaon[Jtr]->Write();
			H_KFP_y_Vz[Jtr]->Write();
			H_KFP_m2_nSigmaTOFKaon[Jtr]->Write();
			H_KFP_Pt_m2[Jtr]->Write();

		}
		// KFP PID end

		fout->cd();

		////////////////////////////////////// Reconstruction-QA //////////////////////////////////////
		folder_ReconsQA->cd();
		for (int i=0;i<PDG2NameSize;i++){
			KFPRecons[i]  = folder_ReconsQA->mkdir(NameList[i]);
			KFPRecons[i] -> cd();
			H_ALL_NO_CUT[i]->Write();
			// H_DaughterDCA[i]->Write();
			H_Hyperon_Rap[i]->Write();

		}
		fout->cd();
	}
	cout<<"T_T:"<<endl;
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
	int Recorded_Hyperon = 0;
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
	// if ((!mEvent->isTrigger(530003))
	//   &&(!mEvent->isTrigger(530002))
	//   &&(!mEvent->isTrigger(530806))
	//   &&(!mEvent->isTrigger(530101))
	//   &&(!mEvent->isTrigger(530102))
	//   &&(!mEvent->isTrigger(530201))
	//   &&(!mEvent->isTrigger(530202))
	//   &&(!mEvent->isTrigger(530213))
	//   &&(!mEvent->isTrigger(530851))// bbc
	//   &&(!mEvent->isTrigger(530852))// zdce
	//   &&(!mEvent->isTrigger(530853))// vpd-10
	//   &&(!mEvent->isTrigger(530854))// vpd
	//   &&(!mEvent->isTrigger(530855))// zdc
	//   &&(!mEvent->isTrigger(530856))// bbcnotac
	//   &&(!mEvent->isTrigger(530857))// zdc-notac
	//   &&(!mEvent->isTrigger(530861))// bbc
	//   &&(!mEvent->isTrigger(530866))// bbcnotac
	//   &&(!mEvent->isTrigger(2)) //17132063
	//   &&(!mEvent->isTrigger(3))
	//   &&(!mEvent->isTrigger(4)) //VPD-5
	//   &&(!mEvent->isTrigger(6)) //highMult-VPD-5
	//   &&(!mEvent->isTrigger(15))
	//   &&(!mEvent->isTrigger(16)) //BHT2-VPD-30
	//   &&(!mEvent->isTrigger(17))
	//   &&(!mEvent->isTrigger(54))
	//   &&(!mEvent->isTrigger(55))
	//   &&(!mEvent->isTrigger(56))
	//   &&(!mEvent->isTrigger(57))
	//   &&(!mEvent->isTrigger(58))
	//   &&(!mEvent->isTrigger(59))
	//   &&(!mEvent->isTrigger(61)) //vpd-30
	//   )return kStOK;
	TriggerList Trigger_List_Data(DataName);
	std::vector<int> Trigger_List = Trigger_List_Data.GetTriggerList();
	int TriggerListLength = Trigger_List.size();
	int TriggerID_in_TriggerList = 0;
	bool IfTriggerMatch = false;
	for (int i = 0;i<TriggerListLength;i++){
		if (mEvent->isTrigger(Trigger_List[i])) {
			IfTriggerMatch = true;
			TriggerID_in_TriggerList = i;
			TriggerID = Trigger_List[i];
			break;
		}
	}
	if (IfTriggerMatch == false) {return kStOK;}
	hEventNum -> Fill(1);

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

	// float DiffDVZCenter[2] = { -30.0 , 50.5 };// d+Au@200GeV RUN 16
	// if (!(6.5<fabs(vpdVz-VertexZ) && fabs(vpdVz-VertexZ)<9.5 && fabs(vpdVz)<20)) return kStOK; // band test
	// if (!(fabs(vpdVz-VertexZ)<3 && 20<fabs(vpdVz) && (DiffDVZCenter[0] <= vpdVz) && (vpdVz <= DiffDVZCenter[1]))) return kStOK; // center test

	if (IfQAMode){
		hVertex2D ->Fill(VertexZ,vpdVz);
		hDiffVz   ->Fill(VertexZ-vpdVz); 
	}

	const double DVz = VertexZ-vpdVz;

	// d+Au@200GeV RUN21 https://drupal.star.bnl.gov/STAR/system/files/pwg5.pdf
	if(fabs(VertexZ - 5) > 50) return kStOK; // AuAu27 80 ; dAu@39 25
	if(sqrt(pow(VertexX + 0.4,2.)+pow(VertexY,2.))>2.0) return kStOK; 
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
	// if( centrality<0||centrality>=(nCent-1)) return kStOK;
	int cent = centrality+1;  

	double mWght = refmultWght;
	double mult_corr = refmultCorr;
	bool IfHelix = true;

	if (IfQAMode){
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
	}

	// ======= KFParticle ======= //
	std::vector<StLambdaDecayPair> KFParticleLambdaDecayPair;
	std::vector<KFParticle> ParticleVec , OmegaVec , LambdaVec;
	ParticleVec.resize(0);OmegaVec.resize(0);LambdaVec.resize(0);

	float dcatoPV_hi = 3.0; // Upper limit of DCA to PVs
	float pT_trig_lo = 0.2;
	float pT_trig_hi = 2.0;
	float eta_trig_cut = 1.0;
  
	std::vector<vector<int> > Correlatted_ID_List_T;
	CrefMult = refMult;CgrefMult = grefMult;
	PDG            .resize(0);
	px             .resize(0);
	py             .resize(0);
	pz             .resize(0);
	QA_eta         .resize(0);
	QA_dEdx        .resize(0);
	QA_m2          .resize(0);
	QA_DCA_V0_PV   .resize(0);
	QA_nSigmaProton.resize(0);
	QA_nSigmaPion  .resize(0);
	QA_nSigmaKaon  .resize(0);
	InvariantMass  .resize(0);
	QA_Decay_Length.resize(0);
	QA_Chi2        .resize(0);

	Recorded_KFP_ID.resize(0);
	ParentA.resize(0);
	ParentB.resize(0);
	ParentC.resize(0);
	ParentD.resize(0);
	ParentE.resize(0);
	Int_t nTracks = mPicoDst->numberOfTracks();
	// Calculating Nch
	int NumCharge = 0;
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
		StPicoTrack *track = mPicoDst->track(iTrack);
    	if (! track)            continue;
    	if (! track->charge())  continue;
    	if (  track->nHitsFit() < 15) continue;
		// if (  track->nHitsDedx() < 15) continue;
		// if (  track->nHitsFit()*1.0 / track->nHitsMax() < 0.52 || track->nHitsFit()*1.0 / track->nHitsMax() > 1.05) continue;
		if (  track->nHitsFit()*1.0 / track->nHitsMax() < 0.52) continue;
		// if (  track->dEdxError() < 0.04 || track->dEdxError() > 0.12) continue; // same as kfp
		// if (! track->isPrimary()) continue;
		if (track->gMom().Perp() < 0.2 || track->gMom().Perp() > 2.0) continue;
		if (track->gDCA(Vertex3D).Mag() > 3) continue;
		if (fabs(track->gMom().Eta()) > 1.5) continue;
		NumCharge++;
	}
	Nch = NumCharge;


	std::vector<int> DaughterParticle,MatherPartiecle;DaughterParticle.resize(0);MatherPartiecle.resize(0);
	if (!(DataName == "pp_200_15")){
		SetupKFParticle();
		// cout<<"2b"<<endl;
		// if (InterfaceCantProcessEvent) return;
		// cout<<"InterfaceCantProcessEvent OK"<<endl;


		// // Check If one Particle is used to reconstructed to two more particles
		// ReCons_TrackID.resize(0);std::vector<int> DaughterParticle,MatherPartiecle,MultyReconMather;DaughterParticle.resize(0);MatherPartiecle.resize(0);MultyReconMather.resize(0);
		// // cout<<"KFParticlePerformanceInterface->GetNReconstructedParticles() = "<<KFParticlePerformanceInterface->GetNReconstructedParticles()<<endl;
		// for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++){
		// 	KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
		// 	if (particle.GetPDG() == -1){continue;}
		// 	if ( (fabs(particle.GetPDG()) == OmegaPdg) || (fabs(particle.GetPDG()) == XiPdg)  || (fabs(particle.GetPDG()) == LambdaPdg) ) {
		// 		// cout<<"###############################################"<<endl;
		// 		// cout<<"iKFParticle = "<<iKFParticle<<endl;
		// 		// cout<<"particle.GetPDG() = "<<particle.GetPDG()<<endl;
		// 		for (int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++){
		// 			const int daughterId = particle.DaughterIds()[iDaughter];
		// 			const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
		// 			// cout<<"daughterId = "<<daughterId<<endl;
		// 			// cout<<"daughter.GetPDG() = "<<daughter.GetPDG()<<endl;
		// 			if (daughter.GetPDG() != -1) {
		// 				DaughterParticle.push_back(daughterId);
		// 				MatherPartiecle.push_back(iKFParticle);
		// 			}
		// 			else {
		// 				for (int jDaughter=0; jDaughter < daughter.NDaughters(); jDaughter++){
		// 					const int GdaughterId = daughter.DaughterIds()[jDaughter];
		// 					DaughterParticle.push_back(GdaughterId);
		// 					MatherPartiecle.push_back(iKFParticle);
		// 					// cout<<"GrandDaughterId = "<<GdaughterId<<endl;
		// 					// cout<<"GrandDaughter.GetPDG() = "<<(KFParticleInterface->GetParticles()[GdaughterId]).GetPDG()<<endl;
		// 				}
		// 			}
		// 		}
		// 		// cout<<"MatherPartiecle = [";
		// 		// for (int Itr = 0;Itr < DaughterParticle.size();Itr++) {
		// 		// 	cout<<" "<<MatherPartiecle[Itr]<<" ";
		// 		// 	if (Itr == DaughterParticle.size()-1) {cout<<"]"<<endl;}
		// 		// 	else{cout<<",";}
		// 		// }
		// 		// cout<<"DaughterParticle = [";
		// 		// for (int Itr = 0;Itr < DaughterParticle.size();Itr++) {
		// 		// 	cout<<" "<<DaughterParticle[Itr]<<" ";
		// 		// 	if (Itr == DaughterParticle.size()-1) {cout<<"]"<<endl;}
		// 		// 	else{cout<<",";}
		// 		// }
		// 	}
		// 	// DaughterParticle.resize(0);MatherPartiecle.resize(0);
		// }
		// for (int Itr = 0;Itr < DaughterParticle.size();Itr++){
		// 	for (int Jtr = Itr + 1;Jtr < DaughterParticle.size();Jtr++){
		// 		if (DaughterParticle[Itr] == DaughterParticle[Jtr]) {
		// 			MultyReconMather.push_back(MatherPartiecle[Itr]);
		// 			MultyReconMather.push_back(MatherPartiecle[Jtr]);
		// 		}
		// 	}
		// }

		Omega_Omegab_Num = 0;
		for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++){ 
			KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];

			bool IfWellConstrcuted = true;

			#ifdef DEBUGGING
			std::cout << "Parsing refMult : " << refMult <<std::endl;
			std::cout << "Parsed CrefMult : " << CrefMult <<std::endl;
			#endif
			
			if (IfQAMode){
				// Fill track from KFP
				if ((abs(particle.GetPDG()) == 2212) || (abs(particle.GetPDG()) == 211) || (abs(particle.GetPDG()) == 321)) { // Proton , Pion or Kaon
					for (int Itr = PDG2NameSize;Itr < PDG2NameSize + PDG2NameSize2;Itr++){
						int Jtr = Itr - PDG2NameSize;
						if (PDGList[Itr] != particle.GetPDG()) {continue;}

						// Find in PicoDST
						int iTrack = 0;
						const int globalTrackId = particle.DaughterIds()[0];
						Int_t iTrackStart = globalTrackId - 1;
						if (globalTrackId >= nTracks) {iTrackStart = nTracks - 1;}
						for (Int_t jTrack = iTrackStart;jTrack >= 0;jTrack--){
							StPicoTrack *track = mPicoDst->track(jTrack);
							if (track->id() == globalTrackId){
								iTrack = jTrack;
								break;
							}
						}

						StPicoTrack *track = mPicoDst->track(iTrack);
						float p = track->gMom().Mag();
						float pt = track->gMom().Perp();
						float phi = track->gMom().Phi();
						float eta = track->gMom().Eta();
						double track_px = track->gMom().X();
						double track_py = track->gMom().Y();
						double track_pz = track->gMom().Z();
						if (pt < 0.4 || pt > 1.5) continue;
						if (p < 0.4  || p > 2   ) continue;
						H_KFP_Pt[Jtr] -> Fill(pt);
						H_KFP_P[Jtr] -> Fill(p);
						float tEnergy = pow(pow(track->gMom().Mag(),2) + pow(StKFParticleAnalysisMaker::massList(particle.GetPDG()),2),0.5);
						float rap = 0.5*log((tEnergy+track_pz)/(tEnergy-track_pz));
						H_KFP_rapidity[Jtr]->Fill(rap);
						H_KFP_y_Pt[Jtr]->Fill(rap,pt);
						H_KFP_y_nSigmaKaon[Jtr]->Fill(rap,track->nSigmaKaon());
						H_KFP_y_nSigmaPion[Jtr]->Fill(rap,track->nSigmaPion());
						H_KFP_y_nSigmaElectron[Jtr]->Fill(rap,track->nSigmaElectron());
						H_KFP_y_nHitsFit[Jtr]->Fill(rap,track->nHitsFit());
						H_KFP_y_nHitsDedx[Jtr]->Fill(rap,track->nHitsDedx());
						H_KFP_y_nHitsFit2nHitsMax[Jtr]->Fill(rap,track->nHitsFit()*1.0 / track->nHitsMax());
						H_KFP_y_eta[Jtr]->Fill(rap,eta);
						H_KFP_y_Pz[Jtr]->Fill(rap,track_pz);
						H_KFP_Pxy[Jtr]->Fill(track_px,track_py);
						H_KFP_Pz[Jtr]->Fill(track_pz);
						H_KFP_y_Vz[Jtr]->Fill(rap,VertexZ);
						H_KFP_y_Nch[Jtr]->Fill(rap,NumCharge);
						H_KFP_Pz_Nch[Jtr]->Fill(track_pz,NumCharge);
						for (int Ntr=0;Ntr<PDG2NameSize3;Ntr++){
							// cout<<"CNameList["<<Ktr<<"] = "<<CNameList[Ktr]<<endl;
							if ( CNameList[Ntr] == "Kaon"){H_KFP_Pt_nSigma[Jtr][Ntr]->Fill(track->nSigmaKaon(),track->gMom().Perp());}
							if ( CNameList[Ntr] == "Pion"){H_KFP_Pt_nSigma[Jtr][Ntr]->Fill(track->nSigmaPion(),track->gMom().Perp());}
							if ( CNameList[Ntr] == "Proton"){H_KFP_Pt_nSigma[Jtr][Ntr]->Fill(track->nSigmaProton(),track->gMom().Perp());}
						}

						H_KFP_dEdx_p[Jtr]->Fill(1.0*track->charge()*track->gMom().Mag(),track->dEdx());
						H_KFP_DCAtoPV[Jtr]->Fill(track->gDCA(Vertex3D).Mag());
						H_KFP_eta[Jtr]->Fill(eta);
						H_KFP_nHitsFit_p[Jtr]->Fill(track->nHitsFit(),track->gMom().Mag());
						H_KFP_nHitsFit_nHitsMax[Jtr]->Fill((track->nHitsFit()*1.0)/(track->nHitsMax()*1.0));
						H_KFP_ndEdx[Jtr]->Fill((track->nHitsDedx()));
						bool hasTOF = false;
						int tofindex = track->bTofPidTraitsIndex();
						if (tofindex >= 0) 
						{
							int tofflag = (mPicoDst->btofPidTraits(tofindex))->btofMatchFlag();
							float tof = (mPicoDst->btofPidTraits(tofindex))->btof();
							if((tofflag >= 1) && (tof > 0)) hasTOF = true;
						}
						if (hasTOF) {
							StPicoPhysicalHelix helix = track->helix(magnet);
							TVector3 pkaon = helix.momentum(magnet*kilogauss);
							double beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
							double m2 = pkaon.Mag2()*(1.0 / beta / beta - 1.0);
							if(abs(PDGList[Itr])==321){
								H_KFP_nSigmaTOF_p[Jtr]->Fill((mPicoDst->btofPidTraits(track->bTofPidTraitsIndex()))->nSigmaKaon(),track->gMom().Mag());
								h_KFP_tofYlocal[Jtr]->Fill(rap,(mPicoDst->btofPidTraits(track->bTofPidTraitsIndex()))->btofYLocal());
								H_KFP_y_nSigmaTOFKaon[Jtr]->Fill(rap,(mPicoDst->btofPidTraits(track->bTofPidTraitsIndex()))->nSigmaKaon());
								H_KFP_m2_nSigmaTOFKaon[Jtr]->Fill(m2,(mPicoDst->btofPidTraits(track->bTofPidTraitsIndex()))->nSigmaKaon());
								H_KFP_Pt_m2[Jtr]->Fill(pt,m2);
								H_KFP_y_m2[Jtr]->Fill(rap,m2);
							}
						}

					}

				}
			}



			if (
				(fabs(particle.GetPDG()) != OmegaPdg ) && 
				(fabs(particle.GetPDG()) != XiPdg    ) && 
				(fabs(particle.GetPDG()) != LambdaPdg) &&
				(fabs(particle.GetPDG()) != K0SPdg   ) && 
				(fabs(particle.GetPDG()) != PhiPdg)
			) {continue;}
			Recorded_Hyperon ++;

			// Check if wrong daughters
			// bool IfCorrectDaughter = true;
			// for (int Itr = 0;Itr < MultyReconMather.size();Itr++) {
			// 	if ( iKFParticle == MultyReconMather[Itr]) {
			// 		for (int Jtr = 0;Jtr < PDG2NameSize;Jtr++){
			// 			if (particle.GetPDG() == PDGList[Jtr]){
			// 				H_WrongDaughter[Jtr]->Fill(particle.GetMass());
			// 				IfCorrectDaughter = false;
			// 				break;
			// 			}
			// 		}
			// 		break;
			// 	}
			// 	if (Itr == MultyReconMather.size() - 1) {
			// 		for (int Jtr = 0;Jtr < PDG2NameSize;Jtr++){
			// 			if (particle.GetPDG() == PDGList[Jtr]){
			// 				H_CrectDaughter[Jtr]->Fill(particle.GetMass());
			// 				break;
			// 			}
			// 		}
			// 	}
			// }
			// if (IfCorrectDaughter == false) {continue;}

			//SCHEME 1: reconstruction of V0, the parent particle
			int iTrack,kTrack;
			// cout<<"particle.NDaughters() = "<<particle.NDaughters()<<endl;
			for (int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++){
				const int daughterId = particle.DaughterIds()[iDaughter];
				// cout<<"daughterId = "<<daughterId<<endl;
				const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
				if (particle.GetPDG() == 3334 || particle.GetPDG() == -3334) {
					// cout<<"daughter ID = "<<daughter.GetPDG()<<endl;
					if (daughter.GetPDG() == -1) {
						for (int jDaughter=0; jDaughter < daughter.NDaughters(); jDaughter++){
							// cout<<"Here is good 7"<<endl;
							int DdaughterId = daughter.DaughterIds()[jDaughter];
							KFParticle Ddaughter = KFParticleInterface->GetParticles()[DdaughterId];
							// cout<<"Grand daughter ID = "<<Ddaughter.GetPDG()<<endl;
						}
					}
				}
				// cout<<"Here is good 4"<<endl;
				const int globalTrackId = daughter.DaughterIds()[0];
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
			// cout<<"Here is good 5"<<endl;
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
			if (IfQAMode){
				for (int Itr = 0;Itr < PDG2NameSize;Itr++){
					if (particle.GetPDG() == PDGList[Itr]){
						float PPx = particle.GetPx();
						float PPy = particle.GetPy();
						float PPz = particle.GetPz();
						if (StKFParticleAnalysisMaker::IfGoodDaughterDCA(mPicoDst , iKFParticle , magnet , 0.6 , 0.6)){
							H_ALL_NO_CUT[Itr]->Fill(particle.GetMass());
							H_DaughterDCA[Itr]->Fill(particle.GetMass());
							// QA_IfBadReconstructed.emplace_back(0);
						}
						else
						{
							H_ALL_NO_CUT[Itr]->Fill(particle.GetMass());
							// IfWellConstrcuted = false;
							// QA_IfBadReconstructed.emplace_back(1);
						}
						float PMass = 0.0;
						if      ((abs(PDGList[Itr]) == K0SPdg   ) && (fabs(particle.GetMass() - K0SPdgMass   ) < 3*K0SPdgMassSigma   )){ PMass = K0SPdgMass;}
						else if ((abs(PDGList[Itr]) == LambdaPdg) && (fabs(particle.GetMass() - LambdaPdgMass) < 3*LambdaPdgMassSigma)){ PMass = LambdaPdgMass;}
						else if ((abs(PDGList[Itr]) == XiPdg    ) && (fabs(particle.GetMass() - XiPdgMass    ) < 3*XiPdgMassSigma    )){ PMass = XiPdgMass;}
						else if ((abs(PDGList[Itr]) == OmegaPdg ) && (fabs(particle.GetMass() - OmegaPdgMass ) < 3*OmegaPdgMassSigma )){ PMass = OmegaPdgMass;}
						else if ((abs(PDGList[Itr]) == PhiPdg   ) && (fabs(particle.GetMass() - PhiPdgMass   ) < 3*PhiPdgMassSigma   )){ PMass = PhiPdgMass;}
						else if ((abs(PDGList[Itr]) == LambdaPdg) || (abs(PDGList[Itr]) == XiPdg) || (abs(PDGList[Itr]) == OmegaPdg) || (abs(PDGList[Itr]) == K0SPdg) || (abs(PDGList[Itr]) == PhiPdg)) {break;}
						else {
							continue;
						}
						// 
						bool IfPass = true;
						if ((particle.GetPDG() == 310)) {
							for (int iDaughter = 1;iDaughter<particle.NDaughters();iDaughter++){
								const int daughterId = particle.DaughterIds()[iDaughter];
								// cout<<"daughterId = "<<daughterId<<endl;
								const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
								if (abs(daughter.GetPDG()) != 211) continue;
								int iTrack = 0;
								const int globalTrackId = (KFParticleInterface->GetParticles()[daughterId]).DaughterIds()[0];
								Int_t iTrackStart = globalTrackId - 1;
								if (globalTrackId >= nTracks) {iTrackStart = nTracks - 1;}
								for (Int_t jTrack = iTrackStart;jTrack >= 0;jTrack--){
									StPicoTrack *track = mPicoDst->track(jTrack);
									if (track->id() == globalTrackId){
										iTrack = jTrack;
										break;
									}
								}
								StPicoTrack *track = mPicoDst->track(iTrack);
								std::vector<int> TestPDG;TestPDG.push_back((KFParticleInterface->GetParticles()[daughterId]).GetPDG());
								std::vector<bool> PDGBool = StKFParticleAnalysisMaker::TrackPID(TestPDG , track , Vertex3D);
								if(PDGBool[0] == false) IfPass = false;
							}
						}
						if (IfPass == false) continue;
						//
						float PEnergy = pow(PPx*PPx + PPy*PPy + PPz*PPz + PMass*PMass , 0.5);
						H_Hyperon_Rap[Itr]->Fill(0.5*log((PEnergy+PPz)/(PEnergy-PPz)));
						
						// hHM_ParentDCA->Fill(particle.GetMass(),TrackDCA);
						// QA_DCA_Daughters.emplace_back(-1.0);
						// cout<<"DCA to Parent = "<<DistanceBetween(LTrackI , LTrackK)<<endl;
						break;
					}
				}
			}
			
			// StPicoTrack* mTrackI = (StPicoTrack*)mPicoDst->track(iTrack);
			// StPicoTrack* mTrackK = (StPicoTrack*)mPicoDst->track(kTrack);
			// TVector3 xv0, op1, op2;
			// cout<<"T2";
			// double dca1to2 = closestDistance(mTrackI, mTrackK, magnet, Vertex3D, xv0, op1, op2);
			// cout<<"T3";
			// TVector3 pv0 = op1 + op2;
			// TVector3 xv0toPV = xv0 - Vertex3D;
			// double rdotp = xv0toPV.Dot(pv0);
			// double dcav0toPV = rdotp*rdotp/pv0.Mag2();
			// dcav0toPV = sqrt(xv0toPV.Mag2() - dcav0toPV);
			// double v0decaylength = xv0toPV.Mag();
			// double v0cosrdotp = rdotp/v0decaylength/pv0.Mag();// cout<<"SCHEME 1: DecayLength = "<<v0decaylength<<";  ";
			if (IfTree){
				bool CheckPass = true;
				vector<int> Temp;Temp.resize(0);
				vector<int> TempT;TempT.resize(0);TempT.push_back(iKFParticle);
				if      ((abs(particle.GetPDG()) == PhiPdg)    && (fabs(particle.GetMass() - PhiPdgMass)    > 9*PhiPdgMassSigma))    {CheckPass = false;}
				else if ((abs(particle.GetPDG()) == K0SPdg)    && (fabs(particle.GetMass() - K0SPdgMass)    > 9*K0SPdgMassSigma))    {CheckPass = false;}
				else if ((abs(particle.GetPDG()) == LambdaPdg) && (fabs(particle.GetMass() - LambdaPdgMass) > 9*LambdaPdgMassSigma)) {CheckPass = false;}
				else if ((abs(particle.GetPDG()) == XiPdg)     && (fabs(particle.GetMass() - XiPdgMass)     > 9*XiPdgMassSigma))     {CheckPass = false;}
				else if ((abs(particle.GetPDG()) == OmegaPdg)  && (fabs(particle.GetMass() - OmegaPdgMass)  > 9*OmegaPdgMassSigma))  {CheckPass = false;}
				if (CheckPass == true) {
					for (int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++){
						TempT.push_back(particle.DaughterIds()[iDaughter]);
					}
					int Itr = 1;
					cout<<"_____________________________________________________"<<endl;
					while (Itr < TempT.size()) {
						if ((KFParticleInterface->GetParticles()[TempT[Itr]]).GetPDG() == -1){
							KFParticle daughter = KFParticleInterface->GetParticles()[TempT[Itr]];
							for (int iDaughter=0; iDaughter < daughter.NDaughters(); iDaughter++){
								if (daughter.DaughterIds()[iDaughter] == TempT[Itr]) continue;
								TempT.push_back(daughter.DaughterIds()[iDaughter]);
								if (Itr > 2) cout<<"TempT = ";StKFParticleAnalysisMaker::print(TempT);
							}
						}else{
							Temp.push_back(TempT[Itr]);
							if (Itr > 2) cout<<"Temp = ";StKFParticleAnalysisMaker::print(Temp);
						}
						Itr++;
					}
					for (int i = 0;i<Temp.size();i++){
						for (int j = i+1;j<Temp.size();j++){
							if (Temp[i] == Temp[j]) {
								CheckPass = false;
								break;
							}
						}
					}
					if (Itr > 2) cout<<"====================================================="<<endl;
				}
				if (CheckPass == true){
					Recorded_KFP_ID.push_back(Temp);
				}
				// if (CheckPass == true) { // cuts for Pions used to reconstruct K0S
				// 	KFParticle NKFParticle = KFParticleInterface->GetParticles()[Temp[0]];
				// 	if ((NKFParticle.GetPDG() == 310)) {
				// 		for (int iDaughter = 1;iDaughter<Temp.size();iDaughter++){
				// 			if (abs((KFParticleInterface->GetParticles()[Temp[iDaughter]]).GetPDG()) != 211) continue;
				// 			int iTrack = 0;
				// 			const int globalTrackId = (KFParticleInterface->GetParticles()[Temp[iDaughter]]).DaughterIds()[0];
				// 			Int_t iTrackStart = globalTrackId - 1;
				// 			if (globalTrackId >= nTracks) {iTrackStart = nTracks - 1;}
				// 			for (Int_t jTrack = iTrackStart;jTrack >= 0;jTrack--){
				// 				StPicoTrack *track = mPicoDst->track(jTrack);
				// 				if (track->id() == globalTrackId){
				// 					iTrack = jTrack;
				// 					break;
				// 				}
				// 			}
				// 			StPicoTrack *track = mPicoDst->track(iTrack);
				// 			std::vector<int> TestPDG;TestPDG.push_back((KFParticleInterface->GetParticles()[Temp[iDaughter]]).GetPDG());
				// 			std::vector<bool> PDGBool = StKFParticleAnalysisMaker::TrackPID(TestPDG , track , Vertex3D);
				// 			if(PDGBool[0] == false) CheckPass = false;
				// 		}
				// 	}
				// }
			}
			// cout<<"CrefMult:"<<CrefMult<<endl;
			// cout<<"PDG:"<<particle.GetPDG()<<endl; 

			// cout<<"Here is good 1"<<endl;

			int upQ; // cout<<"Here is good 2"<<endl;
			if (particle.GetPDG() == LambdaPdg) {upQ = 1;} 
			else if (particle.GetPDG() == -1*LambdaPdg) {upQ = -1;} 
			else continue;
		
		
			int eLambda = -(upQ-1)/2; // 0 if Lambda, 1 if AntiLambda

			SetDaughterTrackPointers(iKFParticle);
			if (ProtonTrackIndex == -99999 || PionTrackIndex == -99999) continue; if(!ProtonTrack) continue; if(!PionTrack) continue;

			double dmass = -999; // just a placeholder
			TLorentzVector p4Pair, p4Proton; // just a placeholder
			StLambdaDecayPair TmpLambdaDecayPair(p4Pair, p4Proton, ProtonTrackIndex, PionTrackIndex, (eLambda==0), dmass);
			KFParticleLambdaDecayPair.push_back(TmpLambdaDecayPair);
		} // End loop over KFParticles
		if (IfTree) {
			for (int i=0;i<Recorded_KFP_ID.size();i++) {
				for (int j=1;j<Recorded_KFP_ID[i].size();j++) {
					const KFParticle particle = KFParticleInterface->GetParticles()[Recorded_KFP_ID[i][j]];
					if ( (abs(particle.GetPDG()) == 211)  || 
						 (abs(particle.GetPDG()) == 2212) || 
						 (abs(particle.GetPDG()) == 321) ) {
						int iTrack = 0;
						const int globalTrackId = (particle).DaughterIds()[0];
						Int_t iTrackStart = globalTrackId - 1;
						if (globalTrackId >= nTracks) {iTrackStart = nTracks - 1;}
						for (Int_t jTrack = iTrackStart;jTrack >= 0;jTrack--){
							StPicoTrack *track = mPicoDst->track(jTrack);
							if (track->id() == globalTrackId){
								iTrack = jTrack;
								break;
							}
						}
						Recorded_KFP_ID[i][j] = iTrack;
					}
					else{
						Recorded_KFP_ID[i][j] = -1;
					}
				}
			}
		}// 自此，Recorded_KFP_ID[:][0]是KFP中的位置，其余为DST中的位置或者标识错误的-1

	}
	
	for (int iKFParticle = 0;iKFParticle<Recorded_KFP_ID.size();iKFParticle++){
		if (IfTree){ // Model
			KFParticle particle = KFParticleInterface->GetParticles()[Recorded_KFP_ID[iKFParticle][0]];
			//SCHEME 2:
			KFParticle tempParticle(particle);
			float l,dl;
			KFParticle pv(KFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
			// pv += particle;
			tempParticle.SetProductionVertex(pv);
			tempParticle.GetDecayLength(l, dl);// cout<<"SCHEME 2: DecayLength = "<<l<<";  ";if (fabs(v0decaylength/l)>1.15 || fabs(v0decaylength/l)<0.95){cout<<particle.GetPDG()<<"  "<<particle.GetMass()<<endl;}else{cout<<" "<<endl;}

			if (IfHelix && ((abs(particle.GetPDG()) == OmegaPdg) || (abs(particle.GetPDG()) == XiPdg))) {

				// helix
				TVector3 MomentumOfParticle(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 PositionOfParticle(particle.GetX(), particle.GetY(), particle.GetZ());
				TLorentzVector OmegaLorentz(MomentumOfParticle, particle.GetE());
				StPicoPhysicalHelix heliPositionOfParticle(MomentumOfParticle, PositionOfParticle, magnet*kilogauss, particle.GetQ());
				double pathlength = heliPositionOfParticle.pathLength(Vertex3D, false);
				TVector3 MomentumOfParticle_tb = heliPositionOfParticle.momentumAt(pathlength, magnet*kilogauss); 
				PDG.emplace_back(particle.GetPDG());
				px.emplace_back(MomentumOfParticle_tb.X());
				py.emplace_back(MomentumOfParticle_tb.Y());
				pz.emplace_back(MomentumOfParticle_tb.Z());
				float pt = pow(pow(MomentumOfParticle_tb.X(),2) + pow(MomentumOfParticle_tb.Y(),2),0.5);
				float TanTheta = 1.0*pt/fabs(MomentumOfParticle_tb.Z());
				float Theta = atan(TanTheta);
				QA_eta.emplace_back(-log(tan(0.5*Theta)) * fabs(MomentumOfParticle_tb.Z())/MomentumOfParticle_tb.Z());
				InvariantMass.emplace_back(particle.GetMass());//cout<<"particle.GetAtProductionVertex() = "<<particle.GetAtProductionVertex()<<endl;
				// float DL = 0. , eDL = 0.;particle.GetDecayLength(DL,eDL);
				QA_Decay_Length.emplace_back(l);
				OmegaVec.push_back(particle);ParticleVec.push_back(particle);
				QA_Chi2.emplace_back(particle.GetChi2());
				QA_DCA_V0_PV.emplace_back(-999);
				QA_dEdx.emplace_back(-999);
				QA_nSigmaProton.emplace_back(-999);
				QA_nSigmaPion.emplace_back(-999);
				QA_nSigmaKaon.emplace_back(-999);
				// cout<<"particle.GetPz()="<<particle.GetPz()<<", "<<"MomentumOfParticle_tb.Z()="<<MomentumOfParticle_tb.Z()<<endl; 
				// cout<<"kilogauss = "<<kilogauss<<endl;
				// cout<<"MomentumOfParticle_tb.Mag() = "<<MomentumOfParticle_tb.Mag()<<endl;
				// cout<<"MomentumOfParticle.Mag() = "<<MomentumOfParticle.Mag()<<endl;
				QA_m2.emplace_back(-999);
			}
			else if ((abs(particle.GetPDG()) == LambdaPdg) || (abs(particle.GetPDG()) == K0SPdg) || (abs(particle.GetPDG()) == PhiPdg))
			{
				PDG.emplace_back(particle.GetPDG());
				px.emplace_back(particle.GetPx());
				py.emplace_back(particle.GetPy());
				pz.emplace_back(particle.GetPz());
				QA_eta.emplace_back(particle.GetEta());
				InvariantMass.emplace_back(particle.GetMass());
				QA_Chi2.emplace_back(particle.GetChi2());
				QA_Decay_Length.emplace_back(l);
				QA_DCA_V0_PV.emplace_back(-1);
				QA_dEdx.emplace_back(-999);
				QA_nSigmaProton.emplace_back(-999);
				QA_nSigmaPion.emplace_back(-999);
				QA_nSigmaKaon.emplace_back(-999);
				QA_m2.emplace_back(-999);

			}
		}
	}

	std::vector<int> NeedPDG; NeedPDG.resize(0);
	NeedPDG.push_back( 2212);NeedPDG.push_back( 211);NeedPDG.push_back( 321);
	NeedPDG.push_back(-2212);NeedPDG.push_back(-211);NeedPDG.push_back(-321);
	std::vector<int> track_index;
	double Total_Pz = 0.0, Total_Px = 0.0, Total_Py = 0.0;
	// Filling Track
	if (IfQAMode) {
		hNch_per_VertexZ->Fill(VertexZ,NumCharge);
		H_Nch_triggerBIN[TriggerID_in_TriggerList][2]->Fill(Nch);
	}
	
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
		StPicoTrack *track = mPicoDst->track(iTrack);
		if (IfQAMode) {hdEdx_pQ->Fill(1.0*track->charge()*track->gMom().Mag(),track->dEdx());}
    	if (! track)            continue;
    	if (! track->charge())  continue;
    	if (  track->nHitsFit() < 15) continue;
		if (  track->nHitsDedx() < 15) continue;
		if (  track->nHitsFit()*1.0 / track->nHitsMax() < 0.52 || track->nHitsFit()*1.0 / track->nHitsMax() > 1.05) continue;
		// if (  track->nHitsFit()*1.0 / track->nHitsMax() < 0.52 || track->nHitsFit()*1.0 / track->nHitsMax() > 0.9) continue;
		if (  track->dEdxError() < 0.04 || track->dEdxError() > 0.12) continue; // same as kfp
		if (! track->isPrimary()) continue;
		track_index.push_back(iTrack);
		// H_All_nSigmaKaon_y->Fill();
		Total_Pz += track->gMom().Z();
		Total_Px += track->gMom().X();
		Total_Py += track->gMom().Y();


		// track info
		float p = track->gMom().Mag();
		float pt = track->gMom().Perp();
		float phi = track->gMom().Phi();
		float eta = track->gMom().Eta();
		double track_px = track->gMom().X();
		double track_py = track->gMom().Y();
		double track_pz = track->gMom().Z();
		// float dcatopv = track->gDCA(Vertex3D).Mag();
		// float nSigmaKaon = track->nSigmaKaon();
		// float nSigmaPion = track->nSigmaPion();
		// float nSigmaProton = track->nSigmaProton();
		// float pt_prim = track->pMom().Perp();
		// float phi_prim = track->pMom().Phi();
		// float eta_prim = track->pMom().Eta();

		// TOF Info
		// bool hasTOF = false;
		// int tofindex = track->bTofPidTraitsIndex();
		// float m2 = -999.;
		// float beta = -999.;
		// if (tofindex >= 0) 
		// {
		// 	int tofflag = (mPicoDst->btofPidTraits(tofindex))->btofMatchFlag();
		// 	float tof = (mPicoDst->btofPidTraits(tofindex))->btof();
		// 	float BtofYLocal = (mPicoDst->btofPidTraits(tofindex))->btofYLocal();
		// 	// hgbtofYlocal->Fill(BtofYLocal);
		// 	if((tofflag >= 1) && (tof > 0) && (BtofYLocal > -1.8) && (BtofYLocal < 1.8)) hasTOF = true;
		// }
		// StPicoPhysicalHelix helix = track->helix(magnet);
		// TVector3 pkaon = helix.momentum(magnet*kilogauss);
		// if (hasTOF)
		// {
		// 	beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
		// 	m2 = pkaon.Mag2()*(1.0 / beta / beta - 1.0);

		// 	// some kaon QA
		// 	//if (track->nSigmaKaon() >  6) hgptm2_largenSigmaKaon->Fill(track->gMom().Perp(), m2);
		// 	//if (track->nSigmaKaon() < -6) hgptm2_smallnSigmaKaon->Fill(track->gMom().Perp(), m2);
		// 	zTOF_proton = 1/beta - sqrt(ProtonPdgMass*ProtonPdgMass/pkaon.Mag2()+1);
		// 	zTOF_pion   = 1/beta - sqrt(PionPdgMass*PionPdgMass/pkaon.Mag2()+1);
		// 	zTOF_kaon   = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
		// }

		// Fill tracks
		// bool IfRecordThisTrack = false;
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
		// bool proton_cut = true;
		// if (fabs(nSigmaProton) > 2) proton_cut = false;
		// if (pt < pT_trig_lo || pt > pT_trig_hi) proton_cut = false; 
		// if (fabs(eta_prim) > eta_trig_cut) proton_cut = false;
		// if (!hasTOF && pt > proton_pT_TOFth) proton_cut = false;
		// if (pt > proton_pT_TOFth && (m2 > proton_m2_hi || m2 < proton_m2_lo)) proton_cut = false;
		// ProtonPID proton_pid(0., nSigmaProton, pt); // not using zTOF
		// if (!proton_pid.IsProtonSimple(2., track->charge())) proton_cut = false; // only 0.2 < pt < 2.0!!!
		// if (dcatopv > dcatoPV_hi) proton_cut = false;
		
		// // Test if Pion
		// bool pion_cut = true;
		// if (fabs(nSigmaPion) > 2) pion_cut = false;
		// if (pt < pT_trig_lo || pt > pT_trig_hi) pion_cut = false; // use p < 2
		// if (fabs(eta_prim) > eta_trig_cut) pion_cut = false;
		// PionPID pion_pid(0., nSigmaPion, pt); // not using zTOF
		// if (!pion_pid.IsPionSimple(2., track->charge())) pion_cut = false; // only 0.2 < pt < 2.0!!!
		// if (dcatopv > dcatoPV_hi) pion_cut = false;
		// if (!hasTOF && pt > pion_pT_TOFth) pion_cut = false;
		// if (pt > pion_pT_TOFth && (m2 > pion_m2_hi || m2 < pion_m2_lo)) pion_cut = false;

		// // Test if Kaon
		// bool kaon_cut = true;
		// if (fabs(nSigmaKaon) > 2) kaon_cut = false;
		// if (pt < pT_trig_lo || pt > 1.6) kaon_cut = false; // use p < 1.6
		// if (fabs(eta_prim) > eta_trig_cut) kaon_cut = false;
		// if (!hasTOF && pt > 0.4) kaon_cut = false;
		// if (pt > 0.4 && (m2 > 0.34 || m2 < 0.15)) kaon_cut = false;
		// double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
		// KaonPID kaon_pid(zTOF, nSigmaKaon, pt); // not using zTOF
		// if (!kaon_pid.IsKaonSimple(2., track->charge())) kaon_cut = false; // only 0.2 < pt < 2.0!!!
		// if (dcatopv > 2.0) kaon_cut = false;
		// // for (int i = 0; i < OmegaVec.size(); i++) if (IsKaonOmegaDaughter(OmegaVec[i], track->id())) kaon_cut = false;
		// // for (int i = 0; i < OmegaVec.size(); i++) if (IsTrackParticleDaughter(OmegaVec[i], track->id())) kaon_cut = false;

		// if (proton_cut + pion_cut + kaon_cut == 1) {IfRecordThisTrack = true;QA_IfConfuse.emplace_back(0);}
		// if (proton_cut + pion_cut + kaon_cut > 1){IfRecordThisTrack = true;QA_IfConfuse.emplace_back(1);}

//////////////////////////////////// Used for test //////////////////////////////////////////////////////////////////////////////////////
		// Raw Data bTOF
		bool RawTOF = true;
		bool hasTOF = false;
		bool IfeTof = false;
		float m2 = -999.;
		if (RawTOF){
			int tofindex = track->bTofPidTraitsIndex();
			float beta = -999.;
			if (tofindex >= 0) 
			{
				int tofflag = (mPicoDst->btofPidTraits(tofindex))->btofMatchFlag();
				float tof = (mPicoDst->btofPidTraits(tofindex))->btof();
				float BtofYLocal = (mPicoDst->btofPidTraits(tofindex))->btofYLocal();
				if((tofflag >= 1) && (tof > 0) && (BtofYLocal > -1.8) && (BtofYLocal < 1.8)) hasTOF = true;
				// if((tofflag == 1) && (tof > 0)) hasTOF = true;
				// if ((BtofYLocal < -1.8) && (BtofYLocal > 1.8)) {cout<<"BtofYLocal = "<<BtofYLocal<<endl;}
				// if ((mPicoDst->btofPidTraits(tofindex))->nSigmaKaon() != 0){
				// 	cout<<"nSigmaKaon() = "<<(mPicoDst->btofPidTraits(tofindex))->nSigmaKaon()<<endl;
				// }
			}
			StPicoPhysicalHelix helix = track->helix(magnet);
			TVector3 pkaon = helix.momentum(magnet*kilogauss);
			if (hasTOF)
			{
				beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
				m2 = pkaon.Mag2()*(1.0 / beta / beta - 1.0);
				if (IfQAMode) {
					H_Pt_m2->Fill(track->gMom().Mag(),m2);
					H_Pt_nSigmaKaonTOF->Fill(track->gMom().Mag(),(mPicoDst->btofPidTraits(tofindex))->nSigmaKaon());
					H_m2_nSigmaKaon_Pt->Fill(m2,track->nSigmaKaon(),pt);
					// cout<<"nsigmaTOF = "<<(mPicoDst->btofPidTraits(tofindex))->nSigmaKaon()<<endl;
					if (fabs(1/beta-1)<0.03) {
						hdEdx_pQ_1cut->Fill(1.0*track->charge()*track->gMom().Mag(),track->dEdx());
					}
				}
			}
		}
		// Raw Data eTOF
		if (!hasTOF && RawTOF){
			int tofindex = track->eTofPidTraitsIndex();
			// cout<<"tofindex = "<<tofindex<<endl;
			float beta = -999.;
			if (tofindex >= 0){
				int tofflag = (mPicoDst->etofPidTraits(tofindex))->matchFlag();
				float tof = (mPicoDst->etofPidTraits(tofindex))->tof();
				float EtofYLocal = (mPicoDst->etofPidTraits(tofindex))->deltaY();
				// cout<<"EtofYLocal = "<<EtofYLocal<<endl;
				if((tofflag >= 1) && (tof > 0) && (EtofYLocal > -100) && (EtofYLocal < 100)) hasTOF = true;
				hasTOF = true;
			}
			StPicoPhysicalHelix helix = track->helix(magnet);
			TVector3 pkaon = helix.momentum(magnet*kilogauss);
			if (hasTOF)
			{
				IfeTof = true;
				beta = (mPicoDst->etofPidTraits(tofindex))->beta();
				m2 = pkaon.Mag2()*(1.0 / beta / beta - 1.0);
				if (IfQAMode) {
					H_Pt_m2->Fill(track->gMom().Mag(),m2);
					H_Pt_nSigmaKaonTOF->Fill(track->gMom().Mag(),(mPicoDst->btofPidTraits(tofindex))->nSigmaKaon());
					H_m2_nSigmaKaon_Pt->Fill(m2,track->nSigmaKaon(),pt);
					// cout<<"nsigmaTOF = "<<(mPicoDst->btofPidTraits(tofindex))->nSigmaKaon()<<endl;
					if (fabs(1/beta-1)<0.03) {
						hdEdx_pQ_1cut->Fill(1.0*track->charge()*track->gMom().Mag(),track->dEdx());
					}
				}
			}
		}
		if (IfQAMode) {
			H_Pt_nSigmaKaon->Fill(track->gMom().Mag(),track->nSigmaKaon());
			H_eta_nSigmaKaon  [TriggerID_in_TriggerList][2]->Fill(eta,track->nSigmaKaon());
			H_eta_nSigmaPion  [TriggerID_in_TriggerList][2]->Fill(eta,track->nSigmaPion());
			H_eta_nSigmaProton[TriggerID_in_TriggerList][2]->Fill(eta,track->nSigmaProton());
			H_eta_m2          [TriggerID_in_TriggerList][2]->Fill(eta,m2);
			H_eta_PVz         [TriggerID_in_TriggerList][2]->Fill(eta,VertexZ);
			H_eta_PVr         [TriggerID_in_TriggerList][2]->Fill(eta,VertexR);
			H_eta_DVz         [TriggerID_in_TriggerList][2]->Fill(eta,DVz);
			H_eta_triggerBIN  [TriggerID_in_TriggerList][2]->Fill(eta);
			if (hasTOF) {H_eta_triggerBIN_hasTOF[TriggerID_in_TriggerList][2]->Fill(eta);}
			if (track->charge() > 0) {
				H_eta_nSigmaKaon  [TriggerID_in_TriggerList][0]->Fill(eta,track->nSigmaKaon());
				H_eta_nSigmaPion  [TriggerID_in_TriggerList][0]->Fill(eta,track->nSigmaPion());
				H_eta_nSigmaProton[TriggerID_in_TriggerList][0]->Fill(eta,track->nSigmaProton());
				H_eta_m2          [TriggerID_in_TriggerList][0]->Fill(eta,m2);
				H_eta_PVz         [TriggerID_in_TriggerList][0]->Fill(eta,VertexZ);
				H_eta_PVr         [TriggerID_in_TriggerList][0]->Fill(eta,VertexR);
				H_eta_DVz         [TriggerID_in_TriggerList][0]->Fill(eta,DVz);
				H_eta_triggerBIN  [TriggerID_in_TriggerList][0]->Fill(eta);
				if (hasTOF) {H_eta_triggerBIN_hasTOF[TriggerID_in_TriggerList][0]->Fill(eta);}
			}
			if (track->charge() < 0) {
				H_eta_nSigmaKaon  [TriggerID_in_TriggerList][1]->Fill(eta,track->nSigmaKaon());
				H_eta_nSigmaPion  [TriggerID_in_TriggerList][1]->Fill(eta,track->nSigmaPion());
				H_eta_nSigmaProton[TriggerID_in_TriggerList][1]->Fill(eta,track->nSigmaProton());
				H_eta_m2          [TriggerID_in_TriggerList][1]->Fill(eta,m2);
				H_eta_PVz         [TriggerID_in_TriggerList][1]->Fill(eta,VertexZ);
				H_eta_PVr         [TriggerID_in_TriggerList][1]->Fill(eta,VertexR);
				H_eta_DVz         [TriggerID_in_TriggerList][1]->Fill(eta,DVz);
				H_eta_triggerBIN  [TriggerID_in_TriggerList][1]->Fill(eta);
				if (hasTOF) {H_eta_triggerBIN_hasTOF[TriggerID_in_TriggerList][1]->Fill(eta);}
			}
			H_eta_trigger     ->Fill(eta,TriggerID_in_TriggerList);
		}

		if (pt > 1.4) {continue;}
		std::vector<bool> PDGBool = StKFParticleAnalysisMaker::TrackPID(NeedPDG , track , Vertex3D);
		for (int Ktr = 0;Ktr < PDGBool.size();Ktr++) {
			if (PDGBool[Ktr] == true) {
				for (int Itr = PDG2NameSize;Itr < PDG2NameSize + PDG2NameSize2;Itr++){
					int Jtr = Itr - PDG2NameSize;
					if (NeedPDG[Ktr] != PDGList[Itr]){continue;}
					//// For Kaon
					if (abs(PDGList[Itr])==KaonPdg) {
						// if (
						// 	(
						// 		((0.196 <= m2)&&(m2 <= 0.292)) // Tight 0.5 < $p_t$ < 0.6 GeV
						// 		// ((0.123 <= m2)&&(m2 <= 0.359)) // Soft 1.0 < $p_t$ < 1.1 GeV
						// 		// ((0.0674 <= m2)&&(m2 <= 0.407)) // Soft 1.2 < $p_t$ < 1.3 GeV
						// 		//  || (fabs(track->nSigmaPion())>3&&fabs(track->nSigmaProton())>3)
						// 	) == false
						// )
						// {continue;}
						std::vector<float> m2Zone = TPCandTOF_Gen.KaonTOFm2(pt, DataName);
						if (
							(
								// (0.2 <= pt && pt < 0.3) ||
								((m2Zone[0] <= m2)&&(m2 <= m2Zone[1]))
							) == false
						)
						{continue;}
					}
					//// For Pion
					if (abs(PDGList[Itr])==PionPdg) {
						if (
							(
								true
								// ((-0.012 <= m2)&&(m2 <= 0.049)) // Tight 0.5 < $p_t$ < 0.6 GeV
								//  || (fabs(track->nSigmaPion())>3&&fabs(track->nSigmaProton())>3)
							) == false
						)
						{continue;}
					}
					//// For Proton
					if (abs(PDGList[Itr])==ProtonPdg) {
						if (
							(
								((0.746 <= m2)&&(m2 <= 1.015)) // Tight 0.5 < $p_t$ < 0.6 GeV
								//  || (fabs(track->nSigmaPion())>3&&fabs(track->nSigmaProton())>3)
							) == false
						)
						{continue;}
					}
					//// 
					float tEnergy = pow(pow(track->gMom().Mag(),2) + pow(StKFParticleAnalysisMaker::massList(NeedPDG[Ktr]),2),0.5);
					float rap = 0.5*log((tEnergy+track_pz)/(tEnergy-track_pz));
					if (IfTree) {
						std::vector<int> Temp;Temp.resize(0);
						Temp.push_back(iTrack);Temp.push_back(iTrack);
						Recorded_KFP_ID.push_back(Temp);
						QA_Chi2.emplace_back(-999);
						QA_Decay_Length.emplace_back(-999);
						PDG.emplace_back(NeedPDG[Ktr]);
						px.emplace_back(track_px);
						py.emplace_back(track_py);
						pz.emplace_back(track_pz);
						QA_eta.emplace_back(eta);
						QA_dEdx.emplace_back(track->dEdx());
						QA_nSigmaProton.emplace_back(track->nSigmaProton());
						QA_nSigmaPion.emplace_back(track->nSigmaPion());
						QA_nSigmaKaon.emplace_back(track->nSigmaKaon());
						QA_DCA_V0_PV.emplace_back(track->gDCA(Vertex3D).Mag());
						QA_m2.emplace_back(m2);
						InvariantMass.emplace_back(-999); 
					}
					if (IfQAMode) {
						H_Pt[Jtr] -> Fill(pt);
						H_P[Jtr] -> Fill(p);
						H_rapidity[Jtr]->Fill(rap);
						if (IfeTof) H_rapidity_eTOF[Jtr]->Fill(rap);
						H_y_Pt[Jtr]->Fill(rap,pt);
						H_y_P[Jtr]->Fill(rap,p);
						if (hasTOF) H_y_m2[Jtr]->Fill(rap,m2);
						H_y_nSigmaKaon[Jtr]->Fill(rap,track->nSigmaKaon());
						H_y_nSigmaPion[Jtr]->Fill(rap,track->nSigmaPion());
						H_y_nSigmaElectron[Jtr]->Fill(rap,track->nSigmaElectron());
						H_y_nHitsFit[Jtr]->Fill(rap,track->nHitsFit());
						H_y_nHitsDedx[Jtr]->Fill(rap,track->nHitsDedx());
						H_y_nHitsFit2nHitsMax[Jtr]->Fill(rap,track->nHitsFit()*1.0 / track->nHitsMax());
						H_y_eta[Jtr]->Fill(rap,eta);
						H_y_Pz[Jtr]->Fill(rap,track_pz);
						H_Pxy[Jtr]->Fill(track_px,track_py);
						H_Pz[Jtr]->Fill(track_pz);
						H_y_Vz[Jtr]->Fill(rap,VertexZ);
						H_y_Nch[Jtr]->Fill(rap,NumCharge);
						H_Pz_Nch[Jtr]->Fill(track_pz,NumCharge);
						if(abs(PDGList[Itr])==321 && hasTOF) {
							H_nSigmaTOF_p[Jtr]->Fill((mPicoDst->btofPidTraits(track->bTofPidTraitsIndex()))->nSigmaKaon(),track->gMom().Mag());
							hgbtofYlocal[Jtr]->Fill(rap,(mPicoDst->btofPidTraits(track->bTofPidTraitsIndex()))->btofYLocal());
							H_y_nSigmaTOFKaon[Jtr]->Fill(rap,(mPicoDst->btofPidTraits(track->bTofPidTraitsIndex()))->nSigmaKaon());
							H_m2_nSigmaTOFKaon[Jtr]->Fill(m2,(mPicoDst->btofPidTraits(track->bTofPidTraitsIndex()))->nSigmaKaon());
						}
						H_dEdx_p[Jtr]->Fill(1.0*track->charge()*track->gMom().Mag(),track->dEdx());
						H_DCAtoPV[Jtr]->Fill(track->gDCA(Vertex3D).Mag());
						H_eta[Jtr]->Fill(eta);
						H_nHitsFit_p[Jtr]->Fill(track->nHitsFit(),track->gMom().Mag());
						H_nHitsFit_nHitsMax[Jtr]->Fill((track->nHitsFit()*1.0)/(track->nHitsMax()*1.0));
						H_ndEdx[Jtr]->Fill((track->nHitsDedx()));
						for (int Ntr=0;Ntr<PDG2NameSize3;Ntr++){
							// cout<<"CNameList["<<Ktr<<"] = "<<CNameList[Ktr]<<endl;
							if ( CNameList[Ntr] == "Kaon"){H_Pt_nSigma[Jtr][Ntr]->Fill(track->nSigmaKaon(),track->gMom().Perp());InvariantMass.emplace_back(KaonPdgMass);}
							if ( CNameList[Ntr] == "Pion"){H_Pt_nSigma[Jtr][Ntr]->Fill(track->nSigmaPion(),track->gMom().Perp());InvariantMass.emplace_back(PionPdgMass);}
							if ( CNameList[Ntr] == "Proton"){H_Pt_nSigma[Jtr][Ntr]->Fill(track->nSigmaProton(),track->gMom().Perp());InvariantMass.emplace_back(ProtonPdgMass);}
						}
					}
					
					break;
				}
			}
		}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// if (IfRecordThisTrack == true) {
		// 	hdEdx_pQ_2cut->Fill(1.0*track->charge()*track->gMom().Mag(),track->dEdx());

		// 	px.emplace_back(track->gMom().X());
		// 	py.emplace_back(track->gMom().Y());
		// 	pz.emplace_back(track->gMom().Z());
		// 	if      (proton_cut) {
		// 		IfRecordThisTrack = true;
		// 		if (track->charge() > 0) {PDG.emplace_back( 2212);}
		// 		else                     {PDG.emplace_back(-2212);}
		// 		InvariantMass.emplace_back(ProtonPdgMass);
		// 	}
		// 	else if (pion_cut) {
		// 		IfRecordThisTrack = true;
		// 		if (track->charge() > 0) {PDG.emplace_back( 211);}
		// 		else                     {PDG.emplace_back(-211);}
		// 		InvariantMass.emplace_back(PionPdgMass);
		// 	}
		// 	else if (kaon_cut) {
		// 		IfRecordThisTrack = true;
		// 		if (track->charge() > 0) {PDG.emplace_back( 321);}
		// 		else                     {PDG.emplace_back(-321);}
		// 		InvariantMass.emplace_back(KaonPdgMass);
		// 	}

		// 	// Filling QA
		// 	if (hasTOF) {
		// 		QA_zTOF_proton.emplace_back(zTOF_proton);QA_zTOF_pion.emplace_back(zTOF_pion);QA_zTOF_kaon.emplace_back(zTOF_kaon);
		// 		QA_m2.emplace_back(m2);
		// 	}
		// 	else{
		// 		QA_zTOF_proton.emplace_back(0.);QA_zTOF_pion.emplace_back(0.);QA_zTOF_kaon.emplace_back(0.);
		// 		QA_m2.emplace_back(0.);
		// 	}
		// 	QA_nSigmaProton.emplace_back(nSigmaProton);
		// 	QA_nSigmaPion.emplace_back(nSigmaPion);
		// 	QA_nSigmaKaon.emplace_back(nSigmaKaon);
		// 	QA_dEdx.emplace_back(track->dEdx());QA_DCA_V0_PV.emplace_back(dcatopv);
		// 	QA_Decay_Length.emplace_back(-1.0);
		// 	QA_Chi2.emplace_back(-1.0);
		// 	QA_IfBadReconstructed.emplace_back(-1);
		// 	QA_DCA_Daughters.emplace_back(-1.0);
		// }

	}
	if (IfTree) {
		Correlatted_ID_List_T.resize(0);
		for (int iRecorded_KFP=0;iRecorded_KFP<Recorded_KFP_ID.size();iRecorded_KFP++){
			std::vector<int> Temp;Temp.resize(0);
			Correlatted_ID_List_T.push_back(Temp);
		}
		for (int iRecorded_KFP=0;iRecorded_KFP<Recorded_KFP_ID.size();iRecorded_KFP++){
			for (int jRecorded_KFP=iRecorded_KFP+1;jRecorded_KFP<Recorded_KFP_ID.size();jRecorded_KFP++){
				bool IfCorrelated = false;
				for (int kRecorded_KFP=1;kRecorded_KFP < Recorded_KFP_ID[iRecorded_KFP].size();kRecorded_KFP++){
					if (IfCorrelated == true) break;
					for (int nRecorded_KFP=1;nRecorded_KFP < Recorded_KFP_ID[jRecorded_KFP].size();nRecorded_KFP++){
						if ( Recorded_KFP_ID[iRecorded_KFP][kRecorded_KFP] == Recorded_KFP_ID[jRecorded_KFP][nRecorded_KFP] ){
							Correlatted_ID_List_T[iRecorded_KFP].push_back(jRecorded_KFP);
							Correlatted_ID_List_T[jRecorded_KFP].push_back(iRecorded_KFP);
							IfCorrelated = true;
							break;
						}
					}
				}
			}
		}
		for (int iRecorded_KFP=0;iRecorded_KFP<Correlatted_ID_List_T.size();iRecorded_KFP++){
			switch (Correlatted_ID_List_T[iRecorded_KFP].size()){
				case 0:
					ParentA.emplace_back(-1);
					ParentB.emplace_back(-1);
					ParentC.emplace_back(-1);
					ParentD.emplace_back(-1);
					ParentE.emplace_back(-1);
					break;
				case 1:
					ParentA.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][0]);
					ParentB.emplace_back(-1);
					ParentC.emplace_back(-1);
					ParentD.emplace_back(-1);
					ParentE.emplace_back(-1);
					break;
				case 2:
					ParentA.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][0]);
					ParentB.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][1]);
					ParentC.emplace_back(-1);
					ParentD.emplace_back(-1);
					ParentE.emplace_back(-1);
					break;
				case 3:
					ParentA.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][0]);
					ParentB.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][1]);
					ParentC.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][2]);
					ParentD.emplace_back(-1);
					ParentE.emplace_back(-1);
					break;
				case 4:
					ParentA.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][0]);
					ParentB.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][1]);
					ParentC.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][2]);
					ParentD.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][3]);
					ParentE.emplace_back(-1);
					break;
				default:
					ParentA.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][0]);
					ParentB.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][1]);
					ParentC.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][2]);
					ParentD.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][3]);
					ParentE.emplace_back(Correlatted_ID_List_T[iRecorded_KFP][4]);
					break;
			}
		}
		// cout<<"_____________________________________________"<<endl;
		// cout<<"Recorded_KFP_ID              = {"<<endl;
		// for (int i=0;i<Recorded_KFP_ID.size();i++) {
		// 	cout<<"                                  { ";
		// 	for (int j=0;j<Recorded_KFP_ID[i].size();j++) {cout<<Recorded_KFP_ID[i][j];if (j<(Recorded_KFP_ID[i].size()-1)) cout<<" , ";}
		// 	cout<<" }"<<endl;
		// }
		// cout<<"                                }"<<endl;
		// cout<<"PDG.size()                   = "<<PDG.size()<<endl;
		// cout<<"Correlatted_ID_List_T        = {"<<endl;
		// for (int i=0;i<Correlatted_ID_List_T.size();i++) {
		// 	cout<<"                                  { ";
		// 	for (int j=0;j<Correlatted_ID_List_T[i].size();j++) {cout<<Correlatted_ID_List_T[i][j];if (j<(Correlatted_ID_List_T[i].size()-1)) cout<<" , ";}
		// 	cout<<" }"<<endl;
		// }
		// cout<<"                                }"<<endl;
		// cout<<"Correlatted_ID_List.size()   = "<<Correlatted_ID_List.size()<<endl;
		// cout<<"Correlatted_ID_Sta.size()    = "<<Correlatted_ID_Sta.size()<<endl;
		// cout<<"Correlatted_ID_End.size()    = "<<Correlatted_ID_End.size()<<endl;
		// cout<<"_____________________________________________"<<endl;
	}
	// cout<<"Total_Pz = "<<Total_Pz<<endl;
	if (IfQAMode) {
		H_Total_Pz->Fill(Total_Pz);
		H_Total_Pxy->Fill(Total_Px,Total_Py);
	}

// ======= KFParticle end ======= //

// ======= Lambda loop ======= //
	// for(int j=0; j<KFParticleLambdaDecayPair.size(); j++) {
	// 	int i = KFParticleLambdaDecayPair[j].get_idxProton();
	// 	int k = KFParticleLambdaDecayPair[j].get_idxPion();
	// 	if(k == i) continue;

	// 	StPicoTrack* mTrackI = (StPicoTrack*)mPicoDst->track(i);

	// 	int    mchgI = mTrackI->charge();
	// 	int    mhitI = mTrackI->nHitsFit();
	// 	double mMomI = mTrackI->gMom().Mag();
	// 	double mp0xI = mTrackI->gMom().X();
	// 	double mp0yI = mTrackI->gMom().Y();
	// 	double mp0zI = mTrackI->gMom().Z();
	// 	double mpt0I = mTrackI->gMom().Perp();
	// 	double mphiI = mTrackI->gMom().Phi();
	// 	double metaI = mTrackI->gMom().PseudoRapidity();
	// 	double mdcaI = mTrackI->gDCA(Vertex3D).Mag();

	// 	if(mphiI<0) mphiI += 2*M_PI;
	// 	if(mphiI>=2*M_PI) mphiI -= 2*M_PI;

	// 	StPicoTrack* mTrackK = (StPicoTrack*)mPicoDst->track(k);

	// 	int    mchgK = mTrackK->charge();
	// 	int    mhitK = mTrackK->nHitsFit();
	// 	double mMomK = mTrackK->gMom().Mag();
	// 	double mp0xK = mTrackK->gMom().X();
	// 	double mp0yK = mTrackK->gMom().Y();
	// 	double mp0zK = mTrackK->gMom().Z();
	// 	double mpt0K = mTrackK->gMom().Perp();
	// 	double mphiK = mTrackK->gMom().Phi();
	// 	double metaK = mTrackK->gMom().PseudoRapidity();
	// 	double mdcaK = mTrackK->gDCA(Vertex3D).Mag();

	// 	if(mphiK<0) mphiK += 2*M_PI;
	// 	if(mphiK>=2*M_PI) mphiK -= 2*M_PI;

	// 	bool isPP = mchgI>0 && mchgK>0;
	// 	bool isNN = mchgI<0 && mchgK<0;
	// 	bool isPN = mchgI>0 && mchgK<0;
	// 	bool isNP = mchgI<0 && mchgK>0;
	// 	bool isSelfLambda = isPN;
	// 	bool isAntiLambda = isNP;

	// 	// remove the SS cases
	// 	if(isPP || isNN) continue;

	// 	//reconstruction of V0, the parent particle
	// 	TVector3 xv0, op1, op2;
	// 	double dca1to2 = closestDistance(mTrackI, mTrackK, magnet, Vertex3D, xv0, op1, op2);
	// 	TVector3 pv0 = op1 + op2;
	// 	TVector3 xv0toPV = xv0 - Vertex3D;
	// 	double rdotp = xv0toPV.Dot(pv0);
	// 	double dcav0toPV = rdotp*rdotp/pv0.Mag2();
	// 	dcav0toPV = sqrt(xv0toPV.Mag2() - dcav0toPV);
	// 	double v0decaylength = xv0toPV.Mag();
	// 	double v0cosrdotp = rdotp/v0decaylength/pv0.Mag();

	// 	TLorentzVector p4ProtonI, p4PionK;
	// 	p4ProtonI.SetPxPyPzE(op1.X(), op1.Y(), op1.Z(), sqrt(op1.Mag2() + pmass*pmass));
	// 	p4PionK.SetPxPyPzE(  op2.X(), op2.Y(), op2.Z(), sqrt(op2.Mag2() + pimass*pimass));

	// 	// ProtonI & PionK
	// 	TLorentzVector p4Pair = p4ProtonI + p4PionK;
	// 	double massPair = p4Pair.M();
	// 	double ptPair   = p4Pair.Pt();
	// 	double etaPair  = p4Pair.Eta();
	// 	double phiPair  = p4Pair.Phi();
	// }
// ======= Lambda loop ends ======= //

	if (PDG.size()>0){
		PDGMult = PDG.size(); // This is multiplicity of Recorded Particles
		if (Recorded_Hyperon != 0){
			// cout<<"Found Hyperon"<<endl;
			hadronTree->Fill();
		}
		// hadronTree->Fill();
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

	// cout<<"S4"<<endl;
	if ( !KFParticleInterface->ProcessEvent(PicoDst, triggeredTracks,H_ProcessEventNum) ) InterfaceCantProcessEvent = true; else InterfaceCantProcessEvent = false;

	// cout<<"S5"<<endl;
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

std::vector<bool> StKFParticleAnalysisMaker::TrackPID(std::vector<int>& TestPDG , StPicoTrack *track , TVector3 Vertex3D) {
	float TrackID_pt = track->gMom().Perp();
	float TrackID_dcatopv = track->gDCA(Vertex3D).Mag();
	float TrackID_nSigmaKaon = track->nSigmaKaon();
	float TrackID_nSigmaPion = track->nSigmaPion();
	float TrackID_nSigmaProton = track->nSigmaProton();
	float TrackID_eta_prim = track->pMom().Eta();

	float dcatoPV_hi = 3.0; // Upper limit of DCA to PVs
	float pT_trig_lo = 0.2; // 0.2
	float pT_trig_hi = 1.4; // 2
	float eta_trig_cut = 1.0;

	std::vector<bool> result;result.resize(0);
	for (int Itr = 0;Itr < TestPDG.size();Itr++){
		if (abs(TestPDG[Itr]) == 2212){// Proton
			// Test if Proton
			bool proton_cut = true;
			if (fabs(TrackID_nSigmaProton) > 3) proton_cut = false;
			if (TrackID_pt < pT_trig_lo || TrackID_pt > pT_trig_hi) proton_cut = false; 
			// if (fabs(TrackID_eta_prim) > eta_trig_cut) proton_cut = false;
			if (TrackID_dcatopv > dcatoPV_hi) proton_cut = false;

			if (track->charge() * TestPDG[Itr] > 0) {result.push_back(proton_cut);}
			else {result.push_back(false);}
			
		}
		if (abs(TestPDG[Itr]) == 211){// Pion
			// Test if Pion
			bool pion_cut = true;
			if (fabs(TrackID_nSigmaPion) > 3) pion_cut = false;
			if (TrackID_pt < pT_trig_lo || TrackID_pt > pT_trig_hi) pion_cut = false; 
			// if (fabs(TrackID_eta_prim) > eta_trig_cut) pion_cut = false;
			if (TrackID_dcatopv > dcatoPV_hi) pion_cut = false;

			if (track->charge() * TestPDG[Itr] > 0) {result.push_back(pion_cut);}
			else {result.push_back(false);}
			
		}
		if (abs(TestPDG[Itr]) == 321){// Kaon
			// Test if Kaon
			bool kaon_cut = true;
			if (fabs(TrackID_nSigmaKaon - TPCandTOF_Gen.KaonTPCCenter(TrackID_pt , DataName)) > 2) kaon_cut = false;
			// if ((0.2 <= TrackID_pt) && (TrackID_pt < 0.3) && (1 < track->nSigmaKaon()) && (track->nSigmaKaon() < 3)) kaon_cut = true;
			if (TrackID_pt < pT_trig_lo || TrackID_pt > pT_trig_hi) kaon_cut = false; 
			// if (fabs(TrackID_eta_prim) > eta_trig_cut) kaon_cut = false;
			if (TrackID_dcatopv > dcatoPV_hi) kaon_cut = false;

			if (track->charge() * TestPDG[Itr] > 0) {result.push_back(kaon_cut);}
			else {result.push_back(false);}
			
		}
	}
	return result;
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
			// cout<<"daughter ID = "<<daughter.GetPDG()<<endl;
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
			if (result == false) {return result;}

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
	// cout<<"TrackDCA = "<<TrackDCA<<endl;

	if (TrackDCA < Gen1_DCALim) {return result;}
	else {return false;}
}

Double_t StKFParticleAnalysisMaker::massList(int PID)
{
    Double_t Result;
	switch (abs(PID))
	{
		case 2212:
			Result = 0.938272088;
			break;
		case 321 :
			Result = 0.493677;
			break;
		case 211 :
			Result = 0.13957;
			break;
		case 3334 :// OmegaFitMass
			Result = 1.6725;
			break;
		case 3312 :// XiFitMass
			Result = 1.3217;
			break;
		case 3122 :// LambdaFitMass
			Result = 1.11568;
			break;
		default :
			Result = 0;
	}
    return Result;
}

void StKFParticleAnalysisMaker::print(std::vector<int> Temp)
{
	cout<<"{";
    for (int i = 0;i<Temp.size();i++){
		cout<<" "<<Temp[i];
		if (i != (Temp.size() - 1)) cout<<" ,"; 
	}
	cout<<" }"<<endl;
    return ;
}

void StKFParticleAnalysisMaker::print(std::vector<std::vector<int> > Temp)
{
	cout<<"{"<<endl;
    for (int i = 0;i<Temp.size();i++){
		StKFParticleAnalysisMaker::print(Temp[i]);
	}
	cout<<"}"<<endl;
    return ;
}

// StPicoHelix StKFParticleAnalysisMaker::StPicoTrack2StPicoHelix(StPicoTrack* Track){
// 	StPicoHelix Result;
// 	Result.setParameters(Double_t c, Double_t dip, Double_t phase,const TVector3& o, Int_t h)
// }