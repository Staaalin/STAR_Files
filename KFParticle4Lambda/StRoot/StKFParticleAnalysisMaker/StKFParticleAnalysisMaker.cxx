#include "StKFParticleAnalysisMaker.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "StMessMgr.h"
#include <algorithm>
#include <TMath.h>

#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"

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
#define cen_cut            9
#define num_pt_bin         10
#define num_phi_bin        10
#define num_mult_bin       5
#define num_vz_bin         16
#define num_EP_bin         6
#define shift_order_asso   24
#define shift_order_POI    24
#define shift_order_EP     5
#define TPCAssoEtaCut      0.5


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
	mRunStart = 50000;
	mRunEnd   = 100000;

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

	// correlation centrality
	min_cent = 8;
	max_cent = 9;

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

	//////////////////////////////
    int buffer_size = 5000000;
    TTree *hadronTree = new TTree("hadronTree", "Tree_STAR");
    hadronTree->Branch("buffer_size"   ,&buffer_size      ,"buffer_size/I"                   );
    hadronTree->Branch("refMult"       ,&CrefMult         ,"refMult/I"                       );
    hadronTree->Branch("grefMult"      ,&CgrefMult        ,"grefMult/I"                      );
    hadronTree->Branch("EventID"       ,&evtID            ,"EventID/I"                       );
    hadronTree->Branch("RunID"         ,&runID            ,"RunID/I"                         );
    hadronTree->Branch("PDG"           ,&PDG               ,"PDG[buffer_size]/I"              );
    hadronTree->Branch("mix_px"        ,&px               ,"mix_px[buffer_size]/F"           );
    hadronTree->Branch("mix_py"        ,&py               ,"mix_py[buffer_size]/F"           );
    hadronTree->Branch("mix_pz"        ,&pz               ,"mix_pz[buffer_size]/F"           );
    hadronTree->Branch("InvarentMass"  ,&InvarentMass     ,"InvarentMass[buffer_size]/F"     );

	cout << "-----------------------------------------" << endl;
	cout << "------- histograms & tree claimed -------" << endl;
	cout << "-----------------------------------------" << endl;

	return;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::WriteHistograms() {

	///////////////////
	hNRefMult  ->Write();  
	hNRefMultA ->Write();  
	hNRefMultB ->Write();  
	hVertexXY  ->Write();  
	hVertexZ   ->Write();  
	hVertex2D  ->Write();
	hDiffVz    ->Write();
	hcent      ->Write();  
	hcentw     ->Write();  
 
	hcentRefM  ->Write();
	hcentRefW  ->Write();

	hadronTree ->Write();

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

	const int  runID    = mEvent->runId();
	const int  evtID    = mEvent->eventId();
	const int  refMult  = mEvent->refMult();
	const int grefMult  = mEvent->grefMult();
	const int  ranking  = mEvent->ranking();
	const int tofMatch  = mEvent->nBTOFMatch();

	const double magnet = mEvent->bField();

	if(               (!mEvent->isTrigger(610001))
			&&(!mEvent->isTrigger(610011))
			&&(!mEvent->isTrigger(610021))
			&&(!mEvent->isTrigger(610031))
			&&(!mEvent->isTrigger(610041))
			&&(!mEvent->isTrigger(610051))
	  )return kStOK; 

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
	vector<StLambdaDecayPair> KFParticleLambdaDecayPair;
	std::vector<KFParticle> OmegaVec, OmegaSidebandVec, OmegaVecAll;

	SetupKFParticle();
	if (InterfaceCantProcessEvent) return;
	vector<int> PDG;
	vector<float> px,py,pz,InvarentMass;
	PDG.resize(0);px.resize(0);py.resize(0);pz.resize(0);
	int CrefMult,CgrefMult;

	for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++){ 
		const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle]; 
		int upQ; if (particle.GetPDG() == LambdaPdg) upQ = 1; else if (particle.GetPDG() == -1*LambdaPdg) upQ = -1; else continue;
		int eLambda = -(upQ-1)/2; // 0 if Lambda, 1 if AntiLambda

		SetDaughterTrackPointers(iKFParticle);
		if (ProtonTrackIndex == -99999 || PionTrackIndex == -99999) continue; if(!ProtonTrack) continue; if(!PionTrack) continue;

		double dmass = -999; // just a placeholder
		TLorentzVector p4Pair, p4Proton; // just a placeholder
		StLambdaDecayPair TmpLambdaDecayPair(p4Pair, p4Proton, ProtonTrackIndex, PionTrackIndex, (eLambda==0), dmass);
		KFParticleLambdaDecayPair.push_back(TmpLambdaDecayPair);

		if(particle.GetPDG() == LambdaPdg ||
		   particle.GetPDG() == OmegaPdg )// baryons
		{

		/* if (IsOmega && isGoodOmega(cent, particle)) */OmegaVec.push_back(particle);
			CrefMult  = refMult;
			CgrefMult = grefMult;
			PDG.push_back(particle.GetPDG());
			// helix
			TVector3 MomentumOfParticle(particle.GetPx(), particle.GetPy(), particle.GetPz());
			TVector3 PositionOfParticle(particle.GetX(), particle.GetY(), particle.GetZ());
			TLorentzVector OmegaLorentz(MomentumOfParticle, particle.GetE());
			StPicoPhysicalHelix heliPositionOfParticle(MomentumOfParticle, PositionOfParticle, magnet*kilogauss, particle.GetQ());
			double pathlength = heliPositionOfParticle.pathLength(Vertex3D, false);
			TVector3 MomentumOfParticle_tb = heliPositionOfParticle.momentumAt(pathlength, magnet*kilogauss); 
			px.push_back(MomentumOfParticle_tb.X());
			py.push_back(MomentumOfParticle_tb.Y());
			pz.push_back(MomentumOfParticle_tb.Z());
			InvarentMass.push_back(particle.GetE());
		}


	} // End loop over KFParticles
// ======= KFParticle end ======= //

// ======= Lambda loop ======= //

  	Int_t nTracks = mPicoDst->numberOfTracks();
	std::vector<int> track_index;
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) 
	{
    	StPicoTrack *track = mPicoDst->track(iTrack);
    	if (! track)            continue;
    	if (! track->charge())  continue;
    	if (  track->nHitsFit() < 15) continue;
		if (  track->nHitsDedx() < 15) continue;
		if (  track->nHitsFit()*1.0 / track->nHitsMax() < 0.52 || track->nHitsFit()*1.0 / track->nHitsMax() > 1.05) continue;
		//if (  track->dEdxError() < 0.04 || track->dEdxError() > 0.12) continue; // same as kfp
		if (! track->isPrimary()) continue;
		track_index.push_back(iTrack);

		// track info
		float p = track->gMom().Mag();
		float pt = track->gMom().Perp();
		float phi = track->gMom().Phi();
		float eta = track->gMom().Eta();
		float dcatopv = track->gDCA(Vertex3D).Mag();
		float nSigmaKaon = track->nSigmaKaon();
		float nSigmaPion = track->nSigmaPion();
		float nSigmaProton = track->nSigmaProton();
		// if (hTOFEff[8] != 0) hTOFEff_check->Fill(pt, hTOFEff[8]->GetEfficiency(hTOFEff[8]->FindFixBin(pt)));

		float pt_prim = track->pMom().Perp();
		float phi_prim = track->pMom().Phi();
		float eta_prim = track->pMom().Eta();

		// EP
		bool isGoodAsso = true;
		// isGoodAsso &= track->isPrimary();
		isGoodAsso &= pT_asso_lo < pt && pt < pT_asso_hi;
		// if (isGoodAsso)
		// {
		// // fill phi shift for asso
		// 	if (fabs(eta_prim) > TPCAssoEtaCut)
		// {
		// 	for (int i = 1; i <= shift_order_asso; i++) // fill correction for output
		// 	{
		// 		hTPCAssoShiftOutput_cos->Fill(dayPointer, i, cent-1, cos(i*1.0*phi_prim), mWght);
		// 		hTPCAssoShiftOutput_sin->Fill(dayPointer, i, cent-1, sin(i*1.0*phi_prim), mWght);
		// 	}
		// 	// shift asso phi
		// 		float phi_shifted_asso = ShiftAssoPhi(phi_prim, dayPointer+1, cent);
		// 		hTPCAssoPhi->Fill(phi_prim, mWght);
		// 		hTPCAssoPhi_shifted->Fill(phi_shifted_asso, mWght);
		// 		hTPCAssoPhi_2D->Fill(dayPointer, phi_prim, mWght);
		// 		hTPCAssoPhi_2D_shifted->Fill(dayPointer, phi_shifted_asso, mWght);
		// 		if (eta_prim > TPCAssoEtaCut) // construct east EP
		// 	{
		// 		Qx2e += pt*cos(2*phi_shifted_asso);
		// 		Qy2e += pt*sin(2*phi_shifted_asso);
		// 	}
		// 		else if (eta_prim < -1.0*TPCAssoEtaCut) //construct west EP
		// 	{
		// 		Qx2w += pt*cos(2*phi_shifted_asso);
		// 		Qy2w += pt*sin(2*phi_shifted_asso);
		// 	}
		// }
		// }
		bool isGoodPOI = true;
		// isGoodPOI &= track->isPrimary();
		// isGoodPOI &= pT_trig_lo < pt && pt < pT_trig_hi;
		// isGoodPOI &= fabs(eta_prim) < eta_trig_cut;
		// if (isGoodPOI)
		// {
		// // fill phi shift for POI
		// for (int i = 1; i <= shift_order_asso; i++) // fill correction for output
		// {
		// 		hTPCPOIShiftOutput_cos->Fill(dayPointer, i, cent-1, cos(i*1.0*phi_prim), mWght);
		// 		hTPCPOIShiftOutput_sin->Fill(dayPointer, i, cent-1, sin(i*1.0*phi_prim), mWght);
		// }
		// // shift POI phi
		// 	float phi_shifted_POI = ShiftPOIPhi(phi_prim, dayPointer+1, cent);
		// 	hTPCPOIPhi->Fill(phi_prim, mWght);
		// 	hTPCPOIPhi_shifted->Fill(phi_shifted_POI, mWght);
		// 	hTPCPOIPhi_2D->Fill(dayPointer, phi_prim, mWght);
		// 	hTPCPOIPhi_2D_shifted->Fill(dayPointer, phi_shifted_POI, mWght);
		// }

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
			hgbtofYlocal->Fill(BtofYLocal);
			if((tofflag >= 1) && (tof > 0) && (BtofYLocal > -1.8) && (BtofYLocal < 1.8)) hasTOF = true;
		}

		// fill QA
		StPicoPhysicalHelix helix = track->helix(magnet);
		TVector3 pkaon = helix.momentum(magnet*kilogauss);
		// hgpdEdx   ->Fill(pkaon.Mag(), track->dEdx());
		// hgdEdxErr ->Fill(track->dEdxError());
		// hgp       ->Fill(p);
		// hgpT[cent-1]      ->Fill(pt);
		// hgpTeta[cent-1]   ->Fill(pt, eta);
		// if (hasTOF) 
		// {
		// 	hgpT_TOF[cent-1]->Fill(pt);
		// 	hgpTeta_TOF[cent-1]->Fill(pt, eta);
		// }
		// hgDCAtoPV ->Fill(dcatopv);
		
		// int ptbin = static_cast<int>(floor(pt/0.2));
		// double zTPC = TMath::Log(track->dEdx() / 1e6 / StdEdxPull::EvalPred(pkaon.Mag()/KaonPdgMass,1,1)); 
		// hgnSigmaDiff->Fill(zTPC / track->dEdxError() - nSigmaKaon);
		// hgptnSigmaKaon->Fill(pt, nSigmaKaon);
		// hgptnSigmaPion->Fill(pt, nSigmaPion);
		// hgptnSigmaProton->Fill(pt, nSigmaProton);
		
		// if (ptbin >= 0 && ptbin <= 14) hgzTPC_pt[ptbin]->Fill(track->nSigmaKaon());
		if (hasTOF)
		{
			beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
			m2 = pkaon.Mag2()*(1.0 / beta / beta - 1.0);
			// hgpinvbeta->Fill(pkaon.Mag(), 1./beta);
			// hgm2  ->Fill(m2);
			// hgpm2 ->Fill(pkaon.Mag(), m2);
			// hgptm2->Fill(pt, m2);
			// hgm2nSigmaKaon  ->Fill(nSigmaKaon, m2);
			// hgm2nSigmaPion  ->Fill(nSigmaPion, m2);
			// hgm2nSigmaProton->Fill(nSigmaProton, m2);

			// some kaon QA
			//if (track->nSigmaKaon() >  6) hgptm2_largenSigmaKaon->Fill(track->gMom().Perp(), m2);
			//if (track->nSigmaKaon() < -6) hgptm2_smallnSigmaKaon->Fill(track->gMom().Perp(), m2);
			double zTOF_proton = 1/beta - sqrt(ProtonPdgMass*ProtonPdgMass/pkaon.Mag2()+1);
			double zTOF_pion   = 1/beta - sqrt(PionPdgMass*PionPdgMass/pkaon.Mag2()+1);
			double zTOF_kaon   = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() > 0) hgPID2D_proton_pt    [ptbin]->Fill(nSigmaProton, zTOF_proton);
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() < 0) hgPID2D_antiproton_pt[ptbin]->Fill(nSigmaProton, zTOF_proton);
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() > 0) hgPID2D_piplus_pt    [ptbin]->Fill(nSigmaPion  , zTOF_pion  );
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() < 0) hgPID2D_piminus_pt   [ptbin]->Fill(nSigmaPion  , zTOF_pion  );
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() > 0) hgPID2D_kplus_pt     [ptbin]->Fill(nSigmaKaon  , zTOF_kaon  );
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() < 0) hgPID2D_kminus_pt    [ptbin]->Fill(nSigmaKaon  , zTOF_kaon  );
		}

		// primary proton cut for coalescence test
		bool proton_cut = true;
		// if (pt < pT_trig_lo || pt > pT_trig_hi) proton_cut = false; 
		// if (fabs(eta_prim) > eta_trig_cut) proton_cut = false;
		// ProtonPID proton_pid(0., nSigmaProton, pt); // not using zTOF
		// if ((pt > proton_pT_lo && pt < proton_pT_TOFth) && hasTOF) // test efficacy of ProtonPID.h
		// {
		// 	hm2proton_b->Fill(m2);
		// 	if (proton_pid.IsProtonSimple(2., track->charge())) hm2proton_a->Fill(m2);
		// 	if (fabs(nSigmaProton) < 2.)       hm2proton_r->Fill(m2);
		// }
		// if (!hasTOF && pt > proton_pT_TOFth) proton_cut = false;
		// if (pt > proton_pT_TOFth && (m2 > proton_m2_hi || m2 < proton_m2_lo)) proton_cut = false;
		// if (!proton_pid.IsProtonSimple(2., track->charge())) proton_cut = false; // only 0.2 < pt < 2.0!!!
		// // if (fabs(nSigmaProton) > 3) proton_cut = false;
		// if (dcatopv > dcatoPV_hi) proton_cut = false;
		// if (proton_cut)
		// {	
		// 	TLorentzVector lv_proton; lv_proton.SetVectM(track->gMom(), ProtonPdgMass);
		// 	if (track->charge() > 0) {hProtony    ->Fill(lv_proton.Rapidity()); if (fabs(lv_proton.Rapidity()) < 0.5) pct++; }
		// 	else					 {hAntiProtony->Fill(lv_proton.Rapidity()); if (fabs(lv_proton.Rapidity()) < 0.5) pbct++;}

		// 	proton_tracks.push_back(iTrack);
		// }

		// primary pion cut for coalescence test
		bool pion_cut = true;
		// if (pt < pT_trig_lo || pt > pT_trig_hi) pion_cut = false; // use p < 2
		// if (fabs(eta_prim) > eta_trig_cut) pion_cut = false;
		// PionPID pion_pid(0., nSigmaPion, pt); // not using zTOF
		// if ((pt > pion_pT_lo && pt < pion_pT_TOFth) && hasTOF) // test efficacy of ProtonPID.h
		// {
		// 	hm2pion_b->Fill(m2);
		// 	if (pion_pid.IsPionSimple(2., track->charge())) hm2pion_a->Fill(m2);
		// 	if (fabs(nSigmaPion) < 2.)	   hm2pion_r->Fill(m2);				   
		// }
		// if (!hasTOF && pt > pion_pT_TOFth) pion_cut = false;
		// if (pt > pion_pT_TOFth && (m2 > pion_m2_hi || m2 < pion_m2_lo)) pion_cut = false;
		// if (!pion_pid.IsPionSimple(2., track->charge())) pion_cut = false; // only up to pt < 1.8!!!
		// // if (fabs(nSigmaPion) > 3) pion_cut = false;
		// if (dcatopv > dcatoPV_hi) pion_cut = false;
		// if (pion_cut) pion_tracks.push_back(iTrack);

		// primary kaon cut
		/******** looser cut ********/
		bool kaon_cut = true;
		if (pt < pT_trig_lo || pt > pT_trig_hi) kaon_cut = false; // use p < 1.6
		if (fabs(eta_prim) > eta_trig_cut) kaon_cut = false;
		if (!hasTOF && pt > 0.4) kaon_cut = false;
		if (pt > 0.4 && (m2 > 0.34 || m2 < 0.15)) kaon_cut = false;
		double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);

		KaonPID decider(zTOF, nSigmaKaon, pt);
		if (!decider.IsKaonSimple(2., track->charge())) kaon_cut = false;
		if (dcatopv > 2) kaon_cut = false;
		bool isDaughter = false;
		for (int i = 0; i < OmegaVec.size(); i++) if (IsKaonOmegaDaughter(OmegaVec[i], track->id())) isDaughter = true;
		if (isDaughter) kaon_cut = false;
		// if (kaon_cut)
		// {
		// 	if (track->charge() > 0) kpct++; 
		// 	else                     kmct++;
		// 	kaon_tracks.push_back(iTrack);
		// }
		
		/******** stricter cut ********/
		/*
		if (pt < 0.15 || pt> 1.6) continue; // use p < 1.6
		if (!hasTOF) continue;pt
		double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
		KaonPID decider(zTOF, track->nSigmaKaon(), track->gMom().Perp());
		if (!decider.IsKaon()) continue;
		*/

		// kaon QA and correlation
		if (kaon_cut)
		{	
			float TOFEff = 1.0;
			float TOFEff2D = 1.0;
			// float TPCEff = track->charge() > 0? eff_finder.GetEfficiency1D(pt, cent, "Kp") : eff_finder.GetEfficiency1D(pt, cent, "Km");
			// float TPCEff2D = track->charge() > 0? eff_finder.GetEfficiency2D(pt, eta, cent, "Kp") : eff_finder.GetEfficiency2D(pt, eta, cent, "Km");
			// if (hTOFEff[cent-1] != 0 && pt > 0.4) TOFEff = hTOFEff[cent-1]->GetEfficiency(hTOFEff[cent-1]->FindFixBin(pt));
			// if (hTOFEff_2D[cent-1] != 0 && pt > 0.4) TOFEff2D = hTOFEff_2D[cent-1]->GetEfficiency(hTOFEff_2D[cent-1]->FindFixBin(pt, eta));
			// hgKpdEdx    ->Fill(pkaon.Mag(), track->dEdx());
			// hgKp       ->Fill(p);
			// hgKpT      ->Fill(pt);
			// hgKDCAtoPV ->Fill(dcatopv);
			if (hasTOF)
			{
				float beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
				// hgKpinvbeta->Fill(pkaon.Mag(), 1./beta);
				// hgKm2  ->Fill(m2);
				// hgKpm2 ->Fill(pkaon.Mag(), m2);

				// // check kaons misidentified as pions
				// if (m2 > -0.06 && m2 < 0.1) hgKpionpdEdx->Fill(pkaon.Mag(), track->dEdx());
			}
		
			// Omega loop
			if (cent < min_cent || cent > max_cent) continue;
			const int kaonindex = track->id();
			for (int iOmega=0; iOmega < OmegaVec.size(); iOmega++)
			{ 
				const KFParticle particle = OmegaVec[iOmega]; 
				// if (IsKaonOmegaDaughter(particle, kaonindex)) continue;

				// couting Omega
				// if (!fill_nomega)
				// {
				// 	if (particle.GetPDG() > 0) hOmegaUsed->Fill(0.);
				// 	if (particle.GetPDG() < 0) hOmegaUsed->Fill(1.);
				// 	if ( hasLambdabar && particle.GetPDG() > 0) hOmegaUsed_wlb->Fill(0);
				// 	if ( hasLambdabar && particle.GetPDG() < 0) hOmegaUsed_wlb->Fill(1);
				// 	if (!hasLambdabar && particle.GetPDG() > 0) hOmegaUsed_wolb->Fill(0);
				// 	if (!hasLambdabar && particle.GetPDG() < 0) hOmegaUsed_wolb->Fill(1);
				// 	if ( hasLambda && particle.GetPDG() > 0) hOmegaUsed_wl->Fill(0);
				// 	if ( hasLambda && particle.GetPDG() < 0) hOmegaUsed_wl->Fill(1);
				// 	if (!hasLambda && particle.GetPDG() > 0) hOmegaUsed_wol->Fill(0);
				// 	if (!hasLambda && particle.GetPDG() < 0) hOmegaUsed_wol->Fill(1);
				// }

				// pair-wise should be added after this line
				/* */

				// Omega momentum at DCA to PV
				// TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				// TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				// StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());
				// double pathlength = helixOmega.pathLength(Vertex3D, false);
				// TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet*kilogauss); 

				// k*
				// TLorentzVector lv1; lv1.SetVectM(pOmega_tb, OmegaPdgMass);
				TLorentzVector lv2; lv2.SetVectM(track->gMom(), KaonPdgMass);
				// TLorentzVector lv2_ori; lv2_ori.SetVectM(track->gMom(), KaonPdgMass);
				// double dpt = fabs(lv1.Perp()-lv2.Perp());
				// double dy  = fabs(lv1.Rapidity() - lv2.Rapidity());
				// double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
				// TLorentzVector P = lv1 + lv2;
				// TVector3 pair_beta = P.BoostVector();
				// lv1.Boost((-1)*pair_beta); 	
				// lv2.Boost((-1)*pair_beta); 	
				// double kstar = 0.5*(lv1-lv2).Vect().Mag();

				// pT weight if anti-omega
				// float weight = 1.;
				// if (PtReweighting) weight = GetPtWeight(particle);
				
				// float eff = TPCEff * TOFEff;
				// float eff2D = TPCEff2D * TOFEff2D;
				// if (track->charge() > 0 && particle.GetQ() < 0) 
				// {
				// 	hCorrKplusO   ->Fill(kstar, 1./eff);
				// 	hPtCorrKplusO ->Fill(dpt, 1./eff);
				// 	hyCorrKplusO  ->Fill(dy, 1./eff);
				// 	hphiCorrKplusO->Fill(dphi), 1./eff;

				// 	hKpluspt_omega ->Fill(lv2_ori.Perp(), 1./eff);
				// 	hKpluseta_omega->Fill(lv2_ori.Eta()), 1./eff;
				// 	hKplusphi_omega->Fill(lv2_ori.Phi()), 1./eff;
				// 	hKplusy_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	// hCorrKplusO_y_pT  ->Fill(dpt, dy);
				// 	// hCorrKplusO_y_phi ->Fill(dphi, dy);
				// 	// hCorrKplusO_phi_pT->Fill(dpt, dphi);

				// 	if (hasLambdabar)
				// 	{
				// 		hKpluspt_wlb_omega ->Fill(lv2_ori.Perp(), 1./eff);
				// 		hKplusy_wlb_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	}
				// 	else
				// 	{
				// 		hKpluspt_wolb_omega ->Fill(lv2_ori.Perp(), 1./eff);
				// 		hKplusy_wolb_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	}
				// }
				// if (track->charge() > 0 && particle.GetQ() > 0) 
				// {
				// 	hCorrKplusObar   ->Fill(kstar, weight / eff);
				// 	hPtCorrKplusObar ->Fill(dpt, weight / eff);
				// 	hyCorrKplusObar  ->Fill(dy, weight / eff);
				// 	hphiCorrKplusObar->Fill(dphi, weight/ eff);

				// 	hKpluspt_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
				// 	hKpluseta_omegabar->Fill(lv2_ori.Eta(), 1./eff);
				// 	hKplusphi_omegabar->Fill(lv2_ori.Phi(), 1./eff);
				// 	hKplusy_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	// hCorrKplusObar_y_pT  ->Fill(dpt, dy, weight);
				// 	// hCorrKplusObar_y_phi ->Fill(dphi, dy, weight);
				// 	// hCorrKplusObar_phi_pT->Fill(dpt, dphi, weight);

				// 	if (hasLambdabar)
				// 	{
				// 		hKpluspt_wlb_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
				// 		hKplusy_wlb_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	}
				// 	else
				// 	{
				// 		hKpluspt_wolb_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
				// 		hKplusy_wolb_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	}
				// }
				// if (track->charge() < 0 && particle.GetQ() < 0)
				// {
				// 	hCorrKminusO   ->Fill(kstar, 1./eff);
				// 	hPtCorrKminusO ->Fill(dpt, 1./eff);
				// 	hyCorrKminusO  ->Fill(dy, 1./eff);
				// 	hphiCorrKminusO->Fill(dphi, 1./eff);
				// 	if (dpt < 0.5) hNegPtDiff_dphi_KmO->Fill(dphi);
				// 	if (dpt > 1.0) hPosPtDiff_dphi_KmO->Fill(dphi);
					
				// 	hKminuspt_omega ->Fill(lv2_ori.Perp(), 1./eff);
				// 	hKminuseta_omega->Fill(lv2_ori.Eta(), 1./eff);
				// 	hKminusphi_omega->Fill(lv2_ori.Phi(), 1./eff);
				// 	hKminusy_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	// hCorrKminusO_y_pT  ->Fill(dpt, dy);
				// 	// hCorrKminusO_y_phi ->Fill(dphi, dy);
				// 	// hCorrKminusO_phi_pT->Fill(dpt, dphi);

				// 	if (hasLambda)
				// 	{
				// 		hKminuspt_wl_omega ->Fill(lv2_ori.Perp(), 1./eff);
				// 		hKminusy_wl_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	}
				// 	else
				// 	{
				// 		hKminuspt_wol_omega ->Fill(lv2_ori.Perp(), 1./eff);
				// 		hKminusy_wol_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	}
				// }
				// if (track->charge() < 0 && particle.GetQ() > 0) 
				// {
				// 	hCorrKminusObar   ->Fill(kstar, weight / eff);
				// 	hPtCorrKminusObar ->Fill(dpt, weight / eff);
				// 	hyCorrKminusObar  ->Fill(dy, weight / eff);
				// 	hphiCorrKminusObar->Fill(dphi, weight / eff);
				// 	if (dpt < 0.5) hNegPtDiff_dphi_KmOb->Fill(dphi, weight);
				// 	if (dpt > 1.0) hPosPtDiff_dphi_KmOb->Fill(dphi, weight);

				// 	hKminuspt_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
				// 	hKminuseta_omegabar->Fill(lv2_ori.Eta(), 1./eff);
				// 	hKminusphi_omegabar->Fill(lv2_ori.Phi(), 1./eff);
				// 	hKminusy_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	// hCorrKminusObar_y_pT  ->Fill(dpt, dy, weight);
				// 	// hCorrKminusObar_y_phi ->Fill(dphi, dy, weight);
				// 	// hCorrKminusObar_phi_pT->Fill(dpt, dphi, weight);

				// 	if (hasLambda)
				// 	{
				// 		hKminuspt_wl_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
				// 		hKminusy_wl_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	}
				// 	else
				// 	{
				// 		hKminuspt_wol_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
				// 		hKminusy_wol_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
				// 	}
				// }
			} // End loop over regular Omega
			// fill_nomega = true;

			// Omega sideband loop
			// for (int iOmega=0; iOmega < OmegaSidebandVec.size(); iOmega++)
			// { 
			// 	const KFParticle particle = OmegaSidebandVec[iOmega]; 
			// 	if (IsKaonOmegaDaughter(particle, kaonindex)) continue;

			// 	// couting Omega
			// 	if (!fill_nomega_sideband)
			// 	{
			// 		if (particle.GetPDG() > 0) hOmegaUsed_sideband->Fill(0.);
			// 		if (particle.GetPDG() < 0) hOmegaUsed_sideband->Fill(1.);
			// 		if ( hasLambdabar && particle.GetPDG() > 0) hOmegaUsed_wlb_sideband->Fill(0.);
			// 		if ( hasLambdabar && particle.GetPDG() < 0) hOmegaUsed_wlb_sideband->Fill(1.);
			// 		if (!hasLambdabar && particle.GetPDG() > 0) hOmegaUsed_wolb_sideband->Fill(0.);
			// 		if (!hasLambdabar && particle.GetPDG() < 0) hOmegaUsed_wolb_sideband->Fill(1.);
			// 		if ( hasLambda && particle.GetPDG() > 0) hOmegaUsed_wl_sideband->Fill(0.);
			// 		if ( hasLambda && particle.GetPDG() < 0) hOmegaUsed_wl_sideband->Fill(1.);
			// 		if (!hasLambda && particle.GetPDG() > 0) hOmegaUsed_wol_sideband->Fill(0.);
			// 		if (!hasLambda && particle.GetPDG() < 0) hOmegaUsed_wol_sideband->Fill(1.);
			// 	}

			// 	// pair-wise should be added after this line
			// 	/* */

			// 	// Omega momentum at DCA to PV
			// 	TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
			// 	TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
			// 	StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());
			// 	double pathlength = helixOmega.pathLength(Vertex3D, false);
			// 	TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet*kilogauss); 

			// 	// k*
			// 	TLorentzVector lv1; lv1.SetVectM(pOmega_tb, OmegaPdgMass);
			// 	TLorentzVector lv2; lv2.SetVectM(track->gMom(), KaonPdgMass);
			// 	TLorentzVector lv2_ori; lv2_ori.SetVectM(track->gMom(), KaonPdgMass);
			// 	double dpt = fabs(lv1.Perp()-lv2.Perp());
			// 	double dy  = fabs(lv1.Rapidity() - lv2.Rapidity());
			// 	double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
			// 	TLorentzVector P = lv1 + lv2;
			// 	TVector3 pair_beta = P.BoostVector();
			// 	lv1.Boost((-1)*pair_beta); 	
			// 	lv2.Boost((-1)*pair_beta); 	
			// 	double kstar = 0.5*(lv1-lv2).Vect().Mag();

			// 	// pT weight if anti-omega
			// 	float weight = 1.;
			// 	if (PtReweighting) weight = GetPtWeight(particle);
				
			// 	float eff = TPCEff * TOFEff;
			// 	float eff2D = TPCEff2D * TOFEff2D;
			// 	if (track->charge() > 0 && particle.GetQ() < 0) 
			// 	{
			// 		hCorrKplusO_sideband   ->Fill(kstar, 1./eff);
			// 		hPtCorrKplusO_sideband ->Fill(dpt, 1./eff);
			// 		hyCorrKplusO_sideband  ->Fill(dy, 1./eff);
			// 		hphiCorrKplusO_sideband->Fill(dphi, 1./eff);

			// 		hKpluspt_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 		hKpluseta_omega_sideband->Fill(lv2_ori.Eta(), 1./eff);
			// 		hKplusphi_omega_sideband->Fill(lv2_ori.Phi(), 1./eff);
			// 		hKplusy_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		//hCorrKplusO_y_pT  ->Fill(dpt, dy);
			// 		//hCorrKplusO_y_phi ->Fill(dphi, dy);
			// 		//hCorrKplusO_phi_pT->Fill(dpt, dphi);

			// 		if (hasLambdabar)
			// 		{
			// 			hKpluspt_wlb_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 			hKplusy_wlb_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		}
			// 		else
			// 		{
			// 			hKpluspt_wolb_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 			hKplusy_wolb_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		}
			// 	}
			// 	if (track->charge() > 0 && particle.GetQ() > 0) 
			// 	{
			// 		hCorrKplusObar_sideband   ->Fill(kstar, weight/eff);
			// 		hPtCorrKplusObar_sideband ->Fill(dpt, weight/eff);
			// 		hyCorrKplusObar_sideband  ->Fill(dy, weight/eff);
			// 		hphiCorrKplusObar_sideband->Fill(dphi, weight/eff);

			// 		hKpluspt_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 		hKpluseta_omegabar_sideband->Fill(lv2_ori.Eta(), 1./eff);
			// 		hKplusphi_omegabar_sideband->Fill(lv2_ori.Phi(), 1./eff);
			// 		hKplusy_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		//hCorrKplusObar_y_pT  ->Fill(dpt, dy);
			// 		//hCorrKplusObar_y_phi ->Fill(dphi, dy);
			// 		//hCorrKplusObar_phi_pT->Fill(dpt, dphi);

			// 		if (hasLambdabar)
			// 		{
			// 			hKpluspt_wlb_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 			hKplusy_wlb_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		}
			// 		else
			// 		{
			// 			hKpluspt_wolb_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 			hKplusy_wolb_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		}
			// 	}
			// 	if (track->charge() < 0 && particle.GetQ() < 0)
			// 	{
			// 		hCorrKminusO_sideband   ->Fill(kstar, 1./eff);
			// 		hPtCorrKminusO_sideband ->Fill(dpt, 1./eff);
			// 		hyCorrKminusO_sideband  ->Fill(dy, 1./eff);
			// 		hphiCorrKminusO_sideband->Fill(dphi, 1./eff);

			// 		hKminuspt_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 		hKminuseta_omega_sideband->Fill(lv2_ori.Eta(), 1./eff);
			// 		hKminusphi_omega_sideband->Fill(lv2_ori.Phi(), 1./eff);
			// 		hKminusy_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		//if (dpt < 0.5) hNegPtDiff_dphi_KmO->Fill(dphi);
			// 		//if (dpt > 1.0) hPosPtDiff_dphi_KmO->Fill(dphi);

			// 		//hCorrKminusO_y_pT  ->Fill(dpt, dy);
			// 		//hCorrKminusO_y_phi ->Fill(dphi, dy);
			// 		//hCorrKminusO_phi_pT->Fill(dpt, dphi);

			// 		if (hasLambda)
			// 		{
			// 			hKminuspt_wl_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 			hKminusy_wl_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		}
			// 		else
			// 		{
			// 			hKminuspt_wol_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 			hKminusy_wol_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		}

			// 	}
			// 	if (track->charge() < 0 && particle.GetQ() > 0) 
			// 	{
			// 		hCorrKminusObar_sideband   ->Fill(kstar, weight/eff);
			// 		hPtCorrKminusObar_sideband ->Fill(dpt, weight/eff);
			// 		hyCorrKminusObar_sideband  ->Fill(dy, weight/eff);
			// 		hphiCorrKminusObar_sideband->Fill(dphi, weight/eff);

			// 		hKminuspt_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 		hKminuseta_omegabar_sideband->Fill(lv2_ori.Eta(), 1./eff);
			// 		hKminusphi_omegabar_sideband->Fill(lv2_ori.Phi(), 1./eff);
			// 		hKminusy_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		//if (dpt < 0.5) hNegPtDiff_dphi_KmOb->Fill(dphi);
			// 		//if (dpt > 1.0) hPosPtDiff_dphi_KmOb->Fill(dphi);

			// 		//hCorrKminusObar_y_pT  ->Fill(dpt, dy);
			// 		//hCorrKminusObar_y_phi ->Fill(dphi, dy);
			// 		//hCorrKminusObar_phi_pT->Fill(dpt, dphi);

			// 		if (hasLambda)
			// 		{
			// 			hKminuspt_wl_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 			hKminusy_wl_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		}
			// 		else
			// 		{
			// 			hKminuspt_wol_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
			// 			hKminusy_wol_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
			// 		}
			// 	}
			// } // End loop over sideband Omega
			// fill_nomega_sideband = true;
		}
		
	}
	hadronTree->Fill();

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

	if ( !KFParticleInterface->ProcessEvent(PicoDst, triggeredTracks) ) InterfaceCantProcessEvent = true; else InterfaceCantProcessEvent = false;

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
	ProtonTrackIndex = -99999, PionTrackIndex = -99999;
	for(int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++){ 
		const int daughterId = particle.DaughterIds()[iDaughter]; 
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
		const int globalTrackId = daughter.DaughterIds()[0];
		int trackIndex = trackMap[globalTrackId];    
		if (daughter.GetPDG() == upQ*ProtonPdg) ProtonTrackIndex = trackIndex;
		else if (daughter.GetPDG() == upQ*PionPdg) PionTrackIndex = trackIndex;
	}  // iDaughter
	ProtonTrack = PicoDst->track(ProtonTrackIndex); PionTrack = PicoDst->track(PionTrackIndex);	
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