#include "StV0Maker.h"
#include "StDcaService.h"
#include "StarClassLibrary/SystemOfUnits.h"

#include <iostream>
#include <utility>
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoETofHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoArrays.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoETofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh" 
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StPileupUtil/StPileupUtil.h"
#include "StHbtMaker/Infrastructure/StHbtParticle.hh"
#include "StHbtMaker/Infrastructure/StHbtTrack.hh"
#include "StEventMaker/StEventMaker.h"
#include "StEvent/StEvent.h"
#include "StEvent/StPrimaryVertex.h"

#include "tables/St_vertexSeed_Table.h"
#include "StBTofHeader.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "MeanShift.h"

ClassImp(StV0Maker)                   // Macro for CINT compatibility

StV0Maker::StV0Maker( StPicoDstMaker* maker, const char * name) : StMaker(name)
{ // Initialize and/or zero all public/private data members here.

	//for ( Int_t i = 0 ; i < kMaxNumberOfTH1F ; i++ )  // Zero the histogram pointers, not necessary. it is NULL naturaly.
	//  {
	//    histogram[i] = NULL ;
	//  }

	mPicoDstMaker      = maker ;                    // Pass MuDst pointer to DstAnlysisMaker Class member functions
	mV0Type = kLambda;	//Lambda as default!

	mRotate = false;
	mSameSignPlus = false;
	mSameSignMinus = false;

	mDcaAlgoLong = true;

	mDumpNull = false;

	histogram_output = NULL  ;                    // Zero the Pointer to histogram output file
	v0tree_output = NULL  ;                    // Zero the Pointer to v0 tree output file
	mHistogramOutputFileName = "" ;               // Histogram Output File Name will be set inside the "analysis".C macro
	mV0TreeOutputFileName = "" ;               // V0 Output File Name will be set inside the "analysis".C macro

	mV0Tree = NULL ;

	mEventsProcessed = 0     ;                    // Zero the Number of Events processed by the maker 
	mTestNTrack = 0;
	mTestVZ = 0;

	//mBeamHelix = NULL;

}

StV0Maker::~StV0Maker() 
{ // Destroy and/or zero out all public/private data members here.
}

void StV0Maker::initConst(){
	//initialize the constant for different V0 types.

	if(mV0Type == kLambda || mV0Type == kAntiLambda){
		// for Lambda and AntiLambda
		mMass1      = 0.93827; // mass of proton
		mMass2      = 0.13957; // mass of pion
		mMassV0     = 1.115684;// mass of Lambda

		mMassBachelor = 0.13957; //mass of pion
		mMassXi = 1.32131;//mass of xi

		mCharge1	= 1;
		mCharge2	= -1;
		mChargeBachelor = -1;
		mChargeXi	= -1;

		mCharge1_anti     = -1;
		mCharge2_anti     = 1;
		mChargeBachelor_anti = 1;
		mChargeXi_anti    = 1;

		//parameters for StDcaService.cxx
		kShiftConnect = 0.3;
		kShiftContain = 0.3;
	}

	return;
}

void StV0Maker::initParam(){
	//setup the cut values here. do not hard-code them in ::Make()

	//cutAbsVertexZLeEq  = 8.;
	//cutTriggerIdEq  = 340001;	//auau Run11 19.6GeV

	if(mV0Type == kLambda || mV0Type == kAntiLambda){
		cutNHitsGr = 10;
		cutPtGrEq = 0.05;///0.05

		cutAbsNSigma1Le = 3.2;
		cutAbsNSigma2Le = 3.2;
		//Lambda
		cutDca1GrEq_Lambda  = 0.4; //0.3
		cutDca2GrEq_Lambda  =1.3; //1.0
		cutDca1to2LeEq_Lambda = 1.0; //1.0
		cutV0MassWidthLeEq_Lambda = 0.08;
		cutV0rdotpGr_Lambda  = 0.0;
		cutDcaV0Le_Lambda    = 5.0; //5.0
		cutV0DecLenGrEq_Lambda =3.0;//3.0

		//Xi
		cutDca1GrEq_Xi  = 1.5; 
		cutV0DcaGrEq_Xi = 0.5; 
		cutV0MassWidthLeEq_Xi = 0.006;//0.02
		cutDca1to2LeEq_Xi = 0.8; //0.8
		cutXiMassWidthLeEq_Xi = 0.10;
		cutXirdotpGr_Xi  = 0;
		cutXircrosspLeEq_Xi = 0.5; //0.2
		cutDcaXiLe_Xi    = 0.6;
		cutXiDecLenGrEq_Xi =3.8; //2.0
	}

	return;
}

void StV0Maker::initHisto()
{
	// Create Histograms
	// there is no better way to set QA histograms. do not use histogram arrays or vectors.
	// it is not useful. there are no need to operate all the histograms at the same time.

	const Int_t    nbins    =  100   ;

	//QA for events
	hNPrimVertex  = new TH1F( "PrimVertex", "Number of Primary Vertex", 10, 0.0, 10.0 ) ;
	hVertexZ  = new TH1F( "VertexZ", "Event Vertex Z Position", nbins*4, 190.0, 210.0 ) ; 
	hSelectVertexZ  = new TH1F( "SelectVertexZ", "Event Select Vertex Z Position", nbins*4, 190.0, 210.0 ) ; 
	hvxvy = new TH2F("hvxvy","",200,-10,10,200,-10,10);
	hSelectvxvy = new TH2F("Select_hvxvy","",200,-10,10,200,-10,10);
	hNRefMult  = new TH1F( "RefMult", "Reference Multiplicity", 1000, 0.0, 1000.0 ) ;
	hSelectNRefMult  = new TH1F( "SelectRefMult", "Reference Multiplicity of selected events", 1000, 0.0, 1000.0 ) ;
	hNGRefMult  = new TH1F( "GRefMult", "GReference Multiplicity", 1000, 0.0, 1000.0 ) ;
	hSelectNGRefMult  = new TH1F( "SelectGRefMult", "GReference Multiplicity of selected events", 1000, 0.0, 1000.0 ) ;
	hCent = new TH1F("hCent","hcent",10,-0.5,9.5);
	hTofPID = new TH2F("htofpid","Tof pid;q*p;m2",300,-3,3,300,-0.5,2);
	hVpdVz = new TH2F("VpdVz",";TPC Vz [cm];VPD Vz [cm]",1000,-200,200,1000,-200,200);
	hSelectVpdVz = new TH2F("Select_VpdVz",";TPC Vz [cm];VPD Vz [cm]",1000,-200,200,1000,-200,200);

	hRefTof = new TH2F("Ref_Tof" ," ;RefMult ; TofMult",1000,0,500,1000,0,400);

	hSelectRefTof = new TH2F("Select_Ref_Tof"," ;RefMult ; TofMult",1000,0,400,1000,0,400);
	hSelectGrefTofMatch = new TH2F("Select_Gref_TofMatch"," ;GrefMult ; TofMatch",1000,0,400,1000,0,400);
	//QA for global tracks
	hPtRaw  = new TH1F( "PtRaw", "Transverse Momentum for all particles", nbins*6, 0.0, 30.0 ) ;
	hEtaRaw  = new TH1F( "EtaRaw", "Eta for all particles", nbins*5, -2, 2 ) ;
	hPhiRaw  = new TH1F( "PhiRaw", "Phi for all particles", nbins*10, -TMath::Pi(), TMath::Pi() ) ;
	hPt  = new TH1F( "Pt", "Transverse Momentum for selected particles", nbins*6, 0.0, 30.0 ) ;
	hEta  = new TH1F( "Eta", "Eta for selected particles", nbins*5, -2, 2 ) ;
	hPhi  = new TH1F( "Phi", "Phi for selected particles", nbins*10, -TMath::Pi(), TMath::Pi() ) ;
	hPhiLowPt  = new TH1F( "PhiLowPt", "Phi for selected particles", nbins*10, -TMath::Pi(), TMath::Pi() ) ;
	hPhiHighPt  = new TH1F( "PhiHighPt", "Phi for selected particles", nbins*10, -TMath::Pi(), TMath::Pi() ) ;
	hDedxP  = new TH2F( "DedxP", "dEdx for selected particles;q*p;dEdx", nbins*6, -3, 3, 1000, 0, 50) ;
	honeoverBetaP  = new TH2F( "oneoverBeta", ";q*p;1/Beta", nbins*6, -3, 3, 1000, 0, 2) ;
	hNSigmaPion  = new TH1F( "nSigmaPion", "nSigmaPion for selected particles", nbins*2, -10, 10 ) ;
	hNSigmaProton  = new TH1F( "nSigmaProton", "nSigmaProton for selected particles", nbins*2, -10, 10 ) ;
	hNSigmaKaon  = new TH1F( "nSigmaKaon", "nSigmaKaon for selected particles", nbins*2, -10, 10 ) ;
	hNHitsFit  = new TH1F( "nHitsFit", "nHitsFit for all particles", 80, 0, 80 ) ;
	hNHits  = new TH1F( "nHits", "nHits for all particles", 80, 0, 80 ) ;
	//QA for V0's
	pionplusmass2  = new TH1F( "pionplusmass2", "pionplusmass2", 400, -0.3,1.1 ) ;
	pionminusmass2  = new TH1F( "pionminusmass2", "pionminusmass2", 400, -0.3,1.1 ) ;


	return;
}

void StV0Maker::initTree()
{
	//initialize the TTree for StV0Dst
	mV0Tree = new TTree("V0PicoDst","V0PicoDst from StV0Maker");

	mV0Tree->SetDirectory(v0tree_output);

	mV0Tree->Branch("runnumber",&mV0Dst.runnumber,"runnumber/I");
	mV0Tree->Branch("evtnumber",&mV0Dst.evtnumber,"evtnumber/I");
	mV0Tree->Branch("trgmode",&mV0Dst.trgmode,"trgmode/I");
	mV0Tree->Branch("cent",&mV0Dst.cent,"cent/I");
	mV0Tree->Branch("nrefmult",&mV0Dst.nrefmult,"nrefmult/I");
	mV0Tree->Branch("ntofmult",&mV0Dst.ntofmult,"ntofmult/I");
	mV0Tree->Branch("ntofmatch",&mV0Dst.ntofmatch,"ntofmatch/I");
	mV0Tree->Branch("grefmult",&mV0Dst.grefmult,"grefmult/I");
	mV0Tree->Branch("primvertexX",&mV0Dst.primvertexX,"primvertexX/F");
	mV0Tree->Branch("primvertexY",&mV0Dst.primvertexY,"primvertexY/F");
	mV0Tree->Branch("primvertexZ",&mV0Dst.primvertexZ,"primvertexZ/F");
	mV0Tree->Branch("magn",&mV0Dst.magn,"magn/F");
	mV0Tree->Branch("nPionplus",&mV0Dst.nPionplus,"nPionplus/I");
	mV0Tree->Branch("nPionminus",&mV0Dst.nPionminus,"nPionminus/I");

	//Pionplus
	mV0Tree->Branch("id_Pionplus",mV0Dst.id_Pionplus,"id_Pionplus[nPionplus]/I");
	mV0Tree->Branch("nhits_Pionplus",mV0Dst.nhits_Pionplus,"nhits_Pionplus[nPionplus]/I");
	mV0Tree->Branch("nhitsFit_Pionplus",mV0Dst.nhitsFit_Pionplus,"nhitsFit_Pionplus[nPionplus]/I");
	mV0Tree->Branch("dedx_Pionplus",mV0Dst.dedx_Pionplus,"dedx_Pionplus[nPionplus]/F");
	mV0Tree->Branch("nsigmapion_Pionplus",mV0Dst.nsigmapion_Pionplus,"nsigmapion_Pionplus[nPionplus]/F");
	mV0Tree->Branch("nsigmaproton_Pionplus",mV0Dst.nsigmaproton_Pionplus,"nsigmaproton_Pionplus[nPionplus]/F");
	mV0Tree->Branch("nsigmaelectron_Pionplus",mV0Dst.nsigmaelectron_Pionplus,"nsigmaelectron_Pionplus[nPionplus]/F");
	mV0Tree->Branch("nsigmakaon_Pionplus",mV0Dst.nsigmakaon_Pionplus,"nsigmakaon_Pionplus[nPionplus]/F");
	mV0Tree->Branch("px_Pionplus",mV0Dst.px_Pionplus,"px_Pionplus[nPionplus]/F");
	mV0Tree->Branch("py_Pionplus",mV0Dst.py_Pionplus,"py_Pionplus[nPionplus]/F");
	mV0Tree->Branch("pz_Pionplus",mV0Dst.pz_Pionplus,"pz_Pionplus[nPionplus]/F");
	mV0Tree->Branch("dca_Pionplus",mV0Dst.dca_Pionplus,"dca_Pionplus[nPionplus]/F");
	mV0Tree->Branch("tofflag_Pionplus",mV0Dst.tofflag_Pionplus,"tofflag_Pionplus[nPionplus]/I");
	mV0Tree->Branch("tof_Pionplus",mV0Dst.tof_Pionplus,"tof_Pionplus[nPionplus]/F");
	mV0Tree->Branch("tofpathlen_Pionplus",mV0Dst.tofpathlen_Pionplus,"tofpathlen_Pionplus[nPionplus]/F");
	mV0Tree->Branch("mass2pion_Pionplus",mV0Dst.mass2pion_Pionplus,"mass2pion_Pionplus[nPionplus]/F");
	mV0Tree->Branch("Betapion_Pionplus",mV0Dst.Betapion_Pionplus,"Betapion_Pionplus[nPionplus]/F");
	mV0Tree->Branch("emass2pion_Pionplus",mV0Dst.emass2pion_Pionplus,"emass2pion_Pionplus[nPionplus]/F");
	mV0Tree->Branch("eBetapion_Pionplus",mV0Dst.eBetapion_Pionplus,"eBetapion_Pionplus[nPionplus]/F");
	mV0Tree->Branch("PadRow1_Pionplus",mV0Dst.PadRow1_Pionplus,"PadRow1_Pionplus[nPionplus]/I");
	mV0Tree->Branch("PadRow2_Pionplus",mV0Dst.PadRow2_Pionplus,"PadRow2_Pionplus[nPionplus]/I");
	mV0Tree->Branch("IPadRow_Pionplus",mV0Dst.IPadRow_Pionplus,"IPadRow_Pionplus[nPionplus]/l");

	//Pionminus
	mV0Tree->Branch("id_Pionminus",mV0Dst.id_Pionminus,"id_Pionminus[nPionminus]/I");
	mV0Tree->Branch("nhits_Pionminus",mV0Dst.nhits_Pionminus,"nhits_Pionminus[nPionminus]/I");
	mV0Tree->Branch("nhitsFit_Pionminus",mV0Dst.nhitsFit_Pionminus,"nhitsFit_Pionminus[nPionminus]/I");
	mV0Tree->Branch("dedx_Pionminus",mV0Dst.dedx_Pionminus,"dedx_Pionminus[nPionminus]/F");
	mV0Tree->Branch("nsigmapion_Pionminus",mV0Dst.nsigmapion_Pionminus,"nsigmapion_Pionminus[nPionminus]/F");
	mV0Tree->Branch("nsigmaproton_Pionminus",mV0Dst.nsigmaproton_Pionminus,"nsigmaproton_Pionminus[nPionminus]/F");
	mV0Tree->Branch("nsigmaelectron_Pionminus",mV0Dst.nsigmaelectron_Pionminus,"nsigmaelectron_Pionminus[nPionminus]/F");
	mV0Tree->Branch("nsigmakaon_Pionminus",mV0Dst.nsigmakaon_Pionminus,"nsigmakaon_Pionminus[nPionminus]/F");
	mV0Tree->Branch("px_Pionminus",mV0Dst.px_Pionminus,"px_Pionminus[nPionminus]/F");
	mV0Tree->Branch("py_Pionminus",mV0Dst.py_Pionminus,"py_Pionminus[nPionminus]/F");
	mV0Tree->Branch("pz_Pionminus",mV0Dst.pz_Pionminus,"pz_Pionminus[nPionminus]/F");
	mV0Tree->Branch("dca_Pionminus",mV0Dst.dca_Pionminus,"dca_Pionminus[nPionminus]/F");
	mV0Tree->Branch("tofflag_Pionminus",mV0Dst.tofflag_Pionminus,"tofflag_Pionminus[nPionminus]/I");
	mV0Tree->Branch("tof_Pionminus",mV0Dst.tof_Pionminus,"tof_Pionminus[nPionminus]/F");
	mV0Tree->Branch("tofpathlen_Pionminus",mV0Dst.tofpathlen_Pionminus,"tofpathlen_Pionminus[nPionminus]/F");
	//mV0Tree->Branch("pathlen_Pionminus",mV0Dst.pathlen_Pionminus,"pathlen_Pionminus[nPionminus]/F");
	mV0Tree->Branch("mass2pion_Pionminus",mV0Dst.mass2pion_Pionminus,"mass2pion_Pionminus[nPionminus]/F"); 
	mV0Tree->Branch("Betapion_Pionminus",mV0Dst.Betapion_Pionminus,"Betapion_Pionminus[nPionminus]/F"); 
	mV0Tree->Branch("emass2pion_Pionminus",mV0Dst.emass2pion_Pionminus,"emass2pion_Pionminus[nPionminus]/F");
	mV0Tree->Branch("eBetapion_Pionminus",mV0Dst.eBetapion_Pionminus,"eBetapion_Pionminus[nPionminus]/F");
	mV0Tree->Branch("PadRow1_Pionminus",mV0Dst.PadRow1_Pionminus,"PadRow1_Pionminus[nPionminus]/I");
	mV0Tree->Branch("PadRow2_Pionminus",mV0Dst.PadRow2_Pionminus,"PadRow2_Pionminus[nPionminus]/I");
	mV0Tree->Branch("IPadRow_Pionminus",mV0Dst.IPadRow_Pionminus,"IPadRow_Pionminus[nPionminus]/l");

	return;
}

Int_t StV0Maker::Init( )
{
	// setup the constants according to mV0Type
	initConst();

	// initialize parameters (cuts)
	initParam();

	// Create Histogram output file
	if(mHistogramOutputFileName == "") { 
		//CAUTION: ALWAYS USE { } HERE!!! LOG_XXX is a if()xxx macro!!!
		LOG_ERROR << "StV0Maker: Please specify the histrogram output file" <<endm;
		exit(-1);
	}
	else {
		histogram_output = new TFile( mHistogramOutputFileName, "recreate" ) ;  
	}
	// Book histograms
	initHisto();

	// Create V0 Tree output file
	if(mV0TreeOutputFileName == "") {
		LOG_WARN << "StV0Maker: The V0 tree output file is not specified! output is smeared!" <<endm;
	}
	else {
		v0tree_output = new TFile( mV0TreeOutputFileName, "recreate" ) ;
		// Create V0 Tree
		initTree();
	}
	return kStOK ; 
}


Int_t StV0Maker::Make( )
{ // Do each event

	//if(GetDebug()) LOG_QA<<"in StV0Maker::Make"<<endm;
	// Do some cleaning here, used for StXiMaker or other subsequent makers
	//mPassEventCut = false;

	// Get 'event' data 
	StPicoEvent* muEvent      =  mPicoDstMaker->picoDst()->event() ;
	if (!muEvent) return kStOK;

	//Remove badrun                                            
	//StRefMultCorr* refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;
	//refmultCorrUtil->init(muEvent->runId());                   
	//if ( refmultCorrUtil->isBadRun(muEvent->runId()) ) return kStOK;                                    

	hNPrimVertex -> Fill( 1 );
	//Select trigger
	if(!muEvent->isTrigger(730000))return kStOK;//run11_mb
	StRefMultCorr *mRefMultCorrUtil = CentralityMaker::instance()->getRefMultCorrFxt();
	mRefMultCorrUtil->init(muEvent->runId());
	//remove baddrun
	if ( mRefMultCorrUtil->isBadRun( muEvent->runId() ) ) return kStOk;
	mV0Dst.trgmode=0; //dummy for run 9
	mV0Dst.runnumber = muEvent->runId();
	mV0Dst.evtnumber = muEvent->eventId();

	if ( fabs(muEvent->primaryVertex().x()) < 1e-5 && fabs(muEvent->primaryVertex().y()) < 1e-5 && fabs(muEvent->primaryVertex().z()) < 1e-5 )  return kStOK ;  

	// possible duplicate events.
	if ( mPicoDstMaker->picoDst()->numberOfTracks() == mTestNTrack && mEventsProcessed !=0 && mTestVZ !=0 &&  muEvent->primaryVertex().z() == mTestVZ ) {
		LOG_WARN << mEventsProcessed <<" "<<"seems a duplicated event!"<<endm;
		return kStOK ;
	}
	mTestVZ = muEvent->primaryVertex().z();
	mTestNTrack = mPicoDstMaker->picoDst()->numberOfTracks();

	StPicoTrack* track ;                                             // Pointer to a track

	TVector3 tofpv(-999.,-999.,-999);
	tofpv = muEvent->primaryVertex();
	// Fill some QA plots
	const int ntracks = mPicoDstMaker->picoDst()->numberOfTracks();  
	double RefMult=muEvent->fxtMult();
	hNRefMult -> Fill(RefMult);		
	hNGRefMult -> Fill( muEvent->grefMult() );
	hVertexZ -> Fill( tofpv.z() ) ;
	hvxvy->Fill( tofpv.x() ,tofpv.y() );
	hVpdVz->Fill(tofpv.z() , muEvent->vzVpd() );
	hRefTof -> Fill(RefMult , muEvent->btofTrayMultiplicity() );

	// cut on vertexZ
	if (tofpv.z() > 202 ) return kStOK ;				///run14 run16
	if (tofpv.z() < 198 ) return kStOK ;				///run14 run16
	Float_t vx_ave;
	Float_t vy_ave;
	Float_t vxc, vyc;
	////vx vy shif
	vx_ave = 0.00;  vy_ave = -2.00;//
	vxc = tofpv.x() - vx_ave; vyc = tofpv.y() - vy_ave;

	//if (sqrt(tofpv.x()*tofpv.x()+tofpv.y()*tofpv.y()) > 2.0 ) return kStOK ;
	if( sqrt( vxc*vxc+vyc*vyc )>2 ) return kStOK;
	//if(fabs( tofpv.z() - muEvent->vzVpd() ) >3.0)return kStOK;//

	double ntofmult = muEvent->btofTrayMultiplicity();
	double ntofmatch = muEvent->nBTOFMatch();
	if (mRefMultCorrUtil->isPileUpEvent( RefMult, ntofmatch, tofpv.z() ) ) return kStOk;
	mRefMultCorrUtil->initEvent(RefMult, tofpv.z(), muEvent->ZDCx());
	int centbin = mRefMultCorrUtil->getCentralityBin9();
	hCent->Fill(centbin);

	mV0Dst.primvertexX = tofpv.x();
	mV0Dst.primvertexY = tofpv.y();
	mV0Dst.primvertexZ = tofpv.z();

	mPassEventCut = true;

	//refmultCorrUtil->initEvent(muEvent->refMult(), tofpv.z());
	//hSelectNGRefMult
	hSelectVertexZ->Fill( tofpv.z() );
	hSelectvxvy->Fill( tofpv.x() ,tofpv.y() ); 
	hSelectVpdVz ->Fill(tofpv.z() , muEvent->vzVpd() );
	hSelectRefTof -> Fill(RefMult , muEvent->btofTrayMultiplicity() );
	hSelectGrefTofMatch -> Fill(RefMult,ntofmatch);
	// cut on centrality or reference multiplicity.
	if ( muEvent->refMult() ) {}   //TODO: need to check whether this is the same as in old code. the old code might ignore the case of pile-up.
	mV0Dst.grefmult = muEvent->grefMult();
	mV0Dst.cent = centbin;
	mV0Dst.nrefmult = RefMult;
	mV0Dst.ntofmult = ntofmult;
	mV0Dst.ntofmatch = ntofmatch;

	Float_t weight = 1;                          
	//weight = refmultCorrUtil->getWeight();

	// Do 'event' analysis based on event data 

	// Record some information...
	hSelectNRefMult -> Fill( RefMult ); //this is an ESSENTIAL histogram to record the total number of events for certain centrality. always make sure it is filled AFTER event selection!

	Double_t magn = muEvent->bField();
	mV0Dst.magn = magn;

	StPicoTrack* track_primary ;
	int PrimaryTrackID_Pos[5000];
	double PrimaryTrackPx_Pos[5000],PrimaryTrackPy_Pos[5000],PrimaryTrackPz_Pos[5000];
	int nPrimary_Pos = 0;
	int PrimaryTrackID_Neg[5000];
	double PrimaryTrackPx_Neg[5000],PrimaryTrackPy_Neg[5000],PrimaryTrackPz_Neg[5000];
	int nPrimary_Neg = 0;

	///record the positive and negtive particles
	Int_t nTracks = mPicoDstMaker->picoDst()->numberOfTracks();  
	for(Int_t i=0; i<nTracks; i++)  
	{
		track_primary = mPicoDstMaker->picoDst()->track(i);
		//short flag = track_primary->flag();
		//if(flag <=0  )continue; //or <=0 ?
		if(!track_primary->isPrimary())continue;
		if(abs(track_primary->charge())!=1) continue;

		if(track_primary->charge()>0)
		{
			TVector3 primary_p=track_primary->pMom();
			PrimaryTrackID_Pos[nPrimary_Pos]=track_primary->id();;
			PrimaryTrackPx_Pos[nPrimary_Pos]=primary_p.x();
			PrimaryTrackPy_Pos[nPrimary_Pos]=primary_p.y();
			PrimaryTrackPz_Pos[nPrimary_Pos]=primary_p.z();
			nPrimary_Pos ++;
		}

		if(track_primary->charge()<0)
		{
			TVector3 primary_p=track_primary->pMom();
			PrimaryTrackID_Neg[nPrimary_Neg]=track_primary->id();;
			PrimaryTrackPx_Neg[nPrimary_Neg]=primary_p.x();
			PrimaryTrackPy_Neg[nPrimary_Neg]=primary_p.y();
			PrimaryTrackPz_Neg[nPrimary_Neg]=primary_p.z();
			nPrimary_Neg ++;
		}
	}
	//////////////-----------Select pion---------//////////////
	int nPionplus = 0;
	int nPionminus = 0;
	//test for FMR cut
	for(int ip1=0;ip1<nTracks;ip1++){
		for(int ip2=1+ip1;ip2<nTracks;ip2++){
			StPicoPhysicalHelix tHelix1=  mPicoDstMaker->picoDst()->track(ip1)->helix(magn);
			StPicoPhysicalHelix tHelix2=  mPicoDstMaker->picoDst()->track(ip2)->helix(magn);
			TVector3  PrimVert(0,0,0);//must be kept (0,0,0) or(0,0,200.7)for FXT
			TVector3  SecVert(0,0,0);//for daughter particle
			TVector3 tmpTpcEntrancePoint(0,0,0);
			TVector3 tmpTpcExitPoint(0,0,0);
			TVector3 tmpPosSample(0,0,0);
			float tmpZ1[45];
			float tmpU1[45];
			int tmpSect1[45];
			float tmpZ2[45];
			float tmpU2[45];
			int tmpSect2[45];
			CalculateTpcExitAndEntrancePoints(&tHelix1,&PrimVert,&SecVert,&tmpTpcEntrancePoint,&tmpTpcExitPoint,&tmpPosSample,&tmpZ1[0],&tmpU1[0],&tmpSect1[0]);
			CalculateTpcExitAndEntrancePoints(&tHelix2,&PrimVert,&SecVert,&tmpTpcEntrancePoint,&tmpTpcExitPoint,&tmpPosSample,&tmpZ2[0],&tmpU2[0],&tmpSect2[0]);
			double FMR=calcMergingPar(&tmpU1[0], &tmpU2[0], &tmpZ1[0], &tmpZ2[0], &tmpSect1[0], &tmpSect2[0]);
			if(FMR>0)cout<<"fffffffffffffffffffffffffffff="<<FMR<<endl;
		}
	}
	for(Int_t i=0; i<nTracks; i++)  
	{
		track = mPicoDstMaker->picoDst()->track(i);
		if(!track->isPrimary())continue; 
		TVector3 p = track->pMom();
		double pt = p.Perp();
		double eta = p.PseudoRapidity();
		double phi = p.Phi();                
		hPtRaw -> Fill( pt ) ; 
		hEtaRaw -> Fill( eta ) ;		  //at dca to PV
		hPhiRaw -> Fill( phi ) ;		  //at dca to PV
		Int_t nHits = track->nHits();	//total # of hits in all available detectors
		Int_t nHitsFit = track->nHitsFit();
		Int_t nHitsDedx = track->nHitsDedx();
		Short_t charge = track->charge();

		double nsigmapion = track->nSigmaPion();
		double nsigmaproton = track->nSigmaProton();
		double nsigmakaon = track->nSigmaKaon();
		double nsigmaelectron = track->nSigmaElectron();
		double dedx = track->dEdx();
		//if(p.Mag()<0.6){
		//    if(fabs(nsigmaproton)<2)continue;
		//    if(fabs(nsigmakaon)<2)continue;
		//    if(fabs(nsigmaelectron)<2)continue;
		//}

		//some checks.
		hNHits -> Fill( nHits ) ;
		hNHitsFit -> Fill( nHitsFit ) ;
		hDedxP  -> Fill(p.Mag()*charge ,dedx);
		int indextof = track->bTofPidTraitsIndex();
		float BETA = -1.;              
		if(indextof >= 0){              
			StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(indextof);
			if(tofPid) {             
				BETA = tofPid->btofBeta();
			}
		}
		honeoverBetaP->Fill(p.Mag()*charge,1.0/BETA);

		int hrot;		//helicity of helix, sign of -charge*magn
		if (-charge*magn > 0) hrot = 1;
		else hrot = -1;
		//if(track->vertexIndex()!=StMuDst::currentVertexIndex())continue;
		//if you want to use track->dca(), turn this on. if it is not turned on, that function crashes.

		//OK. let's cut tracks
		if(nHitsFit<=10)continue;
		if(nHitsDedx<=10)continue;
		if(track->nHitsFit() < 0.52 * track->nHitsMax() ) continue;
		if(abs(charge)!=1) continue;
		hPt -> Fill( pt) ; //at dca to PV, for a global track, this value is useless. anyway, the pt value is supposed to be the same anywhere.
		hEta -> Fill( eta ) ;		  //at dca to PV
		hPhi -> Fill( phi ) ;		  //at dca to PV
		if(pt<0.5)hPhiLowPt->Fill( phi );
		else hPhiHighPt->Fill( phi );
		hNSigmaPion->Fill(nsigmapion);
		hNSigmaProton->Fill(nsigmaproton);
		hNSigmaKaon->Fill(nsigmakaon);

		if(pt<cutPtGrEq)continue; //should be larger. like 0.15 or 0.2
		///select pion 
		if(charge == mCharge1 && fabs(nsigmapion)<5){
			//StThreeVectorF pv = muEvent->primaryVertexPosition();
			StPicoPhysicalHelix helix = track->helix(magn);
			double pathlength = helix.pathLength(tofpv, false); 
			TVector3 dca = helix.at(pathlength)-tofpv;
			TVector3 origin = helix.origin();
			if(dca.Mag()>3)continue;//other pico before 24/05/2016 has pion dca < 3.50 cm 

			//kmi
			int index2tof = track->bTofPidTraitsIndex();
			Int_t    btofMatchFlag = 0;  
			double tof= -999;  
			double tofpathlen = -999;
			float beta = -1.;              
			float btofYLocal = -999;              
			float btofZLocal = -999;              
			Float_t btof =0.0;             
			double_t mass2pion = -999;
			if(index2tof >= 0){              
				StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
				btofMatchFlag = tofPid->btofMatchFlag();
				if(tofPid) {             
					beta = tofPid->btofBeta();
					btof = tofPid->btof();
					tof = tofPid->btof();												
					btofYLocal = tofPid->btofYLocal();
					btofZLocal = tofPid->btofZLocal();
					if(fabs(btofYLocal)>1.6)continue;
					if(fabs(btofZLocal)>2.8)continue;
				}                        
				TVector3 btofHitPos_ = tofPid->btofHitPos();
				const StThreeVectorF *btofHitPos = new StThreeVectorF(btofHitPos_.X(),btofHitPos_.Y(),btofHitPos_.Z());        
				const StThreeVectorF *vertexPos_ = new StThreeVectorF(tofpv.x(), tofpv.y(), tofpv.z());

				tofpathlen=tofPathLength(vertexPos_, btofHitPos, helix.curvature());
				//		mass2pion = (p.x()*p.x()+p.y()*p.y()+p.z()*p.z())*(900.0*(tof*tof)/(tofpathlen*tofpathlen)-1.0);
				mass2pion = pow(p.Mag(),2)*(pow(1/beta,2)-1.);										
				double mom = sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
				hTofPID->Fill(mom*charge, mass2pion);
			}
			//20231224 Bijun Fan
			float deltaX = 99999;
			float deltaY = 99999;
			int etofIndex = track->eTofPidTraitsIndex();
			int etofMatchFlag = 0;
			bool is_etof = track->isETofTrack();
			float ebeta = -999;
			float ebetaCheck = -999;
			float ebetaCorrReal = -999;
			float etof = 0;
			float eTrackLength = 0;
			if (etofIndex >= 0)
			{
				StPicoETofPidTraits *etofPid = mPicoDstMaker->picoDst()->etofPidTraits(etofIndex);
				etofMatchFlag = etofPid->matchFlag();
				if (etofPid)
				{
					// eBeta_newT0_Real(mPicoEvent, mPicoTrack, etofPid, &ebetaCorrReal, &etof, &eTrackLength, &ebetaCheck, eTof_tStart);
					Int_t etofHitIndex = etofPid->hitIndex();
					StPicoETofHit* etofHit = mPicoDstMaker->picoDst()->etofHit(etofHitIndex);
					float clusterSize = etofHit->clusterSize();
					deltaX = etofPid->deltaX();
					deltaY = etofPid->deltaY();
					if(fabs(deltaX)>5||fabs(deltaY)>10||clusterSize>=100)continue;
					ebeta = etofPid->beta();
				}
			}
			float emass2pion = -999;
			bool isGoodeTof = etofMatchFlag > 0 && ebeta > 0;
			if (isGoodeTof)
			{
				emass2pion = pow(p.Mag(),2)*(pow(1/ebeta,2)-1.);
			}

			double Tmom_p = sqrt(p.x()*p.x()+p.y()*p.y());

			unsigned long mapMask0 = 0xFFFFFF00;
			unsigned long mapMask1 = 0x1FFFFF;
			ULong64_t     ImapMask = 0x1FFFFFFFFFE;
			//ULong64_t     ImapMask = 0xFFFFFFFFFF;
			unsigned long padRow1To24Track1  = track->topologyMap(0) & mapMask0;
			unsigned long padRow25To45Track1 = track->topologyMap(1) & mapMask1;
			ULong64_t IpadRowTrack1 = track->iTpcTopologyMap() & ImapMask;
			//recording
			mV0Dst.id_Pionplus[nPionplus] = track->id();
			mV0Dst.nhits_Pionplus[nPionplus] = track->nHits();
			mV0Dst.nhitsFit_Pionplus[nPionplus] = track->nHitsFit();
			mV0Dst.dedx_Pionplus[nPionplus] = track->dEdx();
			mV0Dst.nsigmapion_Pionplus[nPionplus] = nsigmapion;
			mV0Dst.nsigmaproton_Pionplus[nPionplus] = track->nSigmaProton();
			mV0Dst.nsigmaelectron_Pionplus[nPionplus] = track->nSigmaElectron();
			mV0Dst.nsigmakaon_Pionplus[nPionplus] = track->nSigmaKaon();
			mV0Dst.px_Pionplus[nPionplus] = p.x();
			mV0Dst.py_Pionplus[nPionplus] = p.y();
			mV0Dst.pz_Pionplus[nPionplus] = p.z();
			mV0Dst.dca_Pionplus[nPionplus] = dca.Mag();
			mV0Dst.tofflag_Pionplus[nPionplus] = btofMatchFlag;
			mV0Dst.tof_Pionplus[nPionplus] = tof;
			mV0Dst.tofpathlen_Pionplus[nPionplus] = tofpathlen;
			mV0Dst.emass2pion_Pionplus[nPionplus] = emass2pion;
			mV0Dst.eBetapion_Pionplus[nPionplus] = ebeta;
			mV0Dst.mass2pion_Pionplus[nPionplus] = mass2pion;
			mV0Dst.Betapion_Pionplus[nPionplus] = beta;
			mV0Dst.PadRow1_Pionplus[nPionplus] = padRow1To24Track1;
			mV0Dst.PadRow2_Pionplus[nPionplus] = padRow25To45Track1;
			mV0Dst.IPadRow_Pionplus[nPionplus] = IpadRowTrack1;
			pionplusmass2->Fill(mass2pion);
			nPionplus++;

		}

		///Select pionminus 
		if(charge == mCharge1_anti && fabs(nsigmapion)<5){
			StPicoPhysicalHelix helix = track->helix(magn);
			double pathlength = helix.pathLength(tofpv, false);       
			TVector3 dca = helix.at(pathlength)-tofpv;
			TVector3 origin = helix.origin();
			if(dca.Mag()>3.0)continue;//other pico before 24/05/2016 has pion dca < 4.0 cm 

			int index2tof = track->bTofPidTraitsIndex();
			Int_t    btofMatchFlag = 0;
			double tof= -999; 
			float beta = -1.;
			Float_t btof =0.0;
			double tofpathlen = -999;
			float btofYLocal = -999;              
			float btofZLocal = -999;              
			double_t mass2pion = -999;
			if(index2tof >= 0){
				StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
				btofMatchFlag = tofPid->btofMatchFlag();
				if(tofPid) {
					beta = tofPid->btofBeta();
					btof = tofPid->btof();
					tof = tofPid->btof();
					btofYLocal = tofPid->btofYLocal();
					btofZLocal = tofPid->btofZLocal();
					if(fabs(btofYLocal)>1.6)continue;
					if(fabs(btofZLocal)>2.8)continue;
				}
				TVector3 btofHitPos_ = tofPid->btofHitPos();
				const StThreeVectorF *btofHitPos = new StThreeVectorF(btofHitPos_.X(),btofHitPos_.Y(),btofHitPos_.Z());                       
				const StThreeVectorF *vertexPos_ = new StThreeVectorF(tofpv.x(), tofpv.y(), tofpv.z());

				tofpathlen=tofPathLength(vertexPos_, btofHitPos, helix.curvature());
				//mass2pion = (p.x()*p.x()+p.y()*p.y()+p.z()*p.z())*(900.0*(tof*tof)/(tofpathlen*tofpathlen)-1.0);
				mass2pion = pow(p.Mag(),2)*(pow(1/beta,2)-1.);
				double mom = sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
				hTofPID->Fill(mom*charge, mass2pion);
			}
			//20231224 Bijun Fan
			float deltaX = 99999;
			float deltaY = 99999;
			int etofIndex = track->eTofPidTraitsIndex();
			int etofMatchFlag = 0;
			bool is_etof = track->isETofTrack();
			float ebeta = -999;
			float ebetaCheck = -999;
			float ebetaCorrReal = -999;
			float etof = 0;
			float eTrackLength = 0;
			if (etofIndex >= 0)
			{
				StPicoETofPidTraits *etofPid = mPicoDstMaker->picoDst()->etofPidTraits(etofIndex);
				etofMatchFlag = etofPid->matchFlag();
				if (etofPid)
				{
					// eBeta_newT0_Real(mPicoEvent, mPicoTrack, etofPid, &ebetaCorrReal, &etof, &eTrackLength, &ebetaCheck, eTof_tStart);
					Int_t etofHitIndex = etofPid->hitIndex();
					StPicoETofHit* etofHit = mPicoDstMaker->picoDst()->etofHit(etofHitIndex);
					float clusterSize = etofHit->clusterSize();
					deltaX = etofPid->deltaX();
					deltaY = etofPid->deltaY();
					if(fabs(deltaX)>5||fabs(deltaY)>10||clusterSize>=100)continue;
					ebeta = etofPid->beta();
				}
			}
			float emass2pion = -999;
			bool isGoodeTof = etofMatchFlag > 0 && ebeta > 0;
			if (isGoodeTof)
			{
				emass2pion = pow(p.Mag(),2)*(pow(1/ebeta,2)-1.);
			}


			double Tmom_p = sqrt(p.x()*p.x()+p.y()*p.y());
			unsigned long mapMask0 = 0xFFFFFF00;
			unsigned long mapMask1 = 0x1FFFFF;
			ULong64_t     ImapMask = 0xFFFFFFFFFF;
			unsigned long padRow1To24Track1  = track->topologyMap(0) & mapMask0;
			unsigned long padRow25To45Track1 = track->topologyMap(1) & mapMask1;
			ULong64_t IpadRowTrack1 = track->iTpcTopologyMap() & ImapMask;
			mV0Dst.id_Pionminus[nPionminus] = track->id();
			mV0Dst.nhits_Pionminus[nPionminus] = track->nHits();
			mV0Dst.nhitsFit_Pionminus[nPionminus] = track->nHitsFit();
			mV0Dst.dedx_Pionminus[nPionminus] = track->dEdx();
			mV0Dst.nsigmapion_Pionminus[nPionminus] = nsigmapion;
			mV0Dst.nsigmaproton_Pionminus[nPionminus] = track->nSigmaProton();
			mV0Dst.nsigmaelectron_Pionminus[nPionminus] = track->nSigmaElectron();
			mV0Dst.nsigmakaon_Pionminus[nPionminus] = track->nSigmaKaon();
			mV0Dst.px_Pionminus[nPionminus] = p.x();
			mV0Dst.py_Pionminus[nPionminus] = p.y();
			mV0Dst.pz_Pionminus[nPionminus] = p.z();
			mV0Dst.dca_Pionminus[nPionminus] = dca.Mag();
			mV0Dst.tofflag_Pionminus[nPionminus] = btofMatchFlag;
			mV0Dst.tof_Pionminus[nPionminus] = tof;
			mV0Dst.tofpathlen_Pionminus[nPionminus] = tofpathlen;
			mV0Dst.emass2pion_Pionminus[nPionminus] = emass2pion;
			mV0Dst.eBetapion_Pionminus[nPionminus] = ebeta;
			mV0Dst.mass2pion_Pionminus[nPionminus] = mass2pion;
			mV0Dst.Betapion_Pionminus[nPionminus] = beta;
			mV0Dst.PadRow1_Pionminus[nPionminus] = padRow1To24Track1;
			mV0Dst.PadRow2_Pionminus[nPionminus] = padRow25To45Track1;
			mV0Dst.IPadRow_Pionminus[nPionminus] = IpadRowTrack1;
			pionminusmass2->Fill(mass2pion);
			nPionminus++;
		}

	}
	mV0Dst.nPionplus = nPionplus;
	mV0Dst.nPionminus = nPionminus;

	if(mDumpNull && mV0Tree) mV0Tree->Fill();
	if(mV0Tree) mV0Tree->Fill();
	//dump v0 vector into a TTree

	mEventsProcessed++ ;
	return kStOK ;
}
int StV0Maker::TpcLocalTransform(TVector3& aPoint, int& aSector, int& aRow, float& aU, double& aPhi){
	static int tNPadAtRow[45]={
		88,96,104,112,118,126,134,142,150,158,166,174,182,
		98,100,102,104,106,106,108,110,112,112,114,116,118,120,122,122,
		124,126,128,128,130,132,134,136,138,138,140,142,144,144,144,144};
	static double tSectToPhi[24]={2.,1.,0.,11.,10.,9.,8. ,7. ,6.,5.,4.,3.,
		4.,5.,6., 7., 8.,9.,10.,11.,0.,1.,2.,3.};
	//static double tPhiToSect[24]={2.,1.,0.,11.,10.,9.,8. ,7. ,6.,5.,4.,3.,
	//			4.,5.,6., 7., 8.,9.,10.,11.,0.,1.,2.,3.};
	static double tPadWidthInner = 0.335;
	static double tPadWidthOuter = 0.67;

	static double tPi = TMath::Pi();
	// --- find sector number
	aPhi = aPoint.Phi();
	if(aPhi<0.) aPhi+=(2*tPi);
	aPhi += tPi/12.;
	if(aPhi>2*tPi) aPhi-=2*tPi;
	int tiPhi = (int) (aPhi/tPi*6.);
	if(aPoint.z()<0) {
		aSector = (tiPhi<3)? 3-tiPhi : 15-tiPhi;
	}
	else{
		aSector = (tiPhi<4)? 21+tiPhi : 9+tiPhi;
	}
	aPhi = tSectToPhi[aSector-1]*tPi/6.;

	// --- calculate local coordinate
	float tR = aPoint.x()*cos(aPhi)+aPoint.y()*sin(aPhi);
	aU =      -aPoint.x()*sin(aPhi)+aPoint.y()*cos(aPhi);

	// --- find pad row 
	if(tR<57.6) {
		aRow = 0;
		return 1;
	}
	float radmax = 62.4;
	float spacing= 4.8;
	aRow=1;
	while(tR>radmax && aRow<46){
		aRow++;
		if(aRow==8){
			radmax = 96.2;
			spacing = 5.2;
		}
		else{
			if (aRow==13){
				radmax = 126.195; // lots of stuf in row 13!
				spacing = 2.0;
			}
			else{
				radmax+=spacing;
			}
		}
	}
	if(aRow>45){
		//cout << "No pad row " << tR << endl;
		return 2;
	}

	// --- Check if u (=aU) inbound
	double tPadWidth = aRow<14? tPadWidthInner : tPadWidthOuter;
	if(fabs(aU) > tNPadAtRow[aRow-1]*tPadWidth/2.){
		return 3;
	}

	return 0;
}
double StV0Maker::calcMergingPar(float *track1_mU, float * track2_mU, float *track1_mZ, float * track2_mZ, int * track1_Sec, int * track2_Sec) const{
	double mMaxDuInner = .8;
	double mMaxDzInner = 3.;
	double mMaxDuOuter = 1.4;
	double mMaxDzOuter = 3.2;

	double tDu, tDz;
	int tN = 0;
	double mFracOfMergedRow = 0.;
	double mWeightedAvSep =0.;
	double tDist;
	double tDistMax = 200.;
	for(int ti=0 ; ti<45 ; ti++){
		if(track1_Sec[ti]==track2_Sec[ti] && track1_Sec[ti]!=-1){
			tDu = fabs(track1_mU[ti]-track2_mU[ti]);
			tDz = fabs(track1_mZ[ti]-track2_mZ[ti]);
			tN++;
			if(ti<13){
				mFracOfMergedRow += (tDu<mMaxDuInner && tDz<mMaxDzInner);
				tDist = ::sqrt(tDu*tDu/mMaxDuInner/mMaxDuInner+
						tDz*tDz/mMaxDzInner/mMaxDzInner);
				//mFracOfMergedRow += (tDu<mMaxDuInner && tDz<mMaxDzInner);
			}
			else{
				mFracOfMergedRow += (tDu<mMaxDuOuter && tDz<mMaxDzOuter);
				tDist = ::sqrt(tDu*tDu/mMaxDuOuter/mMaxDuOuter+
						tDz*tDz/mMaxDzOuter/mMaxDzOuter);
				//mFracOfMergedRow += (tDu<mMaxDuOuter && tDz<mMaxDzOuter);
			}
			if(tDist<tDistMax){
				//	mClosestRowAtDCA = ti+1;
				tDistMax = tDist;
			}
			mWeightedAvSep += tDist;
		}
	}
	if(tN>0){
		mWeightedAvSep /= tN;
		mFracOfMergedRow /= tN;
	}
	else{
		//  mClosestRowAtDCA = -1;
		mFracOfMergedRow = -1.;
		mWeightedAvSep = -1.;
	}
	return mFracOfMergedRow;
}
void StV0Maker::CalculateTpcExitAndEntrancePoints(StPicoPhysicalHelix* tHelix,
		TVector3*  PrimVert,
		TVector3*  SecVert,
		TVector3* tmpTpcEntrancePoint,
		TVector3* tmpTpcExitPoint,
		TVector3* tmpPosSample,
		float* tmpZ,
		float* tmpU,
		int* tmpSect){
	// this calculates the exit point of a secondary track, 
	// either through the endcap or through the Outer Field Cage
	// We assume the track to start at tHelix.origin-PrimaryVertex
	// it also calculates the entrance point of the secondary track, 
	// which is the point at which it crosses the
	// inner field cage
	//  static TVector3 ZeroVec(0.,0.,0.);
	TVector3 ZeroVec(0.,0.,0.);
	//   ZeroVec.SetX(tHelix->origin().x()-PrimVert->X());
	//   ZeroVec.SetY(tHelix->origin().y()-PrimVert->Y());
	//   ZeroVec.SetZ(tHelix->origin().z()-PrimVert->Z());
	ZeroVec.SetX(SecVert->X()-PrimVert->X());
	ZeroVec.SetY(SecVert->Y()-PrimVert->Y());
	ZeroVec.SetZ(SecVert->Z()-PrimVert->Z());
	double dip, curv, phase;
	int h;
	curv = tHelix->curvature();
	dip  = tHelix->dipAngle();
	phase= tHelix->phase();
	h    = tHelix->h();

	StPicoHelix hel(curv,dip,phase,ZeroVec,h);

	std::pair<double,double> candidates;
	double sideLength;  // this is how much length to go to leave through sides of TPC
	double endLength;  // this is how much length to go to leave through endcap of TPC
	// figure out how far to go to leave through side...
	candidates = hel.pathLength(200.0);  // bugfix MAL jul00 - 200cm NOT 2cm
	sideLength = (candidates.first > 0) ? candidates.first : candidates.second;

	static TVector3 WestEnd(0.,0.,200.);  // bugfix MAL jul00 - 200cm NOT 2cm
	static TVector3 EastEnd(0.,0.,-200.); // bugfix MAL jul00 - 200cm NOT 2cm
	static TVector3 EndCapNormal(0.,0.,1.0);

	endLength = hel.pathLength(WestEnd,EndCapNormal);
	if (endLength < 0.0) endLength = hel.pathLength(EastEnd,EndCapNormal);

	if (endLength < 0.0) cout << 
		"StHbtParticle::CalculateTpcExitAndEntrancePoints(): "
			<< "Hey -- I cannot find an exit point out endcaps" << endl;
	// OK, firstExitLength will be the shortest way out of the detector...
	double firstExitLength = (endLength < sideLength) ? endLength : sideLength;
	// now then, let's return the POSITION at which particle leaves TPC...
	*tmpTpcExitPoint = hel.at(firstExitLength);
	// Finally, calculate the position at which the track crosses the inner field cage
	candidates = hel.pathLength(50.0);  // bugfix MAL jul00 - 200cm NOT 2cm

	sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
	//  cout << "sideLength 2 ="<<sideLength << endl;
	*tmpTpcEntrancePoint = hel.at(sideLength);
	// This is the secure way !  
	if (::isnan(tmpTpcEntrancePoint->X()) || 
			::isnan(tmpTpcEntrancePoint->Y()) || 
			::isnan(tmpTpcEntrancePoint->Z()) ){ 
		cout << "tmpTpcEntrancePoint NAN"<< endl; 
		cout << "tmpNominalTpcEntrancePoint = " <<tmpTpcEntrancePoint<< endl;
		tmpTpcEntrancePoint->SetX(-9999.);
		tmpTpcEntrancePoint->SetY(-9999.);
		tmpTpcEntrancePoint->SetZ(-9999.);
	} 

	if (::isnan(tmpTpcExitPoint->X()) || 
			::isnan(tmpTpcExitPoint->Y()) || 
			::isnan(tmpTpcExitPoint->Z()) ) {
		//     cout << "tmpTpcExitPoint NAN set at (-9999,-9999,-9999)"<< endl; 
		//     cout << "tmpTpcExitPoint X= " <<tmpTpcExitPoint->X()<< endl;
		//     cout << "tmpTpcExitPoint Y= " <<tmpTpcExitPoint->Y()<< endl;
		//     cout << "tmpTpcExitPoint Z= " <<tmpTpcExitPoint->Z()<< endl;
		tmpTpcExitPoint->SetX(-9999.);
		tmpTpcExitPoint->SetY(-9999.);
		tmpTpcExitPoint->SetZ(-9999.);
	}


	//   if (::isnan(tmpTpcExitPoint->X())) *tmpTpcExitPoint = TVector3(-9999.,-9999.,-9999); 
	//   if (::isnan(tmpTpcEntrancetPoint->X())) *tmpTpcEntrancePoint = TVector3(-9999.,-9999.,-9999); 
	//  cout << "tmpTpcEntrancePoint"<<*tmpTpcEntrancePoint << endl;

	// 03Oct00 - mal.  OK, let's try something a little more 
	// along the lines of NA49 and E895 strategy.
	// calculate the "nominal" position at N radii (say N=11) 
	// within the TPC, and for a pair cut
	// use the average separation of these N
	int irad = 0;
	candidates = hel.pathLength(50.0);
	sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
	while (irad<11 && !::isnan(sideLength)){
		float radius = 50.0 + irad*15.0;
		candidates = hel.pathLength(radius);
		sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
		tmpPosSample[irad] = hel.at(sideLength);
		if(::isnan(tmpPosSample[irad].x()) ||
				::isnan(tmpPosSample[irad].y()) ||
				::isnan(tmpPosSample[irad].z()) 
		  ){
			cout << "tmpPosSample for radius=" << radius << " NAN"<< endl; 
			//cout << "tmpPosSample=(" <<tmpPosSample[irad]<<")"<< endl;
			tmpPosSample[irad] =  TVector3(-9999.,-9999.,-9999);
		}
		irad++;
		if (irad<11){
			float radius = 50.0 + irad*15.0;
			candidates = hel.pathLength(radius);
			sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
		}
	}
	for (int i = irad; i<11; i++)
	{
		tmpPosSample[i] =  TVector3(-9999.,-9999.,-9999);   
	}

	static float tRowRadius[45] = {60,64.8,69.6,74.4,79.2,84,88.8,93.6,98.8, 
		104,109.2,114.4,119.6,127.195,129.195,131.195,
		133.195,135.195,137.195,139.195,141.195,
		143.195,145.195,147.195,149.195,151.195,
		153.195,155.195,157.195,159.195,161.195,
		163.195,165.195,167.195,169.195,171.195,
		173.195,175.195,177.195,179.195,181.195,
		183.195,185.195,187.195,189.195};
	int tRow,tSect,tOutOfBound;
	double tLength,tPhi;
	float tU;
	TVector3 tPoint;
	TVector3 tn(0,0,0);
	TVector3 tr(0,0,0);
	int ti =0;
	// test to enter the loop
	candidates =  hel.pathLength(tRowRadius[ti]);
	tLength = (candidates.first > 0) ? candidates.first : candidates.second;
	if (::isnan(tLength)){
		cout <<"tLength Init tmp NAN" << endl;
		cout <<"padrow number= "<<ti << "not reached" << endl;
		cout << "*** DO NOT ENTER THE LOOP***" << endl;
		tmpSect[ti]=-1;//sector
	}
	// end test
	while(ti<45 && !::isnan(tLength)){
		candidates =  hel.pathLength(tRowRadius[ti]);
		tLength = (candidates.first > 0) ? candidates.first : candidates.second;
		if (::isnan(tLength)){
			cout <<"tLength loop 1st NAN" << endl;
			cout <<"padrow number=  " << ti << " not reached" << endl;
			cout << "*** THIS IS AN ERROR SHOULDN'T  LOOP ***" << endl;
			tmpSect[ti]=-1;//sector
		}
		tPoint = hel.at(tLength);
		// Find which sector it is on
		TpcLocalTransform(tPoint,tmpSect[ti],tRow,tU,tPhi);
		if (::isnan(tmpSect[ti])){
			cout <<"***ERROR tmpSect"<< endl; 
		}
		if (::isnan(tRow)){
			cout <<"***ERROR tRow"<< endl;
		}
		if (::isnan(tU)){
			cout <<"***ERROR tU"<< endl;
		}
		if (::isnan(tPhi)){
			cout <<"***ERROR tPhi"<< endl;
		}  
		// calculate crossing plane
		tn.SetX(cos(tPhi));
		tn.SetY(sin(tPhi));       
		tr.SetX(tRowRadius[ti]*cos(tPhi));
		tr.SetY(tRowRadius[ti]*sin(tPhi));
		// find crossing point
		tLength = hel.pathLength(tr,tn); 
		if (::isnan(tLength)){
			cout <<"tLength loop 2nd  NAN" << endl;
			cout <<"padrow number=  " << ti << " not reached" << endl;
			tmpSect[ti]=-2;//sector
		}
		tPoint = hel.at(tLength);
		tmpZ[ti] = tPoint.z();
		tOutOfBound = TpcLocalTransform(tPoint,tSect,tRow,tmpU[ti],tPhi);
		if (::isnan(tSect)){
			cout <<"***ERROR tSect 2"<< endl; 
		}
		if (::isnan(tRow)){
			cout <<"***ERROR tRow 2"<< endl;
		}
		if (::isnan(tmpU[ti])){
			cout <<"***ERROR tmpU[ti] 2"<< endl;
		}
		if (::isnan(tPhi)){
			cout <<"***ERROR tPhi 2 "<< endl;
		}  
		if(tOutOfBound || (tmpSect[ti] == tSect && tRow!=(ti+1))){
			tmpSect[ti]=-2;
			//	  cout << "missed once"<< endl;
		}
		else{
			if(tmpSect[ti] != tSect){
				// Try again on the other sector
				tn.SetX(cos(tPhi));
				tn.SetY(sin(tPhi));       
				tr.SetX(tRowRadius[ti]*cos(tPhi));
				tr.SetY(tRowRadius[ti]*sin(tPhi));
				// find crossing point
				tLength = hel.pathLength(tr,tn);
				tPoint = hel.at(tLength);
				if (::isnan(tLength)){
					cout <<"tLength loop 3rd NAN" << endl;
					cout <<"padrow number=  "<< ti << " not reached" << endl;
					tmpSect[ti]=-1;//sector
				}
				tmpZ[ti] = tPoint.z();
				tmpSect[ti] = tSect;
				tOutOfBound = TpcLocalTransform(tPoint,tSect,tRow,tmpU[ti],tPhi);
				if (::isnan(tSect)){
					cout <<"***ERROR tSect 3"<< endl; 
				}
				if (::isnan(tRow)){
					cout <<"***ERROR tRow 3"<< endl;
				}
				if (::isnan(tmpU[ti])){
					cout <<"***ERROR tmpU[ti] 3"<< endl;
				}
				if (::isnan(tPhi)){
					cout <<"***ERROR tPhi 3 "<< endl;
				}  
				if(tOutOfBound || tSect!= tmpSect[ti] || tRow!=(ti+1)){
					tmpSect[ti]=-1;
				}
			}
		}
		if (::isnan(tmpSect[ti])){
			cout << "*******************ERROR***************************" << endl;
			cout <<"StHbtParticle--Fctn tmpSect=" << tmpSect[ti] << endl;
			cout << "*******************ERROR***************************" << endl;
		}
		if (::isnan(tmpU[ti])){
			cout << "*******************ERROR***************************" << endl;
			cout <<"StHbtParticle--Fctn tmpU=" << tmpU[ti] << endl;
			cout << "*******************ERROR***************************" << endl;
		}
		if (::isnan(tmpZ[ti])){
			cout << "*******************ERROR***************************" << endl;
			cout <<"StHbtParticle--Fctn tmpZ=" << tmpZ[ti] << endl;
			cout << "*******************ERROR***************************" << endl;
		}
		// If padrow ti not reached all other beyond are not reached
		// in this case set sector to -1
		if (tmpSect[ti]==-1){
			for (int tj=ti; tj<45;tj++){
				tmpSect[tj] = -1;
				ti=45;
			}
		}
		ti++;
		if (ti<45){
			candidates =  hel.pathLength(tRowRadius[ti]);
			tLength = (candidates.first > 0) ? candidates.first : candidates.second;}
	}
}
Int_t StV0Maker::Finish( )
{ // Do once at the end the analysis

	std::cout << "pass 2"  << std::endl;
	// Write histograms to disk, output miscellaneous other information
	if(histogram_output!=NULL) histogram_output -> Write() ;   // Write all histograms to disk 
	if(v0tree_output!=NULL) v0tree_output -> Write() ;   // Write all histograms to disk 

	cout << "Total Events Processed in StV0Maker " << mEventsProcessed << endl ;

	return kStOk ;  

}


