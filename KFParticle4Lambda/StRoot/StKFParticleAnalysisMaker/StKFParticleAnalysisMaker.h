#ifndef StKFParticleAnalysisMaker_h
#define StKFParticleAnalysisMaker_h

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"

#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "StMaker.h"
#include "TString.h"
#include "TObject.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include "TF1.h"
#include "StTrackHelix.h"
#include "MyConstant.h" // must include

class StPicoDst;
class StPicoDstMaker;
class TString;
class KFParticle;
class StKFParticleInterface;
class StKFParticlePerformanceInterface;
class TH1F;
class TH2F;
class TH2D;
class TH3F;
class TH3D;
class TF1;
class TProfile;
class TProfile2D;
class TProfile3D;
class CentralityMaker;
class StRefMultCorr;

class StKFParticleAnalysisMaker : public StMaker 
{
public:
	StKFParticleAnalysisMaker(const char *name, const char *outName);
	virtual ~StKFParticleAnalysisMaker();

	virtual Int_t Init();
	virtual Int_t Make();
	virtual void  Clear(Option_t *opt="");
	virtual Int_t Finish();

	void    setRunEnergyAndListDir(int run,double energy,char ListDir[256]);            

private:
	// KFParticle
	StKFParticleInterface *KFParticleInterface;
	StKFParticlePerformanceInterface *KFParticlePerformanceInterface;
	void SetupKFParticle();
	void SetDaughterTrackPointers(int iKFParticle);
	bool IsKaonOmegaDaughter(KFParticle particle, int kaonTrackId);
	bool IsTrackParticleDaughter(KFParticle particle, int TrackId);
	void SetDaughterTrackHits(KFParticle particle);
	int TrackID(StPicoTrack *track , TVector3 Vertex3D , double magnet , bool Track_has_tof , float m2 = -999. , float beta = -999.);
	TVector3 StKFParticleAnalysisMaker::LocAfterTransfer(StPicoPhysicalHelix Track , double Length);
	double StKFParticleAnalysisMaker::DistanceBetween(TVector3 LA , TVector3 LB);
	bool StKFParticleAnalysisMaker::IfGoodDaughterDCA(StPicoDst* mPicoDst , int iKFParticle , double magnet , double Gen1_DCALim , double Gen2_DCALim);
	std::vector<bool> StKFParticleAnalysisMaker::TrackPID(std::vector<int>& TestPDG , StPicoTrack *track , TVector3 Vertex3D);
	Double_t StKFParticleAnalysisMaker::massList(int PID);
	bool InterfaceCantProcessEvent;
	int ProtonTrackIndex, PionTrackIndex, KaonTrackIndex;
	vector<int> trackMap;
	StPicoDst *PicoDst;
	StPicoTrack *ProtonTrack, *PionTrack, *KaonTrack;
	void BookVertexPlots();

	StPicoDstMaker *mPicoDstMaker;
	StRefMultCorr *mRefMultCorr;

	// cut params for coalescence
	float pT_lo, pT_hi;
	float pT_asso_lo, pT_asso_hi;
	float pT_trig_lo, pT_trig_hi;
	float eta_trig_cut;
	float pion_pT_lo, pion_pT_hi;
	float proton_pT_lo, proton_pT_hi;
	float pion_pT_TOFth; // threshold above which TOF becomes required
	float proton_pT_TOFth;
	float pion_m2_lo, pion_m2_hi;
	float proton_m2_lo, proton_m2_hi;
	float dcatoPV_hi;

	int        mRun;            
	double     mEnergy;            
	TString    mListDir;            

	TString    mOutName;
	double     PI;
	double     twoPI;

	int        mJob;
	std::vector<int> Recorded_runID;

	////////////////
	TH1F *hNRefMult;
	TH1F *hNRefMultA;
	TH1F *hNRefMultB;
	TH2F *hVertexXY;
	TH1F *hVertexZ;
	TH2F *hNch_per_VertexZ;
	TH2F *hVertex2D; 
	TH1F *hDiffVz  ; 
	TH1F *hcent;
	TH1F *hcentw;
	TH1F *H_Total_Pz;
	TH2F *H_Total_Pxy;

	TFile *fout;
	TDirectory* folder_EventQA;
	TDirectory* folder_PIDQA;
	TDirectory* folder_ReconsQA;

	#define PDG2NameSize 6 // APDGList.size()
	#define PDG2NameSize2 6 // BPDGList.size()
	#define PDG2NameSize3 3 // CPDGList.size()
	int PDGList[PDG2NameSize + PDG2NameSize2];
	TString NameList[PDG2NameSize + PDG2NameSize2];
	int BPDGListMass[PDG2NameSize2];
	int CPDGList[PDG2NameSize3];
	TString CNameList[PDG2NameSize3];
	// std::map<int, TString> PDG2Name;
	// PDG2Name[ 3122] = "Lambda"; 
	// PDG2Name[-3122] = "Lambdab";
	// PDG2Name[ 3334] = "Omega";  
	// PDG2Name[-3334] = "Omegab"; 
	TH1F *H_ALL_NO_CUT[PDG2NameSize];// NO CUT
	TH1F *H_DaughterDCA[PDG2NameSize];// Cut DaughterDCA
	TH1F *H_WrongDaughter[PDG2NameSize];
	TH1F *H_CrectDaughter[PDG2NameSize];
	TH1F *H_Hyperon_Rap[PDG2NameSize];
	TDirectory* KFPRecons[PDG2NameSize];
	TH1F *H_rapidity[PDG2NameSize2];
	TH1F *H_rapidity_eTOF[PDG2NameSize2];
	TH1F *H_P[PDG2NameSize2];
	TH1F *H_Pt[PDG2NameSize2];
	TH1F *H_DCAtoPV[PDG2NameSize2];
	TH1F *H_eta[PDG2NameSize2];
	TH2F *H_y_Pt[PDG2NameSize2];
	TH2F *H_y_P[PDG2NameSize2];
	TH2F *H_y_m2[PDG2NameSize2];
	TH2F *H_y_nSigmaKaon[PDG2NameSize2];
	TH2F *H_y_nSigmaPion[PDG2NameSize2];
	TH2F *H_y_nSigmaElectron[PDG2NameSize2];
	TH2F *H_y_nHitsFit[PDG2NameSize2];
	TH2F *H_y_nHitsDedx[PDG2NameSize2];
	TH2F *H_y_nHitsFit2nHitsMax[PDG2NameSize2];
	TH2F *H_y_eta[PDG2NameSize2];
	TH2F *H_y_Vz[PDG2NameSize2];
	TH2F *hgbtofYlocal[PDG2NameSize2];
	TH2F *H_y_Pz[PDG2NameSize2];
	TH2F *H_Pxy[PDG2NameSize2];
	TH1F *H_Pz[PDG2NameSize2];
	TH2F *H_y_Nch[PDG2NameSize2];
	TH2F *H_Pz_Nch[PDG2NameSize2];
	TH2F *H_y_nSigmaTOFKaon[PDG2NameSize2];
	TH2F *H_m2_nSigmaTOFKaon[PDG2NameSize2];
	TDirectory* PID_Tracks[PDG2NameSize2];

	TH2F *H_eta_nSigmaKaon  [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_nSigmaPion  [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_nSigmaProton[30][3];// Trigger Num not larger than 30
	TH2F *H_eta_m2          [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_PVz         [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_PVr         [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_DVz         [30][3];// Trigger Num not larger than 30
	TH1F *H_eta_triggerBIN  [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_trigger     ;

	TH2F *H_nHitsFit_p[PDG2NameSize2];
	TH1F *H_nHitsFit_nHitsMax[PDG2NameSize2];
	TH1F *H_ndEdx[PDG2NameSize2];
	TH2F *H_nSigmaTOF_p[PDG2NameSize2];
	TH2F *H_dEdx_p[PDG2NameSize2];
	TH2F *H_Pt_nSigma[PDG2NameSize2][PDG2NameSize3];
	TH2F *H_Pt_m2;
	TH2F *H_Pt_nSigmaKaon;
	TH2F *H_Pt_nSigmaKaonTOF;
	TH3F *H_m2_nSigmaKaon_Pt;
	TH2F *H_m2_KSigma_S;
	TH2F *H_m2_KSigma_L;
	TH2F *H_All_nSigmaKaon_y;
	TH2F *H_All_nSigmaKaon_eta;
	TH2F *H_nSigmaKaon_Pt_AllTOF;
	TH2F *H_nSigmaKaon_Pt_HasTOF;
	TH2F *H_nSigmaKaon_Pt_HasTOF_NoNsigmaPion;

	// KFP PID QA
	TDirectory* KFPPIDQA;
	TDirectory* KFPPID[PDG2NameSize2];
	TH1F *H_KFP_rapidity[PDG2NameSize2];
	TH1F *H_KFP_r[PDG2NameSize2]; // No meaning, just for compile successfuly...
	TH1F *H_KFP_P[PDG2NameSize2];
	TH1F *H_KFP_Pt[PDG2NameSize2];
	TH1F *H_KFP_DCAtoPV[PDG2NameSize2];
	TH1F *H_KFP_eta[PDG2NameSize2];
	TH2F *H_KFP_y_Pt[PDG2NameSize2];
	TH2F *H_KFP_y_m2[PDG2NameSize2];
	TH2F *H_KFP_y_nSigmaKaon[PDG2NameSize2];
	TH2F *H_KFP_y_nSigmaPion[PDG2NameSize2];
	TH2F *H_KFP_y_nSigmaElectron[PDG2NameSize2];
	TH2F *H_KFP_y_nHitsFit[PDG2NameSize2];
	TH2F *H_KFP_y_nHitsDedx[PDG2NameSize2];
	TH2F *H_KFP_y_nHitsFit2nHitsMax[PDG2NameSize2];
	TH2F *H_KFP_y_eta[PDG2NameSize2];
	TH2F *H_KFP_y_Vz[PDG2NameSize2];
	TH2F *h_KFP_tofYlocal[PDG2NameSize2];
	TH2F *H_KFP_y_Pz[PDG2NameSize2];
	TH2F *H_KFP_Pxy[PDG2NameSize2];
	TH1F *H_KFP_Pz[PDG2NameSize2];
	TH2F *H_KFP_y_Nch[PDG2NameSize2];
	TH2F *H_KFP_Pz_Nch[PDG2NameSize2];
	TH2F *H_KFP_y_nSigmaTOFKaon[PDG2NameSize2];
	TH2F *H_KFP_m2_nSigmaTOFKaon[PDG2NameSize2];
	TH2F *H_KFP_nHitsFit_p[PDG2NameSize2];
	TH1F *H_KFP_nHitsFit_nHitsMax[PDG2NameSize2];
	TH1F *H_KFP_ndEdx[PDG2NameSize2];
	TH2F *H_KFP_nSigmaTOF_p[PDG2NameSize2];
	TH2F *H_KFP_dEdx_p[PDG2NameSize2];
	TH2F *H_KFP_Pt_nSigma[PDG2NameSize2][PDG2NameSize3];
	TH2F *H_KFP_Pt_m2[PDG2NameSize2];


	TH2F *hdEdx_pQ;
	TH2F *hdEdx_pQ_1cut;
	TH2F *hdEdx_pQ_2cut;
	TH2D *hLN_M;
	TH2D *hXY;
	TH2D *hHXY;
	TH2D *hHM_Chi2;
	TH2D *hHM_ParentDCA;
	TH1D *hEventNum;

	TProfile *hcentRefM ; 
	TProfile *hcentRefW ; 

	TTree *hadronTree;
	int buffer_size,CrefMult,CgrefMult,evtID,runID,PDGMult , Omega_Omegab_Num;
	std::vector<int> PDG , ReCons_TrackID;
	std::vector<float> px,py,pz,InvariantMass;
	double zTOF_proton,zTOF_pion,zTOF_kaon;
	// Used for QA
	std::vector<float> QA_dEdx,QA_m2,QA_nSigmaProton,QA_nSigmaPion,QA_nSigmaKaon,QA_Chi2;
	std::vector<double> QA_zTOF_proton,QA_zTOF_pion,QA_zTOF_kaon,QA_Decay_Length,QA_DCA_V0_PV,QA_DCA_Daughters;
	std::vector<int> QA_hasTOF,QA_IfConfuse,QA_IfBadReconstructed;// Used as bool

	/////////////////////////////////////
	int mStps;  

	std::vector<Int_t> badList;
	std::vector<Int_t> runList;

	Int_t findCentrality(int mult);
	Int_t CheckrunNumber(int runnumber);            
	bool  readRunList();            
	bool  readBadList();            
	bool  removeBadID(int runnumber);            
	Int_t openFile();

	void  DeclareHistograms();
	void  WriteHistograms();

	bool isGoodObs(double obs);
	
	// For Ping Siyuan
	float H_ProcessEventNum;
		
	ClassDef(StKFParticleAnalysisMaker, 1)
};

#endif


