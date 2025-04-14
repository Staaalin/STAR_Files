#ifndef StKFParticleAnalysisMaker_h
#define StKFParticleAnalysisMaker_h

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"

#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"


#include "StRoot/KFParticle/KFPTrack.h"
#include "StRoot/KFParticle/KFPVertex.h"
#include "StRoot/KFParticle/KFParticle.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
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

class KFPTrack;
class KFPVertex;
class KFParticle;

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
	TVector3 LocAfterTransfer(StPicoPhysicalHelix Track , double Length);
	double DistanceBetween(TVector3 LA , TVector3 LB);
	bool IfGoodDaughterDCA(StPicoDst* mPicoDst , int iKFParticle , double magnet , double Gen1_DCALim , double Gen2_DCALim);
	std::vector<bool> TrackPID(std::vector<int>& TestPDG , StPicoTrack *track , TVector3 Vertex3D);
	Double_t massList(int PID);
	void print(std::vector<int> Temp);
	void print(std::vector<std::vector<int> > Temp);
	double getSL(Int_t padRow1To24Track1 ,Int_t padRow25To45Track1 ,ULong64_t IpadRow1 ,Int_t nhits1 , Int_t padRow1To24Track2 , Int_t padRow25To45Track2 ,ULong64_t IpadRow2 ,Int_t nhits2,bool IfITPC_T);
	double getphistar(float phi1, float phi2, float Pt1, float Pt2, int q1, int q2,double Bz, double tpcR);
	float GetPairMass(float p1x , float p1y , float p1z , float p2x , float p2y , float p2z , float AMass , float BMass);
	KFParticle ChangeMass(KFParticle daughter);
	bool InterfaceCantProcessEvent;
	int ProtonTrackIndex, PionTrackIndex, KaonTrackIndex;
	vector<int> trackMap;
	StPicoDst *PicoDst;
	StPicoTrack *ProtonTrack, *PionTrack, *KaonTrack;
	void BookVertexPlots();

	vector<vector<int> > Recorded_KFP_ID; int SplitNum;
	// Recorded_KFP_ID structure:
	// 0             { { Reconstructed 1 location in KFP , daughter track 1 location in DST , daughter track 2 location in DST , ... } }
	// 1             { { Reconstructed 2 location in KFP , daughter track 1 location in DST , daughter track 2 location in DST , ... } }
	//               ... 
	// SplitNum      { { 1 PDG                           , track location in DST            } }
	// SplitNum + 1  { { 2 PDG                           , track location in DST            } }
	//               ...

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

	float MSta , MEnd;

	int N_Entries;

	int        mRun;            
	double     mEnergy;            
	TString    mListDir;            
	double     B_inTesla;

	TString    mOutName;
	double     PI;
	double     twoPI;

	int        mJob;
	std::vector<int> Recorded_runID;

	bool IfITPC;
	// SL value
	float  slcutmin = -0.5;
	float  slcutmax = 0.6 ;
	double  SL_Value;
	unsigned long mapMask0 = 0xFFFFFF00;
	unsigned long mapMask1 = 0x1FFFFF;
	ULong64_t     ImapMask = 0x1FFFFFFFFFE;
	Int_t padRow1to24TrackA, padRow1to24TrackB;
	Int_t padRow25to45TrackA,padRow25to45TrackB;
	ULong64_t     IpadRowTrackA,IpadRowTrackB;

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

	bool  IfFill_BM;
	TH2F* H_ALL_OR_Lambda  ;
	TH2F* H_ALL_OR_Lambdab ;
	TH2F* H_ALL_Lambda  ;
	TH2F* H_ALL_Lambdab ;
	TH2F* H_ALLr_Lambda ;
	TH2F* H_ALLr_Lambdab;
	TH2F* H_ALLp_Lambda ;
	TH2F* H_ALLp_Lambdab;
	TH2F* H_ALL_OR_Xi  ;
	TH2F* H_ALL_OR_Xib ;
	TH2F* H_ALL_Xi      ;
	TH2F* H_ALL_Xib     ;
	TH2F* H_ALLr_Xi     ;
	TH2F* H_ALLr_Xib    ;
	TH2F* H_ALLp_Xi     ;
	TH2F* H_ALLp_Xib    ;
	TH2F* H_ALL_OR_Omega  ;
	TH2F* H_ALL_OR_Omegab ;
	TH2F* H_ALL_Omega   ;
	TH2F* H_ALL_Omegab  ;
	TH2F* H_ALLr_Omega  ;
	TH2F* H_ALLr_Omegab ;
	TH2F* H_ALLp_Omega  ;
	TH2F* H_ALLp_Omegab ;
	float AMass ;
	float BMass ;
	float APx   ;
	float APy   ;
	float APz   ;
	float BPx   ;
	float BPy   ;
	float BPz   ;
	float MPt   ;
	float MPz   ;
	float MEnergy;
	float MRap  ;
	float MRap_T;
	KFPTrack   KFPtrack_A,KFPtrack_B,KFPtrack_C;
	KFParticle Particle_M,*Particle_N2[2],*Particle_N3[3],*KFP_PV_P,Particle_A,Particle_B,Particle_C;
	KFPVertex  KFP_PV;
	KFVertex   KF_PV;

	//       0          1             2         3         4
	// { { PDGID , KFP Loc index , Track 1 , Track 2 , Track 3 , ... } ,
	//   { PDGID , KFP Loc index , Track 1 , Track 2 , Track 3 , ... } ,
	//   ...
	// }
	std::vector<std::vector<int> > KFParticleList; 
	TH1F* H_OmegaR_XiKPi_Mass; // Reconstructed by KFP
	TH1F* H_OmegabR_XiKPi_Mass;
	TH1F* H_OmegaR_XiK_Mass; // Xi- + K0S
	TH1F* H_OmegabR_XiK_Mass; 
	TH1F* H_OmegaR_OmegaPiPi_Mass; // Omega- + pi+ + pi-
	TH1F* H_OmegabR_OmegaPiPi_Mass; 
	TH1F* H_Omega0R_OmegaPi_Mass; // Omega- + pi+
	TH1F* H_Omega0bR_OmegaPi_Mass; 
	TH1F* H_Omega0R_XiK_Mass; // Xi- + K+
	TH1F* H_Omega0bR_XiK_Mass; 
	bool IfPass;

	TFile *fout;
	TDirectory* folder_EventQA;
	TDirectory* folder_PIDQA;
	TDirectory* folder_ReconsQA;
	TDirectory* folder_LoadHY;
	TDirectory* folder_RecNewP;

	#define PDG2NameSize  10 // APDGList.size()
	#define PDG2NameSize2 6  // BPDGList.size()
	#define PDG2NameSize3 3  // CPDGList.size()
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
	TH1F *H_rapidity_Only_eTOF[PDG2NameSize2];
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

	TH2F *H_eta_nSigmaKaon         [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_nSigmaPion         [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_nSigmaProton       [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_m2                 [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_PVz                [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_PVr                [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_DVz                [30][3];// Trigger Num not larger than 30
	TH1F *H_eta_triggerBIN         [30][3];// Trigger Num not larger than 30
	TH1F *H_eta_triggerBIN_hasTOF  [30][3];// Trigger Num not larger than 30
	TH1F *H_Nch_triggerBIN         [30][3];// Trigger Num not larger than 30
	TH2F *H_eta_trigger;

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

	TDirectory* DTest;
	// Splite & Merge Effect
	TH2F *H_Before_Merge_Phi_Eta_Kaon_LambdaDaughter;
	TH2F *H_After_Merge_Phi_Eta_Kaon_LambdaDaughter;

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
	TH2D *hXY;
	TH2D *hHXY;
	TH2D *hHM_Chi2;
	TH2D *hHM_ParentDCA;
	TH1D *hEventNum;

	TProfile *hcentRefM ; 
	TProfile *hcentRefW ; 

	TTree *hadronTree;
	int buffer_size,CrefMult,CgrefMult,evtID,runID,PDGMult , Omega_Omegab_Num , TriggerID , Nch;
	float TPVz , p , pt , phi , eta , tEnergy , rap;
	double track_px , track_py , track_pz;
	float trackA_pT , trackA_phi , trackA_eta ;int trackA_charge;
	float trackB_pT , trackB_phi , trackB_eta ;int trackB_charge;
	std::vector<int> PDG , ReCons_TrackID , ParentList , ParentSta , ParentEnd , SE_ParentList , SE_ParentSta , SE_ParentEnd , ME_ParentList , ME_ParentSta , ME_ParentEnd;
	std::vector<float> px,py,pz,InvariantMass,QA_eta,QA_nHitsFit,QA_nHitsMax;
	double zTOF_proton,zTOF_pion,zTOF_kaon;
	// Used for QA
	std::vector<float> QA_dEdx,QA_m2,QA_nSigmaProton,QA_nSigmaPion,QA_nSigmaKaon,QA_Chi2;
	std::vector<double> QA_zTOF_proton,QA_zTOF_pion,QA_zTOF_kaon,QA_Decay_Length,QA_DCA_V0_PV,QA_DCA_Daughters;
	std::vector<int> QA_IfConfuse,QA_IfBadReconstructed;// Used as bool

	std::vector<int> DaughtersID;


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


