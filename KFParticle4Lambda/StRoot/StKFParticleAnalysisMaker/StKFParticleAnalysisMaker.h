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
	TH2F *hVertex2D; 
	TH1F *hDiffVz  ; 
	TH1F *hcent;
	TH1F *hcentw;
	TH1F *H_GOOD_Mass;
	TH1F *H_BAD_Mass;

	TH2D *hdEdx_pQ;
	TH2D *hdEdx_pQ_1cut;
	TH2D *hdEdx_pQ_2cut;
	TH2D *hLN_M;
	TH2D *hXY;
	TH2D *hHXY;
	TH2D *hHM_Chi2;

	TProfile *hcentRefM ; 
	TProfile *hcentRefW ; 

	TTree *hadronTree;
	int buffer_size,CrefMult,CgrefMult,evtID,runID,PDGMult , Omega_Omegab_Num , Recorded_events;
	std::vector<int> PDG , ReCons_TrackID;
	std::vector<float> px,py,pz,InvariantMass;
	double zTOF_proton,zTOF_pion,zTOF_kaon;
	// Used for QA
	std::vector<float> QA_dEdx,QA_m2,QA_nSigmaProton,QA_nSigmaPion,QA_nSigmaKaon,QA_Chi2;
	std::vector<double> QA_zTOF_proton,QA_zTOF_pion,QA_zTOF_kaon,QA_Decay_Length,QA_DCA_V0_PV;
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


