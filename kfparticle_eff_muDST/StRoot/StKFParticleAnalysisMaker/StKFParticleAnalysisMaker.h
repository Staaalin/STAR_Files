// $Id: StKFParticleAnalysisMaker.h,v 1.16 2014/08/06 11:43:53 jeromel Exp $
/*!
 * \class  StKFParticleAnalysisMaker
 * \author Maksym Zyzak
 * \date   2017/10/17
 * \brief  class for analysis of PicoDst
 */                                                                      
#ifndef STAR_StKFParticleAnalysisMaker
#define STAR_StKFParticleAnalysisMaker
//#define __DEVT__
#ifndef StMaker_H
#include "StMaker.h"
#endif
#include "TMVA/Reader.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TH2D.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"

class StKFParticleInterface;
class StKFParticlePerformanceInterface;
class KFParticle;
class KFPTrack;
class StDcaGeometry;
class StPicoDst;
class StMuDst;
class StMuEvent;
class TNtuple;
class TFile;
class TH1D;
class TH2D;
class TProfile;
class TVector3;
class TString;
class TEfficiency;
class TChain;
class TDirectory;
class StRefMultCorr;

class StKFParticleAnalysisMaker : public StMaker {
 private:
  Char_t                mBeg[1];        //!
  Char_t                mEnd[1];        //!
  TDirectory           *fDirs[2];       //!
  StMuDst                          *fMuDst;
  StMuEvent                        *fMuEvent;
  StPicoDst                        *fPicoDst;                          //!
  StKFParticleInterface            *fStKFParticleInterface;            //!
  StKFParticlePerformanceInterface *fStKFParticlePerformanceInterface; //!
  bool fIsPicoAnalysis;
  Bool_t fProcessSignal;
  Int_t fCentrality;
  Int_t fMCPID;
  Int_t fPDGPID;
  Int_t fDau1MCPID;
  Int_t fDau2MCPID;
  Int_t fDau1PDGPID;
  Int_t fDau2PDGPID;
  Bool_t fCollectTrackHistograms;
  Bool_t fCollectPIDHistograms;
  StRefMultCorr *fRefmultCorrUtil;
  
  // QA efficiency
  TH1D *hCentrality;
  TH1D *hKFPRecoParPt[9]; // KFP pT of reconstructed lambdas, may include ghosts
  TH1D *hMCRecoParPt[9]; // MC pT of reconstructed lambdas (that are matched)
  TH1D *hMCParPt[9]; // MC pT of embedded lambdas
  TH1D *hMCParPt_2[9]; // MC pT of embedded lambdas, with both daughters' MC tracks present
  TH1D *hVertexDiff;
  TH1D *hLambdaMassVerification;
  TH1D *hMCLambdaMassVerification;
  TH1D *hLambdaPtVerification[4];
  TH1D *hMCDauPxPyPzVerification[3];
  TH1D *hDaughterChi2Pri; 
  TH1D *hDaughterNDF;
  TH1D *hDaughterNDF_2;

  TH2D *hDaughterChi2PriPt_dist;
  TH2D *hDauProtonChi2PriPt_dist; // vs lambda pt
  TH2D *hDauPionChi2PriPt_dist; // vs lambda pt
  TH2D *hDauProtonDCAtoPVPt_dist; // vs lambda pt
  TH2D *hDauPionDCAtoPVPt_dist; // vs lambda pt
  TH2D *hMCDauProtonChi2PriPt_dist; // vs lambda pt, including all MC daughters
  TH2D *hMCDauPionChi2PriPt_dist; // vs lambda pt, including all MC daughters
  TH2D *hMCDauProtonDCAtoPVPt_dist; // vs lambda pt, including all MC daughters
  TH2D *hMCDauPionDCAtoPVPt_dist; // vs lambda pt, including all MC daughters
  TH2D *hPrimaryProtonChi2PriPt_dist; // vs track pt
  TH2D *hPrimaryPionChi2PriPt_dist; // vs track pt
  TH2D *hDauProtonChi2PriSelfPt_dist; // vs track pt
  TH2D *hDauPionChi2PriSelfPt_dist; // vs track pt

  // QA Daughter Efficiency
  TH1D *hRecoDauProtonPt; // MuDst track of reconstructed lambda's daughter
  TH1D *hRecoDauPionPt; // MuDst track of reconstructed lambda's daughter
  TH1D *hMCRecoDauProtonPt; // MuDst track of MC lambda's daughter that is reconstructed
  TH1D *hMCRecoDauPionPt; // MuDst track of MC lambda's daughter that is reconstructed
  TH1D *hMCRecoDauProtonLamPt; // MuDst track of MC lambda's daughter that is reconstructed vs lambda pt
  TH1D *hMCRecoDauPionLamPt; // MuDst track of MC lambda's daughter that is reconstructed vs lambda pt
  TH1D *hMCDauProtonPt; // MC track of MC lambda's daughter
  TH1D *hMCDauPionPt; // MC track of MC lambda's daughter
  TH1D *hMCDauProtonLamPt; // MC track of MC lambda's daughter vs lambda pt
  TH1D *hMCDauPionLamPt; // MC track of MC lambda's daughter vs lambda pt

  TProfile *hDaughterChi2PriDCAtoPV;
  TProfile *hDaughterChi2PriPt;
  TProfile *hDauProtonChi2PriPt;
  TProfile *hDauPionChi2PriPt;
  TProfile *hDaughterChi2PriPt_2;
  TProfile *hDaughterKFPDisPt;
  TProfile *hDaughterDCAPt;
  TProfile *hDauProtonDCAtoPVPt;
  TProfile *hDauPionDCAtoPVPt;
  TProfile *hDauProtonTracknHitsPt;
  TProfile *hDauPionTracknHitsPt;
  TProfile *hDauProtonTracknHitsFitPt;
  TProfile *hDauPionTracknHitsFitPt;
  TProfile *hDauProtonTracknHitsdEdxPt;
  TProfile *hDauPionTracknHitsdEdxPt;
  TProfile *hDauProtonFirstPointDistancePt;
  TProfile *hDauPionFirstPointDistancePt;
  TProfile *hDauProtonLastPointDistancePt;
  TProfile *hDauPionLastPointDistancePt;
  TProfile *hDecayLengthPt;
  TProfile *hDecayVertexPt;
  TProfile *hDecayLengthTraditionalPt;
  TProfile *hDecayVertexTraditionalPt;
  TProfile *hChi2NDFPt;
  TProfile *hChi2NDFPt_primVtx;
  TProfile *hLdLPt;

 public: 
  StKFParticleAnalysisMaker(const char *name="KFParticleAnalysis");
  virtual       ~StKFParticleAnalysisMaker();
  virtual Int_t  Init();
  virtual Int_t  InitRun(Int_t runumber);
  void           BookVertexPlots();
  virtual Int_t  Make();
  virtual Int_t  Finish();
  Bool_t         Check();
  void AnalysePicoDst() { fIsPicoAnalysis = true;  }
  void AnalyseMuDst()   { fIsPicoAnalysis = false; }
  static void    PrintMem(const Char_t *opt = "");
  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $ $Id: StKFParticleAnalysisMaker.h,v 1.0 2017/10/07 11:43:53 mzyzak Exp $ built " __DATE__ " " __TIME__ ; 
    return cvs;
  }
  void ProcessSignal() { fProcessSignal = true; }
  void CollectTrackHistograms() { fCollectTrackHistograms = true; }
  void CollectPIDHistograms() { fCollectPIDHistograms = true; }
  void SetCentrality(Int_t centrality) { fCentrality = centrality; }
  void SetMCPID(Int_t pid);
  void SetPDGPID(Int_t pid);
  void DeclareHistograms();
  void WriteHistograms();
  bool GetTrack(const StDcaGeometry& dcaG, KFPTrack& track, int q, int index);
  TVector3 ConvertStThreeVector(const StThreeVectorF& p) const;
  
  ClassDef(StKFParticleAnalysisMaker,0)   //
};
#endif
// $Log: StKFParticleAnalysisMaker.h,v $
