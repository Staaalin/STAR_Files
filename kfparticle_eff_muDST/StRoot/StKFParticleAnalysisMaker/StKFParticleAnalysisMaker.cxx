//*-- Author : Yuri Fisyak 02/02/2016
#include "StKFParticleAnalysisMaker.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
//--- KF particle classes ---
#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "KFTopoPerformance.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"
//--- Pico classes ---
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
//--- Mu classes ---
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
//---- Star classes ---
#include "StarClassLibrary/SystemOfUnits.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
//--- TMVA classes ---
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
//--- StRefMult class ---
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
ClassImp(StKFParticleAnalysisMaker);

#define ProtonPDG 2212
#define PionPDG 211
#define ProtonMass 0.938272
#define PionMass 0.139570

//________________________________________________________________________________
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name) : StMaker(name), fIsPicoAnalysis(true),
  fProcessSignal(false), fCentrality(0), fCollectTrackHistograms(false), fCollectPIDHistograms(false), 
  fRefmultCorrUtil(0)
{
  memset(mBeg,0,mEnd-mBeg+1);
}
//________________________________________________________________________________
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker() 
{
  std::cout << "***TEST: We are at the beginning of the StKFPAnalysisMaker destructor***" << std::endl;
  SafeDelete(fStKFParticleInterface);
  SafeDelete(fStKFParticlePerformanceInterface);
  std::cout << "***TEST: We are at the end of the StKFPAnalysisMaker destructor***" << std::endl;
}

//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Init()
{
  // SetCentrality(9);
  DeclareHistograms();

  TFile *f = GetTFile();
  if(f) 
  {
    f->cd();
    BookVertexPlots();
    if(fCollectTrackHistograms)
      fStKFParticleInterface->CollectTrackHistograms();
    if(fCollectPIDHistograms)
      fStKFParticleInterface->CollectPIDHistograms();
  }
      
  fRefmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();

  return kStOK;
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::DeclareHistograms()
{
  // events QA
  hCentrality = new TH1D("hCentrality", "hCentrality", 9, 0.5, 9.5);

  // efficiency
  for (int i = 0; i < 9; i++)
  {
    hKFPRecoParPt[i] = new TH1D(Form("hKFPRecoParPt_%d", i), Form("hKFPRecoParPt_%d", i), 100, 0, 10);
    hMCRecoParPt[i] = new TH1D(Form("hMCRecoParPt_%d", i), Form("hMCRecoParPt_%d", i), 100, 0, 10);
    hMCParPt[i] = new TH1D(Form("hMCParPt_%d", i), Form("hMCParPt_%d", i), 100, 0, 10);
    hMCParPt_2[i] = new TH1D(Form("hMCParPt_2_%d", i), Form("hMCParPt_2_%d", i), 100, 0, 10);
  }

  // KFP QA
  hVertexDiff = new TH1D("hVertexDiff", "hVertexDiff", 1000, -5., 5.);
  hLambdaMassVerification = new TH1D("hLambdaMassVerification", "hLambdaMassVerification", 100, 1.08, 1.16);
  hMCLambdaMassVerification = new TH1D("hMCLambdaMassVerification", "hMCLambdaMassVerification", 100, 1.08, 1.16);
  hLambdaPtVerification[0] = new TH1D("hLambdaPtVerification0", "pt_KFP - pt_MC", 100, -10., 10.);
  hLambdaPtVerification[1] = new TH1D("hLambdaPtVerification1", "pt_KFP_muDstVtx - pt_MC", 100, -10., 10.);
  hLambdaPtVerification[2] = new TH1D("hLambdaPtVerification2", "pt_KFP_KFPVtx - pt_MC", 100, -10., 10.);
  hLambdaPtVerification[3] = new TH1D("hLambdaPtVerification3", "pt_trad - pt_MC", 100, -10., 10.);
  for (int i = 0; i < 3; i++) hMCDauPxPyPzVerification[i] = new TH1D(Form("hMCDauPxPyPzVerification%d", i), Form("hMCDauPxPyPzVerification%d", i), 100, -5., 5.);
  hDaughterChi2Pri = new TH1D("hDaughterChi2", "hDaughterChi2", 1000, -0.5, 9999.5);

  hDaughterChi2PriPt_dist = new TH2D("hDaughterChi2Pt_dist", "hDaughterChi2Pt_dist", 100, 0, 10, 1000, 0., 500.);
  hDauProtonChi2PriPt_dist = new TH2D("hDauProtonChi2PriPt_dist", "hDauProtonChi2PriPt_dist", 100, 0, 10, 1000, 0., 500.);
  hDauPionChi2PriPt_dist = new TH2D("hDauPionChi2PriPt_dist", "hDauPionChi2PriPt_dist", 100, 0, 10, 1000, 0., 500.);
  hDauProtonDCAtoPVPt_dist = new TH2D("hDauProtonDCAtoPVPt_dist", "hDauProtonDCAtoPVPt_dist", 100, 0, 10, 1000, 0., 200.);
  hDauPionDCAtoPVPt_dist = new TH2D("hDauPionDCAtoPVPt_dist", "hDauPionDCAtoPVPt_dist", 100, 0, 10, 1000, 0., 200.);
  hMCDauProtonChi2PriPt_dist = new TH2D("hMCDauProtonChi2PriPt_dist", "hMCDauProtonChi2PriPt_dist", 100, 0, 10, 1000, 0., 500.);
  hMCDauPionChi2PriPt_dist = new TH2D("hMCDauPionChi2PriPt_dist", "hMCDauPionChi2PriPt_dist", 100, 0, 10, 1000, 0., 500.);
  hMCDauProtonDCAtoPVPt_dist = new TH2D("hMCDauProtonDCAtoPVPt_dist", "hMCDauProtonDCAtoPVPt_dist", 100, 0, 10, 1000, 0., 200.);
  hMCDauPionDCAtoPVPt_dist = new TH2D("hMCDauPionDCAtoPVPt_dist", "hMCDauPionDCAtoPVPt_dist", 100, 0, 10, 1000, 0., 200.);
  hPrimaryProtonChi2PriPt_dist = new TH2D("hPrimaryProtonChi2PriPt_dist", "hPrimaryProtonChi2PriPt_dist", 100, 0, 10, 1000, 0., 500.);
  hPrimaryPionChi2PriPt_dist = new TH2D("hPrimaryPionChi2PriPt_dist", "hPrimaryPionChi2PriPt_dist", 100, 0, 10, 1000, 0., 500.);
  hDauProtonChi2PriSelfPt_dist = new TH2D("hDauProtonChi2PriSelfPt_dist", "hDauProtonChi2PriSelfPt_dist", 100, 0, 10, 1000, 0., 500.);
  hDauPionChi2PriSelfPt_dist = new TH2D("hDauPionChi2PriSelfPt_dist", "hDauPionChi2PriSelfPt_dist", 100, 0, 10, 1000, 0., 500.);

  hRecoDauProtonPt = new TH1D("hRecoDauProtonPt", "hRecoDauProtonPt", 100, 0, 10);
  hRecoDauPionPt = new TH1D("hRecoDauPionPt", "hRecoDauPionPt", 100, 0, 10);
  hMCRecoDauProtonPt = new TH1D("hMCRecoDauProtonPt", "hMCRecoDauProtonPt", 100, 0, 10);
  hMCRecoDauPionPt = new TH1D("hMCRecoDauPionPt", "hMCRecoDauPionPt", 100, 0, 10);
  hMCRecoDauProtonLamPt = new TH1D("hMCRecoDauProtonLamPt", "hMCRecoDauProtonLamPt", 100, 0, 10);
  hMCRecoDauPionLamPt = new TH1D("hMCRecoDauPionLamPt", "hMCRecoDauPionLamPt", 100, 0, 10);
  hMCDauProtonPt = new TH1D("hMCDauProtonPt", "hMCDauProtonPt", 100, 0, 10);
  hMCDauPionPt = new TH1D("hMCDauPionPt", "hMCDauPionPt", 100, 0, 10);
  hMCDauProtonLamPt = new TH1D("hMCDauProtonLamPt", "hMCDauProtonLamPt", 100, 0, 10);
  hMCDauPionLamPt = new TH1D("hMCDauPionLamPt", "hMCDauPionLamPt", 100, 0, 10);

  hDauProtonChi2PriPt = new TProfile("hDauProtonChi2PriPt", "hDauProtonChi2PriPt", 100, 0, 10, 0., 100000.);
  hDauPionChi2PriPt = new TProfile("hDauPionChi2PriPt", "hDauPionChi2PriPt", 100, 0, 10, 0., 100000.);
  hDaughterChi2PriDCAtoPV = new TProfile("hDaughterChi2PriDCAtoPV", "hDaughterChi2PriDCAtoPV", 100, 0, 100, 0., 100000.);
  hDaughterNDF = new TH1D("hDaughterNDF", "hDaughterNDF", 20, -9.5, 10.5);
  hDaughterNDF_2 = new TH1D("hDaughterNDF_2", "hDaughterNDF_2", 20, -9.5, 10.5); // reconstruct KFParticle from StMuTrack
  hDaughterChi2PriPt = new TProfile("hDaughterChi2PriPt", "hDaughterChi2PriPt", 100, 0, 10, 0., 100000.);
  hDaughterChi2PriPt_2 = new TProfile("hDaughterChi2PriPt_2", "hDaughterChi2PriPt_2", 100, 0, 10, 0., 100000.); // reconstruct KFParticle from StMuTrack
  hDaughterKFPDisPt = new TProfile("hDaughterKFPDisPt", "hDaughterKFPDisPt", 100, 0, 10, 0., 100.);
  hDaughterDCAPt = new TProfile("hDaughterDCAPt", "hDaughterDCAPt", 100, 0, 10, 0., 100.);
  hDauProtonDCAtoPVPt = new TProfile("hDauProtonDCAtoPVPt", "hDauProtonDCAtoPVPt", 100, 0, 10, 0., 100.);
  hDauPionDCAtoPVPt = new TProfile("hDauPionDCAtoPVPt", "hDauPionDCAtoPVPt", 100, 0, 10, 0., 100.);
  hDauProtonTracknHitsPt = new TProfile("hDauProtonTracknHitsPt", "hDauProtonTracknHitsPt", 100, 0, 10, 0., 100.);
  hDauPionTracknHitsPt = new TProfile("hDauPionTracknHitsPt", "hDauPionTracknHitsPt", 100, 0, 10, 0., 100.);
  hDauProtonTracknHitsFitPt = new TProfile("hDauProtonTracknHitsFitPt", "hDauProtonTracknHitsFitPt", 100, 0, 10, 0., 100.);
  hDauPionTracknHitsFitPt = new TProfile("hDauPionTracknHitsFitPt", "hDauPionTracknHitsFitPt", 100, 0, 10, 0., 100.);
  hDauProtonTracknHitsdEdxPt = new TProfile("hDauProtonTracknHitsdEdxPt", "hDauProtonTracknHitsdEdxPt", 100, 0, 10, 0., 100.);
  hDauPionTracknHitsdEdxPt = new TProfile("hDauPionTracknHitsdEdxPt", "hDauPionTracknHitsdEdxPt", 100, 0, 10, 0., 100.);
  hDauProtonFirstPointDistancePt = new TProfile("hDauProtonFirstPointDistancePt", "hDauProtonFirstPointDistancePt", 100, 0, 10, 0., 210.);
  hDauPionFirstPointDistancePt = new TProfile("hDauPionFirstPointDistancePt", "hDauPionFirstPointDistancePt", 100, 0, 10, 0., 210.);
  hDauProtonLastPointDistancePt = new TProfile("hDauProtonLastPointDistancePt", "hDauProtonLastPointDistancePt", 100, 0, 10, 0., 210.);
  hDauPionLastPointDistancePt = new TProfile("hDauPionLastPointDistancePt", "hDauPionLastPointDistancePt", 100, 0, 10, 0., 210.);
  hDecayLengthPt = new TProfile("hDecayLengthPt", "hDecayLengthPt", 100, 0, 10, 0., 100.);
  hDecayVertexPt = new TProfile("hDecayVertexPt", "hDecayVertexPt", 100, 0, 10, 0., 100.);
  hDecayLengthTraditionalPt = new TProfile("hDecayLengthTraditionalPt", "hDecayLengthTraditionalPt", 100, 0, 10, 0., 100.);
  hDecayVertexTraditionalPt = new TProfile("hDecayVertexTraditionalPt", "hDecayVertexTraditionalPt", 100, 0, 10, 0., 100.);
  hChi2NDFPt = new TProfile("hChi2NDFPt", "hChi2NDFPt", 100, 0, 10, 0., 100.);
  hChi2NDFPt_primVtx = new TProfile("hChi2NDFPt_primVtx", "hChi2NDFPt_primVtx", 100, 0, 10, 0., 100.);
  hLdLPt = new TProfile("hLdLPt", "hLdLPt", 100, 0, 10, 0., 100.);
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::WriteHistograms()
{
  hCentrality->Write();

  for (int i = 0; i < 9; i++)
  {
    hKFPRecoParPt[i]->Write();
    hMCRecoParPt[i]->Write();
    hMCParPt[i]->Write();
    hMCParPt_2[i]->Write();
  }

  hVertexDiff->Write();
  hLambdaMassVerification->Write();
  hMCLambdaMassVerification->Write();
  for(Int_t i=0; i<4; i++)
    hLambdaPtVerification[i]->Write();
  for (int i=0; i<3; i++)
    hMCDauPxPyPzVerification[i]->Write();
  hDaughterChi2Pri->Write();

  hDaughterChi2PriPt_dist->Write();
  hDauProtonChi2PriPt_dist->Write();
  hDauPionChi2PriPt_dist->Write();
  hDauProtonDCAtoPVPt_dist->Write();
  hDauPionDCAtoPVPt_dist->Write();
  hMCDauProtonChi2PriPt_dist->Write();
  hMCDauPionChi2PriPt_dist->Write();
  hMCDauProtonDCAtoPVPt_dist->Write();
  hMCDauPionDCAtoPVPt_dist->Write();
  hPrimaryProtonChi2PriPt_dist->Write();
  hPrimaryPionChi2PriPt_dist->Write();
  hDauProtonChi2PriSelfPt_dist->Write();
  hDauPionChi2PriSelfPt_dist->Write();

  hRecoDauProtonPt->Write();
  hRecoDauPionPt->Write();
  hMCRecoDauProtonPt->Write();
  hMCRecoDauPionPt->Write();
  hMCRecoDauProtonLamPt->Write();
  hMCRecoDauPionLamPt->Write();
  hMCDauProtonPt->Write();
  hMCDauPionPt->Write();
  hMCDauProtonLamPt->Write();
  hMCDauPionLamPt->Write();

  hDauProtonChi2PriPt->Write();
  hDauPionChi2PriPt->Write();
  hDaughterChi2PriDCAtoPV->Write();
  hDaughterNDF->Write();
  hDaughterNDF_2->Write();
  hDaughterChi2PriPt->Write();
  hDaughterChi2PriPt_2->Write();
  hDaughterKFPDisPt->Write();
  hDaughterDCAPt->Write();
  hDauProtonDCAtoPVPt->Write();
  hDauPionDCAtoPVPt->Write();
  hDauProtonTracknHitsPt->Write();
  hDauPionTracknHitsPt->Write();
  hDauProtonTracknHitsFitPt->Write();
  hDauPionTracknHitsFitPt->Write();
  hDauProtonTracknHitsdEdxPt->Write();
  hDauPionTracknHitsdEdxPt->Write();
  hDauProtonFirstPointDistancePt->Write();
  hDauPionFirstPointDistancePt->Write();
  hDauProtonLastPointDistancePt->Write();
  hDauPionLastPointDistancePt->Write();
  hDecayLengthPt->Write();
  hDecayVertexPt->Write();
  hDecayLengthTraditionalPt->Write();
  hDecayVertexTraditionalPt->Write();
  hChi2NDFPt->Write();
  hChi2NDFPt_primVtx->Write();
  hLdLPt->Write();
}
//________________________________________________________________________________
Int_t StKFParticleAnalysisMaker::InitRun(Int_t runumber) 
{
  return StMaker::InitRun(runumber);
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::PrintMem(const Char_t *opt)
{
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  cout << opt 
       << "\tMemory : Total = " << info.fMemTotal 
       << "\tUsed = " << info.fMemUsed
       << "\tFree = " << info.fMemFree
       << "\tSwap Total = " << info.fSwapTotal
       << "\tUsed = " << info.fSwapUsed
       << "\tFree = " << info.fSwapFree << endl;
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::BookVertexPlots()
{
  // fDirs[2] = {0};
  fDirs[0] = TDirectory::CurrentDirectory(); assert(fDirs[0]);
  fDirs[0]->cd();
  if (! fDirs[0]->GetDirectory("Particles")) {
    fDirs[0]->mkdir("Particles");
  }
  fDirs[1] = fDirs[0]->GetDirectory("Particles"); assert(fDirs[1]);
  fDirs[1]->cd();
  PrintMem(fDirs[1]->GetPath());
  
  fStKFParticleInterface = new StKFParticleInterface;
  bool storeMCHistograms = false;
  if(!fIsPicoAnalysis && fProcessSignal) storeMCHistograms = true;
  fStKFParticlePerformanceInterface = new StKFParticlePerformanceInterface(fStKFParticleInterface->GetTopoReconstructor(), storeMCHistograms);
  fDirs[0]->cd();
  PrintMem(fDirs[1]->GetPath());
}
//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Make()
{  
  if(fIsPicoAnalysis)
  {
    fPicoDst = StPicoDst::instance();
    if(!fPicoDst) return kStOK;
  }
  else
  {  
    fMuDst = StMuDst::instance();
    if(!fMuDst) return kStOK;
    else { if(StMuDst::instance()->numberOfPrimaryVertices() == 0 ) return kStOK; }
    fMuEvent = fMuDst->event();
    if (!fMuEvent) return kStOK;
  }

  // vertex selection
  int const originalVertexId = fMuDst->currentVertexIndex();
  StMuPrimaryVertex *selectedVertex = nullptr;
  // choose the default vertex, i.e. the first vertex
  fMuDst->setVertexIndex(0);
  selectedVertex = fMuDst->primaryVertex();
  // fall back to default vertex if no vertex is selected in the algorithm above.
  // should skip this event in the event cuts below.
  if (!selectedVertex)
  {
      LOG_INFO << "Vertex is not valid" << endm;
      // cout<<originalVertexId<<endl;
      fMuDst->setVertexIndex(originalVertexId);
  }

  // event cuts
  if (!selectedVertex) return kStOK;
  if (fabs(fMuEvent->primaryVertexPosition().z()) > 70.0) return kStOK;
  if (fMuEvent->primaryVertexPosition().perp() > 2.0) return kStOK;

  // Centrality
  int Run = fMuEvent->runId();
  fRefmultCorrUtil->init(Run);
  if (fRefmultCorrUtil->isBadRun(Run)) return kStOK;
  fRefmultCorrUtil->initEvent(fMuEvent->refMult(), fMuEvent->primaryVertexPosition().z(), 0);
  int Centrality = 1 + fRefmultCorrUtil->getCentralityBin9();
  if (fCentrality && fCentrality != Centrality) return kStOK;

  TClonesArray *MuMcVertices = fMuDst->mcArray(0);
  Int_t NoMuMcVertices = MuMcVertices->GetEntriesFast();
  TClonesArray *MuMcTracks = fMuDst->mcArray(1);
  Int_t NoMuMcTracks = MuMcTracks->GetEntriesFast();
  // LOG_INFO << "# of MC tracks = " << NoMuMcTracks << " # of MC vertices = " << NoMuMcVertices << endm;
  if (!NoMuMcVertices || !NoMuMcTracks) return kStOK;

  //find max global track index
  int maxGBTrackIndex = -1;
  if(fIsPicoAnalysis)
  {
    for(unsigned int iTrack = 0; iTrack < fPicoDst->numberOfTracks(); iTrack++) 
    {
      StPicoTrack *gTrack = fPicoDst->track(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  else
  {
    for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfGlobalTracks(); iTrack++) 
    {
      StMuTrack *gTrack = fMuDst->globalTracks(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  vector<KFMCTrack> mcTracks(0);
  vector<int> mcIndices(maxGBTrackIndex+1);
  for(unsigned int iIndex=0; iIndex<mcIndices.size(); iIndex++)
    mcIndices[iIndex] = -1;
  
//   fStKFParticleInterface->SetTriggerMode();
//   fStKFParticleInterface->SetSoftKaonPIDMode();
//   fStKFParticleInterface->SetSoftTofPidMode();
//   fStKFParticleInterface->SetChiPrimaryCut(10);
//   
//   fStKFParticleInterface->SetPtCutCharm(0.5);
//   fStKFParticleInterface->SetChiPrimaryCutCharm(8);
//   fStKFParticleInterface->SetLdLCutCharmManybodyDecays(3);
//   fStKFParticleInterface->SetChi2TopoCutCharmManybodyDecays(10);
//   fStKFParticleInterface->SetChi2CutCharmManybodyDecays(3);
//   fStKFParticleInterface->SetLdLCutCharm2D(3);
//   fStKFParticleInterface->SetChi2TopoCutCharm2D(10);
//   fStKFParticleInterface->SetChi2CutCharm2D(3);
  
  vector<int> triggeredTracks;
  bool isGoodEvent = false;
  
  //Process the event
  if(maxGBTrackIndex > 0)
    fStKFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex+1);
  if(fIsPicoAnalysis)
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fPicoDst, triggeredTracks);
  else
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fMuDst, mcTracks, mcIndices, fProcessSignal);

//   bool openCharmTrigger = false;
//   if(isGoodEvent) openCharmTrigger =  fStKFParticleInterface->OpenCharmTrigger();
//   fStKFParticleInterface->OpenCharmTriggerCompression(triggeredTracks.size(), fPicoDst->numberOfTracks(), openCharmTrigger);
  //collect histograms
  
  if(isGoodEvent)
  {
    hCentrality->Fill(Centrality*1.0);
    int centralityBin = -1;
    if (fCentrality) centralityBin = fCentrality;
    float centralityWeight = 1.;
    
    fStKFParticlePerformanceInterface->SetMCTracks(mcTracks);
    fStKFParticlePerformanceInterface->SetMCIndexes(mcIndices);    
    fStKFParticlePerformanceInterface->SetCentralityBin(centralityBin);
    fStKFParticlePerformanceInterface->SetCentralityWeight(centralityWeight);
    Int_t nevent = 100000;
    fStKFParticlePerformanceInterface->SetPrintEffFrequency(nevent);
    fStKFParticlePerformanceInterface->PerformanceAnalysis();

    // efficiency QA -- MC
    for (Int_t itrk = 0; itrk < NoMuMcTracks; itrk++)
    {
      StMuMcTrack *mcTrack = (StMuMcTrack *)MuMcTracks->UncheckedAt(itrk);
      if (!mcTrack) continue;

      // // Select only Triggered Mc Vertex, i.e. the MC track should originate from PV (IdVx=1)
      // Int_t IdVx = mcTrack->IdVx();
      // if (IdVx != 1) continue;

      const int Gid = mcTrack->CorrectGePid(mcTrack->GePid()); // need this for correct lambda PID
      if (Gid == fMCPID)
      {

        // add MC cuts
        //		if (mcTrack->Pxyz().perp()<0.2 ||mcTrack->Pxyz().perp()>2  ) continue;
        //         	if (fabs(mcTrack->Pxyz().pseudoRapidity())<1) continue;
        //         	if (mcTrack->NoHits() < 15) continue;
        hMCParPt[Centrality-1]->Fill(mcTrack->Pxyz().perp());

        // below does not work, topoPerf cannot be initialized as a singleton
        // KFTopoPerformance *topoPerf = KFTopoPerformance::instance();
        // bool matched = false;
        // for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
        // {
        //   KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
        //   if (particle.GetPDG() != fPDGPID) continue;
        //   if (topoPerf->ParticlesMatch()[iParticle].GetBestMatch() == itrk)
        //   {
        //     matched = true;
        //     break;
        //   }
        // }
        // hEfficiencyPt->Fill(matched, mcTrack->Pxyz().perp());

      }
    }

    // track ID mapping
    std::vector<int> trackIdtoIndex; // track->id() to index (0 to nTracks-1)
    std::vector<int> trackIndextoId; // index (0 to nTracks-1) to track->id()
    int maxId = 0;
    for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfGlobalTracks(); iTrack++) 
    {
      StMuTrack *gTrack = fMuDst->globalTracks(iTrack);
      if (! gTrack) continue;
      int trackId = gTrack->id();
      if(trackId > maxId)
        maxId = trackId;
    }
    trackIdtoIndex.resize(maxId+1, -1);
    trackIndextoId.resize(fMuDst->numberOfGlobalTracks(), -1);
    for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfGlobalTracks(); iTrack++) 
    {
      StMuTrack *gTrack = fMuDst->globalTracks(iTrack);
      if (! gTrack) continue;
      int trackId = gTrack->id();
      trackIdtoIndex[trackId] = iTrack;
      trackIndextoId[iTrack] = trackId;
    }

    // efficiency QA -- Reco
    const KFParticle KFP_vtx = fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex(0);
    KFPVertex temp;
    temp.SetXYZ(fMuEvent->primaryVertexPosition().x(), fMuEvent->primaryVertexPosition().y(), fMuEvent->primaryVertexPosition().z());
    float cv[6] = {fMuEvent->primaryVertexErrors().x()*fMuEvent->primaryVertexErrors().x(), 0, 
                   fMuEvent->primaryVertexErrors().y()*fMuEvent->primaryVertexErrors().y(), 0, 0,
                   fMuEvent->primaryVertexErrors().z()*fMuEvent->primaryVertexErrors().z()}; // cxx, cxy, cyy, cxz, cyz, czz
    temp.SetCovarianceMatrix(cv);
    temp.SetNContributors(selectedVertex->noTracks());
    temp.SetChi2(selectedVertex->chiSquared());
    KFVertex KFP_vtx_muDst(temp);
    // LOG_INFO << "Number of KFTopoReconstructor Vertices: " << fStKFParticleInterface->GetTopoReconstructor()->NPrimaryVertices() << endm;
    TVector3 vtx_toporeco(KFP_vtx.X(), KFP_vtx.Y(), KFP_vtx.Z());
    TVector3 vtx_muDst = ConvertStThreeVector(fMuEvent->primaryVertexPosition());
    hVertexDiff->Fill((vtx_toporeco - vtx_muDst).Mag());

    // check primary tracks
    for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
    {
      KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
      // primary particles
      if (particle.DaughterIds().size() != 1) continue; // is a track
      int trackId = particle.DaughterIds()[0];
      if (trackIdtoIndex[trackId] == -1) continue; // track not found
      StMuTrack *gTrack = fMuDst->globalTracks(trackIdtoIndex[trackId]);
      if (! gTrack) continue; 
      if (gTrack->dcaGlobal().mag() > 3) continue; // primary
      float pt = gTrack->pt();
      float chi2 = particle.GetDeviationFromVertex(KFP_vtx);
      if (particle.GetPDG() == fDau1PDGPID) hPrimaryProtonChi2PriPt_dist->Fill(pt, chi2);
      if (particle.GetPDG() == fDau2PDGPID) hPrimaryPionChi2PriPt_dist->Fill(pt, chi2);
    }

    // check MC daughter eff
    for (Int_t itrk = 0; itrk < NoMuMcTracks; itrk++)
    {
      StMuMcTrack *mcTrack = (StMuMcTrack *)MuMcTracks->UncheckedAt(itrk);
      if (!mcTrack) continue;
      if (mcTrack->GePid() != fDau1MCPID && mcTrack->GePid() != fDau2MCPID) continue;
      StMuMcVertex *mcVertex = (StMuMcVertex*)MuMcVertices->UncheckedAt(mcTrack->IdVx()-1);
      if (!mcVertex) continue;
      int idParent = mcVertex->IdParTrk();
      if (idParent <= 0 || idParent > NoMuMcTracks) continue;
      StMuMcTrack *parentMcTrack = (StMuMcTrack*)MuMcTracks->UncheckedAt(idParent-1);
      if (!parentMcTrack) continue;
      if (parentMcTrack->CorrectGePid(parentMcTrack->GePid()) != fMCPID) continue;

      float pt_dau_MC = mcTrack->Pxyz().perp();
      float pt_lam_MC = parentMcTrack->Pxyz().perp();
      if (mcTrack->GePid() == fDau1MCPID) {hMCDauProtonPt->Fill(pt_dau_MC); hMCDauProtonLamPt->Fill(pt_lam_MC);}
      if (mcTrack->GePid() == fDau2MCPID) {hMCDauPionPt  ->Fill(pt_dau_MC); hMCDauPionLamPt  ->Fill(pt_lam_MC);}

      // check if MC daughter is reconstructed
      StMuTrack *gTrack = 0;
      for (unsigned int iTrack = 0; iTrack < fMuDst->numberOfGlobalTracks(); iTrack++)
      {
        StMuTrack *gTrack_temp = fMuDst->globalTracks(iTrack);
        if (! gTrack_temp) continue;
        if (gTrack_temp->idTruth() == itrk+1) {gTrack = gTrack_temp; break;}
      }
      if (!gTrack) continue;
      if (mcTrack->GePid() == fDau1MCPID) {hMCRecoDauProtonPt->Fill(pt_dau_MC); hMCRecoDauProtonLamPt->Fill(pt_lam_MC);}
      if (mcTrack->GePid() == fDau2MCPID) {hMCRecoDauPionPt  ->Fill(pt_dau_MC); hMCRecoDauPionLamPt  ->Fill(pt_lam_MC);}
    }

    // check all MC daughters
    for (Int_t itrk = 0; itrk < NoMuMcTracks; itrk++)
    {
      for (Int_t jtrk = itrk+1; jtrk < NoMuMcTracks; jtrk++)
      {
        StMuMcTrack *p1McTrack = (StMuMcTrack *)MuMcTracks->UncheckedAt(itrk);
        StMuMcTrack *p2McTrack = (StMuMcTrack *)MuMcTracks->UncheckedAt(jtrk);
        if (!p1McTrack || !p2McTrack) continue;
        if ((p1McTrack->GePid() != fDau1MCPID || p2McTrack->GePid() != fDau2MCPID) && (p1McTrack->GePid() != fDau2MCPID || p2McTrack->GePid() != fDau1MCPID)) continue;
        StMuMcTrack *protonMcTrack = p1McTrack->GePid() == fDau1MCPID ? p1McTrack : p2McTrack;
        StMuMcTrack *pionMcTrack = p1McTrack->GePid() == fDau2MCPID ? p1McTrack : p2McTrack;
        int protonIndex = p1McTrack->GePid() == fDau1MCPID ? itrk : jtrk;
        int pionIndex = p1McTrack->GePid() == fDau2MCPID ? itrk : jtrk;
        if (protonMcTrack->IdVx() != pionMcTrack->IdVx()) continue; // same vertex, could be a short-lived particle
        StMuMcVertex *mcVertex = (StMuMcVertex*)MuMcVertices->UncheckedAt(protonMcTrack->IdVx()-1);
        if (!mcVertex) continue;
        int idParent = mcVertex->IdParTrk();
        if (idParent <= 0 || idParent > NoMuMcTracks) continue;
        StMuMcTrack *parentMcTrack = (StMuMcTrack*)MuMcTracks->UncheckedAt(idParent-1);
        if (!parentMcTrack) continue;
        if (parentMcTrack->CorrectGePid(parentMcTrack->GePid()) != fMCPID) continue;
        float pt_MC = parentMcTrack->Pxyz().perp();
        hMCParPt_2[Centrality-1]->Fill(pt_MC);
        // hMCDauProtonPt->Fill(protonMcTrack->Pxyz().perp());
        // hMCDauPionPt->Fill(pionMcTrack->Pxyz().perp());      

        // find original track
        StMuTrack *protonTrack = 0, *pionTrack = 0;
        for (Int_t ktrk = 0; ktrk < fMuDst->numberOfGlobalTracks(); ktrk++)
        {
          StMuTrack *gTrack = fMuDst->globalTracks(ktrk);
          if (!gTrack) continue;
          if (gTrack->idTruth() == protonIndex+1) protonTrack = gTrack;
          else if (gTrack->idTruth() == pionIndex+1) pionTrack = gTrack;
        }
        if (!protonTrack || !pionTrack) continue;
        // hMCRecoDauProtonPt->Fill(protonMcTrack->Pxyz().perp()); // or protonTrack->pt()?
        // hMCRecoDauPionPt->Fill(pionMcTrack->Pxyz().perp());

        // Check lambda mass
        TLorentzVector proton4vec, pion4vec;
        StPhysicalHelixD helix1 = protonTrack->helix();
        StPhysicalHelixD helix2 = pionTrack->helix();
        pair<double,double> ss = helix1.pathLengths(helix2);
        TVector3 p1 = ConvertStThreeVector(helix1.momentumAt(ss.first, fMuEvent->magneticField()*kilogauss));
        TVector3 p2 = ConvertStThreeVector(helix2.momentumAt(ss.second, fMuEvent->magneticField()*kilogauss));
        proton4vec.SetVectM(p1, ProtonMass);
        pion4vec.SetVectM(p2, PionMass);
        hMCLambdaMassVerification->Fill((proton4vec+pion4vec).M());

        // find chi2_prim
        const int proton_index = protonTrack->id();
        const int pion_index = pionTrack->id();
        float pt_proton_self = protonTrack->pt();
        float pt_pion_self = pionTrack->pt();
        Int_t proton_dcaGeometryIndex = protonTrack->index2Cov();    
        Int_t pion_dcaGeometryIndex = pionTrack->index2Cov(); 
        if (proton_dcaGeometryIndex < 0 || pion_dcaGeometryIndex < 0) continue; // no dca geometry
        StDcaGeometry *dcaG_proton = fMuDst->covGlobTracks(proton_dcaGeometryIndex);
        StDcaGeometry *dcaG_pion = fMuDst->covGlobTracks(pion_dcaGeometryIndex);
        StThreeVectorD V(fMuDst->primaryVertex()->position());
        if (!dcaG_proton || !dcaG_pion) continue;
        THelixTrack t_proton = dcaG_proton->thelix();
        THelixTrack t_pion = dcaG_pion->thelix();
        Double_t dca3D_proton = t_proton.Dca(V.xyz());
        Double_t dca3D_pion = t_pion.Dca(V.xyz());
        hMCDauProtonDCAtoPVPt_dist->Fill(pt_MC, dca3D_proton);
        hMCDauPionDCAtoPVPt_dist->Fill(pt_MC, dca3D_pion);

        KFPTrack proton_kfpTrack, pion_kfpTrack;
        int proton_charge = protonTrack->charge();
        int pion_charge = pionTrack->charge();
        if (!GetTrack(*dcaG_proton, proton_kfpTrack, proton_charge, proton_index)) continue;
        if (!GetTrack(*dcaG_pion, pion_kfpTrack, pion_charge, pion_index)) continue;
        KFParticle protonKFP(proton_kfpTrack, fDau1PDGPID);
        KFParticle pionKFP(pion_kfpTrack, fDau2PDGPID);
        hMCDauProtonChi2PriPt_dist->Fill(pt_MC, protonKFP.GetDeviationFromVertex(KFP_vtx));
        hMCDauPionChi2PriPt_dist->Fill(pt_MC, pionKFP.GetDeviationFromVertex(KFP_vtx));

        // verify KFParticle
        hMCDauPxPyPzVerification[0]->Fill(protonTrack->momentum().x() - protonKFP.GetPx());
        hMCDauPxPyPzVerification[1]->Fill(protonTrack->momentum().y() - protonKFP.GetPy());
        hMCDauPxPyPzVerification[2]->Fill(protonTrack->momentum().z() - protonKFP.GetPz());
        hMCDauPxPyPzVerification[0]->Fill(pionTrack->momentum().x() - pionKFP.GetPx());
        hMCDauPxPyPzVerification[1]->Fill(pionTrack->momentum().y() - pionKFP.GetPy());
        hMCDauPxPyPzVerification[2]->Fill(pionTrack->momentum().z() - pionKFP.GetPz());
      }
    }

    // check reconstructed MC
    for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
    {
      KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
      if (particle.GetPDG() != fPDGPID) continue;
      if (particle.NDaughters() != 2) continue;
      const float pt_KFP = particle.GetPt();
      hKFPRecoParPt[Centrality-1]->Fill(pt_KFP); // Is KFP Pt correct?

      // finding the corresponding track
      KFParticle dau1, dau2, proton, pion;
      dau1 = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]]; 
      dau2 = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
      if (!(dau1.GetPDG() == fDau1PDGPID && dau2.GetPDG() == fDau2PDGPID) &&
          !(dau1.GetPDG() == fDau2PDGPID && dau2.GetPDG() == fDau1PDGPID)) continue;
      proton = dau1.GetPDG() == fDau1PDGPID ? dau1 : dau2;
      pion   = dau1.GetPDG() == fDau1PDGPID ? dau2 : dau1;
      if (trackIdtoIndex[proton.DaughterIds()[0]] == -1 || trackIdtoIndex[pion.DaughterIds()[0]] == -1) continue;
      StMuTrack *protonTrack = fMuDst->globalTracks(trackIdtoIndex[proton.DaughterIds()[0]]);
      StMuTrack *pionTrack   = fMuDst->globalTracks(trackIdtoIndex[pion.DaughterIds()[0]]);
      if (!protonTrack || !pionTrack) continue;
      
      // check if tracks are matched
      StMuMcTrack *protonMcTrack = (StMuMcTrack*)MuMcTracks->UncheckedAt(protonTrack->idTruth()-1);
      StMuMcTrack *pionMcTrack   = (StMuMcTrack*)MuMcTracks->UncheckedAt(pionTrack->idTruth()-1);
      if (!protonMcTrack || !pionMcTrack) continue;
      if (protonMcTrack->GePid() != fDau1MCPID || pionMcTrack->GePid() != fDau2MCPID) continue;
      if (protonMcTrack->IdVx() != pionMcTrack->IdVx()) continue; // same vertex, could be a short-lived particle
      StMuMcVertex *mcVertex = (StMuMcVertex*)MuMcVertices->UncheckedAt(protonMcTrack->IdVx()-1);
      if (!mcVertex) continue;
      int idParent = mcVertex->IdParTrk();
      if (idParent <= 0 || idParent > NoMuMcTracks) continue;
      StMuMcTrack *parentMcTrack = (StMuMcTrack*)MuMcTracks->UncheckedAt(idParent-1);
      if (!parentMcTrack) continue;
      if (parentMcTrack->CorrectGePid(parentMcTrack->GePid()) != fMCPID) continue;
      const float pt_MC = parentMcTrack->Pxyz().perp();
      hMCRecoParPt[Centrality-1]->Fill(pt_MC); 

      /* All cuts are done */

      // Daughter QA
      hRecoDauProtonPt->Fill(protonTrack->pt());
      hRecoDauPionPt->Fill(pionTrack->pt());
      hDaughterKFPDisPt->Fill(pt_MC, dau1.GetDistanceFromParticle(dau2)); // DCA between daughters
      for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
      {
        KFParticle daughter = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[iDaughter]];
        float KFP_vtx_num[3] = {KFP_vtx.X(), KFP_vtx.Y(), KFP_vtx.Z()};
        daughter.TransportToPoint(KFP_vtx_num); // has no effect
        float chi2_dau = daughter.GetDeviationFromVertex(KFP_vtx);
        hDaughterChi2PriPt->Fill(pt_MC, chi2_dau); // chi2 of daughters to primary vertex
        hDaughterChi2PriPt_dist->Fill(pt_MC, chi2_dau); // as a histogram
        if (daughter.GetPDG() == fDau1PDGPID)
        {
          hDauProtonChi2PriPt->Fill(pt_MC, chi2_dau); // chi2 of proton to primary vertex
          hDauProtonChi2PriPt_dist->Fill(pt_MC, chi2_dau); // as a histogram
        }
        else if (daughter.GetPDG() == fDau2PDGPID)
        {
          hDauPionChi2PriPt->Fill(pt_MC, chi2_dau); // chi2 of pion to primary vertex
          hDauPionChi2PriPt_dist->Fill(pt_MC, chi2_dau); // as a histogram
        }
        hDaughterNDF->Fill(daughter.GetNDF()); // NDF of daughters

        StMuTrack *dauTrack = fMuDst->globalTracks(trackIdtoIndex[daughter.DaughterIds()[0]]);   
        const int index = dauTrack->id();
        float pt_self = dauTrack->pt();
        Int_t dcaGeometryIndex = dauTrack->index2Cov();     
        StDcaGeometry *dcaG = fMuDst->covGlobTracks(dcaGeometryIndex);
        if (dcaG)
        {
          Int_t q = 1; if (dauTrack->charge() < 0) q = -1;
          KFPTrack kfpTrack;
          GetTrack(*dcaG, kfpTrack, q, index);
          KFParticle daughterPDG(kfpTrack, daughter.GetPDG());
          hDaughterChi2PriPt_2->Fill(pt_MC, daughterPDG.GetDeviationFromVertex(KFP_vtx)); 
          hDaughterChi2Pri->Fill(daughterPDG.GetDeviationFromVertex(KFP_vtx)); 
          hDaughterNDF_2->Fill(daughterPDG.GetNDF()); // NDF of daughters

          THelixTrack t = dcaG->thelix();
          StThreeVectorD V(fMuDst->primaryVertex()->position());
          Double_t dca3D = t.Dca(V.xyz());
          if (daughter.GetPDG() == fDau1PDGPID) hDauProtonDCAtoPVPt_dist->Fill(pt_MC, dca3D); 
          else if (daughter.GetPDG() == fDau2PDGPID) hDauPionDCAtoPVPt_dist->Fill(pt_MC, dca3D);
        }

        TVector3 firstPoint = ConvertStThreeVector(dauTrack->firstPoint());
        TVector3 lastPoint = ConvertStThreeVector(dauTrack->lastPoint());
        if      (daughter.GetPDG() == fDau1PDGPID) 
        {
          hDauProtonChi2PriSelfPt_dist->Fill(pt_self, chi2_dau); 
          hDauProtonDCAtoPVPt->Fill(pt_MC, dauTrack->dcaGlobal().mag()); // DCA of daughter to primary vertex
          hDauProtonTracknHitsPt->Fill(pt_MC, dauTrack->nHits()); 
          hDauProtonTracknHitsFitPt->Fill(pt_MC, dauTrack->nHitsFit()); 
          hDauProtonTracknHitsdEdxPt->Fill(pt_MC, dauTrack->nHitsDedx()); 
          hDauProtonFirstPointDistancePt->Fill(pt_MC, (firstPoint-vtx_muDst).Mag()); // distance between first point and primary vertex
          hDauProtonLastPointDistancePt->Fill(pt_MC, (lastPoint-vtx_muDst).Mag()); // distance between last point and primary vertex
        }
        else if (daughter.GetPDG() == fDau2PDGPID) 
        {
          hDauPionChi2PriSelfPt_dist->Fill(pt_self, chi2_dau); 
          hDauPionDCAtoPVPt->Fill(pt_MC, dauTrack->dcaGlobal().mag()); // DCA of daughter to primary vertex
          hDauPionTracknHitsPt->Fill(pt_MC, dauTrack->nHits());
          hDauPionTracknHitsFitPt->Fill(pt_MC, dauTrack->nHitsFit());
          hDauPionTracknHitsdEdxPt->Fill(pt_MC, dauTrack->nHitsDedx());
          hDauPionFirstPointDistancePt->Fill(pt_MC, (firstPoint-vtx_muDst).Mag()); // distance between first point and primary vertex
          hDauPionLastPointDistancePt->Fill(pt_MC, (lastPoint-vtx_muDst).Mag()); // distance between last point and primary vertex
        }
        hDaughterChi2PriDCAtoPV->Fill(dauTrack->dcaGlobal().mag(), chi2_dau); // DCA of daughter to primary vertex vs chi2
      }

      // Lambda QA
      StPhysicalHelixD helix1 = protonTrack->helix();
      StPhysicalHelixD helix2 = pionTrack->helix();
      pair<double,double> ss = helix1.pathLengths(helix2);
      TVector3 x1 = ConvertStThreeVector(helix1.at(ss.first));
      TVector3 x2 = ConvertStThreeVector(helix2.at(ss.second));
      TVector3 p1 = ConvertStThreeVector(helix1.momentumAt(ss.first, fMuEvent->magneticField()*kilogauss));
      TVector3 p2 = ConvertStThreeVector(helix2.momentumAt(ss.second, fMuEvent->magneticField()*kilogauss));
      const float pt_trad = (p1+p2).Perp();
      hDaughterDCAPt->Fill(pt_MC, (x1-x2).Mag()); // DCA between daughters

      TVector3 xv0 = (x1 + x2) * 0.5;
      hDecayLengthTraditionalPt->Fill(pt_MC, (xv0 - vtx_muDst).Mag()); // distance between decay vertex and primary vertex
      hDecayVertexTraditionalPt->Fill(pt_MC, xv0.Mag()); // distance between decay vertex and (0,0,0)
      TLorentzVector lv1, lv2;
      lv1.SetVectM(p1, ProtonMass);
      lv2.SetVectM(p2, PionMass);
      hLambdaMassVerification->Fill((lv1+lv2).M());
      
      hChi2NDFPt->Fill(pt_MC, particle.GetChi2()/particle.GetNDF()*1.0); // chi2/ndf of reconstructed particle
      hDecayVertexPt->Fill(pt_MC, TVector3(particle.GetX(), particle.GetY(), particle.GetZ()).Mag()); // position of reconstructed particle
      particle.SetProductionVertex(KFP_vtx);
      const float pt_KFP_vtx = particle.GetPt();
      hChi2NDFPt_primVtx->Fill(pt_MC, particle.GetChi2()/particle.GetNDF()*1.0); // chi2/ndf of reconstructed particle, after setting primary vertex
      float l, dl;
      particle.GetDecayLength(l, dl);
      hDecayLengthPt->Fill(pt_MC, l); // decay length of reconstructed particle
      hLdLPt->Fill(pt_MC, l/dl); // decay length significance of reconstructed particle

      particle.SetProductionVertex(KFP_vtx_muDst);
      const float pt_KFP_muDst = particle.GetPt();
      
      // fill pt diff histograms
      hLambdaPtVerification[0]->Fill(pt_KFP - pt_MC);
      hLambdaPtVerification[1]->Fill(pt_KFP_muDst - pt_MC);
      hLambdaPtVerification[2]->Fill(pt_KFP_vtx - pt_MC);
      hLambdaPtVerification[3]->Fill(pt_trad - pt_MC);
    }
  }
  
  return kStOK;
}

void StKFParticleAnalysisMaker::SetMCPID(Int_t pid)
{
  fMCPID = pid;
  if (pid == 18) {fDau1MCPID = 14; fDau2MCPID = 9;}
  if (pid == 26) {fDau1MCPID = 15; fDau2MCPID = 8;} 
}

void StKFParticleAnalysisMaker::SetPDGPID(Int_t pid)
{
  fPDGPID = pid;
  if (pid ==  3122) {fDau1PDGPID = 2212; fDau2PDGPID = -211;}
  if (pid == -3122) {fDau1PDGPID = -2212; fDau2PDGPID = 211;}
}

TVector3 StKFParticleAnalysisMaker::ConvertStThreeVector(const StThreeVectorF& p) const
{
  return TVector3(p.x(), p.y(), p.z());
}

bool StKFParticleAnalysisMaker::GetTrack(const StDcaGeometry& dcaG, KFPTrack& track, int q, int index)
{
  Double_t xyzp[6], CovXyzp[21];
  dcaG.GetXYZ(xyzp,CovXyzp);
  
  bool goodTrack=1;
  for(int iPar=0; iPar<6; iPar++)
    goodTrack = goodTrack && finite(xyzp[iPar]);
  for(int iC=0; iC<21; iC++)
    goodTrack = goodTrack && finite(CovXyzp[iC]);
  goodTrack &= goodTrack && CovXyzp[0]  >=0.f && CovXyzp[0]  < 100.f;
  goodTrack &= goodTrack && CovXyzp[2]  >=0.f && CovXyzp[2]  < 100.f;
  goodTrack &= goodTrack && CovXyzp[5]  >=0.f && CovXyzp[5]  < 100.f;
  goodTrack &= goodTrack && CovXyzp[9]  >=0.f && CovXyzp[9]  < 1.f;
  goodTrack &= goodTrack && CovXyzp[14] >=0.f && CovXyzp[14] < 1.f;
  goodTrack &= goodTrack && CovXyzp[20] >=0.f && CovXyzp[20] < 1.f;
  if(!goodTrack) return false;
  
  track.SetParameters(xyzp);
  track.SetCovarianceMatrix(CovXyzp);
  track.SetNDF(1);
  //    track.SetChi2(GlobalTracks_mChiSqXY[k]);
  track.SetID(index);

  track.SetCharge(q);
  return true;
}

Int_t StKFParticleAnalysisMaker::Finish() 
{ 
  fDirs[0]->cd();
	WriteHistograms();
	
  return kStOK;
}
