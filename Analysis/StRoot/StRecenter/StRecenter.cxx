// 2022, 02, 25: Like add function "BadQVec()" to reject Qx=Qy=0, replace function "passTrackEtaNumCut()" 
#include "StRecenter.h"
#include "StConstants.h"
#include "StCut.h"
#include "StProManger.h"
#include "StCorrection.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>
#include <array>
#include "badrun.h"
#include "run.h"
#include "para_pileup.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
ClassImp(StRecenter)

StRecenter::StRecenter(const char* name, StPicoDstMaker *picoMaker, char* jobid)
: StMaker(name)
{
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;

    mOutPut_Shift =Form("%s_shift.root",jobid);
    mOutPut_ReCenter=Form("%s_recenter.root",jobid);
}

//----------------------------------------------------------------------------- 
StRecenter::~StRecenter()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StRecenter::Init() 
{
    mRefMultCorrUtil = new StRefMultCorr("refmult");
    mCut = new StCut();
    mProManger = new StProManger();
    mCorrection = new StCorrection();

    mCorrection->InitReCenterCorrection();
    mFile_Corr_Shift = new TFile(mOutPut_Shift.Data(),"RECREATE");
    mProManger->InitShift();
    mFile_Corr_ReCenter = new TFile(mOutPut_ReCenter.Data(),"RECREATE");
    mFile_Corr_ReCenter->cd();
    mCorrection->InitNtuple();
    return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StRecenter::Finish() 
{

    mFile_Corr_ReCenter->cd();
    mCorrection->writeNtuple();
    mFile_Corr_ReCenter->Close();

    mFile_Corr_Shift->cd();
    mProManger->WriteShift();
    mFile_Corr_Shift->Close();
    return kStOK;
}

//----------------------------------------------------------------------------- 
void StRecenter::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StRecenter::Make() 
{
    if(!mPicoDstMaker) 
    {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) 
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }

    event = (StPicoEvent*)mPicoDst->event();
    if(!event)
    {
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }

    int runId = event->runId();
    int runIndex = GetRunIndex(runId);

    // bad runs have been excluded from input list
    for (int i = 0; i < n_badrun; i++)
    {
        if (runId == badrun[i]) return 0;
    }

    if(!mCut->passEventCut(event)) return 0;

    TVector3 pVtx = event->primaryVertex();
    Int_t refMult = event->refMult();
    Int_t nBTofMatched = event->nBTOFMatch();

    // The 3rd parameter for Vz bin, Not used yet
    if(isPileUp(refMult, nBTofMatched, 1)) return 0;

    //-------------------------Self Centrality ----------------------------------
    int cent9 = Centrality(refMult);
    if(cent9 > 8 || cent9 < 0) return kStOK;

    //-------------------------Official Centrality ----------------------------------
    /*mRefMultCorrUtil->init(runId);
    if ( mRefMultCorrUtil->isBadRun( runId ) ) return 0;

    if(!mCut->passEventCut(event)) return 0;

    // IMPORTANT: vertex position is needed for Au+Au 19.6 GeV 2019
    if (mRefMultCorrUtil->isPileUpEvent( refMult, nBTofMatched, pVtx.Z() ) ) return 0;

    mRefMultCorrUtil->initEvent(refMult, pVtx.Z(), event->ZDCx());

    // In case centrality calculation is failed for the given event it will
    if (mRefMultCorrUtil->getCentralityBin16() < 0 ||
        mRefMultCorrUtil->getCentralityBin9() < 0) return 0;

    Int_t cent9 = mRefMultCorrUtil->getCentralityBin9();       // 0: 70-80%, 1: 60-70%, ..., 6: 10-20%, 7: 5-10%, 8: 0-5%
    Double_t weight = mRefMultCorrUtil->getWeight();*/

    // vz sign
    Float_t vz = pVtx.Z();
    Int_t vz_sign;
    if(vz > 0.0) vz_sign = 0;
    else vz_sign = 1;

    const Int_t nTracks = mPicoDst->numberOfTracks();
    for(Int_t i = 0; i < nTracks; i++) // track loop
    {
        StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
        //StPicoPhysicalHelix helix = track->helix(mField);
        //Float_t dca = helix.geometricSignedDistance(pVtx);
        Float_t dca = track->gDCA(pVtx).Mag();
        if(mCut->passTrackEP(track,dca))
        {
            mCorrection->addTrack_Full(track);
            for(Int_t j = 0; j < 4; j++) // eta_gap loop
            {
                if(mCorrection->passTrackEtaEast(track, j, 0)) // neg eta sub
                {
                    mCorrection->addTrack_East(track, j);
                }
                if(mCorrection->passTrackEtaWest(track, j, 0)) // pos eta sub
                {
                    mCorrection->addTrack_West(track, j);
                }
            }
        }
    }

    if ( mCorrection->IsBadQVec() ) return 0; 
    // Event-wise calibration re-centering
    mCorrection->RecenterQVec(cent9, runIndex, vz_sign);

    mCorrection->fillNtuple(mPicoDst,cent9,nBTofMatched,runIndex);

    // full event shift parameter 
    for(Int_t k = 0; k < 5; k++) // ShiftOrder loop
    {
        // Event Plane method
        TVector2 Psi2Vector_Full_EP = mCorrection->calPsi2_Full_EP(k);
        TVector2 Psi3Vector_Full_EP = mCorrection->calPsi3_Full_EP(k);
        mProManger->FillEventFull_EP(Psi2Vector_Full_EP,Psi3Vector_Full_EP,cent9,runIndex,vz_sign,k);// Shift par
    }

    // eta sub event shift parameter
    for(Int_t j = 0; j < 4; j ++)  // eta gap
    {
        for(Int_t k = 0; k < 5; k++) // ShiftOrder loop
        {
            // Event Plane method
            TVector2 Psi2Vector_East_EP = mCorrection->calPsi2_East_EP(j,k);
            TVector2 Psi3Vector_East_EP = mCorrection->calPsi3_East_EP(j,k);
            mProManger->FillEventEast_EP(Psi2Vector_East_EP,Psi3Vector_East_EP,cent9,runIndex,vz_sign,j,k);

            TVector2 Psi2Vector_West_EP = mCorrection->calPsi2_West_EP(j,k);
            TVector2 Psi3Vector_West_EP = mCorrection->calPsi3_West_EP(j,k);
            mProManger->FillEventWest_EP(Psi2Vector_West_EP,Psi3Vector_West_EP,cent9,runIndex,vz_sign,j,k);

        }
    }
    mCorrection->clear();

    return kStOK;
}
//__________________________________________________________________________________

int StRecenter::Centrality(int gRefMult )
{
    int centrality;
    int centFull[9]={6, 11, 21, 38, 61, 95, 141, 205, 249};
    if      (gRefMult>=centFull[8]) centrality=8;
    else if (gRefMult>=centFull[7]) centrality=7;
    else if (gRefMult>=centFull[6]) centrality=6;
    else if (gRefMult>=centFull[5]) centrality=5;
    else if (gRefMult>=centFull[4]) centrality=4;
    else if (gRefMult>=centFull[3]) centrality=3;
    else if (gRefMult>=centFull[2]) centrality=2;
    else if (gRefMult>=centFull[1]) centrality=1;
    else if (gRefMult>=centFull[0]) centrality=0;
    else centrality = 9;

    return centrality;
}
