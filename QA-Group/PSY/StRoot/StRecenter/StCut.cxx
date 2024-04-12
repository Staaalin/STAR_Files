#include "StCut.h"
#include "StConstants.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StMessMgr.h"

ClassImp(StCut)

StCut::StCut()
{
}

StCut::~StCut()
{
}

//---------------------------------------------------------------------------------

bool StCut::passEventCut(StPicoEvent *event)
{
    const Float_t vx = event->primaryVertex().X();
    const Float_t vy = event->primaryVertex().Y();
    const Float_t vz = event->primaryVertex().Z();

    if(fabs(vz) > TriFlow::mVzMax) return kFALSE;
    if(sqrt(vx*vx+vy*vy) > TriFlow::mVrMax) return kFALSE;
    /*if(!(event->isTrigger(640001) || event->isTrigger(640011) || event->isTrigger(640021)
      || event->isTrigger(640031) || event->isTrigger(640041) || event->isTrigger(640051))) return kFALSE;*/

    return kTRUE;
}

bool StCut::passTrackBasic(StPicoTrack *track)
{
    Float_t nHitsFit = track->nHitsFit();
    Float_t nHitsMax = track->nHitsMax();

    if ( nHitsFit < TriFlow::mHitsFitTPCMin ) return kFALSE;
    //if ( nHitsFit < TriFlow::mHitsRatioTPCMin * nHitsMax) return kFALSE; // No Ratio cut

    // eta cut
    Float_t eta = track->pMom().PseudoRapidity();
    if(fabs(eta) > TriFlow::mEtaMax) return kFALSE;

    return kTRUE;
}

//---------------------------------------------------------------------------------

bool StCut::passTrackEP(StPicoTrack *track, float dca)
{
    if(!track) return kFALSE;
    if(!passTrackBasic(track)) return kFALSE;

    if(dca > TriFlow::mDcaEPMax) return kFALSE;
        
    // pt cut 0.2 - 2.0 GeV/c
    Float_t pt = track->pMom().Perp();
    Float_t p  = track->pMom().Mag();
    if(!(pt > TriFlow::mPrimPtMin && pt < TriFlow::mPrimPtMax && p < TriFlow::mPrimMomMax))
    {
        return kFALSE;
    }

    return kTRUE;
}
//---------------------------------------------------------------------------------
