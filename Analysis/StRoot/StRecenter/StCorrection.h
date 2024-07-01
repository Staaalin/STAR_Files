#ifndef StCorrection_h
#define StCorrection_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class TProfile2D;
class TFile;
class TNtuple;

class StCorrection : public TObject
{
  public:
    StCorrection();
    ~StCorrection();

    // ReCenter Correction
    bool passTrackEtaEast(StPicoTrack*,Int_t,Int_t);
    bool passTrackEtaWest(StPicoTrack*,Int_t,Int_t);

    Float_t getWeight(StPicoTrack*);
    TVector2 calq2Vector(StPicoTrack*);
    TVector2 calq3Vector(StPicoTrack*);
    TVector2 getReCenterPar_Full(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign); // order: 0 = 2nd, 1 = 3rd 
    TVector2 getReCenterPar_East(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap); // order: 0 = 2nd, 1 = 3rd
    TVector2 getReCenterPar_West(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap); // order: 0 = 2nd, 1 = 3rd 

    void InitReCenterCorrection();
    void addTrack_Full(StPicoTrack* track);
    void addTrack_East(StPicoTrack* track, Int_t eta_gap);
    void addTrack_West(StPicoTrack* track, Int_t eta_gap);
    void RecenterQVec(Int_t Cent9, Int_t RunIndex, Int_t vz_sign);

    void InitNtuple();
    Int_t fillNtuple(StPicoDst*, Int_t, Int_t, Int_t);
    void writeNtuple();

    bool passTrackEtaNumCut(Int_t); // eta_gap // Not Used
    bool IsBadQVec();

    TVector2 calPsi2_East_EP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 2, 4, 6, 8, 10
    TVector2 calPsi2_West_EP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 2, 4, 6, 8, 10
    TVector2 calPsi2_Full_EP(Int_t); // 1 = ShiftOrder: 2, 4, 6, 8, 10

    TVector2 calPsi3_East_EP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 3, 6, 9, 12, 15
    TVector2 calPsi3_West_EP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 3, 6, 9, 12, 15
    TVector2 calPsi3_Full_EP(Int_t); // 1 = ShiftOrder: 3, 6, 9, 12, 15

    Float_t AngleShift(Float_t Psi_raw, Float_t order);

    void clear();

  private:
    //Event Plane method
    TVector2 mQ2Vector_Full_EP, mQ2Vector_Full_EP_raw, mQ2Vector_East_EP_raw, mQ2Vector_West_EP_raw;
    TVector2 mQ3Vector_Full_EP, mQ3Vector_Full_EP_raw, mQ3Vector_East_EP_raw, mQ3Vector_West_EP_raw;
    TVector2 mQ2Vector_East_EP[4], mQ2Vector_West_EP[4], mQ3Vector_East_EP[4], mQ3Vector_West_EP[4];

    Int_t mQCounter_East[4], mQCounter_West[4], mQCounter_Full;

    TFile *mInPutFile;
    TNtuple *mNtuple;
    Float_t mFillNtuple[69];

    TFile *mInPutFile_Shift;
    TFile *mInPutFile_Res;

    static TString mVStr[2];
    static TString mOrder[2];

  ClassDef(StCorrection,1)
};

#endif
