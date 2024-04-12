#ifndef StProManger_h
#define StProManger_h

#include "TVector2.h"
#include "TString.h"
#include "StMessMgr.h"

class TProfile2D;
class TProfile;
class TH1F;

class StProManger
{
  public:
    StProManger();
    ~StProManger();


    // Shift Correction
    void InitShift();
    void WriteShift();
    // Event Plane method
    void FillEventEast_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k); // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
    void FillEventWest_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k); // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
    void FillEventFull_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k); // i = vertex pos/neg, k = ShiftOrder

  private:

    //raw EventPlane distribution (before recenter and shift)

    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    // Event Plane method
    TProfile2D *p_mq2x_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2y_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2x_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2y_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2x_Full_EP[2];    // 0 = vertex pos/neg
    TProfile2D *p_mq2y_Full_EP[2];    // 0 = vertex pos/neg

    TProfile2D *p_mq3x_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3y_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3x_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3y_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3x_Full_EP[2];    // 0 = vertex pos/neg
    TProfile2D *p_mq3y_Full_EP[2];    // 0 = vertex pos/neg

    // Shift Correction | x axis is RunIndex, y axis is Centrality
    // Event Plane method
    TProfile2D *p_mcos2_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin2_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos2_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin2_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos2_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_msin2_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder

    TProfile2D *p_mcos3_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin3_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos3_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin3_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos3_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_msin3_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder

    static TString mVStr[2];

    ClassDef(StProManger,1)
};

#endif
