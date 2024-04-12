#include "StCorrection.h"
#include "StConstants.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"

ClassImp(StCorrection)

TString StCorrection::mVStr[2] = {"pos","neg"};
TString StCorrection::mOrder[2] = {"2nd","3rd"};
//---------------------------------------------------------------------------------

StCorrection::StCorrection()
{
}

StCorrection::~StCorrection()
{
}

//---------------------------------------------------------------------------------

void StCorrection::InitReCenterCorrection()
{
  TString InPutFile = Form("/star/u/lliu/test/11p5GeV/recenterpar/ReCenterPar.root");

  mInPutFile = TFile::Open(InPutFile.Data());

  for(Int_t i = 0; i < 4; i++)
  {
    mQ2Vector_East_EP[i].Set(0.0,0.0);
    mQ3Vector_East_EP[i].Set(0.0,0.0);
    mQCounter_East[i] = 0;

    mQ2Vector_West_EP[i].Set(0.0,0.0);
    mQ3Vector_West_EP[i].Set(0.0,0.0);
    mQCounter_West[i] = 0;
  }

  mQ2Vector_Full_EP.Set(0.0,0.0);
  mQ3Vector_Full_EP.Set(0.0,0.0);
  mQCounter_Full = 0;

}

//---------------------------------------------------------------------------------


bool StCorrection::passTrackEtaEast(StPicoTrack *track, Int_t i, Int_t Mode) // neg || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
  Float_t eta = track->pMom().PseudoRapidity();
  
  if(Mode == 0) // Event Plane Mode
  {
    // eta_gap between two sub event plane is 2*mEta_Gap[i]
    if(!(eta > -1.0*TriFlow::mEtaMax && eta < -1.0*TriFlow::mEta_Gap[i])) return kFALSE;

    return kTRUE;
  }

  if(Mode == 1) // Flow Mode
  {
    if(!(eta > -1.0*TriFlow::mEtaMax && eta < 0.0)) return kFALSE;

    return kTRUE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StCorrection::passTrackEtaWest(StPicoTrack *track, Int_t i, Int_t Mode) // pos || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
  Float_t eta = track->pMom().PseudoRapidity();

  if(Mode == 0) // Event Plane Mode
  {
    // eta_gap between two sub event plane is 2*mEta_Gap[i]
    if(!(eta > TriFlow::mEta_Gap[i] && eta < TriFlow::mEtaMax)) return kFALSE;

    return kTRUE;
  }

  if(Mode == 1) // Flow Mode
  {
    if(!(eta > 0.0 && eta < TriFlow::mEtaMax)) return kFALSE;

    return kTRUE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------
// Event Plane method and Scalar Product method
TVector2 StCorrection::calq2Vector(StPicoTrack *track)
{
  const Float_t phi = track->pMom().Phi();
  TVector2 q2Vector(0.0,0.0);

  const Float_t q2x = TMath::Cos(2.0*phi);
  const Float_t q2y = TMath::Sin(2.0*phi);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StCorrection::calq3Vector(StPicoTrack *track)
{
  const Float_t phi = track->pMom().Phi();
  TVector2 q3Vector(0.0,0.0);

  const Float_t q3x = TMath::Cos(3.0*phi);
  const Float_t q3y = TMath::Sin(3.0*phi);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
}

Float_t StCorrection::getWeight(StPicoTrack *track)
{
  Float_t pt = track->pMom().Perp();
  Float_t w;
  if(pt <= TriFlow::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > TriFlow::mPrimPtWeight)
  {
    w = TriFlow::mPrimPtWeight;
  }

  return w;
}
//------------------------------------------------------------------------------

TVector2 StCorrection::getReCenterPar_East(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap)
{
  Float_t mean_Qx, mean_Qy;
  TVector2 Qvector(0.0,0.0);

  TString ProName_x = Form("Qx_%s_Vertex_%s_EtaGap_%d_East_EP", mOrder[order].Data(), mVStr[vz_sign].Data(), eta_gap);
  TString ProName_y = Form("Qy_%s_Vertex_%s_EtaGap_%d_East_EP", mOrder[order].Data(), mVStr[vz_sign].Data(), eta_gap);

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_Qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex, (Double_t)Cent9));
  mean_Qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex, (Double_t)Cent9));

  Qvector.Set(mean_Qx,mean_Qy);

  return Qvector;
}

//---------------------------------------------------------------------------------

TVector2 StCorrection::getReCenterPar_West(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap)
{
  Float_t mean_Qx, mean_Qy;
  TVector2 Qvector(0.0,0.0);

  TString ProName_x = Form("Qx_%s_Vertex_%s_EtaGap_%d_West_EP", mOrder[order].Data(), mVStr[vz_sign].Data(), eta_gap);
  TString ProName_y = Form("Qy_%s_Vertex_%s_EtaGap_%d_West_EP", mOrder[order].Data(), mVStr[vz_sign].Data(), eta_gap);

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_Qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex, (Double_t)Cent9));
  mean_Qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex, (Double_t)Cent9));

  Qvector.Set(mean_Qx,mean_Qy);

  return Qvector;
}

//---------------------------------------------------------------------------------

TVector2 StCorrection::getReCenterPar_Full(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign )
{
  Float_t mean_Qx, mean_Qy;
  TVector2 Qvector(0.0,0.0);

  TString ProName_x = Form("Qx_%s_Vertex_%s_Full_EP",mOrder[order].Data(),mVStr[vz_sign].Data() );
  TString ProName_y = Form("Qy_%s_Vertex_%s_Full_EP",mOrder[order].Data(),mVStr[vz_sign].Data() );

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_Qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_Qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  Qvector.Set(mean_Qx,mean_Qy);

  return Qvector;
}

//---------------------------------------------------------------------------------

void StCorrection::addTrack_East(StPicoTrack *track, Int_t j)
{
  Float_t w = getWeight(track);
  mQ2Vector_East_EP[j] += w * calq2Vector(track);
  mQ3Vector_East_EP[j] += w * calq3Vector(track);

  mQCounter_East[j]++;
}

//---------------------------------------------------------------------------------

void StCorrection::addTrack_West(StPicoTrack *track, Int_t j)
{
  Float_t w = getWeight(track);
  mQ2Vector_West_EP[j] += w * calq2Vector(track);
  mQ3Vector_West_EP[j] += w * calq3Vector(track);

  mQCounter_West[j]++;
}

//---------------------------------------------------------------------------------

void StCorrection::addTrack_Full(StPicoTrack *track)
{
  Float_t w = getWeight(track);
  mQ2Vector_Full_EP += w * (calq2Vector(track));
  mQ3Vector_Full_EP += w * (calq3Vector(track));

  mQCounter_Full++;
}

void StCorrection::RecenterQVec(Int_t Cent9, Int_t RunIndex, Int_t vz_sign)
{
  // Only for Save Raw Event Plane
  mQ2Vector_Full_EP_raw = mQ2Vector_Full_EP;
  mQ3Vector_Full_EP_raw = mQ3Vector_Full_EP;
  // For Raw East and West, only save the case which 0.05 < |eta| < 1.5
  mQ2Vector_East_EP_raw = mQ2Vector_East_EP[0];
  mQ3Vector_East_EP_raw = mQ3Vector_East_EP[0];
  mQ2Vector_West_EP_raw = mQ2Vector_West_EP[0];
  mQ3Vector_West_EP_raw = mQ3Vector_West_EP[0];


  mQ2Vector_Full_EP -= getReCenterPar_Full(0, Cent9, RunIndex, vz_sign); 
  mQ3Vector_Full_EP -= getReCenterPar_Full(1, Cent9, RunIndex, vz_sign);

  for (int j = 0; j < 4; j++)
  {
    mQ2Vector_East_EP[j] -= getReCenterPar_East(0, Cent9, RunIndex, vz_sign, j);
    mQ3Vector_East_EP[j] -= getReCenterPar_East(1, Cent9, RunIndex, vz_sign, j);

    mQ2Vector_West_EP[j] -= getReCenterPar_West(0, Cent9, RunIndex, vz_sign, j);
    mQ3Vector_West_EP[j] -= getReCenterPar_West(1, Cent9, RunIndex, vz_sign, j);
  }
}
//---------------------------------------------------------------------------------

void StCorrection::clear()
{
  for(Int_t i = 0; i < 4; i++)
  {
    mQ2Vector_East_EP[i].Set(0.0,0.0);
    mQ3Vector_East_EP[i].Set(0.0,0.0);
    mQCounter_East[i] = 0;

    mQ2Vector_West_EP[i].Set(0.0,0.0);
    mQ3Vector_West_EP[i].Set(0.0,0.0);
    mQCounter_West[i] = 0;
  }
  
  mQ2Vector_Full_EP.Set(0.0,0.0);
  mQ3Vector_Full_EP.Set(0.0,0.0);
  mQCounter_Full = 0;

  mQ2Vector_East_EP_raw.Set(0.0,0.0);
  mQ2Vector_West_EP_raw.Set(0.0,0.0);
  mQ2Vector_Full_EP_raw.Set(0.0,0.0);
  mQ3Vector_East_EP_raw.Set(0.0,0.0);
  mQ3Vector_West_EP_raw.Set(0.0,0.0);
  mQ3Vector_Full_EP_raw.Set(0.0,0.0);

}

//---------------------------------------------------------------------------------

void StCorrection::InitNtuple()
{
  mNtuple = new TNtuple("Ntuple","Ntuple","runId:eventId:RefMult:ZDCx:BBCx:vzVpd:Centrality9:vx:vy:vz:nToFMatched:mQ2X_East_EP_0:mQ2X_East_EP_1:mQ2X_East_EP_2:mQ2X_East_EP_3:mQ2Y_East_EP_0:mQ2Y_East_EP_1:mQ2Y_East_EP_2:mQ2Y_East_EP_3:mQ3X_East_EP_0:mQ3X_East_EP_1:mQ3X_East_EP_2:mQ3X_East_EP_3:mQ3Y_East_EP_0:mQ3Y_East_EP_1:mQ3Y_East_EP_2:mQ3Y_East_EP_3:mQ2X_West_EP_0:mQ2X_West_EP_1:mQ2X_West_EP_2:mQ2X_West_EP_3:mQ2Y_West_EP_0:mQ2Y_West_EP_1:mQ2Y_West_EP_2:mQ2Y_West_EP_3:mQ3X_West_EP_0:mQ3X_West_EP_1:mQ3X_West_EP_2:mQ3X_West_EP_3:mQ3Y_West_EP_0:mQ3Y_West_EP_1:mQ3Y_West_EP_2:mQ3Y_West_EP_3:mQ2X_Full_EP:mQ2Y_Full_EP:mQ3X_Full_EP:mQ3Y_Full_EP:mQCounter_East_0:mQCounter_East_1:mQCounter_East_2:mQCounter_East_3:mQCounter_West_0:mQCounter_West_1:mQCounter_West_2:mQCounter_West_3:mQCounter_Full:mQ2X_East_EP_raw:mQ2Y_East_EP_raw:mQ2X_West_EP_raw:mQ2Y_West_EP_raw:mQ2X_Full_EP_raw:mQ2Y_Full_EP_raw:mQ3X_East_EP_raw:mQ3Y_East_EP_raw:mQ3X_West_EP_raw:mQ3Y_West_EP_raw:mQ3X_Full_EP_raw:mQ3Y_Full_EP_raw:runIndex");
  mNtuple->SetAutoSave(50000000);
}

//---------------------------------------------------------------------------------

Int_t StCorrection::fillNtuple(StPicoDst *pico, Int_t Cent9, Int_t nToFMatched, Int_t runIndex)
{
  StPicoEvent *event = pico->event();
  if(!event)
  {
    return kFALSE;
  }

  // event information
  mFillNtuple[0]  = (Float_t)event->runId();
  mFillNtuple[1]  = (Float_t)event->eventId();
  mFillNtuple[2]  = (Float_t)event->refMult();
  mFillNtuple[3]  = (Float_t)event->ZDCx();
  mFillNtuple[4]  = (Float_t)event->BBCx();
  mFillNtuple[5]  = (Float_t)event->vzVpd();
  mFillNtuple[6]  = (Float_t)Cent9;
  mFillNtuple[7]  = (Float_t)event->primaryVertex().X();
  mFillNtuple[8]  = (Float_t)event->primaryVertex().Y();
  mFillNtuple[9]  = (Float_t)event->primaryVertex().Z();
  mFillNtuple[10] = (Float_t)nToFMatched;

  // Q Vector Event Plane method
  // East
  mFillNtuple[11] = (Float_t)mQ2Vector_East_EP[0].X();
  mFillNtuple[12] = (Float_t)mQ2Vector_East_EP[1].X();
  mFillNtuple[13] = (Float_t)mQ2Vector_East_EP[2].X();
  mFillNtuple[14] = (Float_t)mQ2Vector_East_EP[3].X();

  mFillNtuple[15] = (Float_t)mQ2Vector_East_EP[0].Y();
  mFillNtuple[16] = (Float_t)mQ2Vector_East_EP[1].Y();
  mFillNtuple[17] = (Float_t)mQ2Vector_East_EP[2].Y();
  mFillNtuple[18] = (Float_t)mQ2Vector_East_EP[3].Y();

  mFillNtuple[19] = (Float_t)mQ3Vector_East_EP[0].X();
  mFillNtuple[20] = (Float_t)mQ3Vector_East_EP[1].X();
  mFillNtuple[21] = (Float_t)mQ3Vector_East_EP[2].X();
  mFillNtuple[22] = (Float_t)mQ3Vector_East_EP[3].X();

  mFillNtuple[23] = (Float_t)mQ3Vector_East_EP[0].Y();
  mFillNtuple[24] = (Float_t)mQ3Vector_East_EP[1].Y();
  mFillNtuple[25] = (Float_t)mQ3Vector_East_EP[2].Y();
  mFillNtuple[26] = (Float_t)mQ3Vector_East_EP[3].Y();
  // West
  mFillNtuple[27] = (Float_t)mQ2Vector_West_EP[0].X();
  mFillNtuple[28] = (Float_t)mQ2Vector_West_EP[1].X();
  mFillNtuple[29] = (Float_t)mQ2Vector_West_EP[2].X();
  mFillNtuple[30] = (Float_t)mQ2Vector_West_EP[3].X();

  mFillNtuple[31] = (Float_t)mQ2Vector_West_EP[0].Y();
  mFillNtuple[32] = (Float_t)mQ2Vector_West_EP[1].Y();
  mFillNtuple[33] = (Float_t)mQ2Vector_West_EP[2].Y();
  mFillNtuple[34] = (Float_t)mQ2Vector_West_EP[3].Y();

  mFillNtuple[35] = (Float_t)mQ3Vector_West_EP[0].X();
  mFillNtuple[36] = (Float_t)mQ3Vector_West_EP[1].X();
  mFillNtuple[37] = (Float_t)mQ3Vector_West_EP[2].X();
  mFillNtuple[38] = (Float_t)mQ3Vector_West_EP[3].X();

  mFillNtuple[39] = (Float_t)mQ3Vector_West_EP[0].Y();
  mFillNtuple[40] = (Float_t)mQ3Vector_West_EP[1].Y();
  mFillNtuple[41] = (Float_t)mQ3Vector_West_EP[2].Y();
  mFillNtuple[42] = (Float_t)mQ3Vector_West_EP[3].Y();
  // Full
  mFillNtuple[43] = (Float_t)mQ2Vector_Full_EP.X();
  mFillNtuple[44] = (Float_t)mQ2Vector_Full_EP.Y();
  mFillNtuple[45] = (Float_t)mQ3Vector_Full_EP.X();
  mFillNtuple[46] = (Float_t)mQ3Vector_Full_EP.Y();

  // East
  mFillNtuple[47] = (Float_t)mQCounter_East[0];
  mFillNtuple[48] = (Float_t)mQCounter_East[1];
  mFillNtuple[49] = (Float_t)mQCounter_East[2];
  mFillNtuple[50] = (Float_t)mQCounter_East[3];
  // West
  mFillNtuple[51] = (Float_t)mQCounter_West[0];
  mFillNtuple[52] = (Float_t)mQCounter_West[1];
  mFillNtuple[53] = (Float_t)mQCounter_West[2];
  mFillNtuple[54] = (Float_t)mQCounter_West[3];
  // Full
  mFillNtuple[55] = (Float_t)mQCounter_Full;
  //raw 2nd EP
  mFillNtuple[56] = (Float_t)mQ2Vector_East_EP_raw.X();
  mFillNtuple[57] = (Float_t)mQ2Vector_East_EP_raw.Y();
  mFillNtuple[58] = (Float_t)mQ2Vector_West_EP_raw.X();
  mFillNtuple[59] = (Float_t)mQ2Vector_West_EP_raw.Y();
  mFillNtuple[60] = (Float_t)mQ2Vector_Full_EP_raw.X();
  mFillNtuple[61] = (Float_t)mQ2Vector_Full_EP_raw.Y();

  //raw 3rd EP
  mFillNtuple[62] = (Float_t)mQ3Vector_East_EP_raw.X();
  mFillNtuple[63] = (Float_t)mQ3Vector_East_EP_raw.Y();
  mFillNtuple[64] = (Float_t)mQ3Vector_West_EP_raw.X();
  mFillNtuple[65] = (Float_t)mQ3Vector_West_EP_raw.Y();
  mFillNtuple[66] = (Float_t)mQ3Vector_Full_EP_raw.X();
  mFillNtuple[67] = (Float_t)mQ3Vector_Full_EP_raw.Y();

  // runIndex
  mFillNtuple[68] = (Float_t)runIndex;

  // random sub
  mNtuple->Fill(mFillNtuple);

  return kTRUE;
}

//---------------------------------------------------------------------------------

void StCorrection::writeNtuple()
{
  mNtuple->Write();
}

//---------------------------------------------------------------------------------

bool StCorrection::passTrackEtaNumCut(Int_t j)
{
  if(!(mQCounter_East[j] > TriFlow::mTrackMin && mQCounter_West[j] > TriFlow::mTrackMin)) return kFALSE;

  return kTRUE;
}

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
TVector2 StCorrection::calPsi2_East_EP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_East_EP[j].X();
  Float_t Qy = mQ2Vector_East_EP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder2[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder2[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StCorrection::calPsi2_West_EP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_West_EP[j].X();
  Float_t Qy = mQ2Vector_West_EP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder2[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder2[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StCorrection::calPsi2_Full_EP(Int_t k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_Full_EP.X();
  Float_t Qy = mQ2Vector_Full_EP.Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder2[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder2[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

//---------------------------------------------------------------------------------
// 3rd
TVector2 StCorrection::calPsi3_East_EP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ3Vector_East_EP[j].X();
  Float_t Qy = mQ3Vector_East_EP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/3.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder3[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder3[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StCorrection::calPsi3_West_EP(Int_t j, Int_t k) // j = eta_gap, k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ3Vector_West_EP[j].X();
  Float_t Qy = mQ3Vector_West_EP[j].Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/3.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder3[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder3[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StCorrection::calPsi3_Full_EP(Int_t k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ3Vector_Full_EP.X();
  Float_t Qy = mQ3Vector_Full_EP.Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/3.0;
  Float_t Psi_Sin = TMath::Sin(TriFlow::mShiftOrder3[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(TriFlow::mShiftOrder3[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

//---------------------------------------------------------------------------------

Float_t StCorrection::AngleShift(Float_t Psi_raw, Float_t order)
{
  Float_t Psi_Corr = Psi_raw;
  if(Psi_raw > TMath::Pi()/order)
  {
    Psi_Corr = Psi_raw - 2.0*TMath::Pi()/order;
  }
  if(Psi_raw < -1.0*TMath::Pi()/order)
  {
    Psi_Corr = Psi_raw + 2.0*TMath::Pi()/order;
  }

  return Psi_Corr;
}

//---------------------------------------------------------------------------------

// Event Plane method

bool StCorrection::IsBadQVec()
{
  bool mBadQvec = false;
  for(Int_t i = 0; i < 4; i++)
  {
    if ( ( mQ2Vector_East_EP[i].X() == 0.0 && mQ2Vector_East_EP[i].Y() == 0.0 ) || 
         ( mQ2Vector_West_EP[i].X() == 0.0 && mQ2Vector_West_EP[i].Y() == 0.0 ) )
    { 
      //cout << "east Q is X " << mQ2Vector_East_EP[i].X() << " west Q is " << mQ2Vector_West_EP[i].X() << endl;
      //cout << "east Q is Y " << mQ2Vector_East_EP[i].Y() << " west Q is " << mQ2Vector_West_EP[i].Y() << endl;
      mBadQvec=true;
    }
    //if ( ( mQ3Vector_East_EP[i].X() == 0.0 && mQ3Vector_East_EP[i].Y() == 0.0 ) || 
    //     ( mQ3Vector_West_EP[i].X() == 0.0 && mQ3Vector_West_EP[i].Y() == 0.0 ) ) mBadQvec=true;
  }
  return mBadQvec;
}