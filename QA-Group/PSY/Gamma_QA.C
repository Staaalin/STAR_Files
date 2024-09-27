using namespace std;

#include "stdio.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

/// PicoDst headers
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoDstReader.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoDst.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoEvent.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoTrack.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoBTofHit.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoBTowHit.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoEmcTrigger.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoBTofPidTraits.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StEpdUtil/StEpdEpFinder.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StRefMultCorr/StRefMultCorr.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StRefMultCorr/CentralityMaker.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoEpdHit.h"
// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StEpdUtil//StEpdGeom.h"

/// PicoDst headers
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StEpdUtil//StEpdGeom.h"

//class StRefMultCorr;
//class CentralityMaker;

const float PI = TMath::Pi();
const float MM = 2/PI;
const int Phibin = 80;
const int order = 5;
const float mpi= 0.13957;

const int opt_sys = 0;
const int opt_HighOrder = 0;//0 for gamma112, 1 for 132, not 2
const int opt_useEPD = 11;  //0 is TPC EP, 1 is EPD, 11 is 1st-order EPD, 2 is BBC, 3 is ZDC
const int opt_boost = 0;
const int opt_ESE = 0; //0 means default, 1 will select 0.4 < q < 1.6
const int nHar = 2; //2nd order EP

const float Eta_EP_Cut = 0.55;	//particles in the event plane
const float Eta_ESE_Cut = 0.5;  //particles for event shape engineering in using TPC EP

const float pt_asso_up = 2.0;//1;
const float pt_asso_lo = 0.06;//0.15;
const float EtaCut_asso = 1.4999;
const float DcaCut_asso = 3;         //2
const int nFitHits_asso = 10;	//15, number of fit hits

const float pt_trig_up = 2;
const float pt_trig_lo = 0.06;
const float EtaCut_trig = 1.4999;//1.4999
const float DcaCut_trig = 3;
const int nFitHits_trig = 10;
const int nFitRatio_trig_lo = 0.52;  //0.52
const int nFitRatio_trig_up = 1.05;  //1.05

const Float_t Vz_cut = 2;    	//30, Vz
const Float_t Vr_cut = 2;	//2, Vr
const float Vz_diff = 10;        //5, difference between Vz{TPC} and Vz{VPD}
const int cenDef[9] = {6,11,21,38,61,95,141,205,249}; //Updated 19June 2023 for 11GeV v1

//Full production dataset 
////
static Int_t runmin=19158057; // 19158057
static Int_t runmax=19169016;//need to add 1 from the runmax 19161029
static Int_t runbins=runmax-runmin;
const int run_sta = (int)((runmin%1000000)/1);
const int run_end = (int)((runmax%1000000)/1);


//
const float mEPDthresh = 0.3;
const float mEPDMax = 5.0;

//bad runs
const int Nrun_MB = 0;
const int bad_Ref_day3[Nrun_MB] ={};
// {5040,6007,6008,6013,6014,6015,6027,6029,6031,7034 ,10036,11001,11004,12030,12035,14027,15029,17047,17048,19019,19020,21008,21009,21010,21011};//Updated 19June 2023 for 11GeV v1 from Like-Liu and Eugenia’s study: //https://drupal.star.bnl.gov/STAR/system/files/auau11_run20.pdf


//new efficiency
const float PP0[9] ={8.72972e-01,8.76830e-01,8.68846e-01,8.80217e-01,8.73341e-01,8.72515e-01,8.70356e-01,8.68530e-01,8.72790e-01};
const float PP1[9] = {1.24992e-01,1.41588e-01,1.40217e-01,1.42039e-01,1.43380e-01,1.40686e-01,1.42506e-01,1.44801e-01,1.50136e-01};
const float PP2[9] = {4.82547e+00,7.87162e+00,1.72209e+01,1.18083e+01,1.37230e+01,2.53850e+01,9.90126e+00,1.02641e+01,5.92769e+00};

//defining histograms
	struct WeightHists {
		TProfile2D *pTPCmeanPhi_1 ;
                TProfile2D *pTPCmeanPhi_2 ;
                TProfile2D *pTPCmeanPhi_3 ;
		TProfile2D *pTPCmeanPhiAsso_1 ;
                TProfile2D *pTPCmeanPhiAsso_2 ;
                TProfile2D *pTPCmeanPhiAsso_3 ;
	};
        struct WeightHists WeightHist[15];
        struct WeightHists ReadWeightHist[15];

        struct PhiHists {
		TH2D* Hist_Phi_1;
                TH2D* Hist_Phi_2;
                TH2D* Hist_Phi_3;
                TH2D* Hist_Phi_new_1;
                TH2D* Hist_Phi_new_2;
                TH2D* Hist_Phi_new_3;
	};
        struct PhiHists PhiHist[15];

        TProfile *pTemp_track = new TProfile("pTemp_track","pTemp_track",10,0.5,10.5,-50,50); //pT, eta,phi,charge,DCA,
        TProfile *pTemp_track_dedx = new TProfile("pTemp_track_dedx","pTemp_track_dedx",10,0.5,10.5,0,200); //pT, eta,phi,charge,DCA,
        TProfile *pT_Day3 = new TProfile("pT_Day3","pT_Day3",run_end-run_sta,run_sta,run_end,0,10);
        TProfile *Eta_Day3 = new TProfile("Eta_Day3","Eta_Day3",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *Phi_Day3 = new TProfile("Phi_Day3","Phi_Day3",run_end-run_sta,run_sta,run_end,-4,4);
        TProfile *Charge_Day3 = new TProfile("Charge_Day3","Charge_Day3",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *DCA_Day3 = new TProfile("DCA_Day3","DCA_Day3",run_end-run_sta,run_sta,run_end,0,4);
	TH2D *Ref_TOF = new TH2D("Ref_TOF","Ref_TOF",500,0.5,500.5,5000,0.5,5000.5);
        TProfile *Ref_Day3 = new TProfile("Ref_Day3","RefMult vs Run",run_end-run_sta,run_sta,run_end,0,999);
        TProfile *TOF_Day3 = new TProfile("TOF_Day3","TOFMult vs Run",run_end-run_sta,run_sta,run_end,0,5000);
        TProfile *NPT_Day3 = new TProfile("NPT_Day3","NPTracks vs Run",run_end-run_sta,run_sta,run_end,0,5000);
        TProfile *NPA_Day3 = new TProfile("NPA_Day3","NPAsso vs Run",run_end-run_sta,run_sta,run_end,0,5000);
        TProfile *TPC_Day3_cos2 = new TProfile("TPC_Day3_cos2","cos(2*psi) vs Run", run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *TPC_Day3_sin2 = new TProfile("TPC_Day3_sin2","sin(2*psi) vs Run", run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *BBCe_Day3_cos2= new TProfile("BBCe_Day3_cos2","BBCe_Day3_cos2",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *BBCe_Day3_sin2= new TProfile("BBCe_Day3_sin2","BBCe_Day3_sin2",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *BBCw_Day3_cos2= new TProfile("BBCw_Day3_cos2","BBCw_Day3_cos2",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *BBCw_Day3_sin2= new TProfile("BBCw_Day3_sin2","BBCw_Day3_sin2",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPDe_Day3_cos1= new TProfile("EPDe_Day3_cos1","EPDe_Day3_cos1",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPDe_Day3_sin1= new TProfile("EPDe_Day3_sin1","EPDe_Day3_sin1",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPDw_Day3_cos1= new TProfile("EPDw_Day3_cos1","EPDw_Day3_cos1",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPDw_Day3_sin1= new TProfile("EPDw_Day3_sin1","EPDw_Day3_sin1",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPDe_Day3_cos2= new TProfile("EPDe_Day3_cos2","EPDe_Day3_cos2",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPDe_Day3_sin2= new TProfile("EPDe_Day3_sin2","EPDe_Day3_sin2",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPDw_Day3_cos2= new TProfile("EPDw_Day3_cos2","EPDw_Day3_cos2",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPDw_Day3_sin2= new TProfile("EPDw_Day3_sin2","EPDw_Day3_sin2",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *ZDCe_Day3_cos= new TProfile("ZDCe_Day3_cos","ZDCe_Day3_cos",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *ZDCe_Day3_sin= new TProfile("ZDCe_Day3_sin","ZDCe_Day3_sin",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *ZDCw_Day3_cos= new TProfile("ZDCw_Day3_cos","ZDCw_Day3_cos",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *ZDCw_Day3_sin= new TProfile("ZDCw_Day3_sin","ZDCw_Day3_sin",run_end-run_sta,run_sta,run_end,-1,1);
	TProfile *BBC_cor_Day3 = new TProfile("BBC_cor_Day3","BBC_cor_Day3",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPD1_cor_Day3 = new TProfile("EPD1_cor_Day3","EPD1_cor_Day3",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *EPD_cor_Day3 = new TProfile("EPD_cor_Day3","EPD_cor_Day3",run_end-run_sta,run_sta,run_end,-1,1);
        TProfile *ZDC_cor_Day3 = new TProfile("ZDC_cor_Day3","ZDC_cor_Day3",run_end-run_sta,run_sta,run_end,-1,1);

	TH1D *hTally = new TH1D("hTally","hTally",11,-0.5,10.5);
        TH1D *hTall  = new TH1D("hTall ","hTall ",12,-0.5,11.5);
	TH1D *hZDCcoin = new TH1D("hZDCcoin","hZDCcoin",1000,0,100000);
        TH1D *hTrigger = new TH1D("hTrigger","hTrigger",200, 0.5, 200.5);
        TH1D *hCentrality = new TH1D("hCentrality","hCentrality",11,-1.5,9.5);
	TH2D *hVertexXY = new TH2D("hVertexXY","hVertexXY",300,-3,3, 300,-3,3);
        TH1D *hVertexZ = new TH1D("hVertexZ","hVertexZ",500,-250,250);
	TH1D *hVzDiff = new TH1D("hVzDiff","hVzDiff",500,-250,250);
        TH2D *hMult_Vz = new TH2D("hMult_Vz","hMult_Vz",1000,-0.5,999.5,500,-250,250);
        TH2D *hMult_Vz_new = new TH2D("hMult_Vz_new","hMult_Vz_new",1000,-0.5,999.5,500,-250,250);
        TH1D *BBC1= new TH1D("BBC1","BBC1",6,0.5,6.5);
        TH1D *BBC2= new TH1D("BBC2","BBC2",10,0.5,10.5);
        TH1D *BBC7= new TH1D("BBC7","BBC7",10,0.5,10.5);
        TH1D *BBC8= new TH1D("BBC8","BBC8",6,0.5,6.5);
	TH1D *ZDC_e_h = new TH1D("ZDC_e_h","ZDC_e_h",8,0.5,8.5);
        TH1D *ZDC_e_v = new TH1D("ZDC_e_v","ZDC_e_v",8,0.5,8.5);
        TH1D *ZDC_w_h = new TH1D("ZDC_w_h","ZDC_w_h",8,0.5,8.5);
        TH1D *ZDC_w_v = new TH1D("ZDC_w_v","ZDC_w_v",8,0.5,8.5);
	TProfile2D *pZDCcenter = new TProfile2D("pZDCcenter","pZDCcenter",4,0.5,4.5,(run_end-run_sta)/10,run_sta/10,run_end/10, 0,20);
        TProfile2D *pZDCcenter_new = new TProfile2D("pZDCcenter_new","pZDCcenter_new",4,0.5,4.5,(run_end-run_sta)/10,run_sta/10,run_end/10, -20,20);

        TH1D* Hist_positive = new TH1D("Hist_positive","Hist_positive",500,-0.5,499.5);
        TH1D* Hist_negative = new TH1D("Hist_negative","Hist_negative",500,-0.5,499.5);
        TH1D* Hist_Ch       = new TH1D("Hist_Ch","Hist_Ch",1000,-0.5,999.5);
        TH1D* Hist_netCh    = new TH1D("Hist_netCh","Hist_netCh",999,-499.5,499.5);
        TH1D* Hist_netChAsym= new TH1D("Hist_netChAsym","Hist_netChAsym",600,-1.5+0.0025,1.5+0.0025);
        TH2D* Hist_netChAsym_bin    = new TH2D("Hist_netChAsym_bin","Hist_netChAsym_bin",5,0.5,5.5,600,-1.5+0.0025,1.5+0.0025);
        TProfile *p_netChAsym_RefMult = new TProfile("p_netChAsym_RefMult","p_netChAsym_RefMult",300,-1.5+0.0025,1.5+0.0025, 0.,900);
        TProfile *p_netChAsym_cos     = new TProfile("p_netChAsym_cos","p_netChAsym_cos",300,-1.5+0.0025,1.5+0.0025, -1,1);
        TProfile *Hist_cos_Ach = new TProfile("Hist_cos_Ach","Hist_cos_Ach",5,0.5,5.5,-1,1);
	TProfile *p_v2_Ach = new TProfile("p_v2_Ach","p_v2_Ach",300,-1.5+0.0025,1.5+0.0025, -100,100);
        TH2D *Hist_pt_pos_Ach = new TH2D("Hist_pt_pos_Ach","Hist_pt_pos_Ach",5,0.5,5.5,300,0,15);
	TH2D *Hist_pt_neg_Ach = new TH2D("Hist_pt_neg_Ach","Hist_pt_neg_Ach",5,0.5,5.5,300,0,15);
        TProfile2D *p_v2_pt_pos_Ach = new TProfile2D("p_v2_pt_pos_Ach","p_v2_pt_pos_Ach",5,0.5,5.5,300,0,15,-100,100);
        TProfile2D *p_v2_pt_neg_Ach = new TProfile2D("p_v2_pt_neg_Ach","p_v2_pt_neg_Ach",5,0.5,5.5,300,0,15,-100,100);

        TH2D* Hist_TPC_EP_east = new TH2D("Hist_TPC_EP_east","Hist_TPC_EP_east",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_TPC_EP_west = new TH2D("Hist_TPC_EP_west","Hist_TPC_EP_west",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_TPC_EP_full = new TH2D("Hist_TPC_EP_full","Hist_TPC_EP_full",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_TPC_EP_for = new TH2D("Hist_TPC_EP_for","Hist_TPC_EP_for",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_TPC_EP_bac = new TH2D("Hist_TPC_EP_bac","Hist_TPC_EP_bac",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_BBC_EP_east = new TH2D("Hist_BBC_EP_east","Hist_BBC_EP_east",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_BBC_EP_west = new TH2D("Hist_BBC_EP_west","Hist_BBC_EP_west",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_EPD_EP1_east = new TH2D("Hist_EPD_EP1_east","Hist_EPD_EP1_east",72,0,2*PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_EPD_EP1_west = new TH2D("Hist_EPD_EP1_west","Hist_EPD_EP1_west",72,0,2*PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_EPD_EP_east = new TH2D("Hist_EPD_EP_east","Hist_EPD_EP_east",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_EPD_EP_west = new TH2D("Hist_EPD_EP_west","Hist_EPD_EP_west",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
	TH2D* Hist_ZDC_EP_east = new TH2D("Hist_ZDC_EP_east","Hist_ZDC_EP_east",72,-PI,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_ZDC_EP_west = new TH2D("Hist_ZDC_EP_west","Hist_ZDC_EP_west",72,-PI,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);

        TH2D* Hist_TPC_EP_east_flat = new TH2D("Hist_TPC_EP_east_flat","Hist_TPC_EP_east_flat",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_TPC_EP_west_flat = new TH2D("Hist_TPC_EP_west_flat","Hist_TPC_EP_west_flat",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_TPC_EP_full_flat = new TH2D("Hist_TPC_EP_full_flat","Hist_TPC_EP_full_flat",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_TPC_EP_for_flat = new TH2D("Hist_TPC_EP_for_flat","Hist_TPC_EP_for_flat",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_TPC_EP_bac_flat = new TH2D("Hist_TPC_EP_bac_flat","Hist_TPC_EP_bac_flat",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_BBC_EP_east_flat = new TH2D("Hist_BBC_EP_east_flat","Hist_BBC_EP_east_flat",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_BBC_EP_west_flat = new TH2D("Hist_BBC_EP_west_flat","Hist_BBC_EP_west_flat",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_EPD_EP1_east_flat = new TH2D("Hist_EPD_EP1_east_flat","Hist_EPD_EP1_east_flat",72,-PI,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_EPD_EP1_west_flat = new TH2D("Hist_EPD_EP1_west_flat","Hist_EPD_EP1_west_flat",72,-PI,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_EPD_EP_east_flat = new TH2D("Hist_EPD_EP_east_flat","Hist_EPD_EP_east_flat",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_EPD_EP_west_flat = new TH2D("Hist_EPD_EP_west_flat","Hist_EPD_EP_west_flat",36,0,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_ZDC_EP_east_flat = new TH2D("Hist_ZDC_EP_east_flat","Hist_ZDC_EP_east_flat",72,-PI,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);
        TH2D* Hist_ZDC_EP_west_flat = new TH2D("Hist_ZDC_EP_west_flat","Hist_ZDC_EP_west_flat",72,-PI,PI,(run_end-run_sta)/1000,run_sta/1000,run_end/1000);

	TH1D* Hist_TPC_EP_full_m1 = new TH1D("Hist_TPC_EP_full_m1","Hist_TPC_EP_full - phi 1",36,0,PI);
        TH1D* Hist_TPC_EP_full_m2 = new TH1D("Hist_TPC_EP_full_m2","Hist_TPC_EP_full - phi 1 - phi 2",36,0,PI);
        TH1D* Hist_TPC_EP_full_m1_flat = new TH1D("Hist_TPC_EP_full_m1_flat","Hist_TPC_EP_full - phi 1, flat used for v2",36,0,PI);
	TH1D* Hist_TPC_EP_full_m2_flat = new TH1D("Hist_TPC_EP_full_m2_flat","Hist_TPC_EP_full - phi 1 - phi 2, flat used for gamma112",36,0,PI);
        TProfile2D* pTPC_EP_east = new TProfile2D("pTPC_EP_east","pTPC_EP_east",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pTPC_EP_west = new TProfile2D("pTPC_EP_west","pTPC_EP_west",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pTPC_EP_full = new TProfile2D("pTPC_EP_full","pTPC_EP_full",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pTPC_EP_for = new TProfile2D("pTPC_EP_for","pTPC_EP_for",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pTPC_EP_bac = new TProfile2D("pTPC_EP_bac","pTPC_EP_bac",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pBBC_EP_east = new TProfile2D("pBBC_EP_east","pBBC_EP_east",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pBBC_EP_west = new TProfile2D("pBBC_EP_west","pBBC_EP_west",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pEPD_EP1_east = new TProfile2D("pEPD_EP1_east","pEPD_EP1_east",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pEPD_EP1_west = new TProfile2D("pEPD_EP1_west","pEPD_EP1_west",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pEPD_EP_east = new TProfile2D("pEPD_EP_east","pEPD_EP_east",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pEPD_EP_west = new TProfile2D("pEPD_EP_west","pEPD_EP_west",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pZDC_EP_east = new TProfile2D("pZDC_EP_east","pZDC_EP_east",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
        TProfile2D* pZDC_EP_west = new TProfile2D("pZDC_EP_west","pZDC_EP_west",2*order,0.5,2*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");

	TH2D* Hist_BBC_vs_TPC = new TH2D("Hist_BBC_vs_TPC","Hist_BBC_vs_TPC",30,0,PI,30,0,PI);
        TH2D* Hist_EPD_vs_TPC = new TH2D("Hist_EPD_vs_TPC","Hist_EPD_vs_TPC",30,0,PI,30,0,PI);

        TProfile *Hist_cos = new TProfile("Hist_cos","Hist_cos",4,0.5,4.5,-1,1,"");
        TProfile *Hist_cos_BBC = new TProfile("Hist_cos_BBC","Hist_cos_BBC",4,0.5,4.5,-1,1,"");
        TProfile *Hist_cos_EPD = new TProfile("Hist_cos_EPD","Hist_cos_EPD",7,0.5,7.5,-1,1,"");
        TProfile *Hist_cos_ZDC = new TProfile("Hist_cos_ZDC","Hist_cos_ZDC",4,0.5,4.5,-1,1,"");

	TH1D *Hist_DCA = new TH1D("Hist_DCA","Hist_DCA",100,0,10);
        TH2D *hEtaPtDist = new TH2D("EtaPtDist","EtaPtDist",34, -1.7, 1.7,300,0,15);
	TH2D *hEtaPhiDist= new TH2D("hEtaPhiDist","hEtaPhiDist",34, -1.7, 1.7,Phibin,-PI,PI);
        TH2D *hPhiPtDist = new TH2D("PhiPtDist","PhiPtDist",Phibin,-PI,PI,300,0,15);
        TH1D* Hist_Pt = new TH1D("Hist_Pt","Hist_Pt",300,0,15);
        TH1D* Hist_Pt_TOF = new TH1D("Hist_Pt_TOF","Hist_Pt_TOF",300,0,15);
	TH1D* rc;
	TH2D* Hist_nSigma_pi = new TH2D("Hist_nSigma_pi","Hist_nSigma_pi", 40, 0, 2, 200, -10,10);
        TH2D* Hist_nSigma_pi_TOF = new TH2D("Hist_nSigma_pi_TOF","Hist_nSigma_pi_TOF", 40, 0, 2, 200, -10,10);
        TH2D* Hist_nSigma_proton = new TH2D("Hist_nSigma_proton","Hist_nSigma_proton", 40, 0, 2, 200, -10,10);
	TH2D* wt = new TH2D("Order1etaWeight","Order1etaWeight",100,1.5,6.5,9,0,9);
        TH2D* wt2= new TH2D("Order2etaWeight","Order2etaWeight",100,1.5,6.5,9,0,9);
        TH2D* Hist_Phi = new TH2D("Hist_Phi","Hist_Phi",Phibin,-PI,PI,4,0.5,4.5);
        TH1D* hDpt   = new TH1D("hDpt","hDpt",200,0,2);

        TH1D* hCParent = new TH1D("hCParent","hCParent",2,0.5,2.5);

// yu's parameter for QA
// Run by run QA
TProfile * runidvsrefmult = new TProfile("runidvsrefmult","Run Id-RefMult",runbins,runmin,runmax,0,1000);
TProfile * runidvsepdEhits = new TProfile("runidvsepdEhits","Run Id-EPD East",runbins,runmin,runmax,0,5000);
TProfile * runidvsepdWhits = new TProfile("runidvsepdWhits","Run Id-EPD West",runbins,runmin,runmax,0,5000);
//TProfile * runidvszdcand = new TProfile("runidvszdcand","Run Id-ZDC coincidence",runbins,runmin,runmax,0,5000);
TProfile * runidvstofmult = new TProfile("runidvstofmult","Run Id-TOFMult",runbins,runmin,runmax,0,5000);
TProfile * runidvstofmatched = new TProfile("runidvstofmatched","Run Id-TOFMatched",runbins,runmin,runmax,0,5000);
TProfile * runidvsnptrack = new TProfile("runidvsnptrack","Run Id-track number",runbins,runmin,runmax,0,5000);
TProfile * runidvsdedx = new TProfile("runidvsdedx","Run Id-dedx",runbins,runmin,runmax,0,200);

TProfile * runidvsnpasso = new TProfile("runidvsnpasso","Run Id-asso number",runbins,runmin,runmax,0,5000);
TProfile * runidvsavgpt = new TProfile("runidvsavgpt","Run Id-AvgpT",runbins,runmin,runmax,0,10);//event averaged pt
TProfile * runidvsavgeta = new TProfile("runidvsavgeta","Run Id-Avgeta",runbins,runmin,runmax,-2,2);
TProfile * runidvsavgdca = new TProfile("runidvsavgdca","Run Id-Avgdca",runbins,runmin,runmax,0,10);
TProfile * runidvsgdcaxy = new TProfile("runidvsgdcaxy","Run Id-gdca xy",runbins,runmin,runmax,-20,20);
TProfile * runidvsdcaxysigma = new TProfile("runidvsdcaxysigma","Run Id-gdca xy sigma",runbins,runmin,runmax,-20,20);
TProfile * runidvsavgvz = new TProfile("runidvsavgvz","Run Id-Avgvz",runbins,runmin,runmax,-250,250);
TProfile * runidvsavgphi = new TProfile("runidvsavgphi","Run Id-Avgphi",runbins,runmin,runmax,-4,4);
TProfile * runidvsavgQ1x = new TProfile("runidvsavgQ1x","Run Id-AvgQ1x",runbins,runmin,runmax,-100,100);
TProfile * runidvsavgQ1y = new TProfile("runidvsavgQ1y","Run Id-AvgQ1y",runbins,runmin,runmax,-100,100);
TProfile * runidvsavgQ2x = new TProfile("runidvsavgQ2x","Run Id-AvgQ2x",runbins,runmin,runmax,-100,100);
TProfile * runidvsavgQ2y = new TProfile("runidvsavgQ2y","Run Id-AvgQ2y",runbins,runmin,runmax,-100,100);
TProfile * runidvsTPCcos = new TProfile("runidvsTPCcos","Run Id-TPCcos",runbins,runmin,runmax,-1,1);
TProfile * runidvsTPCsin = new TProfile("runidvsTPCsin","Run Id-TPCsin",runbins,runmin,runmax,-1,1);
TProfile * runidvsEPD2cos = new TProfile("runidvsEPDcos","Run Id-2nd EPDcos",runbins,runmin,runmax,-1,1);
TProfile * runidvsEPD2sin = new TProfile("runidvsEPDsin","Run Id-2nd EPDsin",runbins,runmin,runmax,-1,1);
TProfile * runidvsEPD1cos = new TProfile("runidvsEPD1cos","Run Id-1st EPDcos",runbins,runmin,runmax,-1,1);
TProfile * runidvsEPD1sin = new TProfile("runidvsEPD1sin","Run Id-1st EPDsin",runbins,runmin,runmax,-1,1);


TProfile * runidvsavgEpdQ1x = new TProfile("runidvsavgEpdQ1x","Run Id-AvgEpdQ1x",runbins,runmin,runmax,-500,500);
TProfile * runidvsavgEpdQ1y = new TProfile("runidvsavgEpdQ1y","Run Id-AvgEpdQ1y",runbins,runmin,runmax,-500,500);
TProfile * runidvsavgEpdQ2x = new TProfile("runidvsavgEpdQ2x","Run Id-AvgEpdQ2x",runbins,runmin,runmax,-500,500);
TProfile * runidvsavgEpdQ2y = new TProfile("runidvsavgEpdQ2y","Run Id-AvgEpdQ2y",runbins,runmin,runmax,-500,500);



TH1D* Hist_mEpdHits = new TH1D("Hist_mEpdHits", "Hist_mEpdHits", 1000,0,1000); 



//defining variables
float pVx, pVy, pVz, VPDvz, BBCco, ZDCcoin, net_Nch_Asym, mQx, mQy;                   //run, event info
int   Run, Day, Day2, Day3, Trigger, RefMult, TOFMult, Ntofmatch, Centrality, NPTracks, Fcount, Scount, POIcount;   //
UShort_t FxtMult;
int   Charge, Charge2, ChargeAsso;
float ndEdx, nSigma_p, nSigma_pi, DCAGlobal, Eta, Theta, Phi, Pt, eff, TOFflag;             //track info    
float ndEdx2, nSigma_p2, nSigma_pi2, DCAGlobal2, Eta2, Theta2, Phi2, Pt2, eff2, TOFflag2;   //2nd track info
float DCAGlobalAsso, EtaAsso, PhiAsso, PtAsso, TOFflagAsso;
TVector3 pV;
float gDCAxy,dEdx;
float Eweight = 1;
int Weight_Read = 0;
TH1D *TOF_eff=0;
float PhiMean_sin[order]={0},PhiMean_cos[order]={0}, PhiMeanAsso_sin[order]={0}, PhiMeanAsso_cos[order]={0};
vector<float> PhiAsso_new, Phi_new;
vector<int> iCharge, match;
TProfile2D *Read_TPC_EP_full=0, *Read_TPC_EP_east=0, *Read_TPC_EP_west=0, *Read_TPC_EP_for=0, *Read_TPC_EP_bac=0;
float PsiMean_F[2*order]={0},PsiMean_E[2*order]={0},PsiMean_W[2*order]={0},PsiMean_f[2*order]={0},PsiMean_b[2*order]={0};
TProfile2D *Read_BBC_EP_east=0, *Read_BBC_EP_west=0, *Read_EPD_EP_east=0, *Read_EPD_EP_west=0, *Read_EPD_EP1_east=0, *Read_EPD_EP1_west=0;
float Psi_BBC_E[2*order]={0},Psi_BBC_W[2*order]={0},Psi_EPD_E[2*order]={0},Psi_EPD_W[2*order]={0},Psi1_EPD_E[2*order]={0},Psi1_EPD_W[2*order]={0};
float MeanNetChargeAsym, RMSNetChargeAsym;
float BBC_gain_east[16]={0}, BBC_gain_west[16]={0};
float ZDC_gain_east_v[7] = {0}, ZDC_gain_east_h[8] = {0}, ZDC_gain_west_v[7] = {0}, ZDC_gain_west_h[8] = {0};
TProfile2D *Read_ZDCcenter=0, *Read_ZDC_EP_east=0, *Read_ZDC_EP_west=0;
float zdcsmd_x0_e = 0, zdcsmd_y0_e = 0, zdcsmd_x0_w = 0, zdcsmd_y0_w = 0;
float Psi_ZDC_E[2*order]={0},Psi_ZDC_W[2*order]={0};

float TPC_EP_full = 0, TPC_EP_east = 0, TPC_EP_west = 0, TPC_EP_for = 0, TPC_EP_bac = 0;
float TPC_EP_full_new = 0, TPC_EP_east_new = 0, TPC_EP_west_new = 0, TPC_EP_for_new = 0, TPC_EP_bac_new = 0;
float BBC_EP_east = 0, BBC_EP_west = 0, EPD_EP_east = 0, EPD_EP_west = 0, ZDC_EP_east = 0, ZDC_EP_west = 0;
float BBC_EP_east_new=0, BBC_EP_west_new=0, EPD_EP_east_new=0, EPD_EP_west_new=0, ZDC_EP_east_new=0, ZDC_EP_west_new = 0;
float EPD_EP1_east = 0, EPD_EP1_west = 0, EPD_EP1_east_new = 0, EPD_EP1_west_new = 0;
float  Q2 = 0;
int N_L_P1 = 0, N_L_N1 = 0, N_R_P1 = 0, N_R_N1 = 0, N_T_P1 = 0, N_T_N1 = 0, N_B_P1 = 0, N_B_N1 = 0;
int N_L_P2 = 0, N_L_N2 = 0, N_R_P2 = 0, N_R_N2 = 0, N_T_P2 = 0, N_T_N2 = 0, N_B_P2 = 0, N_B_N2 = 0;

//StRefMultCorr* refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;

//sub-functions
void DefineHistogram();						//define histograms 
void WriteHistogram(int c, int o);  			//into result ROOT file
void WriteWeight(char* OutFileName);		   		//into weight file
int  ReadWeight(char* InFileName);			 	//input from weight file
bool IsGoodEvent(int c);		   			//select good events and fill event level histograms
bool IsGoodBBC(StPicoEvent *e);            			//cuts on BBC ADCs
bool IsGoodZDC(StPicoEvent *e);            			//cuts on ZDC ADCs
bool IsGoodAsso(StPicoTrack *p);			//cuts on EP particles			
bool IsGoodPOI(StPicoTrack *p); 			//cuts on particles of interest
bool IsGoodTrack(StPicoTrack *p);	   			//cuts on tracks
bool IsGoodPionNoTOF(StPicoTrack *p);				//cuts on low pT pions without TOF
bool IsGoodPion(StPicoDst *d, StPicoTrack *p, int opt); 	//cuts on pions
bool IsGoodKaon(StPicoDst *d, StPicoTrack *p, int opt);         //cuts on kaons
bool IsGoodProton(StPicoDst *d, StPicoTrack *p, int opt);       //cuts on protons
bool CountCharge(StPicoDst *d);					//count good tracks, charge asymmetry
void MakeTPC_EP(StPicoDst *d, int *iTr);			//reconstruct TPC EPs
void MakeBBC_EP(StPicoEvent *ev);				//reconstruct BBC EPs
void MakeZDC_EP(StPicoEvent *ev);				//reconstruct ZDC EPs
void FillEP_resolution();					//Fill profiles for EP resolution
void FillPhiPOI();  			   			//particles of interest, shift parameters to make phi distribution flat
void FillPhiAsso();   			   			//particles for EP, shift parameters to make phi distribution flat
void ShiftPhiAsso(int tr);   		   			//flatten the phi distribution
void ShiftPhiPOI(int tr);		   			//flatten the phi distribution
void ShiftPsi();			   			//flatten EPs
float GetPhiInBBC(int e_w, int bbcN);      			//input est_wst(0,1),BBC PMT#
float GetXYInZDC(int e_w, int v_h, int zdcN, int opt_raw = 0);	//input ver_hor(0,1),ZDCSMD slat#
void EPD_hits(TClonesArray *mEpdHits);

//void __attribute__((constructor)) LoadLib();  
/////////////////////////////main program starts here/////////////////////////////////////
void Gamma_QA(int cen=1, int opt_weight =1, const Char_t *inFile = "test.list"){	//main_function
	delete gRandom;
	gRandom = new TRandom(0);

	gSystem->Load("StUtilities");
	gSystem->Load("StRefMultCorr");
	gSystem->Load("StPicoEvent");
	gSystem->Load("StPicoDstMaker");

	DefineHistogram();
        char fname_old[200], fname_new[200];
        sprintf(fname_new,"cen%d.weight_1%d%d_QA_new.root",cen,nHar-1,nHar);
        sprintf(fname_old,"cen%d.weight_1%d%d_QA.root",cen,nHar-1,nHar);
	Weight_Read = ReadWeight(fname_old);

        cout<<"Ready to read"<<endl;
  	StPicoDstReader* picoReader = new StPicoDstReader(inFile);
        cout<<"Finish reading"<<endl;
  	picoReader->Init();
        cout<<"Finish initting"<<endl;
  	if( !picoReader->chain() ) { std::cout << "No chain has been found." << std::endl; }
	Int_t nentries = picoReader->chain()->GetEntries();
	
	StRefMultCorr* refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;

  	// This is a way if you want to spead up IO
        std::cout << "Explicit read status for some branches" << std::endl;
        picoReader->SetStatus("*",0);
        picoReader->SetStatus("Event",1);
        picoReader->SetStatus("Track",1);
        picoReader->SetStatus("BTofHit",1);
        picoReader->SetStatus("BTofPidTraits",1);
        cout<<"1st"<<endl;

	// Prepare EPD
	TClonesArray* mEpdHits = new TClonesArray("StPicoEpdHit");
        unsigned int found;
	picoReader->chain()->SetBranchStatus("EpdHit*",1,&found);
	cout << "StPicoDstMaker::SetStatus "<< found <<" to EpdHit" << endl;
	picoReader->chain()->SetBranchAddress("EpdHit",&mEpdHits);
	StEpdEpFinder* mEpFinder = new StEpdEpFinder(9,fname_new,fname_old);
  	mEpFinder->SetnMipThreshold(0.3);    	// recommended by EPD group
  	mEpFinder->SetMaxTileWeight(2.0);     	// recommended by EPD group
  	mEpFinder->SetEpdHitFormat(2);         	// 2=pico   
	mEpFinder->SetEtaWeights(1,*wt);		// eta weight for 1st-order EP
        mEpFinder->SetEtaWeights(2,*wt2);	// eta weight for 2nd-order EP, select different eta range

	//loop through events
	for(int i = 0; i < nentries; i++){

		if((i+1)%1000==0) cout<<"Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
		Bool_t readEvent = picoReader->readPicoEvent(i);
    		if( !readEvent ) {
      			cout << "Something went wrong, Master! Nothing to analyze..." << endl;
  	    		break;
    		}
		// Retrieve picoDst
		StPicoDst *dst = picoReader->picoDst();
		// Retrieve event information
		StPicoEvent *event = dst->event();
    		if( !event ) {
      			cout << "Something went wrong, Master! Event is hiding from me..." << endl;
      			break;
    		}
        //cout<<"2nd"<<endl;
		hTally->Fill(0);
                // if(!event->isTrigger(710000)&& !event->isTrigger(710010) &&!event->isTrigger(710020)) continue; //Updated 19June 2023 for 11GeV v1
                if( !event->isTrigger(630052) ) continue; // AuAu26p5

		Run	= event->runId();
		pV 	= event->primaryVertex();
		pVz	= pV.Z();
		pVx	= pV.X();
		pVy	= pV.Y();
		VPDvz   = event->vzVpd();
		RefMult = event->refMult();
		FxtMult = event->fxtMult();
		TOFMult = event->btofTrayMultiplicity();
		Ntofmatch = event->nBTOFMatch();
		NPTracks= dst->numberOfTracks();
		BBCco   = event->BBCx();
                ZDCcoin = event->ZDCx();
		Day 	= (int)((Run%1000000)/1000); 
		Day2    = (int)((Run%1000000)/10);
                Day3    = (int)((Run%1000000)/1);

                // Siyuan Ping: Reject Bad Run
                if (Run == 19168001) continue;
                if (Run == 19168002) continue;
                if (Run == 19168003) continue;
                if (Run == 19168016) continue;
                if (Run == 19167054) continue;
                if (Run == 19166003) continue;
                if (Run == 19164002) continue;
                if (Run == 19164021) continue;
                if (Run == 19161031) continue;
                if (Run == 19161032) continue;
                if (Run == 19161033) continue;
                if (Run == 19159042) continue;
                if (Run == 19158053) continue;
                if (Run == 19158054) continue;
                if (Run == 19158055) continue;
                if (Run == 19158056) continue;
                if (Run == 19157033) continue;
                if (Run == 19157034) continue;
                if (Run == 19157035) continue;
                if (Run == 19157036) continue;
                if (Run == 19157037) continue;
                if (Run == 19157038) continue;
                if (Run == 19157039) continue;
                if (Run == 19157040) continue;
                if (Run == 19157041) continue;
                if (Run == 19157042) continue;
                if (Run == 19157043) continue;
                if (Run == 19156034) continue;
                if (Run == 19156035) continue;
                if (Run == 19156036) continue;
                if (Run == 19156038) continue;
                if (Run == 19156039) continue;
                if (Run == 19156069) continue;

		// if(Run<21001000) continue;

/*                refmultCorrUtil->init(Run);
                if ( refmultCorrUtil->isBadRun(Run) ) continue;
		hTally->Fill(1);
                refmultCorrUtil->initEvent(RefMult, pVz, ZDCcoin) ;
                if(!refmultCorrUtil->passnTofMatchRefmultCut(1.*RefMult, 1.*Ntofmatch)) continue;
                hTally->Fill(2);
                Centrality  = 1 + refmultCorrUtil->getCentralityBin9() ;
                Eweight = refmultCorrUtil->getWeight();
                if((i+1)%1000==0) cout<<"centraility = "<<Centrality<<" Eweight = "<<Eweight<<endl;
*/
                // PSY: Loop for Nch to check refmult
                Int_t nTracks = dst->numberOfTracks();
                int NumCharge = 0;
                for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
                        StPicoTrack *track = dst->track(iTrack);
                        if (! track)            continue;
                        if (! track->charge())  continue;
                        if (! track->isPrimary()) continue;
                        NumCharge++;
                        if (  track->nHitsFit() < 10) continue;
                        if (  (track->gDCA(pV)).Mag() > 3.0) continue;
                        // if (  fabs(track->gMom().Eta()) > 1.5) continue;
                        runidvsavgeta->Fill(Run,pTemp_track->GetBinContent(2));
                        if (  (0.06 > track->gMom().Perp()) || (track->gMom().Perp() > 2.0)) continue;
                }
                RefMult = NumCharge;
                // RefMult = FxtMult;
                // cout<<"FxtMult = "<<FxtMult<<endl;


///temp add
	Centrality = 0;
	for(int j=0;j<9;j++) if(RefMult>cenDef[j]) Centrality = j+1;
	Eweight = 1;

//temp stop
//cout<<"3rd"<<endl;
		if(!IsGoodEvent(cen)) continue;
		if(!CountCharge(dst)) continue;;
//shuffle tracks for random EPs
          	int iTrack[Fcount];
	  	Scount = Fcount/2 -1;
          	for(int q=0;q<Fcount;q++) iTrack[q] = q;
          	random_shuffle(iTrack,iTrack+Fcount);

//TPC EP reconstruction
	  	MakeTPC_EP(dst,iTrack);
		random_shuffle(iCharge.begin(),iCharge.end()); //shuffle charge for POI
		if(TPC_EP_for<-9 || TPC_EP_bac<-9) continue;
                if(opt_ESE==1 && (sqrt(Q2)<0.4 || sqrt(Q2)>1.6)) continue;
//BBC EP
          	//if(!IsGoodBBC(event)) continue;
	  	MakeBBC_EP(event);

//EPD EP
	  	StEpdEpInfo result = mEpFinder->Results(mEpdHits,pV,(Centrality>0)? Centrality-1:0);
                EPD_EP1_east = result.EastPhiWeightedPsi(1);
                EPD_EP1_west = result.WestPhiWeightedPsi(1);
	  	EPD_EP_east = result.EastPhiWeightedPsi(nHar);
	  	EPD_EP_west = result.WestPhiWeightedPsi(nHar);
//if(EPD_EP1_west > 0.025 && EPD_EP1_west < 0.035) cout<<EPD_EP1_east<<" "<<EPD_EP1_west<<endl;
if((opt_useEPD == 1 || opt_useEPD == 11) && (EPD_EP_east == EPD_EP_west || (EPD_EP1_east > 0.0264 && EPD_EP1_east < 0.0265) || (EPD_EP1_west >0.0264 && EPD_EP1_west < 0.0265))) continue;
                EPDe_Day3_cos1->Fill(Day3, cos(EPD_EP1_east));
                EPDe_Day3_sin1->Fill(Day3, sin(EPD_EP1_east));
                EPDw_Day3_cos1->Fill(Day3, cos(EPD_EP1_west));
                EPDw_Day3_sin1->Fill(Day3, sin(EPD_EP1_west));
	  	EPDe_Day3_cos2->Fill(Day3, cos(nHar*EPD_EP_east));
          	EPDe_Day3_sin2->Fill(Day3, sin(nHar*EPD_EP_east));
	  	EPDw_Day3_cos2->Fill(Day3, cos(nHar*EPD_EP_west));
	  	EPDw_Day3_sin2->Fill(Day3, sin(nHar*EPD_EP_west));

                double Q1xEpdE_in=result.EastPhiWeightedQ(1).X(); //innerEPD cent 0-8
                double Q1yEpdE_in=result.EastPhiWeightedQ(1).Y(); //innerEPD cent 0-8
                double Q1xEpdW_in=result.WestPhiWeightedQ(1).X(); //innerEPD cent 0-8
                double Q1yEpdW_in=result.WestPhiWeightedQ(1).Y(); //innerEPD cent 0-8

                double Q2xEpdE_out=result.EastPhiWeightedQ(2).X(); //outerEPD cent 9-17
                double Q2yEpdE_out=result.EastPhiWeightedQ(2).Y(); //outerEPD cent 9-17
                double Q2xEpdW_out=result.WestPhiWeightedQ(2).X(); //outerEPD cent 9-17
                double Q2yEpdW_out=result.WestPhiWeightedQ(2).Y(); //outerEPD cent 9-17
                //cout<<"EP ";
                //cout << result.EastPhiWeightedPsi(1)<< " "<<result.EastPhiWeightedPsi(nHar)<<endl;
		//cout<<"flow ";
		//cout<<Q1xEpdE_in << " "<<Q1xEpdW_in<<" "<<Q1yEpdE_in <<" "<< Q1yEpdW_in<<endl;
		//cout<< result.EastRawQ(1).Px()<< " "<< result.EastRawQ(1).Py()<< " "<<result.EastRawQ(2).Px()<< " "<< result.EastRawQ(2).Py()<<endl;
                runidvsavgEpdQ1x->Fill(Run,0.5*(Q1xEpdE_in + Q1xEpdW_in));
                runidvsavgEpdQ1y->Fill(Run,0.5*(Q1yEpdE_in + Q1yEpdW_in));
                runidvsavgEpdQ2x->Fill(Run,0.5*(Q2xEpdE_out + Q2xEpdW_out));
                runidvsavgEpdQ2y->Fill(Run,0.5*(Q2yEpdE_out + Q2yEpdW_out));

		runidvsEPD2cos->Fill(Run,0.5*(cos(nHar*EPD_EP_east)+ cos(nHar*EPD_EP_west)));
                runidvsEPD2sin->Fill(Run,0.5*(sin(nHar*EPD_EP_east)+ sin(nHar*EPD_EP_west)));
                runidvsEPD1cos->Fill(Run,0.5*(cos(EPD_EP_east)+ cos(EPD_EP_west)));
                runidvsEPD1sin->Fill(Run,0.5*(sin(EPD_EP_east)+ sin(EPD_EP_west)));

//cout<<"4th"<<endl;
//ZDC EP
	  	if(!IsGoodZDC(event)) continue;
          	MakeZDC_EP(event);
//flatten EPs
 	  	ShiftPsi();
	  	FillEP_resolution();
                EPD_hits(mEpdHits);

//Store the flattened phi for POI
	  	Phi_new.resize(NPTracks);
          	for(int trki = 0; trki < NPTracks; trki++){
			StPicoTrack *picoTrack = dst->track(trki);
			if(!IsGoodTrack(picoTrack)) continue;
			Pt        = picoTrack->pMom().Pt();
			Eta       = picoTrack->pMom().Eta();
			Charge    = picoTrack->charge();
			Phi       = picoTrack->pMom().Phi();
			DCAGlobal = picoTrack->gDCA(pV).Mag();

			if(!IsGoodPOI(picoTrack)) continue;
			if(!IsGoodPion(dst,picoTrack,1)) continue;
			FillPhiPOI();	//Charge is needed here
			ShiftPhiPOI(trki);
	  	}
//cout<<"6th"<<endl;
//////////Real analysis begins here//////////////////////////////////////////////////////////////////////////////////////

	} // Event

	WriteHistogram(cen,opt_weight);
	if(opt_weight==1) {
		mEpFinder->Finish();
		WriteWeight(fname_new);
	}
	return;
}
//////////////////////////////////
bool IsGoodEvent(int c) {
        hZDCcoin->Fill(ZDCcoin);
	hVertexXY->Fill(pVx, pVy);
	if(sqrt((pVx)*(pVx)+(pVy+2)*(pVy+2)) > Vr_cut) return false;
hTally->Fill(3);
        hVzDiff->Fill(pVz-VPDvz);
        //if(TMath::Abs(pVz-VPDvz)> Vz_diff) return false;
        hMult_Vz->Fill(RefMult,pVz);
        hVertexZ->Fill(pVz);
hTally->Fill(4);
        if(TMath::Abs(pVz-200) > Vz_cut) return false;        
//	if(pVz<0) return false;
        hMult_Vz_new->Fill(RefMult,pVz,Eweight);
	Ref_TOF->Fill(RefMult,TOFMult);
        Ref_Day3->Fill(Day3,RefMult);
        TOF_Day3->Fill(Day3,TOFMult);
        NPT_Day3->Fill(Day3,NPTracks);
//YU's  parameter
	runidvsrefmult->Fill(Run,RefMult);
	runidvstofmult->Fill(Run,TOFMult);
	runidvsnptrack->Fill(Run,NPTracks);
	runidvstofmatched->Fill(Run,Ntofmatch);
	runidvsavgvz->Fill(Run,pVz);
hTally->Fill(5);
        int Bad =0;
        for(int jj=0;jj<Nrun_MB;jj++) if(Day3 == bad_Ref_day3[jj]) {Bad = 1;break;}
        if(Bad) return false;  //bad run
hTally->Fill(6);	

        hCentrality->Fill(Centrality, Eweight);
//        if(c && Centrality != c) return false;
//	if(Centrality<1 && Centrality>9) return false;
hTally->Fill(7);

	return true;
}
//////////////////////////////////////////////
bool CountCharge(StPicoDst *d) {
          int Ntof = 0, Npos = 0, Nneg = 0;
          Fcount = 0, POIcount = 0;
          for(int trk = 0; trk < NPTracks; trk++) {
                StPicoTrack *picoTrack = d->track(trk);
                if(!IsGoodTrack(picoTrack)) continue;
                EtaAsso   = picoTrack->pMom().Eta();
                PtAsso    = picoTrack->pMom().Pt();
                ChargeAsso= picoTrack->charge();
                DCAGlobalAsso = picoTrack->gDCA(pV).Mag();
                nSigma_p  = picoTrack->nSigmaProton();

                if(ChargeAsso >0 && fabs(EtaAsso)<1.5 && PtAsso>0.15 && DCAGlobalAsso<1 && !(fabs(nSigma_p)<3 && PtAsso<0.4)) Npos++;
                if(ChargeAsso <0 && fabs(EtaAsso)<1.5 && PtAsso>0.15 && DCAGlobalAsso<1 && !(fabs(nSigma_p)<3 && PtAsso<0.4)) Nneg++;
                if(picoTrack->isTofTrack()) Ntof++;
                Hist_DCA->Fill(DCAGlobalAsso);
		
                if(IsGoodAsso(picoTrack)) Fcount++;
		if(IsGoodPOI(picoTrack))  POIcount++;
          }
          if(Ntof<2) return false;   //at least 2 tracks match TOF
          NPA_Day3->Fill(Day3,Fcount);
//Yu create NPA
	  runidvsnpasso->Fill(Run,Fcount);

hTally->Fill(9);
          int net_Nch= Npos - Nneg;
          net_Nch_Asym =-99;
          if((Npos + Nneg)> 0) net_Nch_Asym = (Npos - Nneg)/float(Npos + Nneg);
          Hist_positive->Fill(Npos);
          Hist_negative->Fill(Nneg);
          Hist_Ch->Fill(Npos+Nneg);
          Hist_netCh->Fill(net_Nch);
          if (net_Nch_Asym >-99) {
                Hist_netChAsym->Fill(net_Nch_Asym);
                p_netChAsym_RefMult->Fill(net_Nch_Asym, RefMult);
                if(net_Nch_Asym < (MeanNetChargeAsym - RMSNetChargeAsym)) Hist_netChAsym_bin->Fill(1,net_Nch_Asym);
                else if(net_Nch_Asym < (MeanNetChargeAsym - (0.3*RMSNetChargeAsym))) Hist_netChAsym_bin->Fill(2,net_Nch_Asym);
                else if(net_Nch_Asym < (MeanNetChargeAsym + (0.3*RMSNetChargeAsym))) Hist_netChAsym_bin->Fill(3,net_Nch_Asym);
                else if(net_Nch_Asym < (MeanNetChargeAsym + RMSNetChargeAsym)) Hist_netChAsym_bin->Fill(4,net_Nch_Asym);
                else Hist_netChAsym_bin->Fill(5,net_Nch_Asym);

          }
	  return true;
}
/////////////////////////////////////////////
void MakeTPC_EP(StPicoDst *d, int *iTr) {
          iCharge.resize(POIcount);
          match.resize(NPTracks);
          PhiAsso_new.resize(NPTracks);

          TVector2 mQ, mQ1, mQ2, mQ3, mQ4;
          mQx=0., mQy=0.;
          float mQQx_1=0., mQQy_1=0.,mQQx=0., mQQy=0., mQx1=0., mQy1=0., mQx2=0., mQy2=0., mQx3=0., mQy3=0., mQx4=0., mQy4=0.;
          Fcount = 0, POIcount = 0;
          int Qcount = 0;
          
          for(int trk = 0; trk < NPTracks; trk++) {
                StPicoTrack *picoTrack = d->track(trk);
                if(!IsGoodTrack(picoTrack)) continue;
                EtaAsso   = picoTrack->pMom().Eta();
                PtAsso    = picoTrack->pMom().Pt();
                PhiAsso   = picoTrack->pMom().Phi();
                DCAGlobalAsso = picoTrack->gDCA(pV).Mag();
                ChargeAsso= picoTrack->charge();
		gDCAxy = picoTrack->gDCAxy(pVx,pVy);
		dEdx = picoTrack->dEdx();
		//if(dEdx>100 || dEdx<0) cout<<dEdx<<endl;
		//to shuffle charge of POI
		match[trk] = 0;		
		if(IsGoodPOI(picoTrack)) {
			match[trk] = POIcount;
			iCharge[POIcount] = ChargeAsso;
			POIcount++;
		}
                if(!IsGoodAsso(picoTrack)) continue;

                pTemp_track->Fill(1,PtAsso);
                pTemp_track->Fill(2,EtaAsso);
                pTemp_track->Fill(3,PhiAsso);
                pTemp_track->Fill(4,ChargeAsso);
                pTemp_track->Fill(5,DCAGlobalAsso);
		pTemp_track->Fill(6,gDCAxy);
                pTemp_track_dedx->Fill(1,dEdx);
                FillPhiAsso();  //ChargeAsso is needed here
                ShiftPhiAsso(trk);
		//this is q2
		if((opt_useEPD==0 && fabs(EtaAsso)<Eta_ESE_Cut) || opt_useEPD>0) {
			if(IsGoodPOI(picoTrack) && IsGoodPion(d,picoTrack,1)) {
                                mQQx_1 += cos(PhiAsso_new[trk]*nHar);
                                mQQy_1 += sin(PhiAsso_new[trk]*nHar);
                                mQQx += cos(PhiAsso_new[trk]*nHar);
                                mQQy += sin(PhiAsso_new[trk]*nHar);
                                Qcount++;
			}
		}
		mQx +=PtAsso*cos(PhiAsso_new[trk]*nHar); mQy +=PtAsso*sin(PhiAsso_new[trk] * nHar);
                if(iTr[Fcount] > Scount) {mQx1 +=PtAsso*cos(PhiAsso_new[trk]*nHar); mQy1 +=PtAsso*sin(PhiAsso_new[trk] * nHar);}
                else {mQx2 += PtAsso * cos(PhiAsso_new[trk] * nHar); mQy2 += PtAsso * sin(PhiAsso_new[trk] * nHar);}
                if(EtaAsso> Eta_EP_Cut) {mQx3 +=PtAsso*cos(PhiAsso_new[trk]*nHar); mQy3 +=PtAsso*sin(PhiAsso_new[trk] * nHar);}
                if(EtaAsso<-Eta_EP_Cut) {mQx4 +=PtAsso*cos(PhiAsso_new[trk]*nHar); mQy4 +=PtAsso*sin(PhiAsso_new[trk] * nHar);}
                Fcount++;
          }
          mQ.Set(mQx, mQy); mQ1.Set(mQx1, mQy1); mQ2.Set(mQx2, mQy2); mQ3.Set(mQx3, mQy3); mQ4.Set(mQx4, mQy4);
          TPC_EP_full = mQ.Phi()/nHar;
          TPC_EP_east = mQ1.Phi()/nHar;
          TPC_EP_west = mQ2.Phi()/nHar;
          TPC_EP_for  = mQ3.Phi()/nHar;
          TPC_EP_bac  = mQ4.Phi()/nHar;
	  if(mQx3==0 || mQy3==0) TPC_EP_for = -10;
          if(mQx4==0 || mQy4==0) TPC_EP_bac = -10;
          Q2 = (mQQx*mQQx+mQQy*mQQy)/float(Qcount);
          TPC_Day3_cos2->Fill(Day3, cos(nHar*TPC_EP_full));
          TPC_Day3_sin2->Fill(Day3, sin(nHar*TPC_EP_full));
          pT_Day3->Fill(Day3,pTemp_track->GetBinContent(1));
          Eta_Day3->Fill(Day3,pTemp_track->GetBinContent(2));
          Phi_Day3->Fill(Day3,pTemp_track->GetBinContent(3));
          Charge_Day3->Fill(Day3,pTemp_track->GetBinContent(4));
          DCA_Day3->Fill(Day3,pTemp_track->GetBinContent(5));
//add Yu's
        runidvsdedx->Fill(Run,pTemp_track_dedx->GetBinContent(1)); 
	runidvsavgpt->Fill(Run,pTemp_track->GetBinContent(1));
	runidvsavgeta->Fill(Run,pTemp_track->GetBinContent(2));
	runidvsavgdca->Fill(Run,pTemp_track->GetBinContent(5));
	runidvsavgphi->Fill(Run,pTemp_track->GetBinContent(3));
	runidvsgdcaxy->Fill(Run,pTemp_track->GetBinContent(6));
	runidvsdcaxysigma->Fill(Run,pTemp_track->GetBinError(6)*sqrt(pTemp_track->GetBinEntries(6)));
	//float asigma= pTemp_track->GetBinError(6)*sqrt(pTemp_track->GetBinEntries(6));
	//if(asigma>10) cout<<asigma <<endl;//pTemp_track->GetBinError(6)*sqrt(pTemp_track->GetBinEntries(6));
	runidvsavgQ2x->Fill(Run,pow(mQQx,2)/float(Qcount));
	runidvsavgQ2y->Fill(Run,pow(mQQy,2)/float(Qcount));
	runidvsavgQ1x->Fill(Run,pow(mQQx_1,2)/float(Qcount));
	runidvsavgQ1y->Fill(Run,pow(mQQy_1,2)/float(Qcount));
	runidvsTPCcos->Fill(Run,cos(nHar*TPC_EP_full));
        runidvsTPCsin->Fill(Run,sin(nHar*TPC_EP_full));
          pTemp_track->Reset();
}
///////////////////////////////////////////
bool IsGoodBBC(StPicoEvent *e) {
        float BBC_sum_east = 0, BBC_sum_west = 0;
        for(int j=0;j<16;j++) BBC_sum_east += e->bbcAdcEast(j);
        for(int j=0;j<16;j++) BBC_sum_west += e->bbcAdcWest(j);
        if(BBC_sum_east<150 || BBC_sum_west<150) return false;
        int saturate = 0;
        for(int j=0;j<16;j++) if(e->bbcAdcEast(j)>4500 || e->bbcAdcWest(j)>4500) {saturate = 1;break;}
        if(saturate==1) return false;
        hTally->Fill(10);

        return true;
}
//////////////////////////////////////////////
void MakeBBC_EP(StPicoEvent *ev) {
          for(int j=0;j<6;j++) {
                BBC1->Fill(j+1, ev->bbcAdcEast(j)/1000.);
                BBC8->Fill(j+1, ev->bbcAdcWest(j)/1000.);
          }
          for(int j=0;j<10;j++) {
                BBC2->Fill(j+1, ev->bbcAdcEast(j+6)/1000.);
                BBC7->Fill(j+1, ev->bbcAdcWest(j+6)/1000.);
          }
	  TVector2 mBe, mBw;
	  float mBe_x = 0, mBe_y = 0, mBw_x = 0, mBw_y = 0; 
          for(int trk = 0; trk < 16; trk++) {
                float sig_bbc = ev->bbcAdcEast(trk);
                float phi_bbc = GetPhiInBBC(0, trk+1);
                mBe_x += cos(nHar*phi_bbc)*sig_bbc*BBC_gain_east[trk];
                mBe_y += sin(nHar*phi_bbc)*sig_bbc*BBC_gain_east[trk];
          }
          for(int trk = 0; trk < 16; trk++) {
                float sig_bbc = ev->bbcAdcWest(trk);
                float phi_bbc = GetPhiInBBC(1, trk+1);
                mBw_x += cos(nHar*phi_bbc)*sig_bbc*BBC_gain_west[trk];
                mBw_y += sin(nHar*phi_bbc)*sig_bbc*BBC_gain_west[trk];
          }
          if(mBe_x==0 || mBe_y==0 || mBw_x==0 || mBw_y==0) return;
          mBe.Set(mBe_x, mBe_y); mBw.Set(mBw_x, mBw_y);
          BBC_EP_east = mBe.Phi()/nHar;
          BBC_EP_west = mBw.Phi()/nHar;
          BBCe_Day3_cos2->Fill(Day3, cos(nHar*BBC_EP_east));
          BBCe_Day3_sin2->Fill(Day3, sin(nHar*BBC_EP_east));
          BBCw_Day3_cos2->Fill(Day3, cos(nHar*BBC_EP_west));
          BBCw_Day3_sin2->Fill(Day3, sin(nHar*BBC_EP_west));
}
//////////////////////////////////////
bool IsGoodZDC(StPicoEvent *e) {

        return true;
}
/////////////////////////////////////////////
void MakeZDC_EP(StPicoEvent *ev) {

          for(int j=0;j<8;j++) {
                ZDC_e_h->Fill(j+1, ev->ZdcSmdEastHorizontal(j)/1000.);
                ZDC_e_v->Fill(j+1, ev->ZdcSmdEastVertical(j)/1000.);
                ZDC_w_h->Fill(j+1, ev->ZdcSmdWestHorizontal(j)/1000.);
                ZDC_w_v->Fill(j+1, ev->ZdcSmdWestVertical(j)/1000.);
          }

          float eXsum=0.,eYsum=0.,eXWgt=0.,eYWgt=0.;
          float wXsum=0.,wYsum=0.,wXWgt=0.,wYWgt=0.;
          for(int v=0;v<=6;v++) {
                eXsum += GetXYInZDC(0,0,v+1,1)*ZDC_gain_east_v[v]*ev->ZdcSmdEastVertical(v);
                eXWgt += ZDC_gain_east_v[v]*ev->ZdcSmdEastVertical(v);
                wXsum += GetXYInZDC(1,0,v+1,1)*ZDC_gain_west_v[v]*ev->ZdcSmdWestVertical(v);
                wXWgt += ZDC_gain_west_v[v]*ev->ZdcSmdWestVertical(v);
          }
          for(int v=0;v<=7;v++) {
                eYsum += GetXYInZDC(0,1,v+1,1)*ZDC_gain_east_h[v]*ev->ZdcSmdEastHorizontal(v);
                eYWgt += ZDC_gain_east_h[v]*ev->ZdcSmdEastHorizontal(v);
                wYsum += GetXYInZDC(1,1,v+1,1)*ZDC_gain_west_h[v]*ev->ZdcSmdWestHorizontal(v);
                wYWgt += ZDC_gain_west_h[v]*ev->ZdcSmdWestHorizontal(v);
          }
          if(eXWgt) pZDCcenter->Fill(1,Day2, eXsum/eXWgt);
          if(eYWgt) pZDCcenter->Fill(2,Day2, eYsum/eYWgt);
          if(wXWgt) pZDCcenter->Fill(3,Day2, wXsum/wXWgt);
          if(wYWgt) pZDCcenter->Fill(4,Day2, wYsum/wYWgt);

          if(Weight_Read && Read_ZDCcenter->ProjectionY("center1",1,1)->GetBinContent(Day2-run_sta/10+1)>0) {
                zdcsmd_x0_e = Read_ZDCcenter->ProjectionY("center1",1,1)->GetBinContent(Day2-run_sta/10+1);
                zdcsmd_y0_e = Read_ZDCcenter->ProjectionY("center2",2,2)->GetBinContent(Day2-run_sta/10+1);
                zdcsmd_x0_w = Read_ZDCcenter->ProjectionY("center3",3,3)->GetBinContent(Day2-run_sta/10+1);
                zdcsmd_y0_w = Read_ZDCcenter->ProjectionY("center4",4,4)->GetBinContent(Day2-run_sta/10+1);
          }
          eXsum=0.,eYsum=0.,eXWgt=0.,eYWgt=0.;
          wXsum=0.,wYsum=0.,wXWgt=0.,wYWgt=0.;
          for(int v=0;v<=6;v++) {
                eXsum += GetXYInZDC(0,0,v+1,0)*ZDC_gain_east_v[v]*ev->ZdcSmdEastVertical(v);
                eXWgt += ZDC_gain_east_v[v]*ev->ZdcSmdEastVertical(v);
                wXsum += GetXYInZDC(1,0,v+1,0)*ZDC_gain_west_v[v]*ev->ZdcSmdWestVertical(v);
                wXWgt += ZDC_gain_west_v[v]*ev->ZdcSmdWestVertical(v);
          }
          for(int v=0;v<=7;v++) {
                eYsum += GetXYInZDC(0,1,v+1,0)*ZDC_gain_east_h[v]*ev->ZdcSmdEastHorizontal(v);
                eYWgt += ZDC_gain_east_h[v]*ev->ZdcSmdEastHorizontal(v);
                wYsum += GetXYInZDC(1,1,v+1,0)*ZDC_gain_west_h[v]*ev->ZdcSmdWestHorizontal(v);
                wYWgt += ZDC_gain_west_h[v]*ev->ZdcSmdWestHorizontal(v);
          }
	  //these 4 profiles should all give zero
          if(eXWgt) pZDCcenter_new->Fill(1,Day2, eXsum/eXWgt);
          if(eYWgt) pZDCcenter_new->Fill(2,Day2, eYsum/eYWgt);
          if(wXWgt) pZDCcenter_new->Fill(3,Day2, wXsum/wXWgt);
          if(wYWgt) pZDCcenter_new->Fill(4,Day2, wYsum/wYWgt);
          ZDC_EP_east = atan2(eYsum/eYWgt, eXsum/eXWgt);
          ZDC_EP_west = atan2(wYsum/wYWgt, wXsum/wXWgt);
          ZDCe_Day3_cos->Fill(Day3, cos(ZDC_EP_east));
          ZDCe_Day3_sin->Fill(Day3, sin(ZDC_EP_east));
          ZDCw_Day3_cos->Fill(Day3, cos(ZDC_EP_west));
          ZDCw_Day3_sin->Fill(Day3, sin(ZDC_EP_west));
}
/////////////////////////////////////////////////
void FillEP_resolution() {
          float cos_ew = cos(nHar*TPC_EP_east_new-nHar*TPC_EP_west_new);
          float cos_fb = cos(nHar*TPC_EP_for_new -nHar*TPC_EP_bac_new);
	  
          Hist_cos->Fill(1,cos_fb, Eweight);
          Hist_cos->Fill(2,cos_ew, Eweight);
          Hist_cos->Fill(3,cos(nHar*TPC_EP_for_new -nHar*TPC_EP_full_new));
          Hist_cos->Fill(4,cos(nHar*TPC_EP_bac_new -nHar*TPC_EP_full_new));
          Hist_cos_BBC->Fill(1,cos(nHar*BBC_EP_east_new-nHar*BBC_EP_west_new), Eweight);
          Hist_cos_BBC->Fill(2,cos(nHar*BBC_EP_east_new-nHar*TPC_EP_full_new), Eweight);
          Hist_cos_BBC->Fill(3,cos(nHar*BBC_EP_west_new-nHar*TPC_EP_full_new), Eweight);
          Hist_cos_EPD->Fill(1,cos(nHar*EPD_EP_east_new-nHar*EPD_EP_west_new), Eweight);
          Hist_cos_EPD->Fill(2,cos(nHar*EPD_EP_east_new-nHar*TPC_EP_full_new), Eweight);
          Hist_cos_EPD->Fill(3,cos(nHar*EPD_EP_west_new-nHar*TPC_EP_full_new), Eweight);
          Hist_cos_EPD->Fill(4,cos(EPD_EP1_east_new-EPD_EP1_west_new), Eweight);
          Hist_cos_EPD->Fill(5,cos(nHar*TPC_EP_full_new-EPD_EP1_east_new-EPD_EP1_west_new), Eweight);
          Hist_cos_EPD->Fill(6,cos(nHar*EPD_EP_east_new-EPD_EP1_east_new-EPD_EP1_west_new), Eweight);
          Hist_cos_EPD->Fill(7,cos(nHar*EPD_EP_west_new-EPD_EP1_east_new-EPD_EP1_west_new), Eweight);
          Hist_cos_ZDC->Fill(1,cos(ZDC_EP_east_new-ZDC_EP_west_new), Eweight);
	  BBC_cor_Day3->Fill(Day3, cos(nHar*BBC_EP_east_new-nHar*BBC_EP_west_new));
	  EPD_cor_Day3->Fill(Day3, cos(nHar*EPD_EP_east_new-nHar*EPD_EP_west_new));
          EPD1_cor_Day3->Fill(Day3, cos(EPD_EP1_east_new-EPD_EP1_west_new));
	  ZDC_cor_Day3->Fill(Day3, cos(ZDC_EP_east_new-ZDC_EP_west_new));
}
/////////////////////////////////////////////////
bool IsGoodAsso(StPicoTrack *p) {
	if(p->nHitsFit() < nFitHits_asso) return false;
	if(p->pMom().Pt() > pt_asso_up || p->pMom().Pt() < pt_asso_lo) return false;
        if(p->gDCA(pV).Mag() > DcaCut_asso) return false;
        if(fabs(p->pMom().Eta())>=EtaCut_asso) return false;
	return true;
}
/////////////////////////////////////////////////////
bool IsGoodPOI(StPicoTrack *p) {
	int nHF = p->nHitsFit();
        if(nHF < nFitHits_trig) return false;
        int nHM = p->nHitsMax();
        //if(nHF/float(nHM) < nFitRatio_trig_lo || nHF/float(nHM) >= nFitRatio_trig_up) return false;
        if(p->pMom().Pt() > pt_trig_up || p->pMom().Pt() < pt_trig_lo) return false;
        if(p->gDCA(pV).Mag() > DcaCut_trig) return false;
        if(fabs(p->pMom().Eta())>=EtaCut_trig) return false;
        return true;
}
/////////////////////////////////////////////////
bool IsGoodTrack(StPicoTrack *p) {
                if(!p->isPrimary()) return false;
		return true;
}
/////////////////////////////////////////////////
bool IsGoodPionNoTOF(StPicoTrack *p) {
	if(p->gDCA(pV).Mag()>1) return false;
	float nSig_i = p->nSigmaPion();
        float ndEdx_i = p->nHitsDedx();
        if(ndEdx_i<10 || nSig_i>2 || nSig_i<-2) return false;
	return true;
}
//////////////////////////////////////////////////
bool IsGoodPion(StPicoDst *d, StPicoTrack *p, int opt) {
//hTall->Fill(1);
	if(p->gDCA(pV).Mag()>1)	return false;
//hTall->Fill(2);
	float eta_i = p->pMom().Eta();
	if(fabs(eta_i)>0.9) return false;
//hTall->Fill(3);
        float pt_i =  p->pMom().Pt();
	float p_i = pt_i*cosh(eta_i);
	if(pt_i<0.2 || p_i>1.6) return false;
//hTall->Fill(4);
	float nSig_i = p->nSigmaPion();
	float ndEdx_i = p->nHitsDedx();
	if(ndEdx_i<15 || nSig_i>2 || nSig_i<-2) return false;
//hTall->Fill(5);
	if(opt==0) return true;
	if(!(p->isTofTrack())) return false;
//hTall->Fill(6);
	StPicoBTofPidTraits *trait = d->btofPidTraits( p->bTofPidTraitsIndex() );
	if(!trait) return false;
//hTall->Fill(7);
        if(trait->btof()<=0) return false;
//hTall->Fill(8);
	if(fabs(trait->btofYLocal())>1.8) return false;
//hTall->Fill(9);
	float beta_i = trait->btofBeta();
	if(beta_i==0) return false;
//hTall->Fill(10);
	float mass2_i = p_i*p_i*(1.0/beta_i/beta_i - 1.0);
	if(mass2_i<-0.01 || mass2_i>0.1) return false;
//hTall->Fill(11);
	return true;
}
///////////////////////////////////////////////////
bool IsGoodKaon(StPicoDst *d, StPicoTrack *p, int opt) {
	if(p->gDCA(pV).Mag()>1) return false;
        float eta_i = p->pMom().Eta();
        if(fabs(eta_i)>0.9) return false;
        float pt_i =  p->pMom().Pt();
        float p_i = pt_i*cosh(eta_i);
        if(pt_i<0.2 || p_i>1.6) return false;

        float nSig_i = p->nSigmaKaon();
        float ndEdx_i = p->nHitsDedx();
        if(ndEdx_i<15 || nSig_i>2 || nSig_i<-2) return false;

        if(opt==0) return true;
        if(!(p->isTofTrack())) return false;
        StPicoBTofPidTraits *trait = d->btofPidTraits( p->bTofPidTraitsIndex() );
        if(!trait) return false;
        if(trait->btof()<=0) return false;
        if(fabs(trait->btofYLocal())>1.8) return false;
        float beta_i = trait->btofBeta();
        if(beta_i==0) return false;
        float mass2_i = p_i*p_i*(1.0/beta_i/beta_i - 1.0);
        if(mass2_i<0.2 || mass2_i>0.35) return false;
        return true;
}
////////////////////////////////////////////////////
bool IsGoodProton(StPicoDst *d, StPicoTrack *p, int opt) {
        if(p->gDCA(pV).Mag()>1) return false;
        float eta_i = p->pMom().Eta();
        if(fabs(eta_i)>0.9) return false;
        float pt_i =  p->pMom().Pt();
        float p_i = pt_i*cosh(eta_i);
        if(pt_i<0.4 || p_i>2) return false;

        float nSig_i = p->nSigmaProton();
        float ndEdx_i = p->nHitsDedx();
        if(ndEdx_i<15 || nSig_i>2 || nSig_i<-2) return false;
        
        if(opt==0) return true;
        if(!(p->isTofTrack())) return false;
        StPicoBTofPidTraits *trait = d->btofPidTraits( p->bTofPidTraitsIndex() );
        if(!trait) return false;
        if(trait->btof()<=0) return false;
        if(fabs(trait->btofYLocal())>1.8) return false;
        float beta_i = trait->btofBeta();
        if(beta_i==0) return false;
        float mass2_i = p_i*p_i*(1.0/beta_i/beta_i - 1.0);
        if(mass2_i<0.8 || mass2_i>1) return false;
        return true;
}
////////////////////////////////////
void ShiftPsi() {
          Hist_TPC_EP_full->Fill(TPC_EP_full,Day);
          Hist_TPC_EP_east->Fill(TPC_EP_east,Day);
          Hist_TPC_EP_west->Fill(TPC_EP_west,Day);
          Hist_TPC_EP_for->Fill(TPC_EP_for,Day);
          Hist_TPC_EP_bac->Fill(TPC_EP_bac,Day);
	  Hist_BBC_EP_east->Fill(BBC_EP_east,Day);
          Hist_BBC_EP_west->Fill(BBC_EP_west,Day);
          Hist_EPD_EP1_east->Fill(EPD_EP1_east,Day);
          Hist_EPD_EP1_west->Fill(EPD_EP1_west,Day);
          Hist_EPD_EP_east->Fill(EPD_EP_east,Day);
          Hist_EPD_EP_west->Fill(EPD_EP_west,Day);
	  Hist_ZDC_EP_east->Fill(ZDC_EP_east,Day);
          Hist_ZDC_EP_west->Fill(ZDC_EP_west,Day);

          for(int kk=0;kk<order;kk++) {
                pTPC_EP_full->Fill(1+2*kk,Day2,cos(nHar*(kk+1)*TPC_EP_full),Eweight);
                pTPC_EP_full->Fill(2+2*kk,Day2,sin(nHar*(kk+1)*TPC_EP_full),Eweight);
                pTPC_EP_east->Fill(1+2*kk,Day2,cos(nHar*(kk+1)*TPC_EP_east),Eweight);
                pTPC_EP_east->Fill(2+2*kk,Day2,sin(nHar*(kk+1)*TPC_EP_east),Eweight);
                pTPC_EP_west->Fill(1+2*kk,Day2,cos(nHar*(kk+1)*TPC_EP_west),Eweight);
                pTPC_EP_west->Fill(2+2*kk,Day2,sin(nHar*(kk+1)*TPC_EP_west),Eweight);
                pTPC_EP_for->Fill(1+2*kk, Day2,cos(nHar*(kk+1)*TPC_EP_for),Eweight);
                pTPC_EP_for->Fill(2+2*kk, Day2,sin(nHar*(kk+1)*TPC_EP_for),Eweight);
                pTPC_EP_bac->Fill(1+2*kk, Day2,cos(nHar*(kk+1)*TPC_EP_bac),Eweight);
                pTPC_EP_bac->Fill(2+2*kk, Day2,sin(nHar*(kk+1)*TPC_EP_bac),Eweight);
		pBBC_EP_east->Fill(1+2*kk,Day2,cos(nHar*(kk+1)*BBC_EP_east),Eweight);
                pBBC_EP_east->Fill(2+2*kk,Day2,sin(nHar*(kk+1)*BBC_EP_east),Eweight);
                pBBC_EP_west->Fill(1+2*kk,Day2,cos(nHar*(kk+1)*BBC_EP_west),Eweight);
                pBBC_EP_west->Fill(2+2*kk,Day2,sin(nHar*(kk+1)*BBC_EP_west),Eweight);
                pEPD_EP1_east->Fill(1+2*kk,Day2,cos((kk+1)*EPD_EP1_east),Eweight);
                pEPD_EP1_east->Fill(2+2*kk,Day2,sin((kk+1)*EPD_EP1_east),Eweight);
                pEPD_EP1_west->Fill(1+2*kk,Day2,cos((kk+1)*EPD_EP1_west),Eweight);
                pEPD_EP1_west->Fill(2+2*kk,Day2,sin((kk+1)*EPD_EP1_west),Eweight);
                pEPD_EP_east->Fill(1+2*kk,Day2,cos(nHar*(kk+1)*EPD_EP_east),Eweight);
                pEPD_EP_east->Fill(2+2*kk,Day2,sin(nHar*(kk+1)*EPD_EP_east),Eweight);
                pEPD_EP_west->Fill(1+2*kk,Day2,cos(nHar*(kk+1)*EPD_EP_west),Eweight);
                pEPD_EP_west->Fill(2+2*kk,Day2,sin(nHar*(kk+1)*EPD_EP_west),Eweight);
                pZDC_EP_east->Fill(1+2*kk,Day2,cos((kk+1)*ZDC_EP_east),Eweight);
                pZDC_EP_east->Fill(2+2*kk,Day2,sin((kk+1)*ZDC_EP_east),Eweight);
                pZDC_EP_west->Fill(1+2*kk,Day2,cos((kk+1)*ZDC_EP_west),Eweight);
                pZDC_EP_west->Fill(2+2*kk,Day2,sin((kk+1)*ZDC_EP_west),Eweight);
          }

          if(Weight_Read && Read_TPC_EP_full->GetEntries()) {
                for(int k=0;k<2*order;k++) {
                        PsiMean_F[k] = Read_TPC_EP_full->GetBinContent(k+1,Day2-run_sta/10+1);
                        PsiMean_E[k] = Read_TPC_EP_east->GetBinContent(k+1,Day2-run_sta/10+1);
                        PsiMean_W[k] = Read_TPC_EP_west->GetBinContent(k+1,Day2-run_sta/10+1);
                        PsiMean_f[k] = Read_TPC_EP_for->GetBinContent(k+1,Day2-run_sta/10+1);
                        PsiMean_b[k] = Read_TPC_EP_bac->GetBinContent(k+1,Day2-run_sta/10+1);
                        Psi_BBC_E[k] = Read_BBC_EP_east->GetBinContent(k+1,Day2-run_sta/10+1);
                        Psi_BBC_W[k] = Read_BBC_EP_west->GetBinContent(k+1,Day2-run_sta/10+1);
                        Psi1_EPD_E[k] = Read_EPD_EP1_east->GetBinContent(k+1,Day2-run_sta/10+1);
                        Psi1_EPD_W[k] = Read_EPD_EP1_west->GetBinContent(k+1,Day2-run_sta/10+1);
                        Psi_EPD_E[k] = Read_EPD_EP_east->GetBinContent(k+1,Day2-run_sta/10+1);
                        Psi_EPD_W[k] = Read_EPD_EP_west->GetBinContent(k+1,Day2-run_sta/10+1);
                        Psi_ZDC_E[k] = Read_ZDC_EP_east->GetBinContent(k+1,Day2-run_sta/10+1);
                        Psi_ZDC_W[k] = Read_ZDC_EP_west->GetBinContent(k+1,Day2-run_sta/10+1);
                }
          }

          TPC_EP_full_new = TPC_EP_full, TPC_EP_east_new = TPC_EP_east, TPC_EP_west_new = TPC_EP_west;
          TPC_EP_for_new = TPC_EP_for, TPC_EP_bac_new = TPC_EP_bac;
	  BBC_EP_east_new = BBC_EP_east, BBC_EP_west_new = BBC_EP_west;
 	  EPD_EP_east_new = EPD_EP_east, EPD_EP_west_new = EPD_EP_west;
          EPD_EP1_east_new = EPD_EP1_east, EPD_EP1_west_new = EPD_EP1_west;
          ZDC_EP_east_new = ZDC_EP_east, ZDC_EP_west_new = ZDC_EP_west; 
          for(int jj=0;jj<order;jj++) {
                TPC_EP_full_new += -2*PsiMean_F[1+2*jj]*cos(nHar*(jj+1)*TPC_EP_full)/nHar/(jj+1) + 2*PsiMean_F[0+2*jj]*sin(nHar*(jj+1)*TPC_EP_full)/nHar/(jj+1);
                TPC_EP_east_new += -2*PsiMean_E[1+2*jj]*cos(nHar*(jj+1)*TPC_EP_east)/nHar/(jj+1) + 2*PsiMean_E[0+2*jj]*sin(nHar*(jj+1)*TPC_EP_east)/nHar/(jj+1);
                TPC_EP_west_new += -2*PsiMean_W[1+2*jj]*cos(nHar*(jj+1)*TPC_EP_west)/nHar/(jj+1) + 2*PsiMean_W[0+2*jj]*sin(nHar*(jj+1)*TPC_EP_west)/nHar/(jj+1);
                TPC_EP_for_new  += -2*PsiMean_f[1+2*jj]*cos(nHar*(jj+1)*TPC_EP_for)/nHar/(jj+1)  + 2*PsiMean_f[0+2*jj]*sin(nHar*(jj+1)*TPC_EP_for)/nHar/(jj+1);
                TPC_EP_bac_new  += -2*PsiMean_b[1+2*jj]*cos(nHar*(jj+1)*TPC_EP_bac)/nHar/(jj+1)  + 2*PsiMean_b[0+2*jj]*sin(nHar*(jj+1)*TPC_EP_bac)/nHar/(jj+1);
		BBC_EP_east_new += -2*Psi_BBC_E[1+2*jj]*cos(nHar*(jj+1)*BBC_EP_east)/nHar/(jj+1) + 2*Psi_BBC_E[0+2*jj]*sin(nHar*(jj+1)*BBC_EP_east)/nHar/(jj+1);
		BBC_EP_west_new += -2*Psi_BBC_W[1+2*jj]*cos(nHar*(jj+1)*BBC_EP_west)/nHar/(jj+1) + 2*Psi_BBC_W[0+2*jj]*sin(nHar*(jj+1)*BBC_EP_west)/nHar/(jj+1);
		EPD_EP1_east_new += -2*Psi1_EPD_E[1+2*jj]*cos((jj+1)*EPD_EP1_east)/(jj+1) + 2*Psi1_EPD_E[0+2*jj]*sin((jj+1)*EPD_EP1_east)/(jj+1);
                EPD_EP1_west_new += -2*Psi1_EPD_W[1+2*jj]*cos((jj+1)*EPD_EP1_west)/(jj+1) + 2*Psi1_EPD_W[0+2*jj]*sin((jj+1)*EPD_EP1_west)/(jj+1);
		EPD_EP_east_new += -2*Psi_EPD_E[1+2*jj]*cos(nHar*(jj+1)*EPD_EP_east)/nHar/(jj+1) + 2*Psi_EPD_E[0+2*jj]*sin(nHar*(jj+1)*EPD_EP_east)/nHar/(jj+1);
                EPD_EP_west_new += -2*Psi_EPD_W[1+2*jj]*cos(nHar*(jj+1)*EPD_EP_west)/nHar/(jj+1) + 2*Psi_EPD_W[0+2*jj]*sin(nHar*(jj+1)*EPD_EP_west)/nHar/(jj+1);
                ZDC_EP_east_new += -2*Psi_ZDC_E[1+2*jj]*cos((jj+1)*ZDC_EP_east)/(jj+1) + 2*Psi_ZDC_E[0+2*jj]*sin((jj+1)*ZDC_EP_east)/(jj+1);
		ZDC_EP_west_new += -2*Psi_ZDC_W[1+2*jj]*cos((jj+1)*ZDC_EP_west)/(jj+1) + 2*Psi_ZDC_W[0+2*jj]*sin((jj+1)*ZDC_EP_west)/(jj+1);
          }
          if(TPC_EP_full_new>2.*PI/nHar) TPC_EP_full_new -= 2.*PI/nHar; if(TPC_EP_full_new< 0) TPC_EP_full_new += 2.*PI/nHar;
          if(TPC_EP_east_new>2.*PI/nHar) TPC_EP_east_new -= 2.*PI/nHar; if(TPC_EP_east_new< 0) TPC_EP_east_new += 2.*PI/nHar;
          if(TPC_EP_west_new>2.*PI/nHar) TPC_EP_west_new -= 2.*PI/nHar; if(TPC_EP_west_new< 0) TPC_EP_west_new += 2.*PI/nHar;
          if(TPC_EP_for_new>2.*PI/nHar)  TPC_EP_for_new  -= 2.*PI/nHar; if(TPC_EP_for_new< 0)  TPC_EP_for_new  += 2.*PI/nHar;
          if(TPC_EP_bac_new>2.*PI/nHar)  TPC_EP_bac_new  -= 2.*PI/nHar; if(TPC_EP_bac_new< 0)  TPC_EP_bac_new  += 2.*PI/nHar;
          if(BBC_EP_east_new>2.*PI/nHar) BBC_EP_east_new -= 2.*PI/nHar; if(BBC_EP_east_new< 0) BBC_EP_east_new += 2.*PI/nHar;
          if(BBC_EP_west_new>2.*PI/nHar) BBC_EP_west_new -= 2.*PI/nHar; if(BBC_EP_west_new< 0) BBC_EP_west_new += 2.*PI/nHar;
	  if(EPD_EP1_east_new>PI) EPD_EP1_east_new -= 2*PI; if(EPD_EP1_east_new<-PI) EPD_EP1_east_new += 2*PI;
          if(EPD_EP1_west_new>PI) EPD_EP1_west_new -= 2*PI; if(EPD_EP1_west_new<-PI) EPD_EP1_west_new += 2*PI;
          if(EPD_EP_east_new>2.*PI/nHar) EPD_EP_east_new -= 2.*PI/nHar; if(EPD_EP_east_new< 0) EPD_EP_east_new += 2.*PI/nHar;
          if(EPD_EP_west_new>2.*PI/nHar) EPD_EP_west_new -= 2.*PI/nHar; if(EPD_EP_west_new< 0) EPD_EP_west_new += 2.*PI/nHar;
	  if(ZDC_EP_east_new>PI) ZDC_EP_east_new -= 2*PI; if(ZDC_EP_east_new<-PI) ZDC_EP_east_new += 2*PI;
          if(ZDC_EP_west_new>PI) ZDC_EP_west_new -= 2*PI; if(ZDC_EP_west_new<-PI) ZDC_EP_west_new += 2*PI;
          Hist_TPC_EP_east_flat->Fill(TPC_EP_east_new,Day);
          Hist_TPC_EP_west_flat->Fill(TPC_EP_west_new,Day);
          Hist_TPC_EP_for_flat->Fill(TPC_EP_for_new,Day);
          Hist_TPC_EP_bac_flat->Fill(TPC_EP_bac_new,Day);
          Hist_TPC_EP_full_flat->Fill(TPC_EP_full_new,Day);
          Hist_BBC_EP_east_flat->Fill(BBC_EP_east_new,Day);
          Hist_BBC_EP_west_flat->Fill(BBC_EP_west_new,Day);
          Hist_EPD_EP1_east_flat->Fill(EPD_EP1_east_new,Day);
          Hist_EPD_EP1_west_flat->Fill(EPD_EP1_west_new,Day);
          Hist_EPD_EP_east_flat->Fill(EPD_EP_east_new,Day);
          Hist_EPD_EP_west_flat->Fill(EPD_EP_west_new,Day);
	  Hist_ZDC_EP_east_flat->Fill(ZDC_EP_east_new,Day);
          Hist_ZDC_EP_west_flat->Fill(ZDC_EP_west_new,Day);

	  float BBC_EP_full = atan2(sin(nHar*BBC_EP_east_new)+sin(nHar*BBC_EP_west_new),cos(nHar*BBC_EP_east_new)+cos(nHar*BBC_EP_west_new))/nHar;
	  if(BBC_EP_full< 0) BBC_EP_full += 2.*PI/nHar;
	  float EPD_EP_full = atan2(sin(nHar*EPD_EP_east_new)+sin(nHar*EPD_EP_west_new),cos(nHar*EPD_EP_east_new)+cos(nHar*EPD_EP_west_new))/nHar;
          if(EPD_EP_full< 0) EPD_EP_full += 2.*PI/nHar;
	  Hist_BBC_vs_TPC->Fill(TPC_EP_full_new,BBC_EP_full);
	  Hist_EPD_vs_TPC->Fill(TPC_EP_full_new,EPD_EP_full);
}       
//////////////////////////////////
void ShiftPhiPOI(int tr) {
        int index = 0, index2 = 0;
        if(pVz>0) index2 = (Charge > 0)? 1:2;
        else index2 = (Charge > 0)? 3:4;
        int Eta_index = int((Eta+1.5)*5);
        if(Weight_Read && (ReadWeightHist[Eta_index].pTPCmeanPhi_1->GetEntries())) {
                for(int kk=0;kk<order;kk++) {
                        if(pVz>0) index = (Charge > 0)? 1+8*kk  :3+8*kk;
                        else      index = (Charge > 0)? 1+4+8*kk:3+4+8*kk;

			if(Pt<0.5) {
				PhiMean_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_1->GetBinContent(index,  Day2-run_sta/10+1);
                                PhiMean_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_1->GetBinContent(index+1,Day2-run_sta/10+1);
			}
			else if(Pt<1) {
				PhiMean_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_2->GetBinContent(index,  Day2-run_sta/10+1);
                                PhiMean_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_2->GetBinContent(index+1,Day2-run_sta/10+1);
			}
			else {
                                PhiMean_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_3->GetBinContent(index,  Day2-run_sta/10+1);
                                PhiMean_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_3->GetBinContent(index+1,Day2-run_sta/10+1);
			}
                }
        }

        Phi_new[tr] = Phi;   //store the shifted angles
        for(int jj=0;jj<order;jj++) Phi_new[tr] += -2*PhiMean_sin[jj]*cos(jj*Phi+Phi)/(jj+1) + 2*PhiMean_cos[jj]*sin(jj*Phi+Phi)/(jj+1);
        if(Phi_new[tr]> PI) Phi_new[tr] -= 2*PI;
        if(Phi_new[tr]<-PI) Phi_new[tr] += 2*PI;

	if(Pt<0.5) PhiHist[Eta_index].Hist_Phi_new_1->Fill(Phi_new[tr],index2,Eweight);
	else if(Pt<1) PhiHist[Eta_index].Hist_Phi_new_2->Fill(Phi_new[tr],index2,Eweight);
	else PhiHist[Eta_index].Hist_Phi_new_3->Fill(Phi_new[tr],index2,Eweight);
}
//////////////////////////////////////////
void ShiftPhiAsso(int tr) {
        int index = 0;
        int Eta_index = int((EtaAsso+1.5)*5);//1
        if(Weight_Read && (ReadWeightHist[Eta_index].pTPCmeanPhiAsso_1->GetEntries())) {
                for(int kk=0;kk<order;kk++) {
                        if(pVz>0) index = (ChargeAsso > 0)? 1+8*kk  :3+8*kk;
                        else      index = (ChargeAsso > 0)? 1+4+8*kk:3+4+8*kk;

			if(PtAsso<0.5) {
				PhiMeanAsso_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_1->GetBinContent(index,  Day2-run_sta/10+1);
				PhiMeanAsso_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_1->GetBinContent(index+1,Day2-run_sta/10+1);
			}
			else if(PtAsso<1) {
				PhiMeanAsso_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_2->GetBinContent(index,  Day2-run_sta/10+1);
                                PhiMeanAsso_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_2->GetBinContent(index+1,Day2-run_sta/10+1);
			}
			else {
				PhiMeanAsso_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_3->GetBinContent(index,  Day2-run_sta/10+1);
                                PhiMeanAsso_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_3->GetBinContent(index+1,Day2-run_sta/10+1);
			}
                }
        }

        PhiAsso_new[tr] = PhiAsso;
        for(int jj=0;jj<order;jj++) PhiAsso_new[tr] += -2*PhiMeanAsso_sin[jj]*cos(jj*PhiAsso+PhiAsso)/(jj+1) + 2*PhiMeanAsso_cos[jj]*sin(jj*PhiAsso+PhiAsso)/(jj+1);
}
/////////////////////////////////////
void FillPhiAsso() {  // shift parameters for Particles of EP
	int index = 0;
	int Eta_index = int((EtaAsso+1.5)*5);
    	for(int kk=0;kk<order;kk++) {
		if(pVz>0) index = (ChargeAsso > 0)? 1+8*kk  :3+8*kk;
		else      index = (ChargeAsso > 0)? 1+4+8*kk:3+4+8*kk;
		
		if(PtAsso<0.5) {
			WeightHist[Eta_index].pTPCmeanPhiAsso_1->Fill(index,  Day2,cos(kk*PhiAsso+PhiAsso),Eweight);
			WeightHist[Eta_index].pTPCmeanPhiAsso_1->Fill(index+1,Day2,sin(kk*PhiAsso+PhiAsso),Eweight);
		}
		else if(PtAsso<1) {
			WeightHist[Eta_index].pTPCmeanPhiAsso_2->Fill(index,  Day2,cos(kk*PhiAsso+PhiAsso),Eweight);
                        WeightHist[Eta_index].pTPCmeanPhiAsso_2->Fill(index+1,Day2,sin(kk*PhiAsso+PhiAsso),Eweight);
		}
		else {
			WeightHist[Eta_index].pTPCmeanPhiAsso_3->Fill(index,  Day2,cos(kk*PhiAsso+PhiAsso),Eweight);
                        WeightHist[Eta_index].pTPCmeanPhiAsso_3->Fill(index+1,Day2,sin(kk*PhiAsso+PhiAsso),Eweight);
		}
    	}
}
///////////////////////////////////
void FillPhiPOI() { // shift parameters for Particles of Interest
        int index = 0, index2 = 0;
	if(pVz>0) index2 = (Charge > 0)? 1:2;
	else index2 = (Charge > 0)? 3:4;
	int Eta_index = int((Eta+1.5)*5);
        for(int kk=0;kk<order;kk++) {
                if(pVz>0) index = (Charge > 0)? 1+8*kk  :3+8*kk;
                else      index = (Charge > 0)? 1+4+8*kk:3+4+8*kk;
                if(Pt<0.5) {
                        WeightHist[Eta_index].pTPCmeanPhi_1->Fill(index,  Day2,cos(kk*Phi+Phi),Eweight);
                        WeightHist[Eta_index].pTPCmeanPhi_1->Fill(index+1,Day2,sin(kk*Phi+Phi),Eweight);
			if(kk==0) PhiHist[Eta_index].Hist_Phi_1->Fill(Phi,index2,Eweight);
                }
	        else if(Pt<1) {
			WeightHist[Eta_index].pTPCmeanPhi_2->Fill(index,  Day2,cos(kk*Phi+Phi),Eweight);
                        WeightHist[Eta_index].pTPCmeanPhi_2->Fill(index+1,Day2,sin(kk*Phi+Phi),Eweight);
                        if(kk==0) PhiHist[Eta_index].Hist_Phi_2->Fill(Phi,index2,Eweight);
		}
		else {
			WeightHist[Eta_index].pTPCmeanPhi_3->Fill(index,  Day2,cos(kk*Phi+Phi),Eweight);
                        WeightHist[Eta_index].pTPCmeanPhi_3->Fill(index+1,Day2,sin(kk*Phi+Phi),Eweight);
                        if(kk==0) PhiHist[Eta_index].Hist_Phi_3->Fill(Phi,index2,Eweight);
		}
	}
}
/////////////////////////////////
void WriteHistogram(int c, int o) {
	char fname_out[200];
	sprintf(fname_out,"cen%d.gamma1%d%d_pipi_EP%d_Boost%d.root",c,nHar-1,nHar,opt_useEPD,opt_boost);
	if(opt_HighOrder==1) sprintf(fname_out,"cen%d.gamma1%d%d_pipi_EP%d_Boost%d.root",c,nHar+1,nHar,opt_useEPD,opt_boost);
	if(opt_HighOrder==2) sprintf(fname_out,"cen%d.gamma422_pipi_EP%d_Boost%d.root",c,opt_useEPD,opt_boost);
	if(o==1) sprintf(fname_out,"cen%d.v%d.root",c,nHar);
        TFile *fout = new TFile(fname_out,"RECREATE");

//        gROOT->GetList()->ls();
        TList* list = gROOT->GetList();//GetListOfKeys();
        TIter next(list);
        TKey* key;
        TObject* obj;
        while ((key = (TKey*)next())) {
		TString tempStr(key->GetName());
		if (tempStr.Contains("Temp")) continue;
		if(o==1 && (tempStr.Contains("Parity")||tempStr.Contains("Delta"))) continue;
        	obj = gROOT->Get(key->GetName());
                if (!obj) continue;
                if(obj->IsA() == TDirectory::Class()){
                        delete obj;
                        obj = NULL;
                        continue;
        	}
        	obj->Write();
   	}
	fout->Write();
	fout->Close();
}
//////////////////////////////////////////////////
void WriteWeight(char* OutFileName) {
        TFile *fWgtNew = new TFile(OutFileName,"UPDATE");
	for(int i=0;i<15;i++) {
		WeightHist[i].pTPCmeanPhi_1->Write();
                WeightHist[i].pTPCmeanPhi_2->Write();
                WeightHist[i].pTPCmeanPhi_3->Write();
                WeightHist[i].pTPCmeanPhiAsso_1->Write();
                WeightHist[i].pTPCmeanPhiAsso_2->Write();
                WeightHist[i].pTPCmeanPhiAsso_3->Write();
	}
	BBC1->Write();
        BBC2->Write();
        BBC7->Write();
        BBC8->Write();
	ZDC_e_h->Write();
        ZDC_e_v->Write();
        ZDC_w_h->Write();
        ZDC_w_v->Write();
	pZDCcenter->Write();
        Hist_netChAsym->Write();
        pTPC_EP_east->Write();
        pTPC_EP_west->Write();
        pTPC_EP_for->Write();
        pTPC_EP_bac->Write();
        pTPC_EP_full->Write();
        pBBC_EP_east->Write();
        pBBC_EP_west->Write();
        pEPD_EP1_east->Write();
        pEPD_EP1_west->Write();
        pEPD_EP_east->Write();
        pEPD_EP_west->Write();
        pZDC_EP_east->Write();
        pZDC_EP_west->Write();
//	rc->Write();
        fWgtNew->Close();
}
/////////////////////////////////////////////////////
int ReadWeight(char* InFileName) {
        TFile *fWgt=new TFile(InFileName,"READ");
        if(!fWgt->IsOpen()) return 0;
        if(fWgt->IsOpen()) {
		TString* histTitle;
        	for(int i=0;i<15;i++) {
                	histTitle = new TString("TPCmeanPhi_1_");
                	*histTitle += i+1;
                	ReadWeightHist[i].pTPCmeanPhi_1 = (TProfile2D*)fWgt->Get(histTitle->Data());
                	delete histTitle;

			histTitle = new TString("TPCmeanPhi_2_");
                        *histTitle += i+1;
                        ReadWeightHist[i].pTPCmeanPhi_2 = (TProfile2D*)fWgt->Get(histTitle->Data());
                        delete histTitle;

                        histTitle = new TString("TPCmeanPhi_3_");
                        *histTitle += i+1;
                        ReadWeightHist[i].pTPCmeanPhi_3 = (TProfile2D*)fWgt->Get(histTitle->Data());
                        delete histTitle;
		
                        histTitle = new TString("TPCmeanPhiAsso_1_");
                        *histTitle += i+1;
                        ReadWeightHist[i].pTPCmeanPhiAsso_1 = (TProfile2D*)fWgt->Get(histTitle->Data());
                        delete histTitle;

                        histTitle = new TString("TPCmeanPhiAsso_2_");
                        *histTitle += i+1;
                        ReadWeightHist[i].pTPCmeanPhiAsso_2 = (TProfile2D*)fWgt->Get(histTitle->Data());
                        delete histTitle;

                        histTitle = new TString("TPCmeanPhiAsso_3_");
                        *histTitle += i+1;
                        ReadWeightHist[i].pTPCmeanPhiAsso_3 = (TProfile2D*)fWgt->Get(histTitle->Data());
                        delete histTitle;
		}

                TOF_eff = (TH1D*)fWgt->Get("rc");
                if(TOF_eff && TOF_eff->GetEntries()) {float cont = TOF_eff->GetBinContent(20); TOF_eff->Scale(1.25/cont);}
                Read_TPC_EP_full = (TProfile2D*)fWgt->Get("pTPC_EP_full");
                Read_TPC_EP_east = (TProfile2D*)fWgt->Get("pTPC_EP_east");
                Read_TPC_EP_west = (TProfile2D*)fWgt->Get("pTPC_EP_west");
                Read_TPC_EP_for = (TProfile2D*)fWgt->Get("pTPC_EP_for");
                Read_TPC_EP_bac = (TProfile2D*)fWgt->Get("pTPC_EP_bac");
                Read_BBC_EP_east = (TProfile2D*)fWgt->Get("pBBC_EP_east");
                Read_BBC_EP_west = (TProfile2D*)fWgt->Get("pBBC_EP_west");
                Read_EPD_EP1_east = (TProfile2D*)fWgt->Get("pEPD_EP1_east");
                Read_EPD_EP1_west = (TProfile2D*)fWgt->Get("pEPD_EP1_west");
                Read_EPD_EP_east = (TProfile2D*)fWgt->Get("pEPD_EP_east");
                Read_EPD_EP_west = (TProfile2D*)fWgt->Get("pEPD_EP_west");
                Read_ZDC_EP_east = (TProfile2D*)fWgt->Get("pZDC_EP_east");
                Read_ZDC_EP_west = (TProfile2D*)fWgt->Get("pZDC_EP_west");
cout<<"Loaded: TPC/BBC/EPD/ZDC EP corrections"<<endl;
		TH1D* Read_netChAsym   = (TH1D*)fWgt->Get("Hist_netChAsym");
		MeanNetChargeAsym= Read_netChAsym->GetMean();
cout<<"Loaded: MeanNetChargeAsym "<<endl;
		RMSNetChargeAsym = Read_netChAsym->GetRMS();
cout<<"Loaded: charge asymmetry "<<endl;
		TH1D* BBC_east_1 = (TH1D*)fWgt->Get("BBC1;1");
        	TH1D* BBC_east_2 = (TH1D*)fWgt->Get("BBC2;1");
        	TH1D* BBC_west_2 = (TH1D*)fWgt->Get("BBC7;1");
        	TH1D* BBC_west_1 = (TH1D*)fWgt->Get("BBC8;1");
        	float east_mean1 = (BBC_east_1->GetSum())/6.;
        	float east_mean2 = (BBC_east_2->GetSum())/10.;
        	float west_mean2 = (BBC_west_2->GetSum())/10.;
        	float west_mean1 = (BBC_west_1->GetSum())/6.;

		float content;
		for(int i=0;i<6;i++) {
			content = BBC_east_1->GetBinContent(i+1);
			BBC_gain_east[i] = (content>0)? east_mean1/content : 1;
			content = BBC_west_1->GetBinContent(i+1);
			BBC_gain_west[i] = (content>0)? west_mean1/content : 1;
		}
		for(int i=0;i<10;i++) {
			content = BBC_east_2->GetBinContent(i+1);
			BBC_gain_east[i+6] = (content>0)? 0.2*east_mean2/content : 1;
			content = BBC_west_2->GetBinContent(i+1);
			BBC_gain_west[i+6] = (content>0)? 0.2*west_mean2/content : 1;
		}
cout<<"Loaded: BBC gains"<<endl;
		float lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
  		float cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
		float par1[9] = {474.651,474.651,474.651,474.651,474.651,3.27243e+02,1.72351,1.72351,1.72351};
		float par2[9] = {3.55515,3.55515,3.55515,3.55515,3.55515,3.56938,-7.69075,-7.69075,-7.69075};
		float par3[9] = {1.80162,1.80162,1.80162,1.80162,1.80162,1.67113,4.72770,4.72770,4.72770};

    		for (int iy=1; iy<=9; iy++){
			for (int ix=1; ix<101; ix++){
      				double eta = wt->GetXaxis()->GetBinCenter(ix);
      				wt->SetBinContent(ix,iy, 1);
//      				wt->SetBinContent(ix,iy, (fabs(eta)>3.2)? lin[iy-1]*eta+cub[iy-1]*pow(eta,3):0);
//				wt2->SetBinContent(ix,iy, sqrt(1-1/par1[iy-1]/par1[iy-1]/cosh(eta)/cosh(eta))
//								/(1+exp((abs(eta)-par2[iy-1])/par3[iy-1])));
				wt2->SetBinContent(ix,iy, 1);
    			}
  		}
cout<<"Loaded: EPD 1st-order EP weights"<<endl;
cout<<"Loaded: EPD 2nd-order EP weights"<<endl;
		TH1D* Read_ZDC_e_h = (TH1D*)fWgt->Get("ZDC_e_h;1");
		TH1D* Read_ZDC_e_v = (TH1D*)fWgt->Get("ZDC_e_v;1");
                TH1D* Read_ZDC_w_h = (TH1D*)fWgt->Get("ZDC_w_h;1");
                TH1D* Read_ZDC_w_v = (TH1D*)fWgt->Get("ZDC_w_v;1");
		east_mean1 = (Read_ZDC_e_v->GetSum())/7.;
                east_mean2 = (Read_ZDC_e_h->GetSum())/8.;
                west_mean1 = (Read_ZDC_w_v->GetSum())/7.;
                west_mean2 = (Read_ZDC_w_h->GetSum())/8.;
		for(int i=0;i<7;i++) {
			content = Read_ZDC_e_v->GetBinContent(i+1);
			ZDC_gain_east_v[i] = (content>0)? east_mean1/content : 1;
                        content = Read_ZDC_w_v->GetBinContent(i+1);
                        ZDC_gain_west_v[i] = (content>0)? west_mean1/content : 1;
		}
                for(int i=0;i<8;i++) {
                        content = Read_ZDC_e_h->GetBinContent(i+1);
                        ZDC_gain_east_h[i] = (content>0)? east_mean2/content : 1;
                        content = Read_ZDC_w_h->GetBinContent(i+1);
                        ZDC_gain_west_h[i] = (content>0)? west_mean2/content : 1;
                }
cout<<"Loaded: ZDCSMD gains"<<endl;
		Read_ZDCcenter = (TProfile2D*)fWgt->Get("pZDCcenter");
cout<<"Loaded: ZDCSMD beam center"<<endl;
 
        }
	return 1;
}


/////////////////////////////
float GetPhiInBBC(int e_w, int bbcN) { //bbcN=1 to 24
const float phi_div=PI/6;
float bbc_phi=phi_div;
switch(bbcN) {
case 1: bbc_phi=3*phi_div;
break;
case 2: bbc_phi=phi_div;
break;
case 3: bbc_phi=-1*phi_div;
break;
case 4: bbc_phi=-3*phi_div;
break;
case 5: bbc_phi=-5*phi_div;
break;
case 6: bbc_phi=5*phi_div;
break;
case 7: bbc_phi= (gRandom->Rndm()>0.5) ? 2*phi_div:4*phi_div;
break;
case 8: bbc_phi=3*phi_div;
break;
case 9: bbc_phi=phi_div;
break;
case 10: bbc_phi=0.;
break;
case 11: bbc_phi=-phi_div;
break;
case 12: bbc_phi=(gRandom->Rndm()>0.5) ? -2*phi_div:-4*phi_div;
break;
case 13: bbc_phi=-3*phi_div;
break;
case 14: bbc_phi=-5*phi_div;
break;
case 15: bbc_phi=PI;
break;
case 16: bbc_phi=5*phi_div;
break;
case 17: bbc_phi=3*phi_div;
break;
case 18: bbc_phi=0.;
break;
case 19: bbc_phi=-3*phi_div;
break;
case 20: bbc_phi= PI;
break;
case 21: bbc_phi=3*phi_div;
break;
case 22: bbc_phi=0.;
break;
case 23: bbc_phi=-3*phi_div;
break;
case 24: bbc_phi= PI;
break;
}
if(e_w==0) {if(bbc_phi > -0.001) {bbc_phi = PI-bbc_phi;}
                else {bbc_phi = -PI-bbc_phi;}}
return bbc_phi;
}
///////////////////////////////////////////
float GetXYInZDC(int e_w, int v_h, int zdcN, int opt_raw ) { //zdcN=1 to 7 or 8
  //get position of each slat;strip starts from 1
	float zdcsmd_x[7] = {0.5,2,3.5,5,6.5,8,9.5};
 	float zdcsmd_y[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};
     	float zdcsmd_xx = zdcsmd_x0_e, zdcsmd_yy = zdcsmd_y0_e;
	if(1==e_w) {zdcsmd_xx = zdcsmd_x0_w; zdcsmd_yy = zdcsmd_y0_w;}
	if(1==opt_raw) {
		if(0==v_h) return zdcsmd_x[zdcN-1];
		else if(1==v_h) return zdcsmd_y[zdcN-1];
	}
	else {
		if(0==v_h) return zdcsmd_x[zdcN-1] - zdcsmd_xx;
		else if(1==v_h) return zdcsmd_y[zdcN-1] - zdcsmd_yy;
	}
	return 0;
}
/////////////////////////////////////////////
void DefineHistogram() {
	TString* histTitle;
	for(int i=0;i<15;i++) {
                histTitle = new TString("TPCmeanPhi_1_");
                *histTitle += i+1;
                WeightHist[i].pTPCmeanPhi_1 = new TProfile2D(histTitle->Data(),histTitle->Data(),8*order,0.5,8*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
                delete histTitle;

		histTitle = new TString("TPCmeanPhi_2_");
                *histTitle += i+1;
                WeightHist[i].pTPCmeanPhi_2 = new TProfile2D(histTitle->Data(),histTitle->Data(),8*order,0.5,8*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
                delete histTitle;

                histTitle = new TString("TPCmeanPhi_3_");
                *histTitle += i+1;
                WeightHist[i].pTPCmeanPhi_3 = new TProfile2D(histTitle->Data(),histTitle->Data(),8*order,0.5,8*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
                delete histTitle;

                histTitle = new TString("TPCmeanPhiAsso_1_");
                *histTitle += i+1;
                WeightHist[i].pTPCmeanPhiAsso_1 = new TProfile2D(histTitle->Data(),histTitle->Data(),8*order,0.5,8*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
                delete histTitle;

                histTitle = new TString("TPCmeanPhiAsso_2_");
                *histTitle += i+1;
                WeightHist[i].pTPCmeanPhiAsso_2 = new TProfile2D(histTitle->Data(),histTitle->Data(),8*order,0.5,8*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
                delete histTitle;

                histTitle = new TString("TPCmeanPhiAsso_3_");
                *histTitle += i+1;
                WeightHist[i].pTPCmeanPhiAsso_3 = new TProfile2D(histTitle->Data(),histTitle->Data(),8*order,0.5,8*order+0.5,(run_end-run_sta)/10,run_sta/10,run_end/10,-1,1,"");
                delete histTitle;

		histTitle = new TString("Hist_Phi_1_");
                *histTitle += i+1;
                PhiHist[i].Hist_Phi_1 = new TH2D(histTitle->Data(),histTitle->Data(),Phibin,-PI,PI,4,0.5,4.5);
                delete histTitle;

                histTitle = new TString("Hist_Phi_2_");
                *histTitle += i+1;
                PhiHist[i].Hist_Phi_2 = new TH2D(histTitle->Data(),histTitle->Data(),Phibin,-PI,PI,4,0.5,4.5);
                delete histTitle;

                histTitle = new TString("Hist_Phi_3_");
                *histTitle += i+1;
                PhiHist[i].Hist_Phi_3 = new TH2D(histTitle->Data(),histTitle->Data(),Phibin,-PI,PI,4,0.5,4.5);
                delete histTitle;

                histTitle = new TString("Hist_Phi_new_1_");
                *histTitle += i+1;
                PhiHist[i].Hist_Phi_new_1 = new TH2D(histTitle->Data(),histTitle->Data(),Phibin,-PI,PI,4,0.5,4.5);
                delete histTitle;

                histTitle = new TString("Hist_Phi_new_2_");
                *histTitle += i+1;
                PhiHist[i].Hist_Phi_new_2 = new TH2D(histTitle->Data(),histTitle->Data(),Phibin,-PI,PI,4,0.5,4.5);
                delete histTitle;

                histTitle = new TString("Hist_Phi_new_3_");
                *histTitle += i+1;
                PhiHist[i].Hist_Phi_new_3 = new TH2D(histTitle->Data(),histTitle->Data(),Phibin,-PI,PI,4,0.5,4.5);
                delete histTitle;
        }

}

void EPD_hits(TClonesArray* mEpdHits){
                long nepdMIPsE = 0.;
                long nepdMIPsW = 0.;
		Hist_mEpdHits->Fill(mEpdHits->GetEntries());
          for(int hit=0;hit<mEpdHits->GetEntries();hit++){
                    int tileId, ring, TT, PP, EW, ADC;
                    float nMip;
                    StPicoEpdHit* epdHit = (StPicoEpdHit*)((*mEpdHits)[hit]);
                    tileId = epdHit->side();//+1 if tile is on West side; -1 if on East side
                    EW = (tileId<0)?0:1; //E 0, W 1
                    nMip = epdHit->nMIP();

                    if(nMip<mEPDthresh) continue;
                    if(EW == 0){ nepdMIPsE += (nMip>mEPDMax ? nMip : 1);}
                    if(EW == 1){ nepdMIPsW += (nMip>mEPDMax ? nMip : 1);}

          }
                  runidvsepdEhits->Fill(Run,nepdMIPsE);
                  runidvsepdWhits->Fill(Run,nepdMIPsW);


}

