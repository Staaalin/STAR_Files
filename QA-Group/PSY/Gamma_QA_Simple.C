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
                }
                RefMult = NumCharge;
                // RefMult = FxtMult;
                // cout<<"FxtMult = "<<FxtMult<<endl;
		if(!IsGoodEvent(cen)) continue;

                for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
                        StPicoTrack *picoTrack = dst->track(iTrack);
                        if (! picoTrack)            continue;
                        if (! picoTrack->charge())  continue;
                        if (! picoTrack->isPrimary()) continue;

                        EtaAsso   = picoTrack->pMom().Eta();
                        PtAsso    = picoTrack->pMom().Pt();
                        PhiAsso   = picoTrack->pMom().Phi();
                        DCAGlobalAsso = picoTrack->gDCA(pV).Mag();
                        ChargeAsso= picoTrack->charge();
                        gDCAxy = picoTrack->gDCAxy(pVx,pVy);
                        dEdx = picoTrack->dEdx();

                        pTemp_track->Fill(1,PtAsso);
                        pTemp_track->Fill(2,EtaAsso);
                        pTemp_track->Fill(3,PhiAsso);
                        pTemp_track->Fill(4,ChargeAsso);
                        pTemp_track->Fill(5,DCAGlobalAsso);
                        pTemp_track->Fill(6,gDCAxy);
                        pTemp_track_dedx->Fill(1,dEdx);
                        
                }
        //add Yu's
                runidvsdedx->Fill(Run,pTemp_track_dedx->GetBinContent(1)); 
                runidvsavgpt->Fill(Run,pTemp_track->GetBinContent(1));
                runidvsavgeta->Fill(Run,pTemp_track->GetBinContent(2));
                runidvsavgdca->Fill(Run,pTemp_track->GetBinContent(5));
                runidvsavgphi->Fill(Run,pTemp_track->GetBinContent(3));
                runidvsgdcaxy->Fill(Run,pTemp_track->GetBinContent(6));
                runidvsdcaxysigma->Fill(Run,pTemp_track->GetBinError(6)*sqrt(pTemp_track->GetBinEntries(6)));
                
        }
}

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
