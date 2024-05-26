#include <stdlib.h>
// #include <sys/types.h>
// #include <sys/stat.h>
// #include <dirent.h>
// #include <random>
#include "math.h"
#include "string.h"
#include <vector>
// #ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TGraph.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TUnixSystem.h"
#include "TRandom3.h"
// #endif
#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
using namespace std;

Double_t massList(int PID , float Energy)
{
    Double_t Result;
    if (Energy == 62.0){
        switch (PID)
        {
            case 321 :
                Result = 0.493677;
                break;
            case -321 :
                Result = 0.493677;
                break;
            case 211 :
                Result = 0.13957;
                break;
            case -211 :
                Result = 0.13957;
                break;
            case 3334 :// OmegaFitMass
                Result = 1.6725;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 1.6727;
                break;
            case 3312 :// XiFitMass
                Result = 1.3223;
                break;
            case -3312 :// XiBarFitMass
                Result = 1.3223;
                break;
            case 3122 :// LambdaFitMass
                Result = 1.1161;
                break;
            case -3122 :// LambdaBarFitMass
                Result = 1.1161;
                break;
            default :
                Result = 0;
        }
    }
    return Result;
}

Double_t massListSigma(int PID , float Energy)
{
    Double_t Result;
    if (Energy == 62.0){
        switch (PID)
        {
            case 3334 :// OmegaFitMass
                Result = 0.0029;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 0.0024;
                break;
            case 3312 :// XiFitMass
                Result = 0.0024;
                break;
            case -3312 :// XiBarFitMass
                Result = 0.0024;
                break;
            case 3122 :// LambdaFitMass
                Result = 0.0020;
                break;
            case -3122 :// LambdaBarFitMass
                Result = 0.0020;
                break;
            default :
                Result = 0;
        }
    }
    return Result;
}

void Subtract(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,
              float Energy)
{

    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) 
        std::vector<Int_t>   *id              = nullptr;
        std::vector<Int_t>   *evtID           = nullptr;
        std::vector<Int_t>   *runID           = nullptr;
        std::vector<Float_t> *px              = nullptr;
        std::vector<Float_t> *py              = nullptr;
        std::vector<Float_t> *pz              = nullptr;
        std::vector<Float_t> *mass            = nullptr;
        std::vector<Float_t> *energy          = nullptr;
        
        TBranch *bid              = nullptr;
        TBranch *bevtID           = nullptr;
        TBranch *brunID           = nullptr;
        TBranch *bpx              = nullptr;
        TBranch *bpy              = nullptr;
        TBranch *bpz              = nullptr;
        TBranch *bmass            = nullptr;
        TBranch *benergy          = nullptr;

        std::vector<Float_t> *dEdx               = nullptr;
        std::vector<Float_t> *m2                 = nullptr;
        std::vector<Float_t> *dcatopv            = nullptr;
        std::vector<int>     *hasTOF             = nullptr;
        std::vector<Float_t> *nSigmaProton       = nullptr;
        std::vector<Float_t> *nSigmaPion         = nullptr;
        std::vector<Float_t> *nSigmaKaon         = nullptr;
        std::vector<double>  *zTOF_proton        = nullptr;
        std::vector<double>  *zTOF_pion          = nullptr;
        std::vector<double>  *zTOF_kaon          = nullptr;
	    std::vector<int>     *IfConfuse          = nullptr;
        std::vector<double>  *Decay_Length       = nullptr;
	    std::vector<int>     *IfBadReconstructed = nullptr;

        TBranch *bdEdx               = nullptr;
        TBranch *bm2                 = nullptr;
        TBranch *bdcatopv            = nullptr;
        TBranch *bhasTOF             = nullptr;
        TBranch *bnSigmaProton       = nullptr;
        TBranch *bnSigmaPion         = nullptr;
        TBranch *bnSigmaKaon         = nullptr;
        TBranch *bzTOF_proton        = nullptr;
        TBranch *bzTOF_pion          = nullptr;
        TBranch *bzTOF_kaon          = nullptr;
	    TBranch *bIfConfuse          = nullptr;
	    TBranch *bDecay_Length       = nullptr;
	    TBranch *bIfBadReconstructed = nullptr;
    
    #else
        #if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
            std::vector<Int_t>   *id              = NULL;
            std::vector<Int_t>   *evtID           = NULL;
            std::vector<Int_t>   *runID           = NULL;
            std::vector<Float_t> *px              = NULL;
            std::vector<Float_t> *py              = NULL;
            std::vector<Float_t> *pz              = NULL;
            std::vector<Float_t> *mass            = NULL;
            std::vector<Float_t> *energy          = NULL;
            
            TBranch *bid              = NULL;
            TBranch *bevtID           = NULL;
            TBranch *brunID           = NULL;
            TBranch *bpx              = NULL;
            TBranch *bpy              = NULL;
            TBranch *bpz              = NULL;
            TBranch *bmass            = NULL;
            TBranch *benergy          = NULL;
        
            std::vector<Float_t> *dEdx               = NULL;
            std::vector<Float_t> *m2                 = NULL;
            std::vector<Float_t> *dcatopv            = NULL;
            std::vector<int>     *hasTOF             = NULL;
            std::vector<Float_t> *nSigmaProton       = NULL;
            std::vector<Float_t> *nSigmaPion         = NULL;
            std::vector<Float_t> *nSigmaKaon         = NULL;
            std::vector<double>  *zTOF_proton        = NULL;
            std::vector<double>  *zTOF_pion          = NULL;
            std::vector<double>  *zTOF_kaon          = NULL;
            std::vector<int>     *IfConfuse          = NULL;
            std::vector<double>  *Decay_Length       = NULL;
	        std::vector<int>     *IfBadReconstructed = NULL;

            TBranch *bdEdx               = NULL;
            TBranch *bm2                 = NULL;
            TBranch *bdcatopv            = NULL;
            TBranch *bhasTOF             = NULL;
            TBranch *bnSigmaProton       = NULL;
            TBranch *bnSigmaPion         = NULL;
            TBranch *bnSigmaKaon         = NULL;
            TBranch *bzTOF_proton        = NULL;
            TBranch *bzTOF_pion          = NULL;
            TBranch *bzTOF_kaon          = NULL;
            TBranch *bIfConfuse          = NULL;
            TBranch *bDecay_Length       = NULL;
	        TBranch *bIfBadReconstructed = NULL;

        #else
            std::vector<Int_t>   *id              = 0;
            std::vector<Int_t>   *evtID           = 0;
            std::vector<Int_t>   *runID           = 0;
            std::vector<Float_t> *px              = 0;
            std::vector<Float_t> *py              = 0;
            std::vector<Float_t> *pz              = 0;
            std::vector<Float_t> *mass            = 0;
            std::vector<Float_t> *energy          = 0;
            
            TBranch *bid              = 0;
            TBranch *bevtID           = 0;
            TBranch *brunID           = 0;
            TBranch *bpx              = 0;
            TBranch *bpy              = 0;
            TBranch *bpz              = 0;
            TBranch *bmass            = 0;
            TBranch *benergy          = 0;

            std::vector<Float_t> *dEdx               = 0;
            std::vector<Float_t> *m2                 = 0;
            std::vector<Float_t> *dcatopv            = 0;
            std::vector<int>     *hasTOF             = 0;
            std::vector<Float_t> *nSigmaProton       = 0;
            std::vector<Float_t> *nSigmaPion         = 0;
            std::vector<Float_t> *nSigmaKaon         = 0;
            std::vector<double>  *zTOF_proton        = 0;
            std::vector<double>  *zTOF_pion          = 0;
            std::vector<double>  *zTOF_kaon          = 0;
            std::vector<int>     *IfConfuse          = 0;
            std::vector<double>  *Decay_Length       = 0;
	        std::vector<int>     *IfBadReconstructed = 0;

            TBranch *bdEdx               = 0;
            TBranch *bm2                 = 0;
            TBranch *bdcatopv            = 0;
            TBranch *bhasTOF             = 0;
            TBranch *bnSigmaProton       = 0;
            TBranch *bnSigmaPion         = 0;
            TBranch *bnSigmaKaon         = 0;
            TBranch *bzTOF_proton        = 0;
            TBranch *bzTOF_pion          = 0;
            TBranch *bzTOF_kaon          = 0;
            TBranch *bIfConfuse          = 0;
            TBranch *bDecay_Length       = 0;
	        TBranch *bIfBadReconstructed = 0;

        #endif
    #endif

    // TString FileName = "output_";
    // cout<<"Start Running"<<endl;
    // cout<<StartFileIndex<<endl;
    // cout<<EndFileIndex<<endl;


    #define HyperonPhaseNum  2
    #define MesonPhaseNum    2
    TString HyperonPhase[HyperonPhaseNum] = {"m1p1Rap","0p1Rap"}
    TString MesonPhase[MesonPhaseNum] = {"Allp","0p5t2p"}
    TH1D *dNdy_Kaonp  [MesonPhaseNum];
    TH1D *dNdy_Kaonm  [MesonPhaseNum];
    TH1D *dNdy_Pionp  [MesonPhaseNum];
    TH1D *dNdy_Pionm  [MesonPhaseNum];
    TH1D *dNdy_Lambda [HyperonPhaseNum];
    TH1D *dNdy_Lambdab[HyperonPhaseNum];
    TH1D *dNdy_Omega  [HyperonPhaseNum];
    TH1D *dNdy_Omegab [HyperonPhaseNum];
    TH1D *dNdy_Xi     [HyperonPhaseNum];
    TH1D *dNdy_Xib    [HyperonPhaseNum];
    TH1D *k_Kp_Omega  [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Omegab [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Omega  [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Omegab [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Xi     [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Xib    [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Xi     [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Xib    [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Lambda [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Lambdab[MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Lambda [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Lambdab[MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Omega  [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Omegab [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Omega  [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Omegab [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Xi     [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Xib    [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Xi     [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Xib    [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Lambda [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Lambdab[MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Lambda [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Lambdab[MesonPhaseNum][HyperonPhaseNum];

    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
        TString HistName = "dNdy_Kaonp_";HistName += MesonPhase[Itr];
        dNdy_Kaonp[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Kaonp[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Kaonp[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Kaonm_";HistName += MesonPhase[Itr];
        dNdy_Kaonm[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Kaonm[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Kaonm[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Pionp_";HistName += MesonPhase[Itr];
        dNdy_Pionp[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Pionp[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Pionp[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Pionm_";HistName += MesonPhase[Itr];
        dNdy_Pionm[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Pionm[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Pionm[Itr]->GetYaxis()->SetTitle("Counts");
    }
    for (int Itr = 0;Itr < HyperonPhaseNum;Itr++) {
        TString HistName = "dNdy_Lambda_";HistName += HyperonPhase[Itr];
        dNdy_Lambda[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Lambda[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Lambda[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Lambdab_";HistName += HyperonPhase[Itr];
        dNdy_Lambdab[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Lambdab[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Lambdab[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Omega_";HistName += HyperonPhase[Itr];
        dNdy_Omega[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Omega[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Omega[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Omegab_";HistName += HyperonPhase[Itr];
        dNdy_Omegab[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Omegab[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Omegab[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Xi_";HistName += HyperonPhase[Itr];
        dNdy_Xi[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Xi[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Xi[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Xib_";HistName += HyperonPhase[Itr];
        dNdy_Xib[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Xib[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Xib[Itr]->GetYaxis()->SetTitle("Counts");
    }
    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
        for (int Jtr = 0;Jtr < HyperonPhaseNum;Jtr++) {
            TString HistName = "k_Kp_Omega_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Omega  [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Omegab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Omegab [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Omega_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Omega  [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Omegab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Omegab [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Xi_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Xi     [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Xib_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Xib    [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Xi_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Xi     [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Xib_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Xib    [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Lambda_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Lambda [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Lambdab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Lambdab[Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Lambda_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Lambda [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Lambdab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Lambdab[Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);

            HistName = "k_Pip_Omega_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Omega  [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Omegab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Omegab [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Omega_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Omega  [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Omegab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Omegab [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Xi_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Xi     [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Xib_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Xib    [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Xi_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Xi     [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Xib_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Xib    [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Lambda_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Lambda [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Lambdab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Lambdab[Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Lambda_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Lambda [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Lambdab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Lambdab[Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);

            k_Kp_Omega   [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Omegab  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Omega   [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Omegab  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Xi      [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Xib     [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Xi      [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Xib     [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Lambda  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Lambdab [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Lambda  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Lambdab [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Omega  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Omegab [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Omega  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Omegab [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Xi     [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Xib    [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Xi     [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Xib    [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Lambda [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Lambdab[Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Lambda [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Lambdab[Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Omega   [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Omegab  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Omega   [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Omegab  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Xi      [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Xib     [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Xi      [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Xib     [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Lambda  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Lambdab [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Lambda  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Lambdab [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Omega  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Omegab [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Omega  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Omegab [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Xi     [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Xib    [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Xi     [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Xib    [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Lambda [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Lambdab[Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Lambda [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Lambdab[Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
        }
    }


    //load data  
    TChain *hadronTree = new TChain("hadronTree");
    for(int i=StartFileIndex;i <= EndFileIndex;i++){
        TString filename = MidName;
        filename+=i;
        filename+=".root";
        hadronTree->Add(filename);
        // cout<<"Add "<<filename<<" Successfully"<<endl;
    }
    
    Int_t mult,npart;

    hadronTree->SetBranchAddress("PDGMult",&mult);
    hadronTree->SetBranchAddress("PDG",&id,&bid);
    hadronTree->SetBranchAddress("mix_px",&px,&bpx);
    hadronTree->SetBranchAddress("mix_py",&py,&bpy);
    hadronTree->SetBranchAddress("mix_pz",&pz,&bpz);
    hadronTree->SetBranchAddress("InvariantMass",&mass,&bmass);
    hadronTree->SetBranchAddress("nSigmaProton",&nSigmaProton,&bnSigmaProton);
    hadronTree->SetBranchAddress("nSigmaPion"  ,&nSigmaPion  ,&bnSigmaPion  );
    hadronTree->SetBranchAddress("nSigmaKaon"  ,&nSigmaKaon  ,&bnSigmaKaon  );

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;
    
    TLorentzVector p0,p3;
    double kstar, rap;
    TVector3 BetaTemp;

    std::vector<float> Kp_px[MesonPhaseNum],Km_px[MesonPhaseNum],Pip_px[MesonPhaseNum],Pim_px[MesonPhaseNum];
    std::vector<float> Kp_py[MesonPhaseNum],Km_py[MesonPhaseNum],Pip_py[MesonPhaseNum],Pim_py[MesonPhaseNum];
    std::vector<float> Kp_pz[MesonPhaseNum],Km_pz[MesonPhaseNum],Pip_pz[MesonPhaseNum],Pim_pz[MesonPhaseNum];
    std::vector<float> Kp_y [MesonPhaseNum],Km_y [MesonPhaseNum],Pip_y [MesonPhaseNum],Pim_y [MesonPhaseNum];
    std::vector<float> Lambda_px[HyperonPhaseNum],Lambdab_px[HyperonPhaseNum],Omega_px[HyperonPhaseNum],Omegab_px[HyperonPhaseNum],Xi_px[HyperonPhaseNum],Xib_px[HyperonPhaseNum];
    std::vector<float> Lambda_py[HyperonPhaseNum],Lambdab_py[HyperonPhaseNum],Omega_py[HyperonPhaseNum],Omegab_py[HyperonPhaseNum],Xi_py[HyperonPhaseNum],Xib_py[HyperonPhaseNum];
    std::vector<float> Lambda_pz[HyperonPhaseNum],Lambdab_pz[HyperonPhaseNum],Omega_pz[HyperonPhaseNum],Omegab_pz[HyperonPhaseNum],Xi_pz[HyperonPhaseNum],Xib_pz[HyperonPhaseNum];
    std::vector<float> Lambda_y [HyperonPhaseNum],Lambdab_y [HyperonPhaseNum],Omega_y [HyperonPhaseNum],Omegab_y [HyperonPhaseNum],Xi_y [HyperonPhaseNum],Xib_y [HyperonPhaseNum];

    int KaonpPID = 321,KaonmPID = -321,PionpPID = 211,PionmPID = -211,LambdaPID = 3122,LambdabPID = -3122,XiPID = 3312,XibPID = -3312,OmegaPID = 3334,OmegabPID = -3334;

    bool MesonStates[MesonPhaseNum],HyperonStates[HyperonPhaseNum];
    TString States;

    for (int i=0;i<nentries;i++){
        hadronTree->GetEntry(i);
        // cout<<mult<<endl;
        // if(b>7){continue;}
        for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
            Kp_px[Itr].resize(0);Km_px[Itr].resize(0);Pip_px[Itr].resize(0);Pim_px[Itr].resize(0);
            Kp_py[Itr].resize(0);Km_py[Itr].resize(0);Pip_py[Itr].resize(0);Pim_py[Itr].resize(0);
            Kp_pz[Itr].resize(0);Km_pz[Itr].resize(0);Pip_pz[Itr].resize(0);Pim_pz[Itr].resize(0);
            Kp_y [Itr].resize(0);Km_y [Itr].resize(0);Pip_y [Itr].resize(0);Pim_y [Itr].resize(0);
        }
        for (int Itr = 0;Itr < HyperonPhaseNum;Itr++) {
            Lambda_px[Itr].resize(0);Lambdab_px[Itr].resize(0);Omega_px[Itr].resize(0);Omegab_px[Itr].resize(0);Xi_px[Itr].resize(0);Xib_px[Itr].resize(0);
            Lambda_py[Itr].resize(0);Lambdab_py[Itr].resize(0);Omega_py[Itr].resize(0);Omegab_py[Itr].resize(0);Xi_py[Itr].resize(0);Xib_py[Itr].resize(0);
            Lambda_pz[Itr].resize(0);Lambdab_pz[Itr].resize(0);Omega_pz[Itr].resize(0);Omegab_pz[Itr].resize(0);Xi_pz[Itr].resize(0);Xib_pz[Itr].resize(0);
            Lambda_y [Itr].resize(0);Lambdab_y [Itr].resize(0);Omega_y [Itr].resize(0);Omegab_y [Itr].resize(0);Xi_y [Itr].resize(0);Xib_y [Itr].resize(0);
        }
        for (int j=0;j < mult;j++){

            // cout<<"id->at(j) = "<<id->at(j)<<endl;
            for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
                MesonStates[Itr] = false;
            }
            for (int Itr = 0;Itr < HyperonPhaseNum;Itr++) {
                HyperonStates[Itr] = false;
            }
            if ((fabs(id->at(j)) == LambdaPID) || (fabs(id->at(j)) == XiPID || (fabs(id->at(j)) == OmegaPID))){
                if ((fabs(mass->at(j) - massList(id->at(j),Energy))) > 3*massListSigma(id->at(j),Energy)) {
                    continue;
                }
                p0.SetXYZM(px->at(j),py->at(j),pz->at(j),massList(id->at(j),Energy));
                rap = p0.Rapidity();
                if (-1 > rap || rap > 1) {
                    continue;
                }

                if (rap > 0) {
                    States = "0p1Rap";
                    for (int Itr = 0;Itr < HyperonPhaseNum;Itr++) {
                        if (HyperonPhase[Itr] == States){
                            HyperonStates[Itr] = true;
                            break;
                        }
                    }
                }
                if (true){
                    States = "m1p1Rap";
                    for (int Itr = 0;Itr < HyperonPhaseNum;Itr++) {
                        if (HyperonPhase[Itr] == States){
                            HyperonStates[Itr] = true;
                            break;
                        }
                    }
                }
                for (int HyperonPhaseIndex = 0;HyperonPhaseIndex < HyperonPhaseNum;HyperonPhaseIndex++){
                    if (HyperonStates[HyperonPhaseIndex] == false) {continue;}
                    if (id->at(j) == LambdaPID) {
                        Lambda_px[HyperonPhaseIndex].push_back(px->at(j));
                        Lambda_py[HyperonPhaseIndex].push_back(py->at(j));
                        Lambda_pz[HyperonPhaseIndex].push_back(pz->at(j));
                        Lambda_y [HyperonPhaseIndex].push_back(rap);
                        dNdy_Lambda[HyperonPhaseIndex]->Fill(rap);
                    }
                    else if (id->at(j) == LambdabPID) {
                        Lambdab_px[HyperonPhaseIndex].push_back(px->at(j));
                        Lambdab_py[HyperonPhaseIndex].push_back(py->at(j));
                        Lambdab_pz[HyperonPhaseIndex].push_back(pz->at(j));
                        Lambdab_y [HyperonPhaseIndex].push_back(rap);
                        dNdy_Lambdab[HyperonPhaseIndex]->Fill(rap);
                    }
                    else if (id->at(j) == XiPID) {
                        Xi_px[HyperonPhaseIndex].push_back(px->at(j));
                        Xi_py[HyperonPhaseIndex].push_back(py->at(j));
                        Xi_pz[HyperonPhaseIndex].push_back(pz->at(j));
                        Xi_y [HyperonPhaseIndex].push_back(rap);
                        dNdy_Xi[HyperonPhaseIndex]->Fill(rap);
                    }
                    else if (id->at(j) == XibPID) {
                        Xib_px[HyperonPhaseIndex].push_back(px->at(j));
                        Xib_py[HyperonPhaseIndex].push_back(py->at(j));
                        Xib_pz[HyperonPhaseIndex].push_back(pz->at(j));
                        Xib_y [HyperonPhaseIndex].push_back(rap);
                        dNdy_Xib[HyperonPhaseIndex]->Fill(rap);
                    }
                    else if (id->at(j) == OmegaPID) {
                        Omega_px[HyperonPhaseIndex].push_back(px->at(j));
                        Omega_py[HyperonPhaseIndex].push_back(py->at(j));
                        Omega_pz[HyperonPhaseIndex].push_back(pz->at(j));
                        Omega_y [HyperonPhaseIndex].push_back(rap);
                        dNdy_Omega[HyperonPhaseIndex]->Fill(rap);
                    }
                    else if (id->at(j) == OmegabPID) {
                        Omegab_px[HyperonPhaseIndex].push_back(px->at(j));
                        Omegab_py[HyperonPhaseIndex].push_back(py->at(j));
                        Omegab_pz[HyperonPhaseIndex].push_back(pz->at(j));
                        Omegab_y [HyperonPhaseIndex].push_back(rap);
                        dNdy_Omegab[HyperonPhaseIndex]->Fill(rap);
                    }
                }
            }

            if ((fabs(id->at(j)) == KaonpPID) || (fabs(id->at(j)) == PionpPID)) {
                if ((fabs(id->at(j)) == KaonpPID) && ((nSigmaKaon->at(j) < -2) || (nSigmaKaon->at(j) > 2))) {
                    continue;
                }
                if ((fabs(id->at(j)) == PionpPID) && ((nSigmaPion->at(j) < -3) || (nSigmaPion->at(j) > 3))) {
                    continue;
                }
                p0.SetXYZM(px->at(j),py->at(j),pz->at(j),massList(id->at(j),Energy));
                rap = p0.Rapidity();
                float Meg = pow(pow(px->at(j),2)+pow(py->at(j),2)+pow(pz->at(j),2),0.5);
                if ((0.5 < Meg) && (Meg < 2)) {
                    States = "0p5t2p";
                    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
                        if (MesonPhase[Itr] == States){
                            MesonStates[Itr] = true;
                            break;
                        }
                    }
                }
                if(true){
                    States = "Allp";
                    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
                        if (MesonPhase[Itr] == States){
                            MesonStates[Itr] = true;
                            break;
                        }
                    }
                }
                for (int MesonPhaseIndex = 0;MesonPhaseIndex < HyperonPhaseNum;MesonPhaseIndex++){
                    if (MesonStates[MesonPhaseIndex] == false) {continue;}
                    if (id->at(j) == PionpPID) {
                        Pip_px[MesonPhaseIndex].push_back(px->at(j));
                        Pip_py[MesonPhaseIndex].push_back(py->at(j));
                        Pip_pz[MesonPhaseIndex].push_back(pz->at(j));
                        Pip_y [MesonPhaseIndex].push_back(rap);
                        dNdy_Pionp[MesonPhaseIndex]->Fill(rap);
                    }
                    else if (id->at(j) == PionmPID) {
                        Pim_px[MesonPhaseIndex].push_back(px->at(j));
                        Pim_py[MesonPhaseIndex].push_back(py->at(j));
                        Pim_pz[MesonPhaseIndex].push_back(pz->at(j));
                        Pim_y [MesonPhaseIndex].push_back(rap);
                        dNdy_Pionm[MesonPhaseIndex]->Fill(rap);
                    }
                    else if (id->at(j) == KaonpPID) {
                        Kp_px[MesonPhaseIndex].push_back(px->at(j));
                        Kp_py[MesonPhaseIndex].push_back(py->at(j));
                        Kp_pz[MesonPhaseIndex].push_back(pz->at(j));
                        Kp_y [MesonPhaseIndex].push_back(rap);
                        dNdy_Kaonp[MesonPhaseIndex]->Fill(rap);
                    }
                    else if (id->at(j) == KaonmPID) {
                        Km_px[MesonPhaseIndex].push_back(px->at(j));
                        Km_py[MesonPhaseIndex].push_back(py->at(j));
                        Km_pz[MesonPhaseIndex].push_back(pz->at(j));
                        Km_y [MesonPhaseIndex].push_back(rap);
                        dNdy_Kaonm[MesonPhaseIndex]->Fill(rap);
                    }
                }
            }
        }
        // cout<<"Here OK"<<endl;

        for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
            for (int Jtr = 0;Jtr < HyperonPhaseNum;Jtr++) {
                // cout<<"KEK0"<<endl;
                for (int Ktr = 0;Ktr < Kp_px.size();Ktr++) {
                    TLorentzVector p1;
                    // float Tpx=Kp_px->at(Ktr),Tpy=Kp_py->at(Ktr),Tpz=Kp_pz->at(Ktr);
                    // p1.SetXYZM(Tpx,Tpy,Tpz,massList(KaonpPID,Energy));
                    p1.SetXYZM(Kp_px->at(Ktr),Kp_py->at(Ktr),Kp_pz->at(Ktr),massList(KaonpPID,Energy));
                    // cout<<p1.Rapidity()<<endl;
                    // cout<<"KEK1"<<endl;
                    // float ALength = pow(pow(A_px->at(j),2)+pow(A_py->at(j),2)+pow(A_pz->at(j),2) ,0.5);
                    for (int Ntr = 0;Ntr < Lambda_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambda_px->at(Ntr),Lambda_py->at(Ntr),Lambda_pz->at(Ntr),massList(LambdaPID,Energy));
                        // cout<<"KEK2"<<endl;

                        p3 = p4 + p2;
                        // float dpt = fabs(p1.Perp() - p2.Perp());
                        // Apt->Fill(dpt);
                        cout<<"2"<<endl;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        cout<<"1"<<endl;
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Lambda[Itr][Jtr]->Fill(kstar);
                        // cout<<"KEK3"<<endl;
                        // Ay->Fill(fabs(A_y[j] - B_y[k]));

                        // float BLength = pow(pow(B_px[k],2)+pow(B_py[k],2)+pow(B_pz[k],2) ,0.5);
                        // float CosAB = (A_px->at(j)*B_px[k] + A_py->at(j)*B_py[k] + A_pz->at(j)*B_pz[k])/(ALength*BLength);
                        // Aphi->Fill(acos(CosAB));
                    }
                    for (int Ntr = 0;Ntr < Lambdab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambdab_px->at(Ntr),Lambdab_py->at(Ntr),Lambdab_pz->at(Ntr),massList(LambdabPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Lambdab[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xi_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xi_px->at(Ntr),Xi_py->at(Ntr),Xi_pz->at(Ntr),massList(XiPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Xi[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xib_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xib_px->at(Ntr),Xib_py->at(Ntr),Xib_pz->at(Ntr),massList(XibPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Xib[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omega_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omega_px->at(Ntr),Omega_py->at(Ntr),Omega_pz->at(Ntr),massList(OmegaPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Omega[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omegab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omegab_px->at(Ntr),Omegab_py->at(Ntr),Omegab_pz->at(Ntr),massList(OmegabPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Omegab[Itr][Jtr]->Fill(kstar);
                    }
                }
                for (int Ktr = 0;Ktr < Km_px.size();Ktr++) {
                    TLorentzVector p1;
                    p1.SetXYZM(Km_px->at(Ktr),Km_py->at(Ktr),Km_pz->at(Ktr),massList(KaonmPID,Energy));
                    // cout<<"KEK4"<<endl;
                    for (int Ntr = 0;Ntr < Lambda_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambda_px->at(Ntr),Lambda_py->at(Ntr),Lambda_pz->at(Ntr),massList(LambdaPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Lambda[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Lambdab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambdab_px->at(Ntr),Lambdab_py->at(Ntr),Lambdab_pz->at(Ntr),massList(LambdabPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Lambdab[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xi_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xi_px->at(Ntr),Xi_py->at(Ntr),Xi_pz->at(Ntr),massList(XiPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Xi[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xib_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xib_px->at(Ntr),Xib_py->at(Ntr),Xib_pz->at(Ntr),massList(XibPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Xib[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omega_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omega_px->at(Ntr),Omega_py->at(Ntr),Omega_pz->at(Ntr),massList(OmegaPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Omega[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omegab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omegab_px->at(Ntr),Omegab_py->at(Ntr),Omegab_pz->at(Ntr),massList(OmegabPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Omegab[Itr][Jtr]->Fill(kstar);
                    }
                }
                for (int Ktr = 0;Ktr < Pip_px.size();Ktr++) {
                    TLorentzVector p1;
                    p1.SetXYZM(Pip_px->at(Ktr),Pip_py->at(Ktr),Pip_pz->at(Ktr),massList(PionpPID,Energy));
                    // cout<<"KEK5"<<endl;
                    for (int Ntr = 0;Ntr < Lambda_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambda_px->at(Ntr),Lambda_py->at(Ntr),Lambda_pz->at(Ntr),massList(LambdaPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Lambda[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Lambdab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambdab_px->at(Ntr),Lambdab_py->at(Ntr),Lambdab_pz->at(Ntr),massList(LambdabPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Lambdab[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xi_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xi_px->at(Ntr),Xi_py->at(Ntr),Xi_pz->at(Ntr),massList(XiPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Xi[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xib_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xib_px->at(Ntr),Xib_py->at(Ntr),Xib_pz->at(Ntr),massList(XibPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Xib[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omega_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omega_px->at(Ntr),Omega_py->at(Ntr),Omega_pz->at(Ntr),massList(OmegaPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Omega[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omegab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omegab_px->at(Ntr),Omegab_py->at(Ntr),Omegab_pz->at(Ntr),massList(OmegabPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Omegab[Itr][Jtr]->Fill(kstar);
                    }
                }
                for (int Ktr = 0;Ktr < Pim_px.size();Ktr++) {
                    TLorentzVector p1;
                    p1.SetXYZM(Pim_px->at(Ktr),Pim_py->at(Ktr),Pim_pz->at(Ktr),massList(PionmPID,Energy));
                    // cout<<"KEK6"<<endl;
                    for (int Ntr = 0;Ntr < Lambda_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambda_px->at(Ntr),Lambda_py->at(Ntr),Lambda_pz->at(Ntr),massList(LambdaPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Lambda[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Lambdab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambdab_px->at(Ntr),Lambdab_py->at(Ntr),Lambdab_pz->at(Ntr),massList(LambdabPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Lambdab[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xi_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xi_px->at(Ntr),Xi_py->at(Ntr),Xi_pz->at(Ntr),massList(XiPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Xi[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xib_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xib_px->at(Ntr),Xib_py->at(Ntr),Xib_pz->at(Ntr),massList(XibPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Xib[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omega_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omega_px->at(Ntr),Omega_py->at(Ntr),Omega_pz->at(Ntr),massList(OmegaPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Omega[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omegab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omegab_px->at(Ntr),Omegab_py->at(Ntr),Omegab_pz->at(Ntr),massList(OmegabPID,Energy));

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Omegab[Itr][Jtr]->Fill(kstar);
                    }
                }
            }
        }
        // cout<<"There OK"<<endl;

    }

    TString OutputFileName = "Cor_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *file = new TFile(OutputFileName, "RECREATE");

    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
        dNdy_Kaonp[Itr]->Write();
        dNdy_Kaonm[Itr]->Write();
        dNdy_Pionp[Itr]->Write();
        dNdy_Pionm[Itr]->Write();
    }
    for (int Itr = 0;Itr < HyperonPhaseNum;Itr++) {
        dNdy_Lambda[Itr]  ->Write();
        dNdy_Lambdab[Itr] ->Write();
        dNdy_Omega[Itr]   ->Write();
        dNdy_Omegab[Itr]  ->Write();
        dNdy_Xi[Itr]      ->Write();
        dNdy_Xib[Itr]     ->Write();
    }
    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
        for (int Jtr = 0;Jtr < HyperonPhaseNum;Jtr++) {

            k_Kp_Omega   [Itr][Jtr]->Write();
            k_Kp_Omegab  [Itr][Jtr]->Write();
            k_Km_Omega   [Itr][Jtr]->Write();
            k_Km_Omegab  [Itr][Jtr]->Write();
            k_Kp_Xi      [Itr][Jtr]->Write();
            k_Kp_Xib     [Itr][Jtr]->Write();
            k_Km_Xi      [Itr][Jtr]->Write();
            k_Km_Xib     [Itr][Jtr]->Write();
            k_Kp_Lambda  [Itr][Jtr]->Write();
            k_Kp_Lambdab [Itr][Jtr]->Write();
            k_Km_Lambda  [Itr][Jtr]->Write();
            k_Km_Lambdab [Itr][Jtr]->Write();
            k_Pip_Omega  [Itr][Jtr]->Write();
            k_Pip_Omegab [Itr][Jtr]->Write();
            k_Pim_Omega  [Itr][Jtr]->Write();
            k_Pim_Omegab [Itr][Jtr]->Write();
            k_Pip_Xi     [Itr][Jtr]->Write();
            k_Pip_Xib    [Itr][Jtr]->Write();
            k_Pim_Xi     [Itr][Jtr]->Write();
            k_Pim_Xib    [Itr][Jtr]->Write();
            k_Pip_Lambda [Itr][Jtr]->Write();
            k_Pip_Lambdab[Itr][Jtr]->Write();
            k_Pim_Lambda [Itr][Jtr]->Write();
            k_Pim_Lambdab[Itr][Jtr]->Write();
        }
    }

    file->Write();

    return;
}
