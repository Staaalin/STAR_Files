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
#include <map>
#include <stdio.h>
using namespace std;

void print(std::vector<int> Temp);
void print(std::vector<float> Temp);
float CenCorr(float Vz,TString Name);
Double_t massList(int PID,TString Name);
Double_t massListSigma(int PID,TString Name);
bool IfInVector(int Num , std::vector<int> V);
bool IfCommonElement(std::vector<int> A , std::vector<int> B);
void DltElement(std::vector<int> &V , int ID);
std::vector<int> GetDaughterPDGLit(int ID);
std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize);

struct Particle{
    int PDG;
    unsigned int EvtID;

    float Px;
    float Py;
    float Pz;

    float Rap;
}

struct ParticlePool{
    unsigned int EvtID;
    vector<Particle> ListA;
    vector<Particle> ListB;
}

// const int CentralityBin[] = {0 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 60 , 70 , 80};// %
const int CentralityBin[] = {0 , 10 , 30 , 50 , 100};// %
const float PVzBin[] = {-45.0 , -35.0 , -25.0 , -15.0 , -5.0 , 5.0 , 15.0 , 25.0 , 35.0 , 45.0 , 55.0}; // Primary Vertex Z (cm) d+Au@200 GeV RUN 21 : -45 ~ 55 cm
const float yBin[]  = {-1.0 , 0.0 , 1.0}; // B_y
const float AyCut[] = {-1.0 , 1.0}; // A_y
int FeedDown[] = { 3334 , -3334};
// int FeedDown[] = {0};
const float EtaCut[] = {-1.0 , 1.0}; // EtaCut for both A and B
const TString DataName = "dAu_200_21";

const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1
const Int_t PVzBinNum = sizeof(PVzBin)/sizeof(PVzBin[0]) - 1; // -1
const Int_t yBinNum = sizeof(yBin)/sizeof(yBin[0]) - 1; // -1
const Int_t FeedDownNum = sizeof(FeedDown)/sizeof(FeedDown[0]);

///////////       Main       ///////////
void NR(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,TString OutMidName,
        int A_PDG,int B_PDG,int Mode = 0)
{
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) 

        std::vector<int>     *PDG                = nullptr;
        std::vector<Float_t> *mix_px             = nullptr;
        std::vector<Float_t> *mix_py             = nullptr;
        std::vector<Float_t> *mix_pz             = nullptr;
        std::vector<Float_t> *QA_eta             = nullptr;
        std::vector<Float_t> *dEdx               = nullptr;
        std::vector<Float_t> *m2                 = nullptr;
        std::vector<Float_t> *dcatopv            = nullptr;
        std::vector<Float_t> *nSigmaProton       = nullptr;
        std::vector<Float_t> *nSigmaPion         = nullptr;
        std::vector<Float_t> *nSigmaKaon         = nullptr;
        std::vector<Float_t> *InvariantMass      = nullptr;
        std::vector<Float_t> *Decay_Length       = nullptr;
        std::vector<Float_t> *Chi2               = nullptr;
        std::vector<int>     *ParentList         = nullptr;
        std::vector<int>     *ParentSta          = nullptr;
        std::vector<int>     *ParentEnd          = nullptr;

        TBranch *bPDG                            = nullptr;
        TBranch *bmix_px                         = nullptr;
        TBranch *bmix_py                         = nullptr;
        TBranch *bmix_pz                         = nullptr;
        TBranch *bQA_eta                         = nullptr;
        TBranch *bdEdx                           = nullptr;
        TBranch *bm2                             = nullptr;
        TBranch *bdcatopv                        = nullptr;
        TBranch *bnSigmaProton                   = nullptr;
        TBranch *bnSigmaPion                     = nullptr;
        TBranch *bnSigmaKaon                     = nullptr;
        TBranch *bInvariantMass                  = nullptr;
        TBranch *bDecay_Length                   = nullptr;
        TBranch *bChi2                           = nullptr;
        TBranch *bParentList                     = nullptr;
        TBranch *bParentSta                      = nullptr;
        TBranch *bParentEnd                      = nullptr;
    
    #else
        #if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)

            std::vector<int>     *PDG                = NULL;
            std::vector<Float_t> *mix_px             = NULL;
            std::vector<Float_t> *mix_py             = NULL;
            std::vector<Float_t> *mix_pz             = NULL;
            std::vector<Float_t> *QA_eta             = NULL;
            std::vector<Float_t> *dEdx               = NULL;
            std::vector<Float_t> *m2                 = NULL;
            std::vector<Float_t> *dcatopv            = NULL;
            std::vector<Float_t> *nSigmaProton       = NULL;
            std::vector<Float_t> *nSigmaPion         = NULL;
            std::vector<Float_t> *nSigmaKaon         = NULL;
            std::vector<Float_t> *InvariantMass      = NULL;
            std::vector<Float_t> *Decay_Length       = NULL;
            std::vector<Float_t> *Chi2               = NULL;
            std::vector<int>     *ParentList         = NULL;
            std::vector<int>     *ParentSta          = NULL;
            std::vector<int>     *ParentEnd          = NULL;

            TBranch *bPDG                            = NULL;
            TBranch *bmix_px                         = NULL;
            TBranch *bmix_py                         = NULL;
            TBranch *bmix_pz                         = NULL;
            TBranch *bQA_eta                         = NULL;
            TBranch *bdEdx                           = NULL;
            TBranch *bm2                             = NULL;
            TBranch *bdcatopv                        = NULL;
            TBranch *bnSigmaProton                   = NULL;
            TBranch *bnSigmaPion                     = NULL;
            TBranch *bnSigmaKaon                     = NULL;
            TBranch *bInvariantMass                  = NULL;
            TBranch *bDecay_Length                   = NULL;
            TBranch *bChi2                           = NULL;
            TBranch *bParentList                     = NULL;
            TBranch *bParentSta                      = NULL;
            TBranch *bParentEnd                      = NULL;

        #else
    
            std::vector<int>     *PDG                = 0;
            std::vector<Float_t> *mix_px             = 0;
            std::vector<Float_t> *mix_py             = 0;
            std::vector<Float_t> *mix_pz             = 0;
            std::vector<Float_t> *QA_eta             = 0;
            std::vector<Float_t> *dEdx               = 0;
            std::vector<Float_t> *m2                 = 0;
            std::vector<Float_t> *dcatopv            = 0;
            std::vector<Float_t> *nSigmaProton       = 0;
            std::vector<Float_t> *nSigmaPion         = 0;
            std::vector<Float_t> *nSigmaKaon         = 0;
            std::vector<Float_t> *InvariantMass      = 0;
            std::vector<Float_t> *Decay_Length       = 0;
            std::vector<Float_t> *Chi2               = 0;
            std::vector<int>     *ParentList         = 0;
            std::vector<int>     *ParentSta          = 0;
            std::vector<int>     *ParentEnd          = 0;
    
            TBranch *bPDG                            = 0;
            TBranch *bmix_px                         = 0;
            TBranch *bmix_py                         = 0;
            TBranch *bmix_pz                         = 0;
            TBranch *bQA_eta                         = 0;
            TBranch *bdEdx                           = 0;
            TBranch *bm2                             = 0;
            TBranch *bdcatopv                        = 0;
            TBranch *bnSigmaProton                   = 0;
            TBranch *bnSigmaPion                     = 0;
            TBranch *bnSigmaKaon                     = 0;
            TBranch *bInvariantMass                  = 0;
            TBranch *bDecay_Length                   = 0;
            TBranch *bChi2                           = 0;
            TBranch *bParentList                     = 0;
            TBranch *bParentSta                      = 0;
            TBranch *bParentEnd                      = 0;

        #endif
    #endif

    std::vector<int> Temp;
    std::vector<float> CMass , CMassSigma;
    bool IfRecord = true , IfRemoveFeedPair = false;
    float BMass = massList(B_PDG)           , AMass = massList(A_PDG);
    float BMassSigma = massListSigma(B_PDG) , AMassSigma = massListSigma(A_PDG);

    std::vector<int> NchList = GetNchList(CentralityBin , CentralityBinNum+1);     // centrality
    cout<<"NchList = ";
    print(NchList);
    cout<<" "<<endl;

    // Hist Parameter
    int kStarBinNum = 400;
    float kStarSta = 0 , kStarEnd = 8;
    
    int dRapBinNum = 300;
    float dRapSta = -3 , dRapEnd = 3;
    
    int dPtBinNum = 200;
    float dPtSta = 0 , dPtEnd = 10;
    
    int MBinNum = 500 , MBinPar = 50;
    float MSta = floor((AMass + BMass)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;

    float NNch , Eta;

    //// Define Histgrams
    //                                     centrality          B_y        PVz
    TH1D* H_Kstar                         [15]                 [15]       [15] ;
    TH1D* H_Mix_Kstar                     [15]                 [15]       [15] ;
    TH1D* H_Res_Kstar                     [15]                 [15]       [15] ;
    TH1D* H_ALL_Kstar                                          [15]            ;
    TH1D* H_ALL_Mix_Kstar                                      [15]            ;
    TH1D* H_ALL_Res_Kstar                                      [15]            ;
    TH1D* H_dRap                          [15]                 [15]       [15] ;
    TH1D* H_Mix_dRap                      [15]                 [15]       [15] ;
    TH1D* H_Res_dRap                      [15]                 [15]       [15] ;
    TH1D* H_ALL_dRap                                           [15]            ;
    TH1D* H_ALL_Mix_dRap                                       [15]            ;
    TH1D* H_ALL_Res_dRap                                       [15]            ;
    TH1D* H_dPt                           [15]                 [15]       [15] ;
    TH1D* H_Mix_dPt                       [15]                 [15]       [15] ;
    TH1D* H_Res_dPt                       [15]                 [15]       [15] ;
    TH1D* H_ALL_dPt                                            [15]            ;
    TH1D* H_ALL_Mix_dPt                                        [15]            ;
    TH1D* H_ALL_Res_dPt                                        [15]            ;
    TH1D* H_Mass                          [15]                 [15]       [15] ;
    TH1D* H_Mix_Mass                      [15]                 [15]       [15] ;
    TH1D* H_Res_Mass                      [15]                 [15]       [15] ;
    TH1D* H_ALL_Mass                                           [15]            ;
    TH1D* H_ALL_Mix_Mass                                       [15]            ;
    TH1D* H_ALL_Res_Mass                                       [15]            ;
    TH1D* H_A_Num                         [15]                 [15]       [15] ;
    TH1D* H_Mix_A_Num                     [15]                 [15]       [15] ;
    TH1D* H_Res_A_Num                     [15]                 [15]       [15] ;
    TH1D* H_ALL_A_Num                                          [15]            ;
    TH1D* H_ALL_Mix_A_Num                                      [15]            ;
    TH1D* H_ALL_Res_A_Num                                      [15]            ;
    TH1D* H_B_Num                         [15]                 [15]       [15] ;
    TH1D* H_Mix_B_Num                     [15]                 [15]       [15] ;
    TH1D* H_Res_B_Num                     [15]                 [15]       [15] ;
    TH1D* H_ALL_B_Num                                          [15]            ;
    TH1D* H_ALL_Mix_B_Num                                      [15]            ;
    TH1D* H_ALL_Res_B_Num                                      [15]            ;
    TH1D* H_Event_Num                     [15]                 [15]       [15] ;
    TH1D* H_Res_Event_Num                 [15]                 [15]       [15] ;
    TH1D* H_ALL_Event_Num                                      [15]            ;
    TH1D* H_ALL_Res_Event_Num                                  [15]            ;

    // A/B d (net)Num / d Dy
    TH1D* H_A_Num_Dy                      [15]                 [15]       [15] ;
    TH1D* H_Mix_A_Num_Dy                  [15]                 [15]       [15] ;
    TH1D* H_Res_A_Num_Dy                  [15]                 [15]       [15] ;
    TH1D* H_ALL_A_Num_Dy                                       [15]            ;
    TH1D* H_ALL_Mix_A_Num_Dy                                   [15]            ;
    TH1D* H_ALL_Res_A_Num_Dy                                   [15]            ;
    TH1D* H_B_Num_Dy                      [15]                 [15]       [15] ;
    TH1D* H_Mix_B_Num_Dy                  [15]                 [15]       [15] ;
    TH1D* H_Res_B_Num_Dy                  [15]                 [15]       [15] ;
    TH1D* H_ALL_B_Num_Dy                                       [15]            ;
    TH1D* H_ALL_Mix_B_Num_Dy                                   [15]            ;
    TH1D* H_ALL_Res_B_Num_Dy                                   [15]            ;

    TString HistNameTemp1 , HistNameTemp2;
    for (int i=0;i<CentralityBinNum;i++){
        for (int k=0;k<PVzBinNum;k++){
            for (int j=0;j<yBinNum;j++){
                HistNameTemp1 = "H_";HistNameTemp1+="Kstar_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Kstar, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Kstar    [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,kStarBinNum,kStarSta,kStarEnd);

                HistNameTemp1 = "H_M_";HistNameTemp1+="Kstar_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Mixed Kstar, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Mix_Kstar[i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,kStarBinNum,kStarSta,kStarEnd);

                HistNameTemp1 = "H_R_";HistNameTemp1+="Kstar_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Resed Kstar, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Res_Kstar[i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,kStarBinNum,kStarSta,kStarEnd);
                
                HistNameTemp1 = "H_";HistNameTemp1+="dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_dRap    [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum,dRapSta,dRapEnd);

                HistNameTemp1 = "H_M_";HistNameTemp1+="dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Mixed dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Mix_dRap[i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum,dRapSta,dRapEnd);

                HistNameTemp1 = "H_R_";HistNameTemp1+="dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Resed dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Res_dRap[i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum,dRapSta,dRapEnd);
            }
        }
    }

    return;
}

void print(std::vector<int> Temp)
{
	cout<<"{";
    for (int i = 0;i<Temp.size();i++){
		cout<<" "<<Temp.at(i);
		if (i != (Temp.size() - 1)) cout<<" ,"; 
	}
	cout<<" }"<<endl;
    return ;
}

void print(std::vector<float> Temp)
{
	cout<<"{";
    for (int i = 0;i<Temp.size();i++){
		cout<<" "<<Temp.at(i);
		if (i != (Temp.size() - 1)) cout<<" ,"; 
	}
	cout<<" }"<<endl;
    return ;
}

float CenCorr(float Vz,TString Name)
{
    if (Name == "dAu_200_21") {// data from https://drupal.star.bnl.gov/STAR/system/files/pwg5.pdf
        if      (Vz < -50.0) {
            return 1.0;
        }
        else if (Vz < -40.0) {
            return 1.13833390;
        }
        else if (Vz < -30.0) {
            return 1.06240111;
        }
        else if (Vz < -20.0) {
            return 1.02187042;
        }
        else if (Vz < -10.0) {
            return 1.00557849;
        }
        else if (Vz < 0.0) {
            return 0.99907267;
        }
        else if (Vz < 10.0) {
            return 0.99731279;
        }
        else if (Vz < 20.0) {
            return 0.99807879;
        }
        else if (Vz < 30.0) {
            return 0.99894410;
        }
        else if (Vz < 40.0) {
            return 0.99543646;
        }
        else if (Vz < 50.0) {
            return 0.99522446;
        }
        else                {
            return 1.0;
        }
    }
    return 1.0;
}

Double_t massList(int PID,TString Name)
{
    Double_t Result;
    if (Name == "dAu_200_21"){
        switch (PID)
        {
            case 321 :
                Result = 0.493677;
                break;
            case -321 :
                Result = 0.493677;
                break;
            case 310 :
                Result = 0.49794;
                break;
            case 211 :
                Result = 0.13957;
                break;
            case -211 :
                Result = 0.13957;
                break;
            case 1003314 :// XiRPdgMass
                Result = 1.6725;
                break;
            case -1003314 :// XiRPdgMass
                Result = 1.6727;
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

Double_t massListSigma(int PID,TString Name)
{
    Double_t Result;
    if (Name == "dAu_200_21"){
        switch (PID)
        {
            case 3334 :// OmegaFitMass
                Result = 0.0029;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 0.0024;
                break;
            case 1003314 :// XiRPdgMass
                Result = 0.0029;
                break;
            case -1003314 :// XiRPdgMass
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
                Result = 100;
        }
    }
    if (Name == "dAu_62_16"){// tbd, used as dAu@200R21
        switch (PID)
        {
            case 3334 :// OmegaFitMass
                Result = 0.0029;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 0.0024;
                break;
            case 1003314 :// XiRPdgMass
                Result = 0.0029;
                break;
            case -1003314 :// XiRPdgMass
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
                Result = 100;
        }
    }
    return Result;
}

bool IfInVector(int Num , std::vector<int> V)
{
    for (int i=0;i<V.size();i++) {
        if (Num == V.at(i)){
            return true;
        }
    }
    return false;
}

bool IfCommonElement(std::vector<int> A , std::vector<int> B)
{
    for (int i=0;i<A.size();i++){
        for (int j=0;j<B.size();j++){
            if (A.at(i) == B.at(j)) return true;
        }
    }
    return false;
}

void DltElement(std::vector<int> &V , int ID)
{
    std::vector<int> V_T;V_T.clear();
    for (int i=0;i<V.size();i++){
        V_T.push_back(V.at(i));
    }
    V.clear();
    for (int i=0;i<V_T.size();i++){
        if (i == ID) continue;
        V.push_back(V_T.at(i));
    }
    return;
}

std::vector<int> GetDaughterPDGLit(int ID)
{
    std::vector<int> V_T;V_T.clear();
    switch (ID)
    {
        case 3334 :// Omega
            V_T.push_back(-321);
            V_T.push_back(3122);
            return V_T;
        case -3334 :// OmegaBar
            V_T.push_back(321);
            V_T.push_back(-3122);
            return V_T;
        case 1003314 :// XiR
            V_T.push_back(-321);
            V_T.push_back(3122);
            return V_T;
        case -1003314 :// XiRBar
            V_T.push_back(321);
            V_T.push_back(-3122);
            return V_T;
        case 3312 :// Xi
            V_T.push_back(-211);
            V_T.push_back(3122);
            return V_T;
        case -3312 :// XiBar
            V_T.push_back(211);
            V_T.push_back(-3122);
            return V_T;
        case 3122 :// Lambda
            V_T.push_back(-211);
            V_T.push_back(2212);
            return V_T;
        case -3122 :// LambdaBar
            V_T.push_back(211);
            V_T.push_back(-2212);
            return V_T;
        default :
            return V_T;
    }
}

std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize)
{
    //This is 329
    std::vector<int> Result;Result.clear();
    // int CentralityListSize = sizeof(CentralityList)/sizeof(CentralityList[0]);
    if (DataName == "dAu_200_21") {
        // data from https://drupal.star.bnl.gov/STAR/system/files/pwg5.pdf
        int NchTable[21] = { 10000 , 55 , 47 , 42 , 38 , 35 , 32 , 29 , 26 , 24 , 21 , 19 , 17 , 15 , 13 , 11 , 9 , 7 , 6 , 4 ,  0};
        int CenTable[21] = {     0 ,  5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 55 , 60 , 65 , 70 , 75 ,80 ,85 ,90 ,95 ,100};
        for (int i=0;i<CentralityListSize;i++) {
            for (int j=0;j<21;j++){
                if (CenTable[j] == CentralityList[i]) {
                    Result.push_back(NchTable[j]);
                    break;
                }
            }
        }
    }
    return Result;
}