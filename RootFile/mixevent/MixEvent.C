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
#include <ctime>
// #include "Math/Vector4D.h"
// #include <Math/PtEtaPhiM4D.h>
// #include <Math/Boost.h>
using namespace std;

// #define DataName           "pAu_200_15"
// #define DataName           "AuAu_27_18"
// #define DataName           "dAu_200_16"
#define DataName           "dAu_200_21"
// #define DataName           "dAu_62_16"
// #define DataName           "dAu_39_16"
// #define DataName           "dAu_20_16"
// #define DataName           "pp_200_15"
// #define DataName           "OO_200_21"

#define SpecialMode false

#define Pi 3.1415926535898

// const int CentralityBin[] = {0 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 60 , 70 , 80};// %
const int CentralityBin[] = {0 , 10 , 30 , 50 , 100};// %
const float PVzBin[] = {-45.0 , -35.0 , -25.0 , -15.0 , -5.0 , 5.0 , 15.0 , 25.0 , 35.0 , 45.0 , 55.0}; // Primary Vertex Z (cm) d+Au@200 GeV RUN 21 : -45 ~ 55 cm
const float yBin[]  = {-1.0 , 0.0 , 1.0}; // B_y
const float AyCut[] = {-1.0 , 1.0}; // A_y
int FeedDown[] = { 3334 , -3334};
// int FeedDown[] = {0};
const float EtaCut[] = {-1.0 , 1.0}; // EtaCut for both A and B

const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1
const Int_t PVzBinNum = sizeof(PVzBin)/sizeof(PVzBin[0]) - 1; // -1
const Int_t yBinNum = sizeof(yBin)/sizeof(yBin[0]) - 1; // -1
const Int_t FeedDownNum = sizeof(FeedDown)/sizeof(FeedDown[0]);

#define A_Num_Per_Event 5
#define B_Num_Per_Event 5
#define HowMuchEventMixing 10

TString KindBin[] = {"Mid","Sid"}
#define KindNum 2
TString PatternBin[] = {"AMBM","AMBS","ASBM"};
#define Pattern 3 // 0:A middle B middle , 1:A middle B sideband , 2:A sideband B middle
// Pattern应当大于KindNum


// 计算质心系速度
std::vector<float> calculateBeta(std::vector<float>& p1, std::vector<float>& p2) {
    std::vector<float> Result;
    double totalPx = p1[0] + p2[0];
    double totalPy = p1[1] + p2[1];
    double totalPz = p1[2] + p2[2];
    double totalE = p1[3] + p2[3];
    Result.push_back(totalPx / totalE);
    Result.push_back(totalPy / totalE);
    Result.push_back(totalPz / totalE);
    return Result;
}

// 计算 Lorentz boost
std::vector<float> boost(std::vector<float>& p, std::vector<float>& beta) {
    double beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    double gamma = 1.0 / std::sqrt(1.0 - beta2);

    double bp = beta[0]*p[0] + beta[1]*p[1] + beta[2]*p[2];
    double gamma2 = (beta2 > 0) ? (gamma - 1.0) / beta2 : 0.0;

    std::vector<float> boosted;
    boosted.push_back(p[0] + gamma2 * bp * beta[0] + gamma * beta[0] * p[3]);
    boosted.push_back(p[1] + gamma2 * bp * beta[1] + gamma * beta[1] * p[3]);
    boosted.push_back(p[2] + gamma2 * bp * beta[2] + gamma * beta[2] * p[3]);
    boosted.push_back(gamma * (p[3] + bp));

    return boosted;
}

float Rho(std::vector<float>& p) {
    float Result = 0;
    for (int i = 0;i < 3;i++){
        Result += p.at(i)*p.at(i);
    }
    Result = pow(Result,0.5);
    return Result;
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

float CenCorr(float Vz)
{
    if (DataName == "dAu_200_21") {// data from https://drupal.star.bnl.gov/STAR/system/files/pwg5.pdf
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

Double_t massList(int PID)
{
    Double_t Result;
    if (DataName == "dAu_200_21"){
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

Double_t massListSigma(int PID)
{
    Double_t Result;
    if (DataName == "dAu_200_21"){
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
    if (DataName == "dAu_62_16"){// tbd, used as dAu@200R21
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

// Get Mass in rest frame
float GetPairMass(float p1x,float p1y,float p1z,float m1,float p2x,float p2y,float p2z,float m2) {
    float E1 = pow(p1x*p1x+p1y*p1y+p1z*p1z+m1*m1,0.5);
    float E2 = pow(p2x*p2x+p2y*p2y+p2z*p2z+m2*m2,0.5);
    float Tot_E = E1+E2;
    float beta[3] = { -(p1x+p2x)/Tot_E , -(p1y+p2y)/Tot_E , -(p1z+p2z)/Tot_E };
    float beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    float gamma = 1.0 / std::sqrt(1.0 - beta2);
    float gamma2 = (beta2 > 0) ? (gamma - 1.0) / beta2 : 0.0;

    float bp1 = beta[0]*p1x + beta[1]*p1y + beta[2]*p1z;
    float bp2 = beta[0]*p2x + beta[1]*p2y + beta[2]*p2z;
    return gamma * (E1 + bp1 + E2 + bp2);
}

float* GetPairMassAndKstar(float p1x , float p1y , float p1z , float p2x , float p2y , float p2z , float AMass , float BMass) {
    float E1 = pow(p1x*p1x+p1y*p1y+p1z*p1z+AMass*AMass,0.5);
    float E2 = pow(p2x*p2x+p2y*p2y+p2z*p2z+BMass*BMass,0.5);
    float Tot_E = E1+E2;
    float beta[3] = { -(p1x+p2x)/Tot_E , -(p1y+p2y)/Tot_E , -(p1z+p2z)/Tot_E };
    float beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    float gamma = 1.0 / std::sqrt(1.0 - beta2);
    float gamma2 = (beta2 > 0) ? (gamma - 1.0) / beta2 : 0.0;

    float bp1 = beta[0]*p1x + beta[1]*p1y + beta[2]*p1z;
    float bp2 = beta[0]*p2x + beta[1]*p2y + beta[2]*p2z;

    // float New_Px = p1x + gamma2 * bp1 * beta[0] + gamma * beta[0] * E1;
    // float New_Py = p1y + gamma2 * bp1 * beta[1] + gamma * beta[1] * E1;
    // float New_Pz = p1z + gamma2 * bp1 * beta[2] + gamma * beta[2] * E1;
    float New_Px = (p1x - p2x) + gamma2 * (bp1-bp2) * beta[0] + gamma * beta[0] * (E1-E2);
    float New_Py = (p1y - p2y) + gamma2 * (bp1-bp2) * beta[1] + gamma * beta[1] * (E1-E2);
    float New_Pz = (p1z - p2z) + gamma2 * (bp1-bp2) * beta[2] + gamma * beta[2] * (E1-E2);

    float* MassAndKstar = new float[2];
    MassAndKstar[0] = (gamma * (E1 + bp1 + E2 + bp2));
    MassAndKstar[1] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    return MassAndKstar;
}

void MixEvent(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,TString OutMidName,
              int A_PDG,int B_PDG,int Mode = 0) // Mode = 0: PDGMult 为vector长度
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

    double kstar, rap;
    TVector3 BetaTemp;
    // ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p1 , p2 , p3 , p4 , p5;
    TLorentzVector p1 , p2 , p3;
    TVector3 BV;
    float tEnergy , APx , APy , APz , BPx , BPy , BPz , PairMass , KS , Pt;
    int i , j , k , l , m , n , Aid , Bid , Cid , Aindex , Bindex , RapIndex;
    int A_Kid , B_Kid , Mix_A_Size , Mix_B_Size , A_EID , AidN , BidN;
    std::vector<int> Temp;
    std::vector<float> CMass , CMassSigma;
    bool IfRecord = true , IfRemoveFeedPair = false;
    float BMass = massList(B_PDG)           , AMass = massList(A_PDG);
    float BMassSigma = massListSigma(B_PDG) , AMassSigma = massListSigma(A_PDG);
    // float MassAndKstar[2];

    std::vector<int> NchList = GetNchList(CentralityBin , CentralityBinNum+1);     // centrality
    cout<<"NchList = ";
    print(NchList);
    cout<<" "<<endl;
    int                                   A_Num                       ;
    float                                 A_Px       [A_Num_Per_Event];
    float                                 A_Py       [A_Num_Per_Event];
    float                                 A_Pz       [A_Num_Per_Event];
    int                                   A_TreID    [A_Num_Per_Event];
    std::vector<std::vector<int> >        A_ParID                     ;
    int                                   A_Kind     [A_Num_Per_Event];
    float                                 A_Rap      [A_Num_Per_Event];
    bool                                  A_IfRecord [A_Num_Per_Event];
    int                                   B_Num                       ;
    float                                 B_Px       [B_Num_Per_Event];
    float                                 B_Py       [B_Num_Per_Event];
    float                                 B_Pz       [B_Num_Per_Event];
    int                                   B_TreID    [B_Num_Per_Event];
    std::vector<std::vector<int> >        B_ParID                     ;
    int                                   B_Kind     [B_Num_Per_Event];
    float                                 B_Rap      [B_Num_Per_Event];
    bool                                  B_IfRecord [B_Num_Per_Event];
    std::vector<std::vector<int> >        C_ParID                     ; // 用于存储Residal Effect
    // used as array
    //                                        centrality          B_y        PVz
    unsigned short int Mix_A_Index_T = 0;
    unsigned short int Mix_A_Index        [15]                 [15]       [15]         [2] [2] ;
    float              Mix_A_Px           [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*A_Num_Per_Event];
    float              Mix_A_Py           [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*A_Num_Per_Event];
    float              Mix_A_Pz           [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*A_Num_Per_Event];
    int                Mix_A_TreID        [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*A_Num_Per_Event];
    unsigned int       Mix_A_EvtID        [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*A_Num_Per_Event];
    float              Mix_A_Rap          [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*A_Num_Per_Event];
    bool               Mix_A_IfMadePair   [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*A_Num_Per_Event];
    unsigned short int Mix_A_ID_Index                          [15]                    [2] [2] ;
    std::vector<int>   Mix_A_ID                                [15]                    [2] [2] ;
    unsigned short int Mix_B_Index_T = 0;
    unsigned short int Mix_B_Index        [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*B_Num_Per_Event];
    float              Mix_B_Px           [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*B_Num_Per_Event];
    float              Mix_B_Py           [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*B_Num_Per_Event];
    float              Mix_B_Pz           [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*B_Num_Per_Event];
    int                Mix_B_TreID        [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*B_Num_Per_Event];
    unsigned int       Mix_B_EvtID        [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*B_Num_Per_Event];
    float              Mix_B_Rap          [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*B_Num_Per_Event];
    bool               Mix_B_IfMadePair   [15]                 [15]       [15]         [2] [2] [(HowMuchEventMixing+1)*B_Num_Per_Event];
    unsigned short int Mix_B_ID_Index                          [15]                    [2] [2] ;
    std::vector<int>   Mix_B_ID                                [15]                    [2] [2] ;
    int                Mix_event_Num      [15]                 [15]       [15]         [2] [2] ;
    int                Mix_event_Num_SUM  [15]                 [15]       [15]         [2] [2] ;
    //        
    TH1D* H_Kstar                         [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Mix_Kstar                     [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_dRap                          [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Mix_dRap                      [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_dPt                           [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Mix_dPt                       [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Mass                          [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Mix_Mass                      [15]                 [15]       [15]         [2] [2] ;

    TH1D* H_Res_Kstar                     [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Res_dRap                      [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Res_dPt                       [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Res_Mass                      [15]                 [15]       [15]         [2] [2] ;

    TH1D* H_ALL_Kstar                                          [15]                    [2] [2] ;
    TH1D* H_ALL_Mix_Kstar                                      [15]                    [2] [2] ;
    TH1D* H_ALL_Res_Kstar                                      [15]                    [2] [2] ;
    TH1D* H_ALL_dPt                                            [15]                    [2] [2] ;
    TH1D* H_ALL_Mix_dPt                                        [15]                    [2] [2] ;
    TH1D* H_ALL_Res_dPt                                        [15]                    [2] [2] ;
    TH1D* H_ALL_dRap                                           [15]                    [2] [2] ;
    TH1D* H_ALL_Mix_dRap                                       [15]                    [2] [2] ;
    TH1D* H_ALL_Res_dRap                                       [15]                    [2] [2] ;
    TH1D* H_ALL_Mass                                                                   [2] [2] ;
    TH1D* H_ALL_Mix_Mass                                                               [2] [2] ;
    TH1D* H_A_Num                         [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_B_Num                         [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Res_A_Num                     [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_Res_B_Num                     [15]                 [15]       [15]         [2] [2] ;
    TH1D* H_ALL_A_Num                                          [15]                    [2] [2] ;
    TH1D* H_ALL_B_Num                                          [15]                    [2] [2] ;
    TH1D* H_ALL_Res_A_Num                                      [15]                    [2] [2] ;
    TH1D* H_ALL_Res_B_Num                                      [15]                    [2] [2] ;
    TProfile* H_Event_Num                 [15]                 [15]       [15]         [2] [2] ;
    TProfile* H_Res_Event_Num             [15]                 [15]       [15]         [2] [2] ;
    TProfile* H_ALL_Event_Num                                  [15]                    [2] [2] ;
    TProfile* H_ALL_Res_Event_Num                              [15]                    [2] [2] ;

    // Rotation
    TH1D* R_S_Kstar                       [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_S_dRap                        [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_S_dPt                         [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_S_Mass                        [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_M_Kstar                       [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_M_dRap                        [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_M_dPt                         [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_M_Mass                        [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_A_Num                         [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_B_Num                         [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_S_Kstar                   [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_S_dRap                    [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_S_dPt                     [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_S_Mass                    [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_M_Kstar                   [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_M_dRap                    [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_M_dPt                     [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_M_Mass                    [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_A_Num                     [15]                 [15]       [15]         [2] [2] ;
    TH1D* R_ALL_B_Num                     [15]                 [15]       [15]         [2] [2] ;

    // Store in test
    TH2F* H_ALL_Kstar_dRap                                     [15]                    [2] [2] ;
    TH2F* H_ALL_Mix_Kstar_dRap                                 [15]                    [2] [2] ;

    int EventPatternMatch                 [15]                 [15]       [15]         [2] [2] ;
    // Used for testing
    int TestSum = 0;
    bool IfFoundOmega = false;

    int kStarBinNum = 400;
    float kStarSta = 0 , kStarEnd = 8;
    
    int dRapBinNum = 300;
    float dRapSta = -3 , dRapEnd = 3;
    
    int dPtBinNum = 200;
    float dPtSta = 0 , dPtEnd = 10;
    
    int MBinNum = 500 , MBinPar = 50;
    float MSta = floor((AMass + BMass)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;

    float NNch , Eta;
    
    TString HistNameI  , HistNameJ  , HistNameK  , HistNameL  , HistNameM ;
    TString HistNameIs , HistNameJs , HistNameKs , HistNameLs , HistNameMr;
    TString HistNameIr , HistNameJr , HistNameKr , HistNameLr , HistNameMs;

    for (i = 0;i < FeedDownNum;i++){
        if (abs(FeedDown[i]) == A_PDG) {
            FeedDown[i] = 0;
            CMass.push_back(-100);
            CMassSigma.push_back(-1);
            continue;
        }
        if (abs(FeedDown[i]) == B_PDG) {
            FeedDown[i] = 0;
            CMass.push_back(-100);
            CMassSigma.push_back(-1);
            continue;
        }
        CMass.push_back(massList(FeedDown[i]));
        CMassSigma.push_back(massListSigma(FeedDown[i]));
    }
    cout<<"CMass = ";print(CMass);
    cout<<"CMassSigma = ";print(CMassSigma);

    for (i=0;i<FeedDownNum;i++) {
        if ( IfInVector(A_PDG , GetDaughterPDGLit(FeedDown[i])) && IfInVector(B_PDG , GetDaughterPDGLit(FeedDown[i])) ) IfRemoveFeedPair = true;
    }

    for (i=0;i<CentralityBinNum;i++){
        for (l=0;l<Pattern;l++){
            for (k=0;k<PVzBinNum;k++){
                for (j=0;j<yBinNum;j++){
                    TString HistName1 = "H_";
                    TString HistName2 = "Cen: [";
                    TString HistName3 = "All Cen , ";
                    TString HistName4 = "H_ALL";
                    TString HistName5 = "H_ALL_Mix";
                    TString HistName6 = "H_ALL_Res";
                    HistName1 += i;HistName1 += "_";
                    HistName2 += CentralityBin[i];HistName2 += "% , ";
                    HistName2 += CentralityBin[i+1];HistName2 += "%], ";
                    HistName3 += yBin[j];HistName3 += " < y";HistName3 += B_PDG;HistName3 += " <  ";
                    HistName3 += yBin[j+1];HistName3 += ", All PVz";
                    HistName1 += j;HistName1 += "_";
                    HistName4 += "_";HistName4 += j;HistName4 += "_";
                    HistName5 += "_";HistName5 += j;HistName5 += "_";
                    HistName6 += "_";HistName6 += j;HistName6 += "_";
                    HistName2 += yBin[j];HistName2 += " < y";HistName2 += B_PDG;HistName2 += " <  ";
                    HistName2 += yBin[j+1];HistName2 += ", ";
                    HistName1 += k;
                    HistName2 += PVzBin[k];HistName2 += " < PVz";HistName2 += B_PDG;HistName2 += " <  ";
                    HistName2 += PVzBin[k+1];
                    TString HistName1s = HistName1;
                    TString HistName2s = HistName2;
                    TString HistName4s = HistName4 + PatternBin[l];
                    TString HistName5s = HistName5 + PatternBin[l];
                    TString HistName6s = HistName6 + PatternBin[l];
                    TString HistName4r = HistName4 + PatternBin[l];
                    TString HistName5r = HistName5 + PatternBin[l];
                    TString HistName6r = HistName6 + PatternBin[l];
                    TString HistName4p = HistName4 + PatternBin[l];
                    TString HistName5p = HistName5 + PatternBin[l];
                    TString HistName6p = HistName6 + PatternBin[l];
                    HistName1s += "_S_";
                    HistName1 += "_M_";
                    HistName1 += PatternBin[l];
                    HistName1s += PatternBin[l];
                    HistName2 += ", Mix, ";
                    HistName2 += PatternBin[l];
                    HistName2s += PatternBin[l];
                    // 区分kstar、dRap……
                    HistNameI  = HistName1  + "_kStar";
                    HistNameIs = HistName1s + "_kStar";
                    HistName4s = HistName4s + "_kStar";
                    HistName5s = HistName5s + "_kStar";
                    HistName6s = HistName6s + "_kStar";
                    HistNameIr = HistName1s + "_Res_kStar";
                    HistNameJ  = HistName1  + "_dRap";
                    HistNameJs = HistName1s + "_dRap";
                    HistName4r = HistName4r + "_dRap";
                    HistName5r = HistName5r + "_dRap";
                    HistName6r = HistName6r + "_dRap";
                    HistNameJr = HistName1s + "_Res_dRap";
                    HistNameK  = HistName1  + "_dPt";
                    HistNameKs = HistName1s + "_dPt";
                    HistName4p = HistName4p + "_dPt";
                    HistName5p = HistName5p + "_dPt";
                    HistName6p = HistName6p + "_dPt";
                    HistNameKr = HistName1s + "_Res_dPt";
                    HistNameL  = HistName1  + "_Mass";
                    HistNameLs = HistName1s + "_Mass";
                    HistNameLr = HistName1s + "_Res_Mass";
                    HistNameMs = HistName1s + "_A_Num";
                    HistNameMr = HistName1s + "_B_Num";
                    if (l == 0) { // AMBM
                        H_Kstar                        [i][j][k][0][0] = new TH1D(HistNameIs,HistName2s,kStarBinNum,kStarSta,kStarEnd);
                        H_Mix_Kstar                    [i][j][k][0][0] = new TH1D(HistNameI,HistName2,kStarBinNum,kStarSta,kStarEnd);
                        H_Res_Kstar                    [i][j][k][0][0] = new TH1D(HistNameIr,HistNameIr,kStarBinNum,kStarSta,kStarEnd);
                        H_dRap                         [i][j][k][0][0] = new TH1D(HistNameJs,HistName2s,dRapBinNum,dRapSta,dRapEnd);
                        H_Mix_dRap                     [i][j][k][0][0] = new TH1D(HistNameJ,HistName2,dRapBinNum,dRapSta,dRapEnd);
                        H_Res_dRap                     [i][j][k][0][0] = new TH1D(HistNameJr,HistNameJr,dRapBinNum,dRapSta,dRapEnd);
                        H_dPt                          [i][j][k][0][0] = new TH1D(HistNameKs,HistName2s,dPtBinNum,dPtSta,dPtEnd);
                        H_Mix_dPt                      [i][j][k][0][0] = new TH1D(HistNameK,HistName2,dPtBinNum,dPtSta,dPtEnd);
                        H_Res_dPt                      [i][j][k][0][0] = new TH1D(HistNameKr,HistNameKr,dPtBinNum,dPtSta,dPtEnd);
                        H_Mass                         [i][j][k][0][0] = new TH1D(HistNameLs,HistName2s,MBinNum,MSta,MEnd);
                        H_Mix_Mass                     [i][j][k][0][0] = new TH1D(HistNameL,HistName2,MBinNum,MSta,MEnd);
                        H_Res_Mass                     [i][j][k][0][0] = new TH1D(HistNameLr,HistNameLr,MBinNum,MSta,MEnd);
                        Mix_event_Num                  [i][j][k][0][0] = 0;
                        Mix_event_Num_SUM              [i][j][k][0][0] = 0;
                        HistNameM = HistNameMs + "_AMBM";
                        H_A_Num                        [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistNameMr + "_AMBM";
                        H_B_Num                        [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Res_A_Num";
                        H_Res_A_Num                    [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Res_B_Num";
                        H_Res_B_Num                    [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Event_Num";
                        H_Event_Num                    [i][j][k][0][0] = new TProfile(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Res_Event_Num";
                        H_Res_Event_Num                [i][j][k][0][0] = new TProfile(HistNameM,HistNameM,1,-1,1);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_S_AMBM_Kstar";
                        // R_S_Kstar                      [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_M_AMBM_Kstar";
                        // R_M_Kstar                      [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_S_AMBM_dRap";
                        // R_S_dRap                       [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_M_AMBM_dRap";
                        // R_M_dRap                       [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_S_AMBM_dPt";
                        // R_S_dPt                        [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_M_AMBM_dPt";
                        // R_M_dPt                        [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_A_AMBM_Num";
                        // R_A_Num                        [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_B_AMBM_Num";
                        // R_B_Num                        [i][j][k][0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        if ((i==0)&&(k==0)){
                            H_ALL_Kstar                   [j]   [0][0] = new TH1D(HistName4s,HistName3,kStarBinNum,kStarSta,kStarEnd);
                            H_ALL_Mix_Kstar               [j]   [0][0] = new TH1D(HistName5s,HistName3,kStarBinNum,kStarSta,kStarEnd);
                            H_ALL_Res_Kstar               [j]   [0][0] = new TH1D(HistName6s,HistName3,kStarBinNum,kStarSta,kStarEnd);
                            H_ALL_dRap                    [j]   [0][0] = new TH1D(HistName4r,HistName3,dRapBinNum,dRapSta,dRapEnd);
                            H_ALL_Mix_dRap                [j]   [0][0] = new TH1D(HistName5r,HistName3,dRapBinNum,dRapSta,dRapEnd);
                            H_ALL_Res_dRap                [j]   [0][0] = new TH1D(HistName6r,HistName3,dRapBinNum,dRapSta,dRapEnd);
                            H_ALL_dPt                     [j]   [0][0] = new TH1D(HistName4p,HistName3,dPtBinNum,dPtSta,dPtEnd);
                            H_ALL_Mix_dPt                 [j]   [0][0] = new TH1D(HistName5p,HistName3,dPtBinNum,dPtSta,dPtEnd);
                            H_ALL_Res_dPt                 [j]   [0][0] = new TH1D(HistName6p,HistName3,dPtBinNum,dPtSta,dPtEnd);
                            HistNameM = "H_ALL_";HistNameM += j;HistNameM += "_AMBM_Kstar_dRap";
                            H_ALL_Kstar_dRap              [j]   [0][0] = new TH2F(HistNameM,HistNameM,80,kStarSta,kStarEnd,60,dRapSta,dRapEnd);
                            HistNameM = "H_ALL_Mix_";HistNameM += j;HistNameM += "_AMBM_Kstar_dRap";
                            H_ALL_Mix_Kstar_dRap          [j]   [0][0] = new TH2F(HistNameM,HistNameM,80,kStarSta,kStarEnd,60,dRapSta,dRapEnd);
                            HistNameM = HistName4 + "A_Num_AMBM";
                            H_ALL_A_Num                   [j]   [0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName6 + "A_Num_AMBM";
                            H_ALL_Res_A_Num               [j]   [0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName4 + "B_Num_AMBM";
                            H_ALL_B_Num                   [j]   [0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName6 + "B_Num_AMBM";
                            H_ALL_Res_B_Num               [j]   [0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName4 + "Event_Num_AMBM";
                            H_ALL_Event_Num               [j]   [0][0] = new TProfile(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName6 + "Event_Num_AMBM";
                            H_ALL_Res_Event_Num           [j]   [0][0] = new TProfile(HistNameM,HistNameM,1,-1,1);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_S_AMBM_Kstar";
                            // R_ALL_S_Kstar                 [j]   [0][0] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_M_AMBM_Kstar";
                            // R_ALL_M_Kstar                 [j]   [0][0] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_S_AMBM_dRap";
                            // R_ALL_S_dRap                  [j]   [0][0] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_M_AMBM_dRap";
                            // R_ALL_M_dRap                  [j]   [0][0] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_S_AMBM_dPt";
                            // R_ALL_S_dPt                   [j]   [0][0] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_M_AMBM_dPt";
                            // R_ALL_M_dPt                   [j]   [0][0] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_A_AMBM_Num";
                            // R_ALL_A_Num                   [j]   [0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_B_AMBM_Num";
                            // R_ALL_B_Num                   [j]   [0][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        }
                    }
                    if (l == 1) { // AMBS
                        H_Kstar                        [i][j][k][0][1] = new TH1D(HistNameIs,HistName2s,kStarBinNum,kStarSta,kStarEnd);
                        H_Mix_Kstar                    [i][j][k][0][1] = new TH1D(HistNameI,HistName2,kStarBinNum,kStarSta,kStarEnd);
                        H_Res_Kstar                    [i][j][k][0][1] = new TH1D(HistNameIr,HistNameIr,kStarBinNum,kStarSta,kStarEnd);
                        H_dRap                         [i][j][k][0][1] = new TH1D(HistNameJs,HistName2s,dRapBinNum,dRapSta,dRapEnd);
                        H_Mix_dRap                     [i][j][k][0][1] = new TH1D(HistNameJ,HistName2,dRapBinNum,dRapSta,dRapEnd);
                        H_Res_dRap                     [i][j][k][0][1] = new TH1D(HistNameJr,HistNameJr,dRapBinNum,dRapSta,dRapEnd);
                        H_dPt                          [i][j][k][0][1] = new TH1D(HistNameKs,HistName2s,dPtBinNum,dPtSta,dPtEnd);
                        H_Mix_dPt                      [i][j][k][0][1] = new TH1D(HistNameK,HistName2,dPtBinNum,dPtSta,dPtEnd);
                        H_Res_dPt                      [i][j][k][0][1] = new TH1D(HistNameKr,HistNameKr,dPtBinNum,dPtSta,dPtEnd);
                        H_Mass                         [i][j][k][0][1] = new TH1D(HistNameLs,HistName2s,MBinNum,MSta,MEnd);
                        H_Mix_Mass                     [i][j][k][0][1] = new TH1D(HistNameL,HistName2,MBinNum,MSta,MEnd);
                        H_Res_Mass                     [i][j][k][0][1] = new TH1D(HistNameLr,HistNameLr,MBinNum,MSta,MEnd);
                        Mix_event_Num                  [i][j][k][0][1] = 0;
                        Mix_event_Num_SUM              [i][j][k][0][1] = 0;
                        HistNameM = HistNameMs + "_AMBS";
                        H_A_Num                        [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistNameMr + "_AMBS";
                        H_B_Num                        [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Res_A_Num";
                        H_Res_A_Num                    [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Res_B_Num";
                        H_Res_B_Num                    [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Event_Num";
                        H_Event_Num                    [i][j][k][0][1] = new TProfile(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Res_Event_Num";
                        H_Res_Event_Num                [i][j][k][0][1] = new TProfile(HistNameM,HistNameM,1,-1,1);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_S_AMBS_Kstar";
                        // R_S_Kstar                      [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_M_AMBS_Kstar";
                        // R_M_Kstar                      [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_S_AMBS_dRap";
                        // R_S_dRap                       [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_M_AMBS_dRap";
                        // R_M_dRap                       [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_S_AMBS_dPt";
                        // R_S_dPt                        [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_M_AMBS_dPt";
                        // R_M_dPt                        [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_A_AMBS_Num";
                        // R_A_Num                        [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_B_AMBS_Num";
                        // R_B_Num                        [i][j][k][0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        if ((i==0)&&(k==0)){
                            H_ALL_Kstar                   [j]   [0][1] = new TH1D(HistName4s,HistName3,kStarBinNum,kStarSta,kStarEnd);
                            H_ALL_Mix_Kstar               [j]   [0][1] = new TH1D(HistName5s,HistName3,kStarBinNum,kStarSta,kStarEnd);
                            H_ALL_Res_Kstar               [j]   [0][1] = new TH1D(HistName6s,HistName3,kStarBinNum,kStarSta,kStarEnd);
                            H_ALL_dRap                    [j]   [0][1] = new TH1D(HistName4r,HistName3,dRapBinNum,dRapSta,dRapEnd);
                            H_ALL_Mix_dRap                [j]   [0][1] = new TH1D(HistName5r,HistName3,dRapBinNum,dRapSta,dRapEnd);
                            H_ALL_Res_dRap                [j]   [0][1] = new TH1D(HistName6r,HistName3,dRapBinNum,dRapSta,dRapEnd);
                            H_ALL_dPt                     [j]   [0][1] = new TH1D(HistName4p,HistName3,dPtBinNum,dPtSta,dPtEnd);
                            H_ALL_Mix_dPt                 [j]   [0][1] = new TH1D(HistName5p,HistName3,dPtBinNum,dPtSta,dPtEnd);
                            H_ALL_Res_dPt                 [j]   [0][1] = new TH1D(HistName6p,HistName3,dPtBinNum,dPtSta,dPtEnd);
                            HistNameM = "H_ALL_";HistNameM += j;HistNameM += "_AMBS_Kstar_dRap";
                            H_ALL_Kstar_dRap              [j]   [0][1] = new TH2F(HistNameM,HistNameM,80,kStarSta,kStarEnd,60,dRapSta,dRapEnd);
                            HistNameM = "H_ALL_Mix_";HistNameM += j;HistNameM += "_AMBS_Kstar_dRap";
                            H_ALL_Mix_Kstar_dRap          [j]   [0][1] = new TH2F(HistNameM,HistNameM,80,kStarSta,kStarEnd,60,dRapSta,dRapEnd);
                            HistNameM = HistName4 + "A_Num_AMBS";
                            H_ALL_A_Num                   [j]   [0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName6 + "A_Num_AMBS";
                            H_ALL_Res_A_Num               [j]   [0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName4 + "B_Num_AMBS";
                            H_ALL_B_Num                   [j]   [0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName6 + "B_Num_AMBS";
                            H_ALL_Res_B_Num               [j]   [0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName4 + "Event_Num_AMBS";
                            H_ALL_Event_Num               [j]   [0][1] = new TProfile(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName6 + "Event_Num_AMBS";
                            H_ALL_Res_Event_Num           [j]   [0][1] = new TProfile(HistNameM,HistNameM,1,-1,1);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_S_AMBS_Kstar";
                            // R_ALL_S_Kstar                 [j]   [0][1] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_M_AMBS_Kstar";
                            // R_ALL_M_Kstar                 [j]   [0][1] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_S_AMBS_dRap";
                            // R_ALL_S_dRap                  [j]   [0][1] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_M_AMBS_dRap";
                            // R_ALL_M_dRap                  [j]   [0][1] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_S_AMBS_dPt";
                            // R_ALL_S_dPt                   [j]   [0][1] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_M_AMBS_dPt";
                            // R_ALL_M_dPt                   [j]   [0][1] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_A_AMBS_Num";
                            // R_ALL_A_Num                   [j]   [0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_B_AMBS_Num";
                            // R_ALL_B_Num                   [j]   [0][1] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        }
                    }
                    if (l == 2) { // ASBM
                        H_Kstar                        [i][j][k][1][0] = new TH1D(HistNameIs,HistName2s,kStarBinNum,kStarSta,kStarEnd);
                        H_Mix_Kstar                    [i][j][k][1][0] = new TH1D(HistNameI,HistName2,kStarBinNum,kStarSta,kStarEnd);
                        H_Res_Kstar                    [i][j][k][1][0] = new TH1D(HistNameIr,HistNameIr,kStarBinNum,kStarSta,kStarEnd);
                        H_dRap                         [i][j][k][1][0] = new TH1D(HistNameJs,HistName2s,dRapBinNum,dRapSta,dRapEnd);
                        H_Mix_dRap                     [i][j][k][1][0] = new TH1D(HistNameJ,HistName2,dRapBinNum,dRapSta,dRapEnd);
                        H_Res_dRap                     [i][j][k][1][0] = new TH1D(HistNameJr,HistNameJr,dRapBinNum,dRapSta,dRapEnd);
                        H_dPt                          [i][j][k][1][0] = new TH1D(HistNameKs,HistName2s,dPtBinNum,dPtSta,dPtEnd);
                        H_Mix_dPt                      [i][j][k][1][0] = new TH1D(HistNameK,HistName2,dPtBinNum,dPtSta,dPtEnd);
                        H_Res_dPt                      [i][j][k][1][0] = new TH1D(HistNameKr,HistNameKr,dPtBinNum,dPtSta,dPtEnd);
                        H_Mass                         [i][j][k][1][0] = new TH1D(HistNameLs,HistName2s,MBinNum,MSta,MEnd);
                        H_Mix_Mass                     [i][j][k][1][0] = new TH1D(HistNameL,HistName2,MBinNum,MSta,MEnd);
                        H_Res_Mass                     [i][j][k][1][0] = new TH1D(HistNameLr,HistNameLr,MBinNum,MSta,MEnd);
                        Mix_event_Num                  [i][j][k][1][0] = 0;
                        Mix_event_Num_SUM              [i][j][k][1][0] = 0;
                        HistNameM = HistNameMs + "_ASBM";
                        H_A_Num                        [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistNameMr + "_ASBM";
                        H_B_Num                        [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Res_A_Num";
                        H_Res_A_Num                    [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Res_B_Num";
                        H_Res_B_Num                    [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Event_Num";
                        H_Event_Num                    [i][j][k][1][0] = new TProfile(HistNameM,HistNameM,1,-1,1);
                        HistNameM = HistName1s + "_Res_Event_Num";
                        H_Res_Event_Num                [i][j][k][1][0] = new TProfile(HistNameM,HistNameM,1,-1,1);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_S_ASBM_Kstar";
                        // R_S_Kstar                      [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_M_ASBM_Kstar";
                        // R_M_Kstar                      [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_S_ASBM_dRap";
                        // R_S_dRap                       [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_M_ASBM_dRap";
                        // R_M_dRap                       [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_S_ASBM_dPt";
                        // R_S_dPt                        [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_M_ASBM_dPt";
                        // R_M_dPt                        [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_A_ASBM_Num";
                        // R_A_Num                        [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        // HistNameM = "R_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_B_ASBM_Num";
                        // R_B_Num                        [i][j][k][1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        if ((i==0)&&(k==0)){
                            H_ALL_Kstar                   [j]   [1][0] = new TH1D(HistName4s,HistName3,kStarBinNum,kStarSta,kStarEnd);
                            H_ALL_Mix_Kstar               [j]   [1][0] = new TH1D(HistName5s,HistName3,kStarBinNum,kStarSta,kStarEnd);
                            H_ALL_Res_Kstar               [j]   [1][0] = new TH1D(HistName6s,HistName3,kStarBinNum,kStarSta,kStarEnd);
                            H_ALL_dRap                    [j]   [1][0] = new TH1D(HistName4r,HistName3,dRapBinNum,dRapSta,dRapEnd);
                            H_ALL_Mix_dRap                [j]   [1][0] = new TH1D(HistName5r,HistName3,dRapBinNum,dRapSta,dRapEnd);
                            H_ALL_Res_dRap                [j]   [1][0] = new TH1D(HistName6r,HistName3,dRapBinNum,dRapSta,dRapEnd);
                            H_ALL_dPt                     [j]   [1][0] = new TH1D(HistName4p,HistName3,dPtBinNum,dPtSta,dPtEnd);
                            H_ALL_Mix_dPt                 [j]   [1][0] = new TH1D(HistName5p,HistName3,dPtBinNum,dPtSta,dPtEnd);
                            H_ALL_Res_dPt                 [j]   [1][0] = new TH1D(HistName6p,HistName3,dPtBinNum,dPtSta,dPtEnd);
                            HistNameM = "H_ALL_";HistNameM += j;HistNameM += "_ASBM_Kstar_dRap";
                            H_ALL_Kstar_dRap              [j]   [1][0] = new TH2F(HistNameM,HistNameM,80,kStarSta,kStarEnd,60,dRapSta,dRapEnd);
                            HistNameM = "H_ALL_Mix_";HistNameM += j;HistNameM += "_ASBM_Kstar_dRap";
                            H_ALL_Mix_Kstar_dRap          [j]   [1][0] = new TH2F(HistNameM,HistNameM,80,kStarSta,kStarEnd,60,dRapSta,dRapEnd);
                            HistNameM = HistName4 + "A_Num_ASBM";
                            H_ALL_A_Num                   [j]   [1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName6 + "A_Num_ASBM";
                            H_ALL_Res_A_Num               [j]   [1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName4 + "B_Num_ASBM";
                            H_ALL_B_Num                   [j]   [1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName6 + "B_Num_ASBM";
                            H_ALL_Res_B_Num               [j]   [1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName4 + "Event_Num_ASBM";
                            H_ALL_Event_Num               [j]   [1][0] = new TProfile(HistNameM,HistNameM,1,-1,1);
                            HistNameM = HistName6 + "Event_Num_ASBM";
                            H_ALL_Res_Event_Num           [j]   [1][0] = new TProfile(HistNameM,HistNameM,1,-1,1);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_S_ASBM_Kstar";
                            // R_ALL_S_Kstar                 [j]   [1][0] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_M_ASBM_Kstar";
                            // R_ALL_M_Kstar                 [j]   [1][0] = new TH1D(HistNameM,HistNameM,kStarBinNum,kStarSta,kStarEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_S_ASBM_dRap";
                            // R_ALL_S_dRap                  [j]   [1][0] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_M_ASBM_dRap";
                            // R_ALL_M_dRap                  [j]   [1][0] = new TH1D(HistNameM,HistNameM,dRapBinNum,dRapSta,dRapEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_S_ASBM_dPt";
                            // R_ALL_S_dPt                   [j]   [1][0] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_M_ASBM_dPt";
                            // R_ALL_M_dPt                   [j]   [1][0] = new TH1D(HistNameM,HistNameM,dPtBinNum,dPtSta,dPtEnd);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_A_ASBM_Num";
                            // R_ALL_A_Num                   [j]   [1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                            // HistNameM = "R_ALL_" + std::to_string(j) + "_B_ASBM_Num";
                            // R_ALL_B_Num                   [j]   [1][0] = new TH1D(HistNameM,HistNameM,1,-1,1);
                        }
                    }
                }
            }
        }
    }
    H_ALL_Mass                     [0][0] = new TH1D("H_Mass_AMBM","H_Mass_AMBM",MBinNum,MSta,MEnd);
    H_ALL_Mix_Mass                 [0][0] = new TH1D("H_Mass_Mix_AMBM","H_Mass_AMBM",MBinNum,MSta,MEnd);
    H_ALL_Mass                     [0][1] = new TH1D("H_Mass_AMBS","H_Mass_AMBS",MBinNum,MSta,MEnd);
    H_ALL_Mix_Mass                 [0][1] = new TH1D("H_Mass_Mix_AMBS","H_Mass_AMBS",MBinNum,MSta,MEnd);
    H_ALL_Mass                     [1][0] = new TH1D("H_Mass_ASBM","H_Mass_ASBM",MBinNum,MSta,MEnd);
    H_ALL_Mix_Mass                 [1][0] = new TH1D("H_Mass_Mix_ASBM","H_Mass_ASBM",MBinNum,MSta,MEnd);

    for (int PatternID = 0;PatternID < Pattern+1;PatternID++) {  
        TString TreeName = "hadronTree";
        if ((Mode == 0) && (PatternID < Pattern)) continue;
        if ((Mode != 0) && (PatternID == Pattern)) break;
        if (Mode != 0) {
            TreeName = PatternBin[PatternID] + "/" + TreeName;
        }

        TChain *hadronTree = new TChain(TreeName);
        for(i=StartFileIndex;i <= EndFileIndex;i++){
            TString filename = MidName;
            filename+=i;
            filename+=".root";
            hadronTree->Add(filename);
            // cout<<"Add "<<filename<<" Successfully"<<endl;
        }
        Int_t PDGMult  ;
        Int_t refMult  ;
        Int_t grefMult ;
        Int_t EventID  ;
        Int_t RunID    ;
        Int_t TriggerID;
        Int_t Nch      ;
        float PVz      ;

        hadronTree->SetBranchAddress("PDGMult"  ,&PDGMult  );
        // hadronTree->SetBranchAddress("refMult"  ,&refMult  );
        // hadronTree->SetBranchAddress("grefMult" ,&grefMult );
        hadronTree->SetBranchAddress("EventID"  ,&EventID  );
        // hadronTree->SetBranchAddress("RunID"    ,&RunID    );
        hadronTree->SetBranchAddress("TriggerID",&TriggerID);
        hadronTree->SetBranchAddress("Nch"      ,&Nch      );
        hadronTree->SetBranchAddress("PVz"      ,&PVz      );
        
        hadronTree->SetBranchAddress("PDG"          ,&PDG          ,&bPDG          );
        hadronTree->SetBranchAddress("mix_px"       ,&mix_px       ,&bmix_px       );
        hadronTree->SetBranchAddress("mix_py"       ,&mix_py       ,&bmix_py       );
        hadronTree->SetBranchAddress("mix_pz"       ,&mix_pz       ,&bmix_pz       );
        // hadronTree->SetBranchAddress("QA_eta"       ,&QA_eta       ,&bQA_eta       );
        // hadronTree->SetBranchAddress("dEdx"         ,&dEdx         ,&bdEdx         );
        // hadronTree->SetBranchAddress("m2"           ,&m2           ,&bm2           );
        // hadronTree->SetBranchAddress("dcatopv"      ,&dcatopv      ,&bdcatopv      );
        // hadronTree->SetBranchAddress("nSigmaProton" ,&nSigmaProton ,&bnSigmaProton );
        // hadronTree->SetBranchAddress("nSigmaPion"   ,&nSigmaPion   ,&bnSigmaPion   );
        // hadronTree->SetBranchAddress("nSigmaKaon"   ,&nSigmaKaon   ,&bnSigmaKaon   );
        hadronTree->SetBranchAddress("InvariantMass",&InvariantMass,&bInvariantMass);
        // hadronTree->SetBranchAddress("Decay_Length" ,&Decay_Length ,&bDecay_Length );
        // hadronTree->SetBranchAddress("Chi2"         ,&Chi2         ,&bChi2         );
        hadronTree->SetBranchAddress("ParentList"   ,&ParentList   ,&bParentList   );
        hadronTree->SetBranchAddress("ParentSta"    ,&ParentSta    ,&bParentSta    );
        hadronTree->SetBranchAddress("ParentEnd"    ,&ParentEnd    ,&bParentEnd    );

        const Int_t nentries=hadronTree->GetEntries();
        cout << "file number: " << nentries << endl;

        time_t time_start;
        time_t time_now;
        time(&time_start);
        clock_t Tstart = clock();
        for (int EntriesID = 0 ; EntriesID < nentries ; EntriesID++){
            cout<<"1"<<endl;
            hadronTree->GetEntry(EntriesID);
            if ((EntriesID+1)%200 == 0) {
                time(&time_now);
                int time_diff = (int)difftime(time_now, time_start);
                cout << time_diff/60 << "min " << time_diff%60 << "s: ";
                long long microseconds = (clock() - Tstart)/10000;
                std::cout << "Microseconds: " << microseconds << "  ";
                cout << "Test/Events = " << 1.0*TestSum/50 << "  ";
                cout<<"Calculating Event "<<(EntriesID+1)<<"/"<<nentries<<endl;
                Tstart = clock();
            }
            cout<<"2"<<endl;

            A_Num = -1;B_Num = -1;
            A_ParID.resize(0);B_ParID.resize(0);
            C_ParID.resize(0);// IfFoundOmega = false;

            for (j=0;j<PDGMult;j++){
                if (PDG->at(j) == A_PDG) {
                    if ( PatternID == Pattern ) {
                        if      (fabs(InvariantMass->at(j) - AMass) <= 3*AMassSigma) {A_Num++;A_Kind[A_Num]=0;}
                        // else if (fabs(InvariantMass->at(j) - AMass) <= 6*AMassSigma) {A_Kind.push_back(1);}
                        else{continue;}
                    }
                    else {
                        if      (fabs(InvariantMass->at(j) - AMass) <= 3*AMassSigma) {
                            if ((PatternID == 0) || (PatternID == 1)) {A_Num++;A_Kind[A_Num]=0;}
                            else {continue;}
                        }
                        else if (fabs(InvariantMass->at(j) - AMass) <= 6*AMassSigma) {
                            if ((PatternID == 2))                     {A_Num++;A_Kind[A_Num]=1;}
                            else {continue;}
                        }
                    }
                    A_Px[A_Num]=mix_px->at(j);
                    A_Py[A_Num]=mix_py->at(j);
                    A_Pz[A_Num]=mix_pz->at(j);
                    A_TreID[A_Num]=j;
                    A_IfRecord[A_Num]=true;
                    Temp.clear();Temp.push_back(j);
                    for (k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                        Temp.push_back(ParentList->at(k));
                    }
                    A_ParID.push_back(Temp);
                    tEnergy = pow(pow(mix_px->at(j),2) + pow(mix_py->at(j),2) + pow(mix_pz->at(j),2) + AMass*AMass,0.5);
                    A_Rap[A_Num]=0.5*log((tEnergy+mix_pz->at(j))/(tEnergy-mix_pz->at(j)));
                }
                else if (PDG->at(j) == B_PDG) {
                    if ( PatternID == Pattern ) {
                        if      (fabs(InvariantMass->at(j) - BMass) <= 3*BMassSigma) {B_Num++;B_Kind[B_Num]=0;}
                        // else if (fabs(InvariantMass->at(j) - BMass) <= 6*BMassSigma) {B_Kind.push_back(1);}
                        else{continue;}
                    }
                    else {
                        if      (fabs(InvariantMass->at(j) - BMass) <= 3*BMassSigma) {
                            if ((PatternID == 0) || (PatternID == 2)) {B_Num++;B_Kind[B_Num]=0;}
                            else {continue;}
                        }
                        else if (fabs(InvariantMass->at(j) - BMass) <= 6*BMassSigma) {
                            if (PatternID == 1)                       {B_Num++;B_Kind[B_Num]=1;}
                            else {continue;}
                        }
                    }
                    B_Px[B_Num]=mix_px->at(j);
                    B_Py[B_Num]=mix_py->at(j);
                    B_Pz[B_Num]=mix_pz->at(j);
                    B_TreID[B_Num]=j;
                    B_IfRecord[B_Num]=true;
                    Temp.clear();Temp.push_back(j);
                    for (k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                        Temp.push_back(ParentList->at(k));
                    }
                    B_ParID.push_back(Temp);
                    tEnergy = pow(pow(mix_px->at(j),2) + pow(mix_py->at(j),2) + pow(mix_pz->at(j),2) + BMass*BMass,0.5);
                    B_Rap[B_Num]=0.5*log((tEnergy+mix_pz->at(j))/(tEnergy-mix_pz->at(j)));
                }
                else{
                    for (l = 0;l < FeedDownNum;l++) {
                        if ( abs(PDG->at(j)) == FeedDown[l] ) {
                            if ((fabs(InvariantMass->at(j) - CMass.at(l)) > 3*CMassSigma.at(l))) continue;
                            Temp.clear();Temp.push_back(j);
                            for (k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                                Temp.push_back(ParentList->at(k));
                            }
                            C_ParID.push_back(Temp);
                            // IfFoundOmega = true;
                            // cout<<"Found Omega"<<endl;
                        }
                    }
                }
            }
            cout<<"3"<<endl;

            // if ((C_ParID.size() != 0)) {continue;}
            if ((A_Num == -1) || (B_Num == -1)) {continue;}
            A_Num++;B_Num++;
            
            // if (IfFoundOmega) {
            //     for (int Aid = 0;Aid < A_Px.size();Aid++) {
            //         cout<<"{ "<<A_PDG<<" } "<<A_TreID.at(Aid)<<" th ";print(A_ParID.at(Aid));
            //     }
            //     for (int Bid = 0;Bid < B_Px.size();Bid++) {
            //         cout<<"{ "<<B_PDG<<" } "<<B_TreID.at(Bid)<<" th ";print(B_ParID.at(Bid));
            //     }
            //     for (int Cid = 0;Cid < C_ParID.size();Cid++) {
            //         cout<<"{ "<<FeedDown[0]<<" } "<<(C_ParID.at(Cid)).at(0)<<" th ";print(C_ParID.at(Cid));
            //     }
            // }

            // A rapidity cut
            for (Aid = 0;Aid < A_Num;Aid++) {
                if ((A_Rap[Aid] < AyCut[0]) || (A_Rap[Aid] > AyCut[1])){
                    A_IfRecord[Aid] = false;
                }
            }

            // A Eta Cut
            for (Aid = 0;Aid < A_Num;Aid++) {
                Eta = -1.0*log(tan(0.5*(acos(A_Pz[Aid]/pow(A_Px[Aid]*A_Px[Aid]+A_Py[Aid]*A_Py[Aid]+A_Pz[Aid]*A_Pz[Aid],0.5)))));
                if ((Eta < EtaCut[0]) || (Eta > EtaCut[1])){
                    A_IfRecord[Aid] = false;
                }
            }

            // B Eta Cut
            for (Bid = 0;Bid < B_Num;Bid++) {
                Eta = -1.0*log(tan(0.5*(acos(B_Pz[Bid]/pow(B_Px[Bid]*B_Px[Bid]+B_Py[Bid]*B_Py[Bid]+B_Pz[Bid]*B_Pz[Bid],0.5)))));
                if ((Eta < EtaCut[0]) || (Eta > EtaCut[1])){
                    B_IfRecord[Bid] = false;
                }
            }
            // 如果A、B有血缘关系，保留B
            for (Bid = 0;Bid < B_Num;Bid++) {
                for (Aid = 0;Aid < A_Num;Aid++) {
                    if (IfInVector(A_TreID[Aid] , B_ParID.at(Bid))){
                        A_IfRecord[Aid] = false;
                    }
                }
            }
            
            // 如果A、B与C有血缘关系，不记录A和B
            for (Aid = 0;Aid < A_Num;Aid++) {
                for (Cid = 0;Cid < C_ParID.size();Cid++) {
                    if (IfInVector(A_TreID[Aid] , C_ParID.at(Cid))) {
                        A_IfRecord[Aid] = false;
                    }
                }
            }
            for (Bid = 0;Bid < B_Num;Bid++) {
                for (Cid = 0;Cid < C_ParID.size();Cid++) {
                    if (IfInVector(B_TreID[Bid] , C_ParID.at(Cid))) {
                        B_IfRecord[Bid] = false;
                    }
                }
            }
            cout<<"4"<<endl;

            // 减除Km-Lambda的不变质量疑似为Omega的pair
            if (IfRemoveFeedPair) {
                for (Bid = 0;Bid < B_Num;Bid++) {
                    if (B_IfRecord[Bid]) {
                        for (Aid = 0;Aid < A_Num;Aid++) {
                            if (A_IfRecord[Aid]) {
                                PairMass = GetPairMass(A_Px[Aid],A_Py[Aid],A_Pz[Aid],AMass,B_Px[Bid],B_Py[Bid],B_Pz[Bid],BMass);
                                for (Cid = 0;Cid < FeedDownNum;Cid++) {
                                    if (fabs(PairMass-CMass.at(Cid))<=3*CMassSigma.at(Cid)) {
                                        A_IfRecord[Aid] = false;
                                        B_IfRecord[Bid] = false;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            cout<<"5"<<endl;

            // Event Index
            int CenIndex = -1;
            for (k=0;k<CentralityBinNum;k++){
                NNch = CenCorr(PVz) * Nch;
                // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                    CenIndex = k;
                    break;
                }
            }
            if (CenIndex == -1) continue;
            
            int PVzIndex = -1;
            for (k=0;k<PVzBinNum;k++){
                // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                if ((PVzBin[k] <= PVz) && (PVz < PVzBin[k+1])) {
                    PVzIndex = k;
                    break;
                }
            }
            if (PVzIndex == -1) continue;

            for (i = 0;i < CentralityBinNum;i++) {
                for (j = 0;j < yBinNum;j++) {
                    for (k = 0;k < PVzBinNum;k++) {
                        for (A_Kid = 0;A_Kid < 2;A_Kid++) {
                            for (B_Kid = 0;B_Kid < 2;B_Kid++) {
                                EventPatternMatch[i]  [j] [k][A_Kid][B_Kid] = 0;
                            }
                        }
                    }
                }
            }
            cout<<"6"<<endl;

            for (Bid = 0;Bid < B_Num;Bid++) {

                if (!(B_IfRecord[Bid])) continue;

                BPx = B_Px[Bid] , BPy = B_Py[Bid] , BPz = B_Pz[Bid];

                // B Index
                rap = B_Rap[Bid];
                int RapIndex = -1;
                for (k=0;k<yBinNum;k++){
                    if ((yBin[k] <= rap) && (rap < yBin[k+1])) {
                        RapIndex = k;
                        break;
                    }
                }

                if ((RapIndex == -1)) {
                    continue;
                }

                B_Kid = B_Kind[Bid];
                for (Aid = 0;Aid < A_Num;Aid++) {
                    if (!(A_IfRecord[Aid])) continue;

                    A_Kid = A_Kind[Aid];

                    if (!IfInVector(Aid , Mix_A_ID[RapIndex] [A_Kid][B_Kid])) {Mix_A_ID[RapIndex] [A_Kid][B_Kid].push_back(Aid);Mix_A_ID_Index[RapIndex] [A_Kid][B_Kid]++;}
                    if (!IfInVector(Bid , Mix_B_ID[RapIndex] [A_Kid][B_Kid])) {Mix_B_ID[RapIndex] [A_Kid][B_Kid].push_back(Bid);Mix_B_ID_Index[RapIndex] [A_Kid][B_Kid]++;}

                    TestSum++;
                }
            }
            cout<<"7"<<endl;

            for (RapIndex = 0;RapIndex < yBinNum;RapIndex++) {
                cout<<"71"<<endl;
                for (A_Kid = 0;A_Kid < 2;A_Kid++) {
                    cout<<"72"<<endl;
                    for (B_Kid = 0;B_Kid < 2;B_Kid++) {
                        cout<<"73"<<endl;
                        if ((Mix_A_ID_Index[RapIndex] [A_Kid][B_Kid] != 0) && (Mix_B_ID_Index[RapIndex] [A_Kid][B_Kid] != 0)) {
                            cout<<"74"<<endl;
                            for (i = 0;i < Mix_A_ID_Index[RapIndex] [A_Kid][B_Kid];i++) {
                                cout<<"741"<<endl;
                                Mix_A_Index_T = Mix_A_Index[CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid];
                                cout<<"742"<<endl;
                                cout<<"Mix_A_ID["+RapIndex+"] ["+A_Kid+"]["+B_Kid+"].size() = "<<Mix_A_ID[RapIndex] [A_Kid][B_Kid].size()<<endl;
                                cout<<"Mix_A_ID_Index["+RapIndex+"] ["+A_Kid+"]["+B_Kid+"] = "<<Mix_A_ID_Index[RapIndex] [A_Kid][B_Kid]<<endl;
                                AidN = Mix_A_ID[RapIndex] [A_Kid][B_Kid][i];
                                cout<<"743"<<endl;
                                Mix_A_Px        [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_A_Index_T] = A_Px [AidN];
                                Mix_A_Py        [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_A_Index_T] = A_Py [AidN];
                                Mix_A_Pz        [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_A_Index_T] = A_Pz [AidN];
                                Mix_A_Rap       [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_A_Index_T] = A_Rap[AidN];
                                cout<<"744"<<endl;
                                Mix_A_EvtID     [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_A_Index_T] = EntriesID;
                                cout<<"745"<<endl;
                                Mix_A_IfMadePair[CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_A_Index_T] = false;
                                cout<<"746"<<endl;
                                Mix_A_Index     [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid]++;
                                cout<<"747"<<endl;
                            }
                            cout<<"75"<<endl;
                            for (i = 0;i < Mix_B_ID_Index[RapIndex] [A_Kid][B_Kid];i++) {
                                Mix_B_Index_T = Mix_B_Index[CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid];
                                BidN = Mix_B_ID[RapIndex] [A_Kid][B_Kid][i];
                                Mix_B_Px        [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_B_Index_T] = B_Px [BidN];
                                Mix_B_Py        [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_B_Index_T] = B_Py [BidN];
                                Mix_B_Pz        [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_B_Index_T] = B_Pz [BidN];
                                Mix_B_Rap       [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_B_Index_T] = B_Rap[BidN];
                                Mix_B_EvtID     [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_B_Index_T] = EntriesID;
                                Mix_B_IfMadePair[CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid][Mix_B_Index_T] = false;
                                Mix_B_Index     [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid]++;
                            }
                            EventPatternMatch[CenIndex][RapIndex][PVzIndex][A_Kid][B_Kid]++;
                        }
                        cout<<"76"<<endl;
                        Mix_A_ID      [RapIndex] [A_Kid][B_Kid].clear();
                        Mix_B_ID      [RapIndex] [A_Kid][B_Kid].clear();
                        Mix_A_ID_Index[RapIndex] [A_Kid][B_Kid] = 0;
                        Mix_B_ID_Index[RapIndex] [A_Kid][B_Kid] = 0;
                    }
                }
            }
            cout<<"8"<<endl;

            for (i = 0;i < yBinNum;i++) {
                for (j = 0;j < PVzBinNum;j++) {
                    for (Aid = 0;Aid < 2;Aid++) {
                        for (Bid = 0;Bid < 2;Bid++) {
                            if (EventPatternMatch[CenIndex][i][j][Aid][Bid] != 0) {
                                Mix_event_Num[CenIndex][i][j][Aid][Bid]++;

                                if (Mix_event_Num[CenIndex][i][j][Aid][Bid] == HowMuchEventMixing+1) {
                                    Mix_A_Size = Mix_A_Index[CenIndex][i][j][Aid][Bid];
                                    Mix_B_Size = Mix_B_Index[CenIndex][i][j][Aid][Bid];
                                    for (Aindex = 0;Aindex < Mix_A_Size;Aindex++) {
                                        A_EID = Mix_A_EvtID[CenIndex][i][j][Aid][Bid] [Aindex];
                                        APx   = Mix_A_Px   [CenIndex][i][j][Aid][Bid] [Aindex];
                                        APy   = Mix_A_Py   [CenIndex][i][j][Aid][Bid] [Aindex];
                                        APz   = Mix_A_Pz   [CenIndex][i][j][Aid][Bid] [Aindex];
                                        for (Bindex = 0;Bindex < Mix_B_Size;Bindex++) {
                                            BPx = Mix_B_Px[CenIndex][i][j][Aid][Bid] [Bindex];
                                            BPy = Mix_B_Py[CenIndex][i][j][Aid][Bid] [Bindex];
                                            BPz = Mix_B_Pz[CenIndex][i][j][Aid][Bid] [Bindex];

                                            // p2.SetXYZM(BPx,BPy,BPz,BMass);
                                            // p1.SetXYZM(APx,APy,APz,AMass);
                                            // p3 = p1 + p2;
                                            // BV = -p3.BoostVector();
                                            // p1.Boost( BV);p2.Boost( BV);
                                            // PairMass = p1.Energy()+p2.Energy();

                                            float* MassAndKstar = GetPairMassAndKstar(APx , APy , APz , BPx , BPy , BPz , AMass , BMass);
                                            PairMass = MassAndKstar[0];

                                            if (IfRemoveFeedPair) {
                                                IfRecord = true;
                                                for (Cid = 0;Cid < FeedDownNum;Cid++) {
                                                    if (fabs(PairMass-CMass.at(Cid))<=3*CMassSigma.at(Cid)) {
                                                        IfRecord = false;
                                                        break;
                                                    }
                                                }
                                            }
                                            if (A_EID != Mix_B_EvtID[CenIndex][i][j][Aid][Bid] [Bindex]) {
                                                // KS = 0.5 * (p2 - p1).Rho();
                                                KS = MassAndKstar[1];
                                                rap = Mix_A_Rap[CenIndex][i][j][Aid][Bid] [Aindex] - Mix_B_Rap[CenIndex][i][j][Aid][Bid] [Bindex];
                                                Pt = fabs(pow(APx*APx + APy*APy , 0.5) - pow(BPx*BPx + BPy*BPy , 0.5));
                                                if (IfRecord) {
                                                    H_Mix_Kstar     [CenIndex][i][j][Aid][Bid]->Fill(KS);
                                                    H_ALL_Mix_Kstar           [i]   [Aid][Bid]->Fill(KS);
                                                    H_Mix_Mass      [CenIndex][i][j][Aid][Bid]->Fill(PairMass);
                                                    H_ALL_Mix_Mass                  [Aid][Bid]->Fill(PairMass);
                                                }
                                                H_Mix_dRap          [CenIndex][i][j][Aid][Bid]->Fill(rap);
                                                H_ALL_Mix_dRap                [i]   [Aid][Bid]->Fill(rap);
                                                H_Mix_dPt           [CenIndex][i][j][Aid][Bid]->Fill(Pt);
                                                H_ALL_Mix_dPt                 [i]   [Aid][Bid]->Fill(Pt);
                                                if (SpecialMode) {
                                                    H_ALL_Mix_Kstar_dRap  [i]   [Aid][Bid]->Fill(KS,rap);
                                                }
                                                Mix_A_IfMadePair[CenIndex][i][j][Aid][Bid] [Aindex] = true;
                                                Mix_B_IfMadePair[CenIndex][i][j][Aid][Bid] [Bindex] = true;

                                                // R_S_Kstar       [CenIndex][i][j][Aid][Bid]->Fill(KS);
                                                // R_ALL_S_Kstar             [i]   [Aid][Bid]->Fill(KS);
                                                // R_S_dRap        [CenIndex][i][j][Aid][Bid]->Fill(rap);
                                                // R_ALL_S_dRap              [i]   [Aid][Bid]->Fill(rap);
                                            }
                                            else{
                                                KS = MassAndKstar[1];
                                                rap = Mix_A_Rap[CenIndex][i][j][Aid][Bid] [Aindex] - Mix_B_Rap[CenIndex][i][j][Aid][Bid] [Bindex];
                                                Pt = fabs(pow(APx*APx + APy*APy , 0.5) - pow(BPx*BPx + BPy*BPy , 0.5));
                                                H_Kstar         [CenIndex][i][j][Aid][Bid]->Fill(KS);
                                                H_Res_Kstar     [CenIndex][i][j][Aid][Bid]->Fill(KS);
                                                H_ALL_Kstar               [i]   [Aid][Bid]->Fill(KS);
                                                H_ALL_Res_Kstar           [i]   [Aid][Bid]->Fill(KS);
                                                H_dRap          [CenIndex][i][j][Aid][Bid]->Fill(rap);
                                                H_Res_dRap      [CenIndex][i][j][Aid][Bid]->Fill(rap);
                                                H_ALL_dRap                [i]   [Aid][Bid]->Fill(rap);
                                                H_ALL_Res_dRap            [i]   [Aid][Bid]->Fill(rap);
                                                H_dPt           [CenIndex][i][j][Aid][Bid]->Fill(Pt);
                                                H_Res_dPt       [CenIndex][i][j][Aid][Bid]->Fill(Pt);
                                                H_ALL_dPt                 [i]   [Aid][Bid]->Fill(Pt);
                                                H_ALL_Res_dPt             [i]   [Aid][Bid]->Fill(Pt);
                                                H_Mass          [CenIndex][i][j][Aid][Bid]->Fill(PairMass);
                                                H_ALL_Mass                      [Aid][Bid]->Fill(PairMass);
                                                if (SpecialMode) {
                                                    H_ALL_Kstar_dRap      [i]   [Aid][Bid]->Fill(KS,rap);
                                                }
                                                Mix_A_IfMadePair[CenIndex][i][j][Aid][Bid] [Aindex] = true;
                                                Mix_B_IfMadePair[CenIndex][i][j][Aid][Bid] [Bindex] = true;
                                            }
                                            delete[] MassAndKstar;
                                        }
                                    }
                                    for (Aindex = 0;Aindex < Mix_A_Size;Aindex++) {
                                        if (Mix_A_IfMadePair[CenIndex][i][j][Aid][Bid] [Aindex] == true) {
                                            H_A_Num        [CenIndex][i][j][Aid][Bid]->Fill(0);
                                            H_Res_A_Num    [CenIndex][i][j][Aid][Bid]->Fill(0);
                                            H_ALL_A_Num              [i]   [Aid][Bid]->Fill(0);
                                            H_ALL_Res_A_Num          [i]   [Aid][Bid]->Fill(0);
                                        }
                                    }
                                    for (Bindex = 0;Bindex < Mix_B_Size;Bindex++) {
                                        if (Mix_B_IfMadePair[CenIndex][i][j][Aid][Bid] [Bindex] == true) {
                                            H_B_Num        [CenIndex][i][j][Aid][Bid]->Fill(0);
                                            H_Res_B_Num    [CenIndex][i][j][Aid][Bid]->Fill(0);
                                            H_ALL_B_Num              [i]   [Aid][Bid]->Fill(0);
                                            H_ALL_Res_B_Num          [i]   [Aid][Bid]->Fill(0);
                                        }
                                    }
                                    Mix_event_Num    [CenIndex][i][j][Aid][Bid] = 0;
                                    Mix_event_Num_SUM[CenIndex][i][j][Aid][Bid]++;
                                    H_Event_Num      [CenIndex][i][j][Aid][Bid]->Fill(0,HowMuchEventMixing+1);
                                    H_ALL_Event_Num            [i]   [Aid][Bid]->Fill(0,HowMuchEventMixing+1);
                                    // Mix_A_Px[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_B_Px[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_A_Py[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_B_Py[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_A_Pz[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_B_Pz[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_A_EvtID[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_B_EvtID[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_A_Rap[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_B_Rap[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_A_IfMadePair[CenIndex][i][j][Aid][Bid].clear();
                                    // Mix_B_IfMadePair[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_A_Index      [CenIndex][i][j][Aid][Bid] = 0;
                                    Mix_B_Index      [CenIndex][i][j][Aid][Bid] = 0;
                                }
                            }
                        }
                    }
                }

            }
            cout<<"9"<<endl;

        }
    }
    
    TString OutputFileName = OutMidName;
    OutputFileName += "H_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileA = new TFile(OutputFileName, "RECREATE");
    folder_kStar = fileA->mkdir("kStar");
    folder_dRap  = fileA->mkdir("dRap");
    folder_dPt   = fileA->mkdir("dPt");
    folder_Mass  = fileA->mkdir("Mass");
    folder_Test  = fileA->mkdir("Test");

    cout << "#######################" << endl;
    cout << "# Calculating Summary #" << endl;
    cout << "#######################" << endl;
    fileA->cd();
    H_ALL_Mass    [0][0]->Write();
    H_ALL_Mix_Mass[0][0]->Write();
    H_ALL_Mass    [0][1]->Write();
    H_ALL_Mix_Mass[0][1]->Write();
    H_ALL_Mass    [1][0]->Write();
    H_ALL_Mix_Mass[1][0]->Write();
    for (i=0;i<CentralityBinNum;i++){
        for (j=0;j<yBinNum;j++){
            for (k=0;k<PVzBinNum;k++){
                for (A_Kid=0;A_Kid<2;A_Kid++){
                    for (B_Kid=0;B_Kid<2;B_Kid++) {
                        if (A_Kid == 1 && B_Kid == 1) continue;
                        if ((Mix_event_Num[i][j][k][A_Kid][B_Kid] != 0) || (Mix_event_Num_SUM[i][j][k][A_Kid][B_Kid] != 0)) {
                        // if (true) {
                            TString Name;
                            if (A_Kid == 0 && B_Kid == 0) Name = "AMBM";
                            if (A_Kid == 0 && B_Kid == 1) Name = "AMBS";
                            if (A_Kid == 1 && B_Kid == 0) Name = "ASBM";
                            cout<<"["<<i<<","<<j<<","<<k<<","<<Name<<"] Filled " << Mix_event_Num_SUM[i][j][k][A_Kid][B_Kid] * HowMuchEventMixing << " events, and remain "<<Mix_event_Num[i][j][k][A_Kid][B_Kid]<<" events, "<<endl; 
                            if ((Mix_event_Num[i][j][k][A_Kid][B_Kid] != 0)) { // 存取剩余池子里的events中粒子，填进*_Res_*
                                Mix_A_Size = Mix_A_Index[i][j][k][A_Kid][B_Kid];
                                Mix_B_Size = Mix_B_Index[i][j][k][A_Kid][B_Kid];
                                for (Aindex = 0;Aindex < Mix_A_Size;Aindex++) {
                                    A_EID = Mix_A_EvtID[i][j][k][A_Kid][B_Kid] [Aindex];
                                    APx   = Mix_A_Px   [i][j][k][A_Kid][B_Kid] [Aindex];
                                    APy   = Mix_A_Py   [i][j][k][A_Kid][B_Kid] [Aindex];
                                    APz   = Mix_A_Pz   [i][j][k][A_Kid][B_Kid] [Aindex];
                                    for (Bindex = 0;Bindex < Mix_B_Size;Bindex++) {
                                        BPx = Mix_B_Px[i][j][k][A_Kid][B_Kid] [Bindex];
                                        BPy = Mix_B_Py[i][j][k][A_Kid][B_Kid] [Bindex];
                                        BPz = Mix_B_Pz[i][j][k][A_Kid][B_Kid] [Bindex];

                                        if (SpecialMode) {
                                            if (APz+BPz < 0) continue;
                                        }

                                        // p2.SetXYZM(BPx,BPy,BPz,BMass);
                                        // p1.SetXYZM(APx,APy,APz,AMass);
                                        // p3 = p1 + p2;
                                        // BV = -p3.BoostVector();
                                        // p1.Boost( BV);p2.Boost( BV);
                                        // PairMass = p1.Energy()+p2.Energy();

                                        float* MassAndKstar = GetPairMassAndKstar(APx , APy , APz , BPx , BPy , BPz , AMass , BMass);
                                        PairMass = MassAndKstar[0];

                                        if (IfRemoveFeedPair) {
                                            IfRecord = true;
                                            for (Cid = 0;Cid < FeedDownNum;Cid++) {
                                                if (fabs(PairMass-CMass.at(Cid))<=3*CMassSigma.at(Cid)) {
                                                    IfRecord = false;
                                                    break;
                                                }
                                            }
                                            if (!IfRecord) continue;
                                        }
                                        if (A_EID != Mix_B_EvtID[i][j][k][A_Kid][B_Kid] [Bindex]) {
                                            continue;
                                        }
                                        else{
                                            KS = MassAndKstar[1];
                                            rap = Mix_A_Rap[i][j][k][A_Kid][B_Kid] [Aindex] - Mix_B_Rap[i][j][k][A_Kid][B_Kid] [Bindex];
                                            Pt = fabs(pow(APx*APx + APy*APy , 0.5) - pow(BPx*BPx + BPy*BPy , 0.5));
                                            H_Res_Kstar     [i][j][k][A_Kid][B_Kid]->Fill(KS);
                                            H_Res_dRap      [i][j][k][A_Kid][B_Kid]->Fill(rap);
                                            H_Res_dPt       [i][j][k][A_Kid][B_Kid]->Fill(Pt);
                                            H_ALL_Res_Kstar    [j]   [A_Kid][B_Kid]->Fill(KS);
                                            H_ALL_Res_dRap     [j]   [A_Kid][B_Kid]->Fill(rap);
                                            H_ALL_Res_dPt      [j]   [A_Kid][B_Kid]->Fill(Pt);
                                            H_Res_Mass      [i][j][k][A_Kid][B_Kid]->Fill(PairMass);
                                            Mix_A_IfMadePair[i][j][k][A_Kid][B_Kid] [Aindex] = true;
                                            Mix_B_IfMadePair[i][j][k][A_Kid][B_Kid] [Bindex] = true;
                                        }
                                        delete[] MassAndKstar;
                                    }
                                }
                                for (Aindex = 0;Aindex < Mix_A_Size;Aindex++) {
                                    if (Mix_A_IfMadePair[i][j][k][A_Kid][B_Kid] [Aindex] == true) {
                                        H_Res_A_Num     [i][j][k][A_Kid][B_Kid]->Fill(0);
                                        H_ALL_Res_A_Num    [j]   [A_Kid][B_Kid]->Fill(0);
                                    }
                                }
                                for (Bindex = 0;Bindex < Mix_B_Size;Bindex++) {
                                    if (Mix_B_IfMadePair[i][j][k][A_Kid][B_Kid] [Bindex] == true) {
                                        H_Res_B_Num     [i][j][k][A_Kid][B_Kid]->Fill(0);
                                        H_ALL_Res_B_Num    [j]   [A_Kid][B_Kid]->Fill(0);
                                    }
                                }
                                H_Res_Event_Num    [i][j][k][A_Kid][B_Kid]->Fill(0,Mix_event_Num[i][j][k][A_Kid][B_Kid]);
                                H_ALL_Res_Event_Num   [j]   [A_Kid][B_Kid]->Fill(0,Mix_event_Num[i][j][k][A_Kid][B_Kid]);
                            }
                        }
                        fileA->cd();
                        H_A_Num        [i][j][k][A_Kid][B_Kid]->Write();
                        H_B_Num        [i][j][k][A_Kid][B_Kid]->Write();
                        H_Res_A_Num    [i][j][k][A_Kid][B_Kid]->Write();
                        H_Res_B_Num    [i][j][k][A_Kid][B_Kid]->Write();
                        H_Event_Num    [i][j][k][A_Kid][B_Kid]->Write();
                        H_ALL_Event_Num   [j]   [A_Kid][B_Kid]->Write();
                        H_Res_Event_Num[i][j][k][A_Kid][B_Kid]->Write();
                        folder_kStar->cd();
                        if(H_Kstar             [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Kstar             [i][j][k][A_Kid][B_Kid]->Write();
                        if(H_Mix_Kstar         [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Mix_Kstar         [i][j][k][A_Kid][B_Kid]->Write();
                        if(H_Res_Kstar         [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Res_Kstar         [i][j][k][A_Kid][B_Kid]->Write();
                        folder_dRap->cd();
                        if(H_dRap              [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_dRap              [i][j][k][A_Kid][B_Kid]->Write();
                        if(H_Mix_dRap          [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Mix_dRap          [i][j][k][A_Kid][B_Kid]->Write();
                        if(H_Res_dRap          [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Res_dRap          [i][j][k][A_Kid][B_Kid]->Write();
                        folder_dPt->cd();
                        if(H_dPt               [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_dPt               [i][j][k][A_Kid][B_Kid]->Write();
                        if(H_Mix_dPt           [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Mix_dPt           [i][j][k][A_Kid][B_Kid]->Write();
                        if(H_Res_dPt           [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Res_dPt           [i][j][k][A_Kid][B_Kid]->Write();
                        folder_Mass->cd();
                        if(H_Mass              [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Mass              [i][j][k][A_Kid][B_Kid]->Write();
                        if(H_Mix_Mass          [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Mix_Mass          [i][j][k][A_Kid][B_Kid]->Write();
                        if(H_Res_Mass          [i][j][k][A_Kid][B_Kid]->GetEntries() != 0) H_Res_Mass          [i][j][k][A_Kid][B_Kid]->Write();
                        folder_Test->cd();
                        if(H_ALL_Mix_Kstar_dRap   [j]   [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Mix_Kstar_dRap   [j]   [A_Kid][B_Kid]->Write();
                        if(H_ALL_Kstar_dRap       [j]   [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Kstar_dRap       [j]   [A_Kid][B_Kid]->Write();
                    }
                }
            }
        }
    }
    fileA->cd();
    for (j=0;j<yBinNum;j++){
        for (A_Kid=0;A_Kid<2;A_Kid++){
            for (B_Kid=0;B_Kid<2;B_Kid++) {
                if ((A_Kid == 1)&&(B_Kid == 1)) continue;
                if(H_ALL_Kstar      [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Kstar      [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_Mix_Kstar  [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Mix_Kstar  [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_Res_Kstar  [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Res_Kstar  [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_dRap       [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_dRap       [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_Mix_dRap   [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Mix_dRap   [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_Res_dRap   [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Res_dRap   [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_dPt        [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_dPt        [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_Mix_dPt    [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Mix_dPt    [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_Res_dPt    [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Res_dPt    [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_A_Num      [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_A_Num      [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_B_Num      [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_B_Num      [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_Res_A_Num  [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Res_A_Num  [j]  [A_Kid][B_Kid]->Write();
                if(H_ALL_Res_B_Num  [j]  [A_Kid][B_Kid]->GetEntries() != 0) H_ALL_Res_B_Num  [j]  [A_Kid][B_Kid]->Write();
                // cout<<"J = "<<j<<" Stored."<<endl;
            }
        }
    }
    fileA->Close();

    cout << "#######################" << endl;
    cout << "# Finish storing Hist #" << endl;
    cout << "#######################" << endl;

    // Writing remaining pool

    OutputFileName = OutMidName;
    OutputFileName += "T_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileB = new TFile(OutputFileName, "RECREATE");
    folder_AMBM = fileB->mkdir("AMBM");
    folder_AMBS = fileB->mkdir("AMBS");
    folder_ASBM = fileB->mkdir("ASBM");
    fileB->cd();
    if (Mode == 0){
        int A_Kid , B_Kid;
        for (int WriteTreeIndex = 0;WriteTreeIndex < Pattern;WriteTreeIndex++) {
            if (WriteTreeIndex == 0) {folder_AMBM->cd();A_Kid = 0;B_Kid = 0;}
            if (WriteTreeIndex == 1) {folder_AMBS->cd();A_Kid = 0;B_Kid = 1;}
            if (WriteTreeIndex == 2) {folder_ASBM->cd();A_Kid = 1;B_Kid = 0;}
            int buffer_size = 5000000;
            int BPDGMult  ;
            int BCrefMult ;
            int BCgrefMult;
            int BevtID    ;
            int BrunID    ;
            int BTriggerID;
            int BNch      ;
            float BPVz    ;
            std::vector<int> BPDG               ;BPDG            .clear();
            std::vector<float> Bpx              ;Bpx             .clear();
            std::vector<float> Bpy              ;Bpy             .clear();
            std::vector<float> Bpz              ;Bpz             .clear();
            std::vector<float> BQA_eta          ;BQA_eta         .clear();
            std::vector<float> BQA_dEdx         ;BQA_dEdx        .clear();
            std::vector<float> BQA_m2           ;BQA_m2          .clear();
            std::vector<float> BQA_DCA_V0_PV    ;BQA_DCA_V0_PV   .clear();
            std::vector<float> BQA_nSigmaProton ;BQA_nSigmaProton.clear();
            std::vector<float> BQA_nSigmaPion   ;BQA_nSigmaPion  .clear();
            std::vector<float> BQA_nSigmaKaon   ;BQA_nSigmaKaon  .clear();
            std::vector<float> BInvariantMass   ;BInvariantMass  .clear();
            std::vector<float> BQA_Decay_Length ;BQA_Decay_Length.clear();
            std::vector<float> BQA_Chi2         ;BQA_Chi2        .clear();
            std::vector<int> BParentList        ;BParentList     .clear();
            std::vector<int> BParentSta         ;BParentSta      .clear();
            std::vector<int> BParentEnd         ;BParentEnd      .clear();
            BhadronTree = new TTree("hadronTree", "Tree_STAR");
            BhadronTree->Branch("PDGMult"            ,&BPDGMult             ,"PDGMult/I"                           );
            // BhadronTree->Branch("refMult"            ,&BCrefMult            ,"refMult/I"                           );
            // BhadronTree->Branch("grefMult"           ,&BCgrefMult           ,"grefMult/I"                          );
            BhadronTree->Branch("EventID"            ,&BevtID               ,"EventID/I"                           );
            // BhadronTree->Branch("RunID"              ,&BrunID               ,"RunID/I"                             );
            BhadronTree->Branch("TriggerID"          ,&BTriggerID           ,"TriggerID/I"                         );
            BhadronTree->Branch("Nch"                ,&BNch                 ,"Nch/I"                               );
            BhadronTree->Branch("PVz"                ,&BPVz                 ,"PVz/F"                               );
            BhadronTree->Branch("PDG"                ,&BPDG                 );
            BhadronTree->Branch("mix_px"             ,&Bpx                  );
            BhadronTree->Branch("mix_py"             ,&Bpy                  );
            BhadronTree->Branch("mix_pz"             ,&Bpz                  );
            // BhadronTree->Branch("QA_eta"             ,&BQA_eta              );

            // Used for PID QA
            // BhadronTree->Branch("dEdx"               ,&BQA_dEdx              );
            // BhadronTree->Branch("m2"                 ,&BQA_m2                );
            // BhadronTree->Branch("dcatopv"            ,&BQA_DCA_V0_PV         );
            // BhadronTree->Branch("nSigmaProton"       ,&BQA_nSigmaProton      );
            // BhadronTree->Branch("nSigmaPion"         ,&BQA_nSigmaPion        );
            // BhadronTree->Branch("nSigmaKaon"         ,&BQA_nSigmaKaon        );
            
            // Used for Reconstruction QA
            BhadronTree->Branch("InvariantMass"      ,&BInvariantMass        );
            // BhadronTree->Branch("Decay_Length"       ,&BQA_Decay_Length      );
            // BhadronTree->Branch("Chi2"               ,&BQA_Chi2              );
            
            // Used for restore corralated information
            BhadronTree->Branch("ParentList"         ,&BParentList     );
            BhadronTree->Branch("ParentSta"          ,&BParentSta      );
            BhadronTree->Branch("ParentEnd"          ,&BParentEnd      );

            std::vector<Int_t> Mix_EvtID;
            for (i=0;i<CentralityBinNum;i++){
                for (j=0;j<yBinNum;j++){
                    for (k=0;k<PVzBinNum;k++){
                        for (m=0;m<Mix_B_EvtID[i][j][k][A_Kid][B_Kid].size();m++){
                            int nIndex = -1;
                            for (n=0;n<Mix_EvtID.size();n++){
                                if (Mix_B_EvtID[i][j][k][A_Kid][B_Kid].at(m) == Mix_EvtID.at(n)){
                                    nIndex = n;
                                    break;
                                }
                            }
                            if (nIndex == -1){
                                Mix_EvtID.push_back(Mix_B_EvtID[i][j][k][A_Kid][B_Kid].at(m));
                                nIndex = Mix_EvtID.size() - 1;
                            }
                        }
                        for (m=0;m<Mix_A_EvtID[i][j][k][A_Kid][B_Kid].size();m++){
                            int nIndex = -1;
                            for (n=0;n<Mix_EvtID.size();n++){
                                if (Mix_A_EvtID[i][j][k][A_Kid][B_Kid].at(m) == Mix_EvtID.at(n)){
                                    nIndex = n;
                                    break;
                                }
                            }
                            if (nIndex == -1){
                                Mix_EvtID.push_back(Mix_A_EvtID[i][j][k][A_Kid][B_Kid].at(m));
                                nIndex = Mix_EvtID.size() - 1;
                            }
                        }
                        
                    }
                }
            }
            for (i=0;i<Mix_EvtID.size();i++){
                hadronTree->GetEntry(Mix_EvtID.at(i));
                BPDGMult   = PDGMult  ;
                // BCrefMult  = refMult  ;
                // BCgrefMult = grefMult ;
                BevtID     = EventID  ;
                // BrunID     = RunID    ;
                BTriggerID = TriggerID;
                BNch       = Nch      ;
                BPVz       = PVz      ;
                for (j=0;j<PDGMult;j++){
                    BPDG            .push_back(PDG          ->at(j));
                    Bpx             .push_back(mix_px       ->at(j));
                    Bpy             .push_back(mix_py       ->at(j));
                    Bpz             .push_back(mix_pz       ->at(j));
                    // BQA_eta         .push_back(QA_eta       ->at(j));
                    // BQA_dEdx        .push_back(dEdx         ->at(j));
                    // BQA_m2          .push_back(m2           ->at(j));
                    // BQA_DCA_V0_PV   .push_back(dcatopv      ->at(j));
                    // BQA_nSigmaProton.push_back(nSigmaProton ->at(j));
                    // BQA_nSigmaPion  .push_back(nSigmaPion   ->at(j));
                    // BQA_nSigmaKaon  .push_back(nSigmaKaon   ->at(j));
                    BInvariantMass  .push_back(InvariantMass->at(j));
                    // BQA_Decay_Length.push_back(Decay_Length ->at(j));
                    // BQA_Chi2        .push_back(Chi2         ->at(j));
                    BParentSta      .push_back(ParentSta    ->at(j));
                    BParentEnd      .push_back(ParentEnd    ->at(j));
                }
                for (j=0;j<ParentList.size();j++){
                    BParentList     .push_back(ParentList   ->at(j));
                }
                BhadronTree->Fill();
                BPDG            .clear();
                Bpx             .clear();
                Bpy             .clear();
                Bpz             .clear();
                // BQA_eta         .clear();
                // BQA_dEdx        .clear();
                // BQA_m2          .clear();
                // BQA_DCA_V0_PV   .clear();
                // BQA_nSigmaProton.clear();
                // BQA_nSigmaPion  .clear();
                // BQA_nSigmaKaon  .clear();
                BInvariantMass  .clear();
                // BQA_Decay_Length.clear();
                // BQA_Chi2        .clear();
                BParentList     .clear();
                BParentSta      .clear();
                BParentEnd      .clear();
            }

            BhadronTree->Write();
        }
    }
    fileB->Write();
    fileB->Close();

    return;
}
