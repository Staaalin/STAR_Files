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

#define Pi 3.1415926535898

const int CentralityBin[] = {0 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 60 , 70 , 80};// %
const float PVzBin[] = {-45.0 , -35.0 , -25.0 , -15.0 , -5.0 , 5.0 , 15.0 , 25.0 , 35.0 , 45.0 , 55.0}; // Primary Vertex Z (cm) d+Au@200 GeV RUN 21 : -45 ~ 55 cm
const float yBin[] = {-8.0 , 0.0 , 0.5 , 8.0}; // B_y
int FeedDown[] = { 3334 , 3312 , 1003314 };

const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1
const Int_t PVzBinNum = sizeof(PVzBin)/sizeof(PVzBin[0]) - 1; // -1
const Int_t yBinNum = sizeof(yBin)/sizeof(yBin[0]) - 1; // -1
const Int_t FeedDownNum = sizeof(FeedDown)/sizeof(FeedDown[0]);

TString KindBin[] = {"Mid","Sid"}
#define KindNum 2
TString PatternBin[] = {"AMBM","AMBS","ASBM"};
#define Pattern 3 // 0:A middle B middle , 1:A middle B sideband , 2:A sideband B middle
// Pattern应当大于KindNum

int HowMuchEventMixing = 10;

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
    float tEnergy , APx , APy , APz , BPx , BPy , BPz;
    int A_Kid , B_Kid , Mix_A_Size , Mix_B_Size , A_EID , AidN , BidN;
    std::vector<int> Temp;
    std::vector<float> CMass , CMassSigma;
    bool IfRecord = true;
    float BMass = massList(B_PDG)           , AMass = massList(A_PDG);
    float BMassSigma = massListSigma(B_PDG) , AMassSigma = massListSigma(A_PDG);

    std::vector<int> NchList = GetNchList(CentralityBin , CentralityBinNum+1);     // centrality
    cout<<"NchList = ";
    print(NchList);
    cout<<" "<<endl;
    //                                        centrality          B_y        PVz
    std::vector<Float_t>                  A_Px           ;
    std::vector<Float_t>                  A_Py           ;
    std::vector<Float_t>                  A_Pz           ;
    std::vector<Int_t>                    A_TreID        ;
    std::vector<std::vector<int> >        A_ParID        ;
    std::vector<int>                      A_Kind         ;
    std::vector<Float_t>                  A_Rap          ;
    std::vector<Int_t>                    A_IfRecord     ;
    std::vector<Float_t>                  B_Px           ;
    std::vector<Float_t>                  B_Py           ;
    std::vector<Float_t>                  B_Pz           ;
    std::vector<Int_t>                    B_TreID        ;
    std::vector<std::vector<int> >        B_ParID        ;
    std::vector<int>                      B_Kind         ;
    std::vector<Float_t>                  B_Rap          ;
    std::vector<Int_t>                    B_IfRecord     ;
    std::vector<std::vector<int> >        C_ParID        ; // 用于存储Residal Effect
    // used as array
    std::vector<float> Mix_A_Px           [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<float> Mix_A_Py           [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<float> Mix_A_Pz           [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<int>   Mix_A_TreID        [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<int>   Mix_A_EvtID        [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<int>   Mix_A_ID                                [yBinNum]               [2] [2] ;
    std::vector<float> Mix_A_Rap          [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<float> Mix_B_Px           [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<float> Mix_B_Py           [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<float> Mix_B_Pz           [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<int>   Mix_B_TreID        [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<int>   Mix_B_EvtID        [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
    std::vector<int>   Mix_B_ID                                [yBinNum]               [2] [2] ;
    std::vector<float> Mix_B_Rap          [CentralityBinNum]   [yBinNum]  [PVzBinNum]  [2] [2] ;
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
    int EventPatternMatch                 [15]                 [15]       [15]         [2] [2] ;
    // Used for testing
    int TestSum = 0;

    int kStarBinNum = 500;
    float kStarSta = 0 , kStarEnd = 10;
    
    int dRapBinNum = 100;
    float dRapSta = -5 , dRapEnd = 5;
    
    int dPtBinNum = 200;
    float dPtSta = 0 , dPtEnd = 10;
    
    int MBinNum = 500 , MBinPar = 50;
    float MSta = floor((AMass + BMass)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;
    
    TString HistNameI  , HistNameJ  , HistNameK  , HistNameL;
    TString HistNameIs , HistNameJs , HistNameKs , HistNameLs;

    for (int i = 0;i < FeedDownNum;i++){
        CMass.push_back(massList(FeedDown[i]));
        CMassSigma.push_back(massListSigma(FeedDown[i]));
        if (abs(FeedDown[i]) == A_PDG) {
            FeedDown[i] = 0;
        }
        if (abs(FeedDown[i]) == B_PDG) {
            FeedDown[i] = 0;
        }
    }

    for (int i=0;i<CentralityBinNum;i++){
        for (int l=0;l<Pattern;l++){
            for (int k=0;k<PVzBinNum;k++){
                for (int j=0;j<yBinNum;j++){
                    TString HistName1 = "H_";
                    TString HistName2 = "Cen: [";
                    HistName1 += i;HistName1 += "_";
                    HistName2 += CentralityBin[i];HistName2 += "% , ";
                    HistName2 += CentralityBin[i+1];HistName2 += "%], ";
                    HistName1 += j;HistName1 += "_";
                    HistName2 += yBin[j];HistName2 += " < y";HistName2 += B_PDG;HistName2 += " <  ";
                    HistName2 += yBin[j+1];HistName2 += ", ";
                    HistName1 += k;
                    HistName2 += PVzBin[k];HistName2 += " < PVz";HistName2 += B_PDG;HistName2 += " <  ";
                    HistName2 += PVzBin[k+1];
                    TString HistName1s = HistName1;
                    TString HistName2s = HistName2;
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
                    HistNameJ  = HistName1  + "_dRap";
                    HistNameJs = HistName1s + "_dRap";
                    HistNameK  = HistName1  + "_dPt";
                    HistNameKs = HistName1s + "_dPt";
                    HistNameL  = HistName1  + "_Mass";
                    HistNameLs = HistName1s + "_Mass";
                    if (l == 0) { // AMBM
                        H_Kstar                        [i][j][k][0][0] = new TH1D(HistNameIs,HistName2s,kStarBinNum,kStarSta,kStarEnd);
                        H_Mix_Kstar                    [i][j][k][0][0] = new TH1D(HistNameI,HistName2,kStarBinNum,kStarSta,kStarEnd);
                        H_dRap                         [i][j][k][0][0] = new TH1D(HistNameJs,HistName2s,dRapBinNum,dRapSta,dRapEnd);
                        H_Mix_dRap                     [i][j][k][0][0] = new TH1D(HistNameJ,HistName2,dRapBinNum,dRapSta,dRapEnd);
                        H_dPt                          [i][j][k][0][0] = new TH1D(HistNameKs,HistName2s,dPtBinNum,dPtSta,dPtEnd);
                        H_Mix_dPt                      [i][j][k][0][0] = new TH1D(HistNameK,HistName2,dPtBinNum,dPtSta,dPtEnd);
                        H_Mass                         [i][j][k][0][0] = new TH1D(HistNameLs,HistName2s,MBinNum,MSta,MEnd);
                        H_Mix_Mass                     [i][j][k][0][0] = new TH1D(HistNameL,HistName2,MBinNum,MSta,MEnd);
                        Mix_event_Num                  [i][j][k][0][0] = 0;
                        Mix_event_Num_SUM              [i][j][k][0][0] = 0;
                    }
                    if (l == 1) { // AMBS
                        H_Kstar                        [i][j][k][0][1] = new TH1D(HistNameIs,HistName2s,kStarBinNum,kStarSta,kStarEnd);
                        H_Mix_Kstar                    [i][j][k][0][1] = new TH1D(HistNameI,HistName2,kStarBinNum,kStarSta,kStarEnd);
                        H_dRap                         [i][j][k][0][1] = new TH1D(HistNameJs,HistName2s,dRapBinNum,dRapSta,dRapEnd);
                        H_Mix_dRap                     [i][j][k][0][1] = new TH1D(HistNameJ,HistName2,dRapBinNum,dRapSta,dRapEnd);
                        H_dPt                          [i][j][k][0][1] = new TH1D(HistNameKs,HistName2s,dPtBinNum,dPtSta,dPtEnd);
                        H_Mix_dPt                      [i][j][k][0][1] = new TH1D(HistNameK,HistName2,dPtBinNum,dPtSta,dPtEnd);
                        H_Mass                         [i][j][k][0][1] = new TH1D(HistNameLs,HistName2s,MBinNum,MSta,MEnd);
                        H_Mix_Mass                     [i][j][k][0][1] = new TH1D(HistNameL,HistName2,MBinNum,MSta,MEnd);
                        Mix_event_Num                  [i][j][k][0][1] = 0;
                        Mix_event_Num_SUM              [i][j][k][0][1] = 0;
                    }
                    if (l == 2) { // ASBM
                        H_Kstar                        [i][j][k][1][0] = new TH1D(HistNameIs,HistName2s,kStarBinNum,kStarSta,kStarEnd);
                        H_Mix_Kstar                    [i][j][k][1][0] = new TH1D(HistNameI,HistName2,kStarBinNum,kStarSta,kStarEnd);
                        H_dRap                         [i][j][k][1][0] = new TH1D(HistNameJs,HistName2s,dRapBinNum,dRapSta,dRapEnd);
                        H_Mix_dRap                     [i][j][k][1][0] = new TH1D(HistNameJ,HistName2,dRapBinNum,dRapSta,dRapEnd);
                        H_dPt                          [i][j][k][1][0] = new TH1D(HistNameKs,HistName2s,dPtBinNum,dPtSta,dPtEnd);
                        H_Mix_dPt                      [i][j][k][1][0] = new TH1D(HistNameK,HistName2,dPtBinNum,dPtSta,dPtEnd);
                        H_Mass                         [i][j][k][1][0] = new TH1D(HistNameLs,HistName2s,MBinNum,MSta,MEnd);
                        H_Mix_Mass                     [i][j][k][1][0] = new TH1D(HistNameL,HistName2,MBinNum,MSta,MEnd);
                        Mix_event_Num                  [i][j][k][1][0] = 0;
                        Mix_event_Num_SUM              [i][j][k][1][0] = 0;
                    }
                }
            }
        }
    }

    for (int PatternID = 0;PatternID < Pattern+1;PatternID++) {  
        TString TreeName = "hadronTree";
        if ((Mode == 0) && (PatternID < Pattern)) continue;
        if ((Mode != 0) && (PatternID == Pattern)) break;
        if (Mode != 0) {
            TreeName = PatternBin[PatternID] + "/" + TreeName;
        }

        TChain *hadronTree = new TChain(TreeName);
        for(int i=StartFileIndex;i <= EndFileIndex;i++){
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

            A_Px   .resize(0);   B_Px.resize(0);
            A_Py   .resize(0);   B_Py.resize(0);
            A_Pz   .resize(0);   B_Pz.resize(0);
            A_ParID.resize(0);B_ParID.resize(0);
            A_TreID.resize(0);B_TreID.resize(0);
            A_Kind .resize(0); B_Kind.resize(0);
            A_Rap  .resize(0);  B_Rap.resize(0);
            A_IfRecord.resize(0);B_IfRecord.resize(0);
            C_ParID.resize(0);

            for (int j=0;j<PDGMult;j++){
                if (PDG->at(j) == A_PDG) {
                    if ( PatternID == Pattern ) {
                        if      (fabs(InvariantMass->at(j) - AMass) <= 3*AMassSigma) {A_Kind.push_back(0);}
                        // else if (fabs(InvariantMass->at(j) - AMass) <= 6*AMassSigma) {A_Kind.push_back(1);}
                        else{continue;}
                    }
                    else {
                        if      (fabs(InvariantMass->at(j) - AMass) <= 3*AMassSigma) {
                            if ((PatternID == 0) || (PatternID == 1)) {A_Kind.push_back(0);}
                            else {continue;}
                        }
                        else if (fabs(InvariantMass->at(j) - AMass) <= 6*AMassSigma) {
                            if ((PatternID == 2))                     {A_Kind.push_back(1);}
                            else {continue;}
                        }
                    }
                    A_Px.push_back(mix_px->at(j));
                    A_Py.push_back(mix_py->at(j));
                    A_Pz.push_back(mix_pz->at(j));
                    A_TreID.push_back(j);
                    A_IfRecord.push_back(1);
                    Temp.clear();Temp.push_back(j);
                    for (int k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                        Temp.push_back(ParentList->at(k));
                    }
                    A_ParID.push_back(Temp);
                    tEnergy = pow(pow(mix_px->at(j),2) + pow(mix_py->at(j),2) + pow(mix_pz->at(j),2) + pow(AMass,2),0.5);
                    A_Rap.push_back(0.5*log((tEnergy+mix_pz->at(j))/(tEnergy-mix_pz->at(j))));
                }
                else if (PDG->at(j) == B_PDG) {
                    if ( PatternID == Pattern ) {
                        if      (fabs(InvariantMass->at(j) - BMass) <= 3*BMassSigma) {B_Kind.push_back(0);}
                        // else if (fabs(InvariantMass->at(j) - BMass) <= 6*BMassSigma) {B_Kind.push_back(1);}
                        else{continue;}
                    }
                    else {
                        if      (fabs(InvariantMass->at(j) - BMass) <= 3*BMassSigma) {
                            if ((PatternID == 0) || (PatternID == 2)) {B_Kind.push_back(0);}
                            else {continue;}
                        }
                        else if (fabs(InvariantMass->at(j) - BMass) <= 6*BMassSigma) {
                            if (PatternID == 1)                       {B_Kind.push_back(1);}
                            else {continue;}
                        }
                    }
                    B_Px.push_back(mix_px->at(j));
                    B_Py.push_back(mix_py->at(j));
                    B_Pz.push_back(mix_pz->at(j));
                    B_TreID.push_back(j);
                    B_IfRecord.push_back(1);
                    Temp.clear();Temp.push_back(j);
                    for (int k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                        Temp.push_back(ParentList->at(k));
                    }
                    B_ParID.push_back(Temp);
                    tEnergy = pow(pow(mix_px->at(j),2) + pow(mix_py->at(j),2) + pow(mix_pz->at(j),2) + pow(BMass,2),0.5);
                    B_Rap.push_back(0.5*log((tEnergy+mix_pz->at(j))/(tEnergy-mix_pz->at(j))));
                }
                else{
                    for (int l = 0;l < FeedDownNum;l++) {
                        if ( abs(PDG->at(j)) == FeedDown[l] ) {
                            if ((fabs(InvariantMass->at(j) - CMass.at(l)) > 3*CMassSigma.at(l))) continue;
                            Temp.clear();Temp.push_back(j);
                            for (int k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                                Temp.push_back(ParentList->at(k));
                            }
                            C_ParID.push_back(Temp);
                        }
                    }
                }
            }

            if ((A_Px.size() == 0) || (B_Px.size() == 0)) {continue;}
            
            // 如果A、B有血缘关系，保留B
            for (int Bid = 0;Bid < B_Px.size();Bid++) {
                for (int Aid = 0;Aid < A_Px.size();Aid++) {
                    if (IfInVector(A_TreID.at(Aid) , B_ParID.at(Bid))){
                        A_IfRecord.at(Aid) = 0;
                        cout<<"Meet1!"<<endl;
                    }
                }
            }
            
            // 如果A、B与C有血缘关系，不记录A和B
            for (int Aid = 0;Aid < A_Px.size();Aid++) {
                for (int Cid = 0;Cid < C_ParID.size();Cid++) {
                    if (IfInVector(A_TreID.at(Aid) , C_ParID.at(Cid))) {
                        A_IfRecord.at(Aid) = 0;
                        cout<<"Meet2!"<<endl;
                    }
                }
            }
            for (int Bid = 0;Bid < B_Px.size();Bid++) {
                for (int Cid = 0;Cid < C_ParID.size();Cid++) {
                    if (IfInVector(B_TreID.at(Bid) , C_ParID.at(Cid))) {
                        B_IfRecord.at(Bid) = 0;
                        cout<<"Meet3!"<<endl;
                    }
                }
            }

            // Event Index
            int CenIndex = -1;
            for (int k=0;k<CentralityBinNum;k++){
                // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                if ((NchList.at(k) >= Nch) && (Nch > NchList.at(k+1))) {
                    CenIndex = k;
                    break;
                }
            }
            if (CenIndex == -1) continue;
            
            int PVzIndex = -1;
            for (int k=0;k<PVzBinNum;k++){
                // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                if ((PVzBin[k] <= PVz) && (PVz < PVzBin[k+1])) {
                    PVzIndex = k;
                    break;
                }
            }
            if (PVzIndex == -1) continue;

            for (int i = 0;i < CentralityBinNum;i++) {
                for (int j = 0;j < yBinNum;j++) {
                    for (int k = 0;k < PVzBinNum;k++) {
                        for (int A_Kid = 0;A_Kid < 2;A_Kid++) {
                            for (int B_Kid = 0;B_Kid < 2;B_Kid++) {
                                EventPatternMatch[i]  [j] [k][A_Kid][B_Kid] = 0;
                            }
                        }
                    }
                }
            }

            for (int Bid = 0;Bid < B_Px.size();Bid++) {

                if (B_IfRecord.at(Bid) == 0) continue;

                BPx = B_Px.at(Bid) , BPy = B_Py.at(Bid) , BPz = B_Pz.at(Bid);

                // B Index
                rap = B_Rap.at(Bid);
                int RapIndex = -1;
                for (int k=0;k<yBinNum;k++){
                    if ((yBin[k] <= rap) && (rap < yBin[k+1])) {
                        RapIndex = k;
                        break;
                    }
                }

                if ((RapIndex == -1)) {
                    continue;
                }

                B_Kid = B_Kind.at(Bid);
                for (int Aid = 0;Aid < A_Px.size();Aid++) {
                    if (A_IfRecord.at(Aid) == 0) continue;

                    A_Kid = A_Kind.at(Aid);

                    if (!IfInVector(Aid , Mix_A_ID[RapIndex] [A_Kid][B_Kid])) Mix_A_ID[RapIndex] [A_Kid][B_Kid].push_back(Aid);
                    if (!IfInVector(Bid , Mix_B_ID[RapIndex] [A_Kid][B_Kid])) Mix_B_ID[RapIndex] [A_Kid][B_Kid].push_back(Bid);

                    TestSum++;
                    // if (TestSum%10 == 0) {
                    //     TLorentzVector Np1 , Np2;
                    //     Np2.SetXYZM(BPx,BPy,BPz,BMass);
                    //     int A_Kid = A_Kind.at(Aid);
                    //     float APx = A_Px.at(Aid) , APy = A_Py.at(Aid) , APz = A_Pz.at(Aid);
                    //     Np1.SetXYZM(APx,APy,APz,AMass);
                    //     TLorentzVector p3 = Np1 + Np2;
                    //     // p3.SetXYZM(APx + BPx,APy + BPy,APz + BPz,AMass + BMass);
                    //     TVector3 BV = -p3.BoostVector();
                    //     Np1.Boost( BV);Np2.Boost( BV);
                    //     if (fabs(Np1.Px() - boostedP1[0]) > 0.0001) {
                    //         cout<<"This is BAD:"<<endl;
                    //         cout<<"Np1 = ( "<<Np1.Px()<<" , "<<Np1.Py()<<" , "<<Np1.Pz()<<" , "<<Np1.E()<<" )"<<endl;
                    //         cout<<"p1 = ";print(boostedP1);
                    //         cout<<Np1.Px()<<" , "<<boostedP1[0]<<endl;
                    //         cout<<"____________________________"<<endl;
                    //     }
                    //     else{
                    //         cout<<"This is GOOD:"<<endl;
                    //         cout<<"Np1 = ( "<<Np1.Px()<<" , "<<Np1.Py()<<" , "<<Np1.Pz()<<" , "<<Np1.E()<<" )"<<endl;
                    //         cout<<"p1 = ";print(boostedP1);
                    //         cout<<Np1.Px()<<" , "<<boostedP1[0]<<endl;
                    //         cout<<"____________________________"<<endl;
                    //     }
                    // }
                }
            }

            for (int RapIndex = 0;RapIndex < yBinNum;RapIndex++) {
                for (int A_Kid = 0;A_Kid < 2;A_Kid++) {
                    for (int B_Kid = 0;B_Kid < 2;B_Kid++) {
                        if ((Mix_A_ID[RapIndex] [A_Kid][B_Kid].size() != 0) && (Mix_B_ID[RapIndex] [A_Kid][B_Kid].size() != 0)) {
                            for (int i = 0;i < Mix_A_ID[RapIndex] [A_Kid][B_Kid].size();i++) {
                                AidN = Mix_A_ID[RapIndex] [A_Kid][B_Kid] .at(i);
                                Mix_A_Px   [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(A_Px.at(AidN));
                                Mix_A_Py   [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(A_Py.at(AidN));
                                Mix_A_Pz   [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(A_Pz.at(AidN));
                                Mix_A_Rap  [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(A_Rap.at(AidN));
                                Mix_A_EvtID[CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(EntriesID);
                            }
                            for (int i = 0;i < Mix_B_ID[RapIndex] [A_Kid][B_Kid].size();i++) {
                                BidN = Mix_B_ID[RapIndex] [A_Kid][B_Kid] .at(i);
                                Mix_B_Px   [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(B_Px.at(BidN));
                                Mix_B_Py   [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(B_Py.at(BidN));
                                Mix_B_Pz   [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(B_Pz.at(BidN));
                                Mix_B_Rap  [CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(B_Rap.at(BidN));
                                Mix_B_EvtID[CenIndex][RapIndex][PVzIndex] [A_Kid][B_Kid].push_back(EntriesID);
                            }
                            EventPatternMatch[CenIndex][RapIndex][PVzIndex][A_Kid][B_Kid]++;
                            // Test imfomation
                            // Mix_A_Size = Mix_A_Px[CenIndex][RapIndex][PVzIndex][A_Kid][B_Kid].size();
                            // Mix_B_Size = Mix_B_Px[CenIndex][RapIndex][PVzIndex][A_Kid][B_Kid].size();
                            // cout<<"####  2rd  #####   "<<EventPatternMatch[CenIndex][RapIndex][PVzIndex][A_Kid][B_Kid]<<"  ["<<CenIndex<<","<<RapIndex<<","<<PVzIndex<<","<<A_Kid<<","<<B_Kid<<"] "<<Mix_event_Num[CenIndex][RapIndex][PVzIndex][A_Kid][B_Kid]<<"   ################"<<endl;
                            // int LoopSize; Mix_B_Size < Mix_A_Size ? LoopSize = Mix_A_Size : LoopSize = Mix_B_Size ;
                            // for (int Bindex = 0;Bindex < LoopSize;Bindex++) {
                            //     cout<<Bindex<<"  ";
                            //     if (Bindex < Mix_A_Size) {cout<<Mix_A_EvtID[CenIndex][RapIndex][PVzIndex][A_Kid][B_Kid].at(Bindex);}
                            //     else {cout<<"       ";}
                            //     cout<<"  ";
                            //     if (Bindex < Mix_B_Size) {cout<<Mix_B_EvtID[CenIndex][RapIndex][PVzIndex][A_Kid][B_Kid].at(Bindex);}
                            //     else {cout<<"       ";}
                            //     cout<<" "<<endl;
                            // }
                        }
                        Mix_A_ID[RapIndex] [A_Kid][B_Kid].resize(0);
                        Mix_B_ID[RapIndex] [A_Kid][B_Kid].resize(0);
                    }
                }
            }

            for (int i = 0;i < yBinNum;i++) {
                for (int j = 0;j < PVzBinNum;j++) {
                    for (int Aid = 0;Aid < 2;Aid++) {
                        for (int Bid = 0;Bid < 2;Bid++) {
                            if (EventPatternMatch[CenIndex][i][j][Aid][Bid] != 0) {
                                Mix_event_Num[CenIndex][i][j][Aid][Bid]++;
                                // Test imfomation
                                // Mix_A_Size = Mix_A_Px[CenIndex][i][j][Aid][Bid].size();
                                // Mix_B_Size = Mix_B_Px[CenIndex][i][j][Aid][Bid].size();
                                // cout<<"####  3rd  #####   "<<EventPatternMatch[CenIndex][i][j][Aid][Bid]<<"  ["<<CenIndex<<","<<i<<","<<j<<","<<Aid<<","<<Bid<<"] "<<Mix_event_Num[CenIndex][i][j][Aid][Bid]<<"   ################"<<endl;
                                // int LoopSize; Mix_B_Size < Mix_A_Size ? LoopSize = Mix_A_Size : LoopSize = Mix_B_Size ;
                                // for (int Bindex = 0;Bindex < LoopSize;Bindex++) {
                                //     cout<<Bindex<<"  ";
                                //     if (Bindex < Mix_A_Size) {cout<<Mix_A_EvtID[CenIndex][i][j][Aid][Bid].at(Bindex);}
                                //     else {cout<<"       ";}
                                //     cout<<"  ";
                                //     if (Bindex < Mix_B_Size) {cout<<Mix_B_EvtID[CenIndex][i][j][Aid][Bid].at(Bindex);}
                                //     else {cout<<"       ";}
                                //     cout<<" "<<endl;
                                // }

                                if (Mix_event_Num[CenIndex][i][j][Aid][Bid] == HowMuchEventMixing+1) {
                                    Mix_A_Size = Mix_A_Px[CenIndex][i][j][Aid][Bid].size();
                                    Mix_B_Size = Mix_B_Px[CenIndex][i][j][Aid][Bid].size();
                                    for (int Aindex = 0;Aindex < Mix_A_Size;Aindex++) {
                                        A_EID = Mix_A_EvtID[CenIndex][i][j][Aid][Bid].at(Aindex);
                                        APx   = Mix_A_Px   [CenIndex][i][j][Aid][Bid].at(Aindex);
                                        APy   = Mix_A_Py   [CenIndex][i][j][Aid][Bid].at(Aindex);
                                        APz   = Mix_A_Pz   [CenIndex][i][j][Aid][Bid].at(Aindex);
                                        for (int Bindex = 0;Bindex < Mix_B_Size;Bindex++) {
                                            BPx = Mix_B_Px[CenIndex][i][j][Aid][Bid].at(Bindex);
                                            BPy = Mix_B_Py[CenIndex][i][j][Aid][Bid].at(Bindex);
                                            BPz = Mix_B_Pz[CenIndex][i][j][Aid][Bid].at(Bindex);

                                            p2.SetXYZM(BPx,BPy,BPz,BMass);
                                            p1.SetXYZM(APx,APy,APz,AMass);
                                            p3 = p1 + p2;
                                            BV = -p3.BoostVector();
                                            p1.Boost( BV);p2.Boost( BV);

                                            if (A_EID != Mix_B_EvtID[CenIndex][i][j][Aid][Bid].at(Bindex)) {
                                                H_Mix_Kstar[CenIndex][i][j][Aid][Bid]->Fill(0.5 * (p2 - p1).Rho());
                                                H_Mix_dRap [CenIndex][i][j][Aid][Bid]->Fill(Mix_A_Rap[CenIndex][i][j][Aid][Bid].at(Aindex) - Mix_B_Rap[CenIndex][i][j][Aid][Bid].at(Bindex));
                                                H_Mix_dPt  [CenIndex][i][j][Aid][Bid]->Fill(fabs(pow(APx*APx + APy*APy , 0.5) - pow(BPx*BPx + BPy*BPy , 0.5)));
                                                H_Mix_Mass [CenIndex][i][j][Aid][Bid]->Fill(p1.Energy()+p2.Energy());
                                            }
                                            else{
                                                H_Kstar[CenIndex][i][j][Aid][Bid]->Fill(0.5 * (p2 - p1).Rho());
                                                H_dRap [CenIndex][i][j][Aid][Bid]->Fill(Mix_A_Rap[CenIndex][i][j][Aid][Bid].at(Aindex) - Mix_B_Rap[CenIndex][i][j][Aid][Bid].at(Bindex));
                                                H_dPt  [CenIndex][i][j][Aid][Bid]->Fill(fabs(pow(APx*APx + APy*APy , 0.5) - pow(BPx*BPx + BPy*BPy , 0.5)));
                                                H_Mass [CenIndex][i][j][Aid][Bid]->Fill(p1.Energy()+p2.Energy());
                                            }
                                        }
                                    }
                                    Mix_event_Num    [CenIndex][i][j][Aid][Bid] = 0;
                                    Mix_event_Num_SUM[CenIndex][i][j][Aid][Bid]++;
                                    Mix_A_Px[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_B_Px[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_A_Py[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_B_Py[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_A_Pz[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_B_Pz[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_A_EvtID[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_B_EvtID[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_A_Rap[CenIndex][i][j][Aid][Bid].clear();
                                    Mix_B_Rap[CenIndex][i][j][Aid][Bid].clear();
                                }
                            }
                        }
                    }
                }

            }

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

    cout << "#######################" << endl;
    cout << "# Calculating Summary #" << endl;
    cout << "#######################" << endl;
    fileA->cd();
    for (int i=0;i<CentralityBinNum;i++){
        for (int j=0;j<yBinNum;j++){
            for (int k=0;k<PVzBinNum;k++){
                for (int A_Kid=0;A_Kid<2;A_Kid++){
                    for (int B_Kid=0;B_Kid<2;B_Kid++) {
                        if (A_Kid == 1 && B_Kid == 1) continue;
                        if ((Mix_event_Num[i][j][k][A_Kid][B_Kid] != 0) || (Mix_event_Num_SUM[i][j][k][A_Kid][B_Kid] != 0)) {
                        // if (true) {
                            TString Name;
                            if (A_Kid == 0 && B_Kid == 0) Name = "AMBM";
                            if (A_Kid == 0 && B_Kid == 1) Name = "AMBS";
                            if (A_Kid == 1 && B_Kid == 0) Name = "ASBM";
                            cout<<"["<<i<<","<<j<<","<<k<<","<<Name<<"] Filled " << Mix_event_Num_SUM[i][j][k][A_Kid][B_Kid] * HowMuchEventMixing << " events, and remain "<<Mix_event_Num[i][j][k][A_Kid][B_Kid]<<" events, "<<endl; 
                        }
                        folder_kStar->cd();
                        H_Kstar[i][j][k][A_Kid][B_Kid]->Write();
                        H_Mix_Kstar[i][j][k][A_Kid][B_Kid]->Write();
                        folder_dRap->cd();
                        H_dRap[i][j][k][A_Kid][B_Kid]->Write();
                        H_Mix_dRap[i][j][k][A_Kid][B_Kid]->Write();
                        folder_dPt->cd();
                        H_dPt[i][j][k][A_Kid][B_Kid]->Write();
                        H_Mix_dPt[i][j][k][A_Kid][B_Kid]->Write();
                        folder_Mass->cd();
                        H_Mass[i][j][k][A_Kid][B_Kid]->Write();
                        H_Mix_Mass[i][j][k][A_Kid][B_Kid]->Write();
                    }
                }
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
            for (int i=0;i<CentralityBinNum;i++){
                for (int j=0;j<yBinNum;j++){
                    for (int k=0;k<PVzBinNum;k++){
                        for (int m=0;m<Mix_B_EvtID[i][j][k][A_Kid][B_Kid].size();m++){
                            int nIndex = -1;
                            for (int n=0;n<Mix_EvtID.size();n++){
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
                        for (int m=0;m<Mix_A_EvtID[i][j][k][A_Kid][B_Kid].size();m++){
                            int nIndex = -1;
                            for (int n=0;n<Mix_EvtID.size();n++){
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
            for (int i=0;i<Mix_EvtID.size();i++){
                hadronTree->GetEntry(Mix_EvtID.at(i));
                BPDGMult   = PDGMult  ;
                // BCrefMult  = refMult  ;
                // BCgrefMult = grefMult ;
                BevtID     = EventID  ;
                // BrunID     = RunID    ;
                BTriggerID = TriggerID;
                BNch       = Nch      ;
                BPVz       = PVz      ;
                for (int j=0;j<PDGMult;j++){
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
                for (int j=0;j<ParentList.size();j++){
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
