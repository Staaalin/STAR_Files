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

const int CentralityBin[] = {0 , 80};// %
const float PVzBin[] = {-45.0 , 55.0}; // Primary Vertex Z (cm) d+Au@200 GeV RUN 21 : -45 ~ 55 cm
const float yBin[] = {-8.0 , 8.0}; // B_y
int FeedDown[] = { 3334 , 3312 , 1003314 };

const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1
const Int_t PVzBinNum = sizeof(PVzBin)/sizeof(PVzBin[0]) - 1; // -1
const Int_t yBinNum = sizeof(yBin)/sizeof(yBin[0]) - 1; // -1
const Int_t FeedDownNum = sizeof(FeedDown)/sizeof(FeedDown[0]);

TString KindBin[] = {"Mid","Sid"};
#define KindNum 2
TString PatternBin[] = {"AMBM","AMBS","ASBM"};
#define Pattern 3 // 0:A middle B middle , 1:A middle B sideband , 2:A sideband B middle
// Pattern应当大于KindNum

int HowMuchEventMixing = 10;


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

std::vector<float> CAB(float APx , float APy , float APz , float BPx , float BPy , float BPz , float fct) // 0 < fct < 1
{
    std::vector<float> Result;Result.clear();
    float CPx = BPx - APx;
    float CPy = BPy - APy;
    float CPz = BPz - APz;
    float CP = pow(CPx*CPx + CPy*CPy + CPz*CPz,0.5);
    float CPL;
    if (CP > 0.2) {CPL = 0.2 + 0.08*(fct - 0.5);}
    else {CPL = CP*fct;}
    Result.push_back(BPx - CPx*(CPL/CP));
    Result.push_back(BPy - CPy*(CPL/CP));
    Result.push_back(BPz - CPz*(CPL/CP));
    return Result;
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

void RandomGenerator(int OutputFileIndex , int RandomSeed , int nentries = 50000) {

    TRandom3 randGen;
    UInt_t A_Num , B_Num;
    // 设置对数正态分布的均值和标准差
    double mean = 1.0;  // 正态分布的均值
    double sigma = 1.0; // 正态分布的标准差
    int A_PDG = 321;
    int B_PDG = 3122;
    int ListIndex;
    TString OutMidName = "Random_";
    // int OutputFileIndex = 0;
    int PatternID = Pattern;

    std::vector<int>     PDG          ;
    std::vector<Float_t> mix_px       ;
    std::vector<Float_t> mix_py       ;
    std::vector<Float_t> mix_pz       ;
    std::vector<Float_t> QA_eta       ;
    std::vector<Float_t> dEdx         ;
    std::vector<Float_t> m2           ;
    std::vector<Float_t> dcatopv      ;
    std::vector<Float_t> nSigmaProton ;
    std::vector<Float_t> nSigmaPion   ;
    std::vector<Float_t> nSigmaKaon   ;
    std::vector<Float_t> InvariantMass;
    std::vector<Float_t> Decay_Length ;
    std::vector<Float_t> Chi2         ;
    std::vector<int>     ParentList   ;
    std::vector<int>     ParentSta    ;
    std::vector<int>     ParentEnd    ;
    TBranch *bPDG                      ;
    TBranch *bmix_px                   ;
    TBranch *bmix_py                   ;
    TBranch *bmix_pz                   ;
    TBranch *bQA_eta                   ;
    TBranch *bdEdx                     ;
    TBranch *bm2                       ;
    TBranch *bdcatopv                  ;
    TBranch *bnSigmaProton             ;
    TBranch *bnSigmaPion               ;
    TBranch *bnSigmaKaon               ;
    TBranch *bInvariantMass            ;
    TBranch *bDecay_Length             ;
    TBranch *bChi2                     ;
    TBranch *bParentList               ;
    TBranch *bParentSta                ;
    TBranch *bParentEnd                ;

    double kstar, rap;
    TVector3 BetaTemp;
    // ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p1 , p2 , p3 , p4 , p5;
    TLorentzVector p1 , p2 , p3;
    TVector3 BV;
    float tEnergy , APx , APy , APz , ARap , APt , AMon , AEta , ATan , ATheta , BPx , BPy , BPz , BRap , BPt , BEta , BMon , BTan , BTheta;
    int A_Kid , B_Kid , Mix_A_Size , Mix_B_Size , A_EID , AidN , BidN;
    std::vector<int> Temp;
    std::vector<float> CMass , CMassSigma , TempF;
    bool IfRecord = true;
    float BMass = massList(B_PDG)           , AMass = massList(A_PDG);
    float BMassSigma = massListSigma(B_PDG) , AMassSigma = massListSigma(A_PDG);
    long long microseconds;

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
    int EventPatternMatch                 [15]                 [15]       [15]         [2] [2] ;
    // Used for testing
    int TestSum = 0;

    int kStarBinNum = 500;
    float kStarSta = 0 , kStarEnd = 10;
    
    int dRapBinNum = 100;
    float dRapSta = -5 , dRapEnd = 5;
    
    int dPtBinNum = 200;
    float dPtSta = 0 , dPtEnd = 10;
    
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
                    if (l == 0) { // AMBM
                        H_Kstar                        [i][j][k][0][0] = new TH1D(HistNameIs,HistName2s,kStarBinNum,kStarSta,kStarEnd);
                        H_Mix_Kstar                    [i][j][k][0][0] = new TH1D(HistNameI,HistName2,kStarBinNum,kStarSta,kStarEnd);
                        H_dRap                         [i][j][k][0][0] = new TH1D(HistNameJs,HistName2s,dRapBinNum,dRapSta,dRapEnd);
                        H_Mix_dRap                     [i][j][k][0][0] = new TH1D(HistNameJ,HistName2,dRapBinNum,dRapSta,dRapEnd);
                        H_dPt                          [i][j][k][0][0] = new TH1D(HistNameKs,HistName2s,dPtBinNum,dPtSta,dPtEnd);
                        H_Mix_dPt                      [i][j][k][0][0] = new TH1D(HistNameK,HistName2,dPtBinNum,dPtSta,dPtEnd);
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
                        Mix_event_Num                  [i][j][k][1][0] = 0;
                        Mix_event_Num_SUM              [i][j][k][1][0] = 0;
                    }
                }
            }
        }
    }

    TH1D *H_A_Pt   = new TH1D("A_Pt","Pt of A",100,0,10);
    TH1D *H_B_Pt   = new TH1D("B_Pt","Pt of B",100,0,10);
    TH1D *H_A_P    = new TH1D("A_P","P of A",100,0,10);
    TH1D *H_B_P    = new TH1D("B_P","P of B",100,0,10);
    TH1D *H_A_Rap  = new TH1D("A_Rap","Rapidity of A",100,-10,10);
    TH1D *H_B_Rap  = new TH1D("B_Rap","Rapidity of B",100,-10,10);

    // int nentries = 50000;

    Int_t PDGMult  ;
    Int_t refMult  ;
    Int_t grefMult ;
    Int_t EventID  ;
    Int_t RunID    ;
    Int_t TriggerID;
    Int_t Nch      ;
    float PVz      ;

    time_t time_start;
    time_t time_now;
    time(&time_start);
    clock_t Tstart = clock();
    for (int EntriesID = 0 ; EntriesID < nentries ; EntriesID++){

        // 随机产生事件
        PDG          .clear();
        mix_px       .clear();
        mix_py       .clear();
        mix_pz       .clear();
        InvariantMass.clear();
        ParentSta    .clear();
        ParentEnd    .clear();
        ParentList   .clear();
        if ((EntriesID)%2000 == 0) {
            microseconds = (clock() - Tstart)/100;
            time(&time_now);
            Tstart = clock();
            int time_diff = (int)difftime(time_now, time_start);
            randGen.SetSeed(microseconds + time_diff + RandomSeed*100);
        }
        ListIndex = 0;
        A_Num = randGen.Integer(3) + 1;
        for (UInt_t UIntI = 0;UIntI < A_Num;UIntI++){
            PDG          .push_back(A_PDG);
            mix_px       .push_back(std::exp(randGen.Gaus(mean, sigma)) * pow(-1.0,randGen.Integer(2)) * 0.12);
            mix_py       .push_back(std::exp(randGen.Gaus(mean, sigma)) * pow(-1.0,randGen.Integer(2)) * 0.12);
            mix_pz       .push_back(std::exp(randGen.Gaus(mean, sigma)) * pow(-1.0,randGen.Integer(2)) * 0.12);
            InvariantMass.push_back(AMass);
            ParentSta    .push_back(ListIndex);
            ParentEnd    .push_back(ListIndex);
            ParentList   .push_back(ListIndex);
            ListIndex++;
            // if (pow(mix_px.at(ListIndex-1)*mix_px.at(ListIndex-1) + mix_py.at(ListIndex-1)*mix_py.at(ListIndex-1),0.5) < 0.2) InvariantMass.at(ListIndex-1) = 5000;
        }
        B_Num = 1;
        for (UInt_t UIntI = 0;UIntI < B_Num;UIntI++){
            PDG          .push_back(B_PDG);
            mix_px       .push_back(std::exp(randGen.Gaus(mean, sigma)) * pow(-1.0,randGen.Integer(2)) * 0.08);
            mix_py       .push_back(std::exp(randGen.Gaus(mean, sigma)) * pow(-1.0,randGen.Integer(2)) * 0.08);
            mix_pz       .push_back(std::exp(randGen.Gaus(mean, sigma)) * pow(-1.0,randGen.Integer(2)) * 0.08);
            InvariantMass.push_back(BMass);
            ParentSta    .push_back(ListIndex);
            ParentEnd    .push_back(ListIndex);
            ParentList   .push_back(ListIndex);
            ListIndex++;
        }
        PDGMult = PDG.size();
        Nch = randGen.Integer(90) + 9;
        PVz = -45.0 + ((randGen.Integer(10000)*1.0)/10000) * 100;

        // 制造相互作用
        // BPx = mix_px.at(A_Num);// in case B_Num = 1;
        // BPy = mix_py.at(A_Num);// in case B_Num = 1;
        // BPz = mix_pz.at(A_Num);// in case B_Num = 1;
        // for (UInt_t UIntI = 0;UIntI < A_Num;UIntI++){
        //     if (randGen.Uniform(0.0, 100.0) < 97.0) continue;
        //     APx = mix_px.at(UIntI);
        //     APy = mix_py.at(UIntI);
        //     APz = mix_pz.at(UIntI);
        //     TempF.clear();TempF = CAB(APx , APy , APz , BPx , BPy , BPz , randGen.Uniform(0.0, 1.0));
        //     mix_px.at(UIntI) = TempF[0];
        //     mix_py.at(UIntI) = TempF[1];
        //     mix_pz.at(UIntI) = TempF[2];
        // }

        // 产生截断
        // for (UInt_t UIntI = 0;UIntI < A_Num;UIntI++){
        //     APx = mix_px.at(UIntI);
        //     APy = mix_py.at(UIntI);
        //     APz = mix_pz.at(UIntI);
        //     tEnergy = pow(APx*APx + APy*APy + APz*APz + AMass*AMass,0.5);
        //     ARap = 0.5*log((tEnergy+APz)/(tEnergy-APz));
        //     if (fabs(randGen.Gaus(0, fabs(1/ARap) + 0.01)) < 1) continue;
        //     InvariantMass.at(ListIndex-1) = 5000;
        // }
        for (UInt_t UIntI = 0;UIntI < A_Num;UIntI++){
            APx = mix_px.at(UIntI);
            APy = mix_py.at(UIntI);
            APz = mix_pz.at(UIntI);
            APt = pow(APx*APx + APy*APy,0.5);
            AMon = pow(APt*APt + APz*APz,0.5);
            ATheta = acos(APz/AMon);
            ATan = tan(0.5*ATheta);
            AEta = (ATan > 0) ? -1.0*log(ATan) : log(-1.0*ATan);
            if (fabs(randGen.Gaus(0 , 10)) < (fabs(AEta)+0.01)) continue;
            InvariantMass.at(ListIndex-1) = 5000;
        }

        // 随机产生事件结束

        if ((EntriesID+1)%200 == 0) {
            cout<<EntriesID<<"/"<<nentries<<" Finished, used ";
            time(&time_now);
            int time_diff = (int)difftime(time_now, time_start);
            cout << time_diff/60 << "min " << time_diff%60 << "s"<<endl;
        }

        A_Px   .resize(0);   B_Px.resize(0);
        A_Py   .resize(0);   B_Py.resize(0);
        A_Pz   .resize(0);   B_Pz.resize(0);
        A_ParID.resize(0);B_ParID.resize(0);
        A_Kind .resize(0); B_Kind.resize(0);
        A_Rap  .resize(0);  B_Rap.resize(0);
        A_IfRecord.resize(0);B_IfRecord.resize(0);
        C_ParID.resize(0);

        for (int j=0;j<PDGMult;j++){
            if (PDG.at(j) == A_PDG) {
                if ( PatternID == Pattern ) {
                    if      (fabs(InvariantMass.at(j) - AMass) <= 3*AMassSigma) {A_Kind.push_back(0);}
                    // else if (fabs(InvariantMass.at(j) - AMass) <= 6*AMassSigma) {A_Kind.push_back(1);}
                    else{continue;}
                }
                else {
                    if      (fabs(InvariantMass.at(j) - AMass) <= 3*AMassSigma) {
                        if ((PatternID == 0) || (PatternID == 1)) {A_Kind.push_back(0);}
                        else {continue;}
                    }
                    else if (fabs(InvariantMass.at(j) - AMass) <= 6*AMassSigma) {
                        if ((PatternID == 2))                     {A_Kind.push_back(1);}
                        else {continue;}
                    }
                }
                A_Px.push_back(mix_px.at(j));
                A_Py.push_back(mix_py.at(j));
                A_Pz.push_back(mix_pz.at(j));
                A_IfRecord.push_back(1);
                Temp.clear();Temp.push_back(j);
                for (int k=ParentSta.at(j);k<=ParentEnd.at(j);k++){
                    Temp.push_back(ParentList.at(k));
                }
                A_ParID.push_back(Temp);
                tEnergy = pow(pow(mix_px.at(j),2) + pow(mix_py.at(j),2) + pow(mix_pz.at(j),2) + pow(AMass,2),0.5);
                A_Rap.push_back(0.5*log((tEnergy+mix_pz.at(j))/(tEnergy-mix_pz.at(j))));
                H_A_Pt ->Fill(pow(pow(mix_px.at(j),2) + pow(mix_py.at(j),2),0.5));
                H_A_P  ->Fill(pow(pow(mix_px.at(j),2) + pow(mix_py.at(j),2) + pow(mix_pz.at(j),2),0.5));
                H_A_Rap->Fill(0.5*log((tEnergy+mix_pz.at(j))/(tEnergy-mix_pz.at(j))));
            }
            else if (PDG.at(j) == B_PDG) {
                if ( PatternID == Pattern ) {
                    if      (fabs(InvariantMass.at(j) - BMass) <= 3*BMassSigma) {B_Kind.push_back(0);}
                    // else if (fabs(InvariantMass.at(j) - BMass) <= 6*BMassSigma) {B_Kind.push_back(1);}
                    else{continue;}
                }
                else {
                    if      (fabs(InvariantMass.at(j) - BMass) <= 3*BMassSigma) {
                        if ((PatternID == 0) || (PatternID == 2)) {B_Kind.push_back(0);}
                        else {continue;}
                    }
                    else if (fabs(InvariantMass.at(j) - BMass) <= 6*BMassSigma) {
                        if (PatternID == 1)                       {B_Kind.push_back(1);}
                        else {continue;}
                    }
                }
                B_Px.push_back(mix_px.at(j));
                B_Py.push_back(mix_py.at(j));
                B_Pz.push_back(mix_pz.at(j));
                B_IfRecord.push_back(1);
                Temp.clear();Temp.push_back(j);
                for (int k=ParentSta.at(j);k<=ParentEnd.at(j);k++){
                    Temp.push_back(ParentList.at(k));
                }
                B_ParID.push_back(Temp);
                tEnergy = pow(pow(mix_px.at(j),2) + pow(mix_py.at(j),2) + pow(mix_pz.at(j),2) + pow(BMass,2),0.5);
                B_Rap.push_back(0.5*log((tEnergy+mix_pz.at(j))/(tEnergy-mix_pz.at(j))));
                H_B_Pt ->Fill(pow(pow(mix_px.at(j),2) + pow(mix_py.at(j),2),0.5));
                H_B_P  ->Fill(pow(pow(mix_px.at(j),2) + pow(mix_py.at(j),2) + pow(mix_pz.at(j),2),0.5));
                H_B_Rap->Fill(0.5*log((tEnergy+mix_pz.at(j))/(tEnergy-mix_pz.at(j))));
            }
            else{
                for (int l = 0;l < FeedDownNum;l++) {
                    if ( abs(PDG.at(j)) == FeedDown[l] ) {
                        if ((fabs(InvariantMass.at(j) - CMass.at(l)) > 3*CMassSigma.at(l))) continue;
                        Temp.clear();Temp.push_back(j);
                        for (int k=ParentSta.at(j);k<=ParentEnd.at(j);k++){
                            Temp.push_back(ParentList.at(k));
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
                if (IfCommonElement(A_ParID.at(Aid) , B_ParID.at(Bid))){
                    A_IfRecord.at(Aid) = 0;
                }
            }
        }
        
        // 如果A、B与C有血缘关系，不记录A和B
        for (int Aid = 0;Aid < A_Px.size();Aid++) {
            for (int Cid = 0;Cid < C_ParID.size();Cid++) {
                if (IfCommonElement(A_ParID.at(Aid) , C_ParID.at(Cid))) {
                    A_IfRecord.at(Aid) = 0;
                }
            }
        }
        for (int Bid = 0;Bid < B_Px.size();Bid++) {
            for (int Cid = 0;Cid < C_ParID.size();Cid++) {
                if (IfCommonElement(B_ParID.at(Bid) , C_ParID.at(Cid))) {
                    B_IfRecord.at(Bid) = 0;
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
                                        }
                                        else{
                                            H_Kstar[CenIndex][i][j][Aid][Bid]->Fill(0.5 * (p2 - p1).Rho());
                                            H_dRap [CenIndex][i][j][Aid][Bid]->Fill(Mix_A_Rap[CenIndex][i][j][Aid][Bid].at(Aindex) - Mix_B_Rap[CenIndex][i][j][Aid][Bid].at(Bindex));
                                            H_dPt  [CenIndex][i][j][Aid][Bid]->Fill(fabs(pow(APx*APx + APy*APy , 0.5) - pow(BPx*BPx + BPy*BPy , 0.5)));
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
    
    TString OutputFileName = OutMidName;
    OutputFileName += "H_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileA = new TFile(OutputFileName, "RECREATE");
    folder_kStar = fileA->mkdir("kStar");
    folder_dRap  = fileA->mkdir("dRap");
    folder_dPt   = fileA->mkdir("dPt");

    cout << "#######################" << endl;
    cout << "# Calculating Summary #" << endl;
    cout << "#######################" << endl;
    H_A_Pt ->Write();
    H_B_Pt ->Write();
    H_A_P  ->Write();
    H_B_P  ->Write();
    H_A_Rap->Write();
    H_B_Rap->Write();
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
                    }
                }
            }
        }
    }
    fileA->Close();
}