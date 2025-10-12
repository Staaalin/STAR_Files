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
#include "TString.h"
// #endif
#include <iostream>
#include <map>
#include <stdio.h>
using namespace std;


// 定义粒子结构体
struct Particle {
    float px;       // x方向动量
    float py;       // y方向动量
    float pz;       // z方向动量
    float mass;     // 质量
    float eta;      // 赝快度
    float y;        // 快度
    float pt;       // 横向动量
    bool  IsRecord; // 是否被记录
    int   TreeID;   // ID in one event
    std::vector<int>   ParentID; // Parent Particle ID in one event
    
    // 构造函数
    Particle(float _px, float _py, float _pz, float _mass, int _TreeID) 
        : px(_px), py(_py), pz(_pz), mass(_mass), TreeID(_TreeID) {
        // 计算赝快度、快度和横向动量
        pt = sqrt(px*px + py*py);
        float p = sqrt(px*px + py*py + pz*pz);
        float E = sqrt(p*p+mass*mass);
        eta = -1.0*log(tan(0.5*(acos(pz/p))));
        y = 0.5 * log((E + pz) / (E - pz));
        IsRecord = false;
        TreeID = 0;
    }
    
    // 计算能量
    float energy() const {
        return sqrt(px*px + py*py + pz*pz + mass*mass);
    }
    
    // 转换为四动量
    TLorentzVector lorentzVector() const {
        return TLorentzVector(px, py, pz, energy());
    }
};

#if defined(__CINT__) || defined(__CLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Particle+;
#pragma link C++ class std::vector<Particle>+;
#endif


#define Pi 3.1415926535898
#define HowMuchEventMixing 10

// int CentralityBin[] = {0 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 55 , 60 , 65 , 70 , 75 , 80 , 85 , 90 , 95 , 100};// %
int CentralityBin[] = {0 , 10 , 20 , 30 , 40 , 50 , 60 , 70 , 80 , 90 , 100};// %
const float PVzBin[] = {-45.0 , -35.0 , -25.0 , -15.0 , -5.0 , 5.0 , 15.0 , 25.0 , 35.0 , 45.0 , 55.0}; // Primary Vertex Z (cm) d+Au@200 GeV RUN 21 : -45 ~ 55 cm
const float yBin[]  = {-1000.0 , 0.0 , 1000.0}; // B_y
const float AyCut[] = {-1000.0   ,     1000.0}; // A_y
int FeedDown[] = { 0 };
const float EtaCut[] = {-1 , 1}; // EtaCut for both A and B
const float MassSigmaWidth = 3.0;
// int MultRemain = 1+197;// p+Au
int MultRemain = 2+197;// d+Au
// int MultRemain = 197+197;// Au+Au

const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1
const Int_t PVzBinNum = sizeof(PVzBin)/sizeof(PVzBin[0]) - 1; // -1
const Int_t yBinNum = sizeof(yBin)/sizeof(yBin[0]) - 1; // -1
const Int_t FeedDownNum = sizeof(FeedDown)/sizeof(FeedDown[0]);

void print(std::vector<int> Temp);
void print(std::vector<float> Temp);
std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize, TString DataName);
bool IfInVector(int Num , std::vector<int> V);
std::vector<int> GetDaughterPDGLit(int ID);
Double_t massList(int PID, TString DataName);
Double_t massListSigma(int PID, TString DataName);
float* GetPairMassAndKstar(float p1x , float p1y , float p1z , float p2x , float p2y , float p2z , float AMass , float BMass);
float* GetPairMassAndKstar(float p1x , float p1y , float p1z , float p2x , float p2y , float p2z , float p3x , float p3y , float p3z , float AMass , float BMass , float CMass);
float* GetPairMassAndKstar(float p1x , float p1y , float p1z , float p2x , float p2y , float p2z , float p3x , float p3y , float p3z , float p4x , float p4y , float p4z , float AMass , float BMass , float CMass , float DMass);
float CenCorr(float Vz, TString DataName);


// 定义事件结构体
struct Event {
    int eventID;                    // 事件ID
    std::vector<Particle> A_particles;  // A类粒子 主粒子
    std::vector<Particle> B_particles;  // B类粒子
    
    // 构造函数
    Event(int _eventID) 
        : eventID(_eventID) {}
};

void print(Event Temp);

void MM(TString MidName,TString DataName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,TString OutMidName,
        int A_PDG,int B_PDG,int Mode = 0,int SP_ME = 0, // Mode = 0: PDGMult 为vector长度 ; SP_Me : if turn on cut of Splite & Merge Effect ; Purity_MC : if turn on 
        int CutID = 0) // 0: default ; 1: nHit ; 2: PVz ; 3: TPC_nSigma ; 4: DCA
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
        std::vector<Float_t> *nHitsFit           = nullptr;
        std::vector<Float_t> *nHitsMax           = nullptr;
        std::vector<int>     *ParentList         = nullptr;
        std::vector<int>     *ParentSta          = nullptr;
        std::vector<int>     *ParentEnd          = nullptr;
        std::vector<int>     *SE_ParentList      = nullptr;
        std::vector<int>     *SE_ParentSta       = nullptr;
        std::vector<int>     *SE_ParentEnd       = nullptr;
        std::vector<int>     *ME_ParentList      = nullptr;
        std::vector<int>     *ME_ParentSta       = nullptr;
        std::vector<int>     *ME_ParentEnd       = nullptr;

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
        TBranch *bnHitsFit                       = nullptr;
        TBranch *bnHitsMax                       = nullptr;
        TBranch *bParentList                     = nullptr;
        TBranch *bParentSta                      = nullptr;
        TBranch *bParentEnd                      = nullptr;
        TBranch *bSE_ParentList                  = nullptr;
        TBranch *bSE_ParentSta                   = nullptr;
        TBranch *bSE_ParentEnd                   = nullptr;
        TBranch *bME_ParentList                  = nullptr;
        TBranch *bME_ParentSta                   = nullptr;
        TBranch *bME_ParentEnd                   = nullptr;
    
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
            std::vector<Float_t> *nHitsFit           = NULL;
            std::vector<Float_t> *nHitsMax           = NULL;
            std::vector<int>     *ParentList         = NULL;
            std::vector<int>     *ParentSta          = NULL;
            std::vector<int>     *ParentEnd          = NULL;
            std::vector<int>     *SE_ParentList      = NULL;
            std::vector<int>     *SE_ParentSta       = NULL;
            std::vector<int>     *SE_ParentEnd       = NULL;
            std::vector<int>     *ME_ParentList      = NULL;
            std::vector<int>     *ME_ParentSta       = NULL;
            std::vector<int>     *ME_ParentEnd       = NULL;

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
            TBranch *bnHitsFit                       = NULL;
            TBranch *bnHitsMax                       = NULL;
            TBranch *bParentList                     = NULL;
            TBranch *bParentSta                      = NULL;
            TBranch *bParentEnd                      = NULL;
            TBranch *bSE_ParentList                  = NULL;
            TBranch *bSE_ParentSta                   = NULL;
            TBranch *bSE_ParentEnd                   = NULL;
            TBranch *bME_ParentList                  = NULL;
            TBranch *bME_ParentSta                   = NULL;
            TBranch *bME_ParentEnd                   = NULL;

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
            std::vector<Float_t> *nHitsFit           = 0;
            std::vector<Float_t> *nHitsMax           = 0;
            std::vector<int>     *ParentList         = 0;
            std::vector<int>     *ParentSta          = 0;
            std::vector<int>     *ParentEnd          = 0;
            std::vector<int>     *SE_ParentList      = 0;
            std::vector<int>     *SE_ParentSta       = 0;
            std::vector<int>     *SE_ParentEnd       = 0;
            std::vector<int>     *ME_ParentList      = 0;
            std::vector<int>     *ME_ParentSta       = 0;
            std::vector<int>     *ME_ParentEnd       = 0;
    
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
            TBranch *bnHitsFit                       = 0;
            TBranch *bnHitsMax                       = 0;
            TBranch *bParentList                     = 0;
            TBranch *bParentSta                      = 0;
            TBranch *bParentEnd                      = 0;
            TBranch *bSE_ParentList                  = 0;
            TBranch *bSE_ParentSta                   = 0;
            TBranch *bSE_ParentEnd                   = 0;
            TBranch *bME_ParentList                  = 0;
            TBranch *bME_ParentSta                   = 0;
            TBranch *bME_ParentEnd                   = 0;

        #endif
    #endif

    bool IfRecord = true , IfRemoveFeedPair = false , IfRemoveSpliteMerge = false , IfRemoveLownHits = false , IfRemoveHighPVz = false , IfRemoveHighTPCsigma = false , IfCutHighDCA = false;
    float kstar, drap , dpt;
    std::vector<float> kstar_Store , drap_Store , dpt_Store , mass_Store;
    std::vector<int>   IfRecorded;
    bool Is2Body = true;
    float NNch , Eta;
    TString TreeName = "hadronTree";

    TVector3 BetaTemp;
    // ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p1 , p2 , p3 , p4 , p5;
    TLorentzVector p1 , p2 , p3;
    TVector3 BV;
    std::vector<int> Temp;
    std::vector<float> MotherMass , MotherMassSigma;
    int AccumSameNum;
    TRandom3 rng(0);
    int i , j , k , l , m , n;
    int RapIndex , CenIndex , PVzIndex;
    std::vector<int> MatchedRap;
    int Aid , Bid , Cid;
    float APx  , BPx ;
    float APy  , BPy ;
    float APz  , BPz ;
    float APt  , BPt ;
    float ARap , BRap;
    float BMass = massList(B_PDG, DataName)           , AMass = massList(A_PDG, DataName)          ;
    float BMassSigma = massListSigma(B_PDG, DataName) , AMassSigma = massListSigma(A_PDG, DataName);
    std::vector<std::vector<int> > C_ParID;

    //                                    centrality    A_Rapidity   PrimaryVertex
    std::vector<Event>    EventPool         [50]           [50]          [50];
    std::vector<Particle> A_Array                          [50]              , B_Array;
    std::vector<Particle> A_List                           [50]              , B_List ;
    TH1F                 *H_Kstar           [50]           [50]          [50];
    TH1F                 *H_ALL_Kstar                      [50]     ;
    TH1F                 *H_Mix_Kstar       [50]           [50]          [50];
    TH1F                 *H_ALL_Mix_Kstar                  [50]     ;
    TH1F                 *H_Tra_Kstar       [50]           [50]          [50];
    TH1F                 *H_ALL_Tra_Kstar                  [50]     ;
    TH1F                 *H_dRap            [50]           [50]          [50];
    TH1F                 *H_ALL_dRap                       [50]     ;
    TH1F                 *H_Mix_dRap        [50]           [50]          [50];
    TH1F                 *H_ALL_Mix_dRap                   [50]     ;
    TH1F                 *H_Tra_dRap        [50]           [50]          [50];
    TH1F                 *H_ALL_Tra_dRap                   [50]     ;
    TH1F                 *H_dPt             [50]           [50]          [50];
    TH1F                 *H_ALL_dPt                        [50]     ;
    TH1F                 *H_Mix_dPt         [50]           [50]          [50];
    TH1F                 *H_ALL_Mix_dPt                    [50]     ;
    TH1F                 *H_Tra_dPt         [50]           [50]          [50];
    TH1F                 *H_ALL_Tra_dPt                    [50]     ;
    TH1F                 *H_ALL_Mass                       [50]     ;
    TH1F                 *H_ALL_Mix_Mass                   [50]     ;
    TH1F                 *H_ALL_Tra_Mass                   [50]     ;
    TH1F                 *H_Rap_A           [50]           [50]          [50];
    TH1F                 *H_ALL_Rap_A                      [50]     ;
    TH1F                 *H_Rap_K_A         [50]           [50]          [50];
    TH1F                 *H_ALL_Rap_K_A                    [50]     ;
    TH1F                 *H_Rap_B           [50]           [50]          [50];
    TH1F                 *H_ALL_Rap_B                      [50]     ;
    TH1F                 *H_Rap_K_B         [50]                         [50];
    TH1F                 *H_ALL_Rap_K_B                             ;
    Particle              A(0,0,0,0,0), B(0,0,0,0,0), C(0,0,0,0,0), D(0,0,0,0,0);
    Event                 TempEvent(0);

    int kStarBinNum = 400;
    float kStarSta = 0 , kStarEnd = 8;
    
    int dRapBinNum = 1000;
    float dRapSta = -10 , dRapEnd = 10;
    
    int SRapBinNum = 1000;
    float SRapSta = -10 , SRapEnd = 10;
    
    int dPtBinNum = 200;
    float dPtSta = 0 , dPtEnd = 10;
    
    int MBinNum = 1000 , MBinPar = 100;
    float MSta = floor((AMass + BMass)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;
    cout<<"Mass Region: [ "<<MSta<<" , "<<MEnd<<" ], BinNum = "<<MBinNum<<". "<<endl;

    for (RapIndex=0;RapIndex<yBinNum;RapIndex++) {
        for (CenIndex=0;CenIndex<CentralityBinNum;CenIndex++) {
            for (PVzIndex=0;PVzIndex<PVzBinNum;PVzIndex++) {
                H_Kstar           [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Kstar_%d_%d_%d"       ,CenIndex,RapIndex,PVzIndex),Form("Kstar, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),kStarBinNum,kStarSta,kStarEnd);
                H_Mix_Kstar       [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Mix_Kstar_%d_%d_%d"   ,CenIndex,RapIndex,PVzIndex),Form("Mix Kstar, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),kStarBinNum,kStarSta,kStarEnd);
                H_Tra_Kstar       [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Tra_Kstar_%d_%d_%d"   ,CenIndex,RapIndex,PVzIndex),Form("Tra Kstar, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),kStarBinNum,kStarSta,kStarEnd);
                H_dRap            [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_dRap_%d_%d_%d"        ,CenIndex,RapIndex,PVzIndex),Form("dRap, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"        ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dRapBinNum,dRapSta,dRapEnd);
                H_Mix_dRap        [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Mix_dRap_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("Mix dRap, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"    ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dRapBinNum,dRapSta,dRapEnd);
                H_Tra_dRap        [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Tra_dRap_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("Tra dRap, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"    ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dRapBinNum,dRapSta,dRapEnd);
                H_dPt             [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_dPt_%d_%d_%d"         ,CenIndex,RapIndex,PVzIndex),Form("dPt, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"         ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dPtBinNum,dPtSta,dPtEnd);
                H_Mix_dPt         [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Mix_dPt_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("Mix dPt, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"     ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dPtBinNum,dPtSta,dPtEnd);
                H_Tra_dPt         [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Tra_dPt_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("Tra dPt, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"     ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dPtBinNum,dPtSta,dPtEnd);
                H_Rap_A           [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Rap_A_%d_%d_%d"       ,CenIndex,RapIndex,PVzIndex),Form("A dN/dy, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"     ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dRapBinNum,dRapSta,dRapEnd);
                H_Rap_K_A         [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Rap_K_A_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("A dN/dy, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"     ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dRapBinNum,dRapSta,dRapEnd);
                H_Rap_B           [CenIndex] [RapIndex] [PVzIndex] = new TH1F(Form("H_Rap_B_%d_%d_%d"       ,CenIndex,RapIndex,PVzIndex),Form("B dN/dy, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"     ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dRapBinNum,dRapSta,dRapEnd);
                if ((RapIndex == 0)) {
                    H_Rap_K_B     [CenIndex]            [PVzIndex] = new TH1F(Form("H_Rap_K_B_%d "       ,CenIndex         ),Form("B dN/dy, [%d,%d]/100"                ,CentralityBin[CenIndex],CentralityBin[CenIndex+1]                                ),dRapBinNum,dRapSta,dRapEnd);
                }
            }
        }
        H_ALL_Kstar                      [RapIndex] = new TH1F(Form("H_ALL_Kstar_%d"      ,         RapIndex),Form("ALL Kstar, %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),kStarBinNum,kStarSta,kStarEnd);
        H_ALL_Mix_Kstar                  [RapIndex] = new TH1F(Form("H_ALL_Mix_Kstar_%d"  ,         RapIndex),Form("ALL Mix_Kstar, %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),kStarBinNum,kStarSta,kStarEnd);
        H_ALL_Tra_Kstar                  [RapIndex] = new TH1F(Form("H_ALL_Tra_Kstar_%d"  ,         RapIndex),Form("ALL Tra_Kstar, %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),kStarBinNum,kStarSta,kStarEnd);
        H_ALL_dRap                       [RapIndex] = new TH1F(Form("H_ALL_dRap_%d"       ,         RapIndex),Form("ALL dRap, %f<A_y<%f"       ,yBin[RapIndex],yBin[RapIndex+1]),dRapBinNum,dRapSta,dRapEnd);
        H_ALL_Mix_dRap                   [RapIndex] = new TH1F(Form("H_ALL_Mix_dRap_%d"   ,         RapIndex),Form("ALL Mix_dRap, %f<A_y<%f"   ,yBin[RapIndex],yBin[RapIndex+1]),dRapBinNum,dRapSta,dRapEnd);
        H_ALL_Tra_dRap                   [RapIndex] = new TH1F(Form("H_ALL_Tra_dRap_%d"   ,         RapIndex),Form("ALL Tra_dRap, %f<A_y<%f"   ,yBin[RapIndex],yBin[RapIndex+1]),dRapBinNum,dRapSta,dRapEnd);
        H_ALL_dPt                        [RapIndex] = new TH1F(Form("H_ALL_dPt_%d"        ,         RapIndex),Form("ALL dPt, %f<A_y<%f"        ,yBin[RapIndex],yBin[RapIndex+1]),dPtBinNum,dPtSta,dPtEnd);
        H_ALL_Mix_dPt                    [RapIndex] = new TH1F(Form("H_ALL_Mix_dPt_%d"    ,         RapIndex),Form("ALL Mix_dPt, %f<A_y<%f"    ,yBin[RapIndex],yBin[RapIndex+1]),dPtBinNum,dPtSta,dPtEnd);
        H_ALL_Tra_dPt                    [RapIndex] = new TH1F(Form("H_ALL_Tra_dPt_%d"    ,         RapIndex),Form("ALL Tra_dPt, %f<A_y<%f"    ,yBin[RapIndex],yBin[RapIndex+1]),dPtBinNum,dPtSta,dPtEnd);
        H_ALL_Mass                       [RapIndex] = new TH1F(Form("H_ALL_Mass_%d"       ,         RapIndex),Form("ALL Mass, %f<A_y<%f"       ,yBin[RapIndex],yBin[RapIndex+1]),MBinNum,MSta,MEnd);
        H_ALL_Mix_Mass                   [RapIndex] = new TH1F(Form("H_ALL_Mix_Mass_%d"   ,         RapIndex),Form("ALL Mix_Mass, %f<A_y<%f"   ,yBin[RapIndex],yBin[RapIndex+1]),MBinNum,MSta,MEnd);
        H_ALL_Tra_Mass                   [RapIndex] = new TH1F(Form("H_ALL_Tra_Mass_%d"   ,         RapIndex),Form("ALL Tra_Mass, %f<A_y<%f"   ,yBin[RapIndex],yBin[RapIndex+1]),MBinNum,MSta,MEnd);
        H_ALL_Rap_A                      [RapIndex] = new TH1F(Form("H_ALL_Rap_A_%d"      ,         RapIndex),Form("A dN/dy, %f<A_y<%f"        ,yBin[RapIndex],yBin[RapIndex+1]),dRapBinNum,dRapSta,dRapEnd);
        H_ALL_Rap_K_A                    [RapIndex] = new TH1F(Form("H_ALL_Rap_K_A_%d"    ,         RapIndex),Form("A dN/dy, %f<A_y<%f"        ,yBin[RapIndex],yBin[RapIndex+1]),dRapBinNum,dRapSta,dRapEnd);
        H_ALL_Rap_B                      [RapIndex] = new TH1F(Form("H_ALL_Rap_B_%d"      ,         RapIndex),Form("B dN/dy, %f<A_y<%f"        ,yBin[RapIndex],yBin[RapIndex+1]),dRapBinNum,dRapSta,dRapEnd);
        if (RapIndex == 0) {
            H_ALL_Rap_K_B                           = new TH1F(     "H_ALL_Rap_K_B"                          ,     "B dN/dy"                                                    ,dRapBinNum,dRapSta,dRapEnd);
        }
    }
    
    std::vector<int> NchList = GetNchList(CentralityBin , CentralityBinNum+1, DataName);     // centrality
    cout<<"NchList = ";
    print(NchList);
    cout<<" "<<endl;

    
    for (i = 0;i < FeedDownNum;i++){
        if (abs(FeedDown[i]) == A_PDG) {
            FeedDown[i] = 0;
            MotherMass.push_back(-100);
            MotherMassSigma.push_back(-1);
            continue;
        }
        if (abs(FeedDown[i]) == B_PDG) {
            FeedDown[i] = 0;
            MotherMass.push_back(-100);
            MotherMassSigma.push_back(-1);
            continue;
        }
        MotherMass.push_back(massList(FeedDown[i], DataName));
        MotherMassSigma.push_back(massListSigma(FeedDown[i], DataName));
    }
    cout<<"MotherMass = ";print(MotherMass);
    cout<<"MotherMassSigma = ";print(MotherMassSigma);

    for (i=0;i<FeedDownNum;i++) {
        if ( IfInVector(A_PDG , GetDaughterPDGLit(FeedDown[i])) && IfInVector(B_PDG , GetDaughterPDGLit(FeedDown[i])) ) IfRemoveFeedPair = true;
    }
    
    if (SP_ME == 1) IfRemoveSpliteMerge = true;
    if (CutID == 1) IfRemoveLownHits = true;
    if (CutID == 2) IfRemoveHighPVz = true;
    if (CutID == 3) IfRemoveHighTPCsigma = true;
    if (CutID == 4) IfCutHighDCA = true;

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
    if(IfCutHighDCA) hadronTree->SetBranchAddress("dcatopv"      ,&dcatopv      ,&bdcatopv      );
    // hadronTree->SetBranchAddress("nSigmaProton" ,&nSigmaProton ,&bnSigmaProton );
    // hadronTree->SetBranchAddress("nSigmaPion"   ,&nSigmaPion   ,&bnSigmaPion   );
    if (IfRemoveHighTPCsigma && ((abs(A_PDG) == 321) || (abs(B_PDG) == 321))){
        hadronTree->SetBranchAddress("nSigmaKaon"   ,&nSigmaKaon   ,&bnSigmaKaon   );
    }
    hadronTree->SetBranchAddress("InvariantMass",&InvariantMass,&bInvariantMass);
    // hadronTree->SetBranchAddress("Decay_Length" ,&Decay_Length ,&bDecay_Length );
    // hadronTree->SetBranchAddress("Chi2"         ,&Chi2         ,&bChi2         );
    if (IfRemoveLownHits) {
        hadronTree->SetBranchAddress("nHitsFit"     ,&nHitsFit     ,&bnHitsFit     );
        hadronTree->SetBranchAddress("nHitsMax"     ,&nHitsMax     ,&bnHitsMax     );
    }
    hadronTree->SetBranchAddress("ParentList"   ,&ParentList   ,&bParentList   );
    hadronTree->SetBranchAddress("ParentSta"    ,&ParentSta    ,&bParentSta    );
    hadronTree->SetBranchAddress("ParentEnd"    ,&ParentEnd    ,&bParentEnd    );
    if (IfRemoveSpliteMerge) {
        hadronTree->SetBranchAddress("SE_ParentList",&SE_ParentList,&bSE_ParentList   );
        hadronTree->SetBranchAddress("SE_ParentSta" ,&SE_ParentSta ,&bSE_ParentSta    );
        hadronTree->SetBranchAddress("SE_ParentEnd" ,&SE_ParentEnd ,&bSE_ParentEnd    );
        hadronTree->SetBranchAddress("ME_ParentList",&ME_ParentList,&bME_ParentList   );
        hadronTree->SetBranchAddress("ME_ParentSta" ,&ME_ParentSta ,&bME_ParentSta    );
        hadronTree->SetBranchAddress("ME_ParentEnd" ,&ME_ParentEnd ,&bME_ParentEnd    );
    }

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;

    time_t time_start;
    time_t time_now;
    time(&time_start);
    clock_t Tstart = clock();
    for (int EntriesID = 0 ; EntriesID < nentries ; EntriesID++){
        hadronTree->GetEntry(EntriesID);
        if ((EntriesID+1)%20000 == 0) {
            time(&time_now);
            int time_diff = (int)difftime(time_now, time_start);
            cout << time_diff/60 << "min " << time_diff%60 << "s: ";
            long long microseconds = (clock() - Tstart)/10000;
            std::cout << "Microseconds: " << microseconds << "  ";
            cout<<"Calculating Event "<<(EntriesID+1)<<"/"<<nentries<<endl;
            Tstart = clock();
        }

        if (IfRemoveHighPVz) {
            if (!((-25.0 <= PVz) && (PVz < 25.0))) continue;
        }

        C_ParID.clear();
        B_List.clear();
        TempEvent.eventID = EntriesID;
        TempEvent.A_particles.clear();
        TempEvent.B_particles.clear();
        for (i=0;i<MatchedRap.size();i++) {
            A_Array[MatchedRap[i]].clear();
            A_List [MatchedRap[i]].clear();
        }
        MatchedRap.clear();
        // 定Centrality
        CenIndex = -1;
        for (k=0;k<CentralityBinNum;k++){
            NNch = CenCorr(PVz, DataName) * Nch;
            // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
            if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                CenIndex = k;
                break;
            }
        }
        if (CenIndex == -1) continue;
        // 定Primary Vertex Z
        PVzIndex = -1;
        for (k=0;k<PVzBinNum;k++){
            // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
            if ((PVzBin[k] <= PVz) && (PVz < PVzBin[k+1])) {
                PVzIndex = k;
                break;
            }
        }
        if (PVzIndex == -1) continue;
        // 遍历粒子，筛选A、B、C、D
        for (i=0;i<PDGMult;i++){
            if (PDG->at(i) == A_PDG) {
                if (fabs(InvariantMass->at(i) - AMass) <= MassSigmaWidth*AMassSigma) {

                    if (IfRemoveHighTPCsigma) {
                        if (abs(A_PDG) == 321) {
                            if (fabs(nSigmaKaon->at(i))>1) continue;
                        }
                    }
                    if (IfRemoveLownHits) {
                        if ((abs(A_PDG) == 321) || (abs(A_PDG) == 211) || (abs(A_PDG) == 2212)) {
                            if (nHitsFit->at(i) < 20) continue;
                        }
                    }
                    if (IfCutHighDCA) {
                        if ((abs(A_PDG) == 321) || (abs(A_PDG) == 211) || (abs(A_PDG) == 2212)) {
                            if ( (0 > dcatopv->at(i)) || (dcatopv->at(i) > 0.5)) continue;
                        }
                    }

                    A = Particle(mix_px->at(i),mix_py->at(i),mix_pz->at(i),AMass,i);
                    A.ParentID.clear();A.ParentID.push_back(i);
                    for (k=ParentSta->at(i);k<=ParentEnd->at(i);k++){
                        A.ParentID.push_back(ParentList->at(k));
                    }
                    if (IfRemoveSpliteMerge) {
                        for (k=SE_ParentSta->at(i);k<=SE_ParentEnd->at(i);k++){
                            A.ParentID.push_back(SE_ParentList->at(k));
                        }
                        for (k=ME_ParentSta->at(i);k<=ME_ParentEnd->at(i);k++){
                            A.ParentID.push_back(ME_ParentList->at(k));
                        }
                    }
                    if ((A.eta < EtaCut[0]) || (A.eta > EtaCut[1])) continue;
                    RapIndex = -1;
                    for (k=0;k<yBinNum;k++){
                        if ((yBin[k] <= A.y) && (A.y < yBin[k+1])) {
                            RapIndex = k;
                            break;
                        }
                    }
                    if (RapIndex == -1) continue;
                    A_List[RapIndex].push_back(A);
                    if (!IfInVector(RapIndex , MatchedRap)) MatchedRap.push_back(RapIndex);
                    H_Rap_K_A[CenIndex][RapIndex][PVzIndex]->Fill(A.y);
                    H_ALL_Rap_K_A      [RapIndex]          ->Fill(A.y);
                    continue;
                }
            }
            else if (PDG->at(i) == B_PDG) {
                if (fabs(InvariantMass->at(i) - BMass) <= MassSigmaWidth*BMassSigma) {

                    if (IfRemoveHighTPCsigma) {
                        if (abs(B_PDG) == 321) {
                            if (fabs(nSigmaKaon->at(i))>1) continue;
                        }
                    }
                    if (IfRemoveLownHits) {
                        if ((abs(B_PDG) == 321) || (abs(B_PDG) == 211) || (abs(B_PDG) == 2212)) {
                            if (nHitsFit->at(i) < 20) continue;
                        }
                    }
                    if (IfCutHighDCA) {
                        if ((abs(B_PDG) == 321) || (abs(B_PDG) == 211) || (abs(B_PDG) == 2212)) {
                            if ( (0 > dcatopv->at(i)) || (dcatopv->at(i) > 0.5)) continue;
                        }
                    }

                    B = Particle(mix_px->at(i),mix_py->at(i),mix_pz->at(i),BMass,i);
                    B.ParentID.clear();B.ParentID.push_back(i);
                    for (k=ParentSta->at(i);k<=ParentEnd->at(i);k++){
                        B.ParentID.push_back(ParentList->at(k));
                    }
                    if (IfRemoveSpliteMerge) {
                        for (k=SE_ParentSta->at(i);k<=SE_ParentEnd->at(i);k++){
                            B.ParentID.push_back(SE_ParentList->at(k));
                        }
                        for (k=ME_ParentSta->at(i);k<=ME_ParentEnd->at(i);k++){
                            B.ParentID.push_back(ME_ParentList->at(k));
                        }
                    }
                    if ((B.eta < EtaCut[0]) || (B.eta > EtaCut[1])) continue;
                    // TempEvent.B_particles.push_back(B);
                    B_List.push_back(B);
                    H_Rap_K_B[CenIndex][PVzIndex]->Fill(B.y);
                    H_ALL_Rap_K_B                ->Fill(B.y);
                    continue;
                }
            }
            else{
                for (l = 0;l < FeedDownNum;l++) {
                    if ( abs(PDG->at(i)) == FeedDown[l] ) {
                        if ((fabs(InvariantMass->at(i) - MotherMass.at(l)) > 3*MotherMassSigma.at(l))) continue;
                        Temp.clear();Temp.push_back(j);
                        for (k=ParentSta->at(i);k<=ParentEnd->at(i);k++){
                            Temp.push_back(ParentList->at(k));
                        }
                        C_ParID.push_back(Temp);
                        // IfFoundOmega = true;
                        // cout<<"Found Omega"<<endl;
                    }
                }
            }
        }
        // 筛选A、B粒子
        for (Bid=0;Bid<B_List.size();Bid++) {
            IfRecord = true;
            // 如果A、B有血缘关系，保留A
            for (i=0;i<MatchedRap.size();i++) {
                for (Aid=0;Aid<A_List[MatchedRap[i]].size();Aid++) {
                    if (IfInVector(B_List[Bid].TreeID , A_List[MatchedRap[i]][Aid].ParentID)) IfRecord = false;
                }
            }
            // 如果A、B与C有血缘关系，不记录A和B
            for (Cid = 0;Cid < C_ParID.size();Cid++) {
                if (IfInVector(B_List[Bid].TreeID , C_ParID.at(Cid))) IfRecord = false;
            }
            if (IfRecord) TempEvent.B_particles.push_back(B_List[Bid]);
        }
        for (i=0;i<MatchedRap.size();i++) {
            for (Aid=0;Aid<A_List[MatchedRap[i]].size();Aid++) {
                IfRecord = true;
                // 如果A、B与C有血缘关系，不记录A和B
                for (Cid = 0;Cid < C_ParID.size();Cid++) {
                    if (IfInVector(B_List[Bid].TreeID , C_ParID.at(Cid))) IfRecord = false;
                }
                if (IfRecord) A_Array[MatchedRap[i]].push_back(A_List[MatchedRap[i]][Aid]);
            }
        }
        if (TempEvent.B_particles.size() >= HowMuchEventMixing+1) continue;
        // 确保同时记录到A、B、...粒子
        if (MatchedRap.size() == 0) continue;                                                        // 有A粒子
        if (TempEvent.B_particles.size() == 0) continue;                                             // 有B粒子
        if (true)
        {
            // 填进池子 & 计算
            if ((CenIndex != -1)&&(PVzIndex != -1)) {
                for (i=0;i<MatchedRap.size();i++) {
                    TempEvent.A_particles.clear();
                    for (j=0;j<A_Array[MatchedRap[i]].size();j++){
                        TempEvent.A_particles.push_back(A_Array[MatchedRap[i]][j]);
                    }
                    RapIndex = MatchedRap[i];
                    EventPool[CenIndex][RapIndex][PVzIndex].push_back(TempEvent);
                    if (EventPool[CenIndex][RapIndex][PVzIndex].size() == HowMuchEventMixing+1) {// 池子填满，开始计算
                        // cout<<"<======= New Round =======>"<<endl;
                        // for(Aid=0;Aid<HowMuchEventMixing+1;Aid++){
                        //     cout<<"This is CenIndex = "<<CenIndex<<", RapIndex = "<<RapIndex<<endl;
                        //     print(EventPool[CenIndex][RapIndex][Aid]);
                        // }
                        if (Is2Body) {
                            for (Aid=0;Aid<HowMuchEventMixing+1;Aid++) {
                                for (Bid=0;Bid<HowMuchEventMixing+1;Bid++) {
                                    for (j=0;j<EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles.size();j++) {
                                        APx  = EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles[j].px;
                                        APy  = EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles[j].py;
                                        APz  = EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles[j].pz;
                                        ARap = EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles[j].y ;
                                        APt  = EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles[j].pt;
                                        for (k=0;k<EventPool[CenIndex][RapIndex][PVzIndex][Bid].B_particles.size();k++) {
                                            // 填入QA，参与Correlation的B粒子的dN/dy
                                            BRap = EventPool[CenIndex][RapIndex][PVzIndex][Bid].B_particles[k].y;

                                            float* MassAndKstar = GetPairMassAndKstar(APx                                                            , APy                                                            , APz                                                            , 
                                                                                      EventPool[CenIndex][RapIndex][PVzIndex][Bid].B_particles[k].px , EventPool[CenIndex][RapIndex][PVzIndex][Bid].B_particles[k].py , EventPool[CenIndex][RapIndex][PVzIndex][Bid].B_particles[k].pz , 
                                                                                      AMass , BMass);
                                            kstar = MassAndKstar[1];
                                            drap  = -ARap+BRap;
                                            dpt   = -APt +EventPool[CenIndex][RapIndex][PVzIndex][Bid].B_particles[k].pt;

                                            // Test
                                            // p2.SetXYZM(EventPool[CenIndex][RapIndex][Bid].B_particles[k].px,EventPool[CenIndex][RapIndex][Bid].B_particles[k].py,EventPool[CenIndex][RapIndex][Bid].B_particles[k].pz,BMass);
                                            // p1.SetXYZM(APx,APy,APz,AMass);
                                            // p3 = p1 + p2;
                                            // BV = -p3.BoostVector();
                                            // p1.Boost( BV);p2.Boost( BV);
                                            // cout<<"##########################################"<<endl;
                                            // cout<<"Self M = "<<MassAndKstar[0]<<", Kstar = "<<MassAndKstar[1]<<endl;
                                            // cout<<"Lonz M = "<<p1.Energy()+p2.Energy()<<", Kstar = "<<0.5 * (p2 - p1).Rho()<<endl;

                                            if (Aid == Bid) {
                                                H_Kstar          [CenIndex] [RapIndex] [PVzIndex] -> Fill(kstar);
                                                H_ALL_Kstar                 [RapIndex]            -> Fill(kstar);
                                                H_dRap           [CenIndex] [RapIndex] [PVzIndex] -> Fill(drap);
                                                H_ALL_dRap                  [RapIndex]            -> Fill(drap);
                                                H_dPt            [CenIndex] [RapIndex] [PVzIndex] -> Fill(dpt);
                                                H_ALL_dPt                   [RapIndex]            -> Fill(dpt);
                                                H_ALL_Mass                  [RapIndex]            -> Fill(MassAndKstar[0]);
                                                AccumSameNum++;
                                            }
                                            else {
                                                H_Mix_Kstar      [CenIndex] [RapIndex] [PVzIndex] -> Fill(kstar);
                                                H_ALL_Mix_Kstar             [RapIndex]            -> Fill(kstar);
                                                H_Mix_dRap       [CenIndex] [RapIndex] [PVzIndex] -> Fill(drap);
                                                H_ALL_Mix_dRap              [RapIndex]            -> Fill(drap);
                                                H_Mix_dPt        [CenIndex] [RapIndex] [PVzIndex] -> Fill(dpt);
                                                H_ALL_Mix_dPt               [RapIndex]            -> Fill(dpt);
                                                H_ALL_Mix_Mass              [RapIndex]            -> Fill(MassAndKstar[0]);
                                            }
                                            delete[] MassAndKstar;
                                            EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles[j].IsRecord = true;
                                            EventPool[CenIndex][RapIndex][PVzIndex][Bid].B_particles[k].IsRecord = true;
                                        }
                                    }
                                }
                                // if (AccumSameNum >= kstar_Store.size()) {
                                //     cout<<"kstar_Store.size() = "<<kstar_Store.size()<<", AccumSameNum = "<<AccumSameNum<<endl;
                                //     cout<<"########################################################"<<endl;
                                //     for(j=0;j<HowMuchEventMixing+1;j++) {
                                //         print(EventPool[CenIndex][RapIndex][j]);
                                //         cout<<"########################################################"<<endl;
                                //     }
                                // }
                            }
                            for (Aid=0;Aid<HowMuchEventMixing+1;Aid++) {
                                for (j=0;j<EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles.size();j++) {
                                    if (   EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles[j].IsRecord) {
                                        H_Rap_A     [CenIndex][RapIndex][PVzIndex]->Fill(EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles[j].y);
                                        H_ALL_Rap_A           [RapIndex]          ->Fill(EventPool[CenIndex][RapIndex][PVzIndex][Aid].A_particles[j].y);
                                    }
                                }
                                for (j=0;j<EventPool[CenIndex][RapIndex][PVzIndex][Aid].B_particles.size();j++) {
                                    if (   EventPool[CenIndex][RapIndex][PVzIndex][Aid].B_particles[j].IsRecord) {
                                        H_Rap_B     [CenIndex][RapIndex][PVzIndex]->Fill(EventPool[CenIndex][RapIndex][PVzIndex][Aid].B_particles[j].y);
                                        H_ALL_Rap_B           [RapIndex]          ->Fill(EventPool[CenIndex][RapIndex][PVzIndex][Aid].B_particles[j].y);
                                    }
                                }
                            }
                        }
                        EventPool[CenIndex][RapIndex][PVzIndex].clear();
                    }
                    TempEvent.A_particles.clear();
                }
            }
            // 清空buffer
            for (i=0;i<MatchedRap.size();i++) {
                A_Array[MatchedRap[i]].clear();
            }
            MatchedRap.clear();
        }
    }
    // 保存.root文件
    
    TString OutputFileName = OutMidName;
    OutputFileName += A_PDG;
    OutputFileName += "_";
    OutputFileName += B_PDG;
    OutputFileName += "_H_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileA = new TFile(OutputFileName, "RECREATE");
    TDirectory *folder_kStar     = fileA->mkdir("kStar");
    TDirectory *folder_dRap      = fileA->mkdir("dRap");
    TDirectory *folder_dPt       = fileA->mkdir("dPt");
    TDirectory *folder_Mass      = fileA->mkdir("Mass");
    TDirectory *folder_A_Num     = fileA->mkdir("A_Num");
    TDirectory *folder_B_Num     = fileA->mkdir("B_Num");
    TDirectory *folder_ALL_A_Num = fileA->mkdir("ALL_A_Num");
    TDirectory *folder_ALL_B_Num = fileA->mkdir("ALL_B_Num");
    TDirectory *ALL_kStar        = folder_kStar->mkdir("ALL");
    TDirectory *ALL_dRap         = folder_dRap ->mkdir("ALL");
    TDirectory *ALL_dPt          = folder_dPt  ->mkdir("ALL");
    TDirectory *ALL_Mass         = folder_Mass ->mkdir("ALL");
    TDirectory *ALL_A_Num        = folder_A_Num->mkdir("ALL");
    TDirectory *ALL_B_Num        = folder_B_Num->mkdir("ALL");
    TDirectory *ALL_ALL_A_Num    = folder_ALL_A_Num->mkdir("ALL");
    TDirectory *ALL_ALL_B_Num    = folder_ALL_B_Num->mkdir("ALL");
    TDirectory *Sep_kStar        = folder_kStar->mkdir("Sep");
    TDirectory *Sep_dRap         = folder_dRap ->mkdir("Sep");
    TDirectory *Sep_dPt          = folder_dPt  ->mkdir("Sep");
    TDirectory *Sep_A_Num        = folder_A_Num->mkdir("Sep");
    TDirectory *Sep_B_Num        = folder_B_Num->mkdir("Sep");
    TDirectory *Sep_ALL_A_Num    = folder_ALL_A_Num->mkdir("Sep");
    TDirectory *Sep_ALL_B_Num    = folder_ALL_B_Num->mkdir("Sep");
    for (RapIndex=0;RapIndex<yBinNum;RapIndex++) {
        for (CenIndex=0;CenIndex<CentralityBinNum;CenIndex++) {
            Sep_kStar->cd();
            H_Kstar                [CenIndex] [RapIndex] [PVzIndex] ->Write();
            H_Mix_Kstar            [CenIndex] [RapIndex] [PVzIndex] ->Write();
            H_Tra_Kstar            [CenIndex] [RapIndex] [PVzIndex] ->Write();
            Sep_dRap->cd();
            H_dRap                 [CenIndex] [RapIndex] [PVzIndex] ->Write();
            H_Mix_dRap             [CenIndex] [RapIndex] [PVzIndex] ->Write();
            H_Tra_dRap             [CenIndex] [RapIndex] [PVzIndex] ->Write();
            Sep_dPt->cd();
            H_dPt                  [CenIndex] [RapIndex] [PVzIndex] ->Write();
            H_Mix_dPt              [CenIndex] [RapIndex] [PVzIndex] ->Write();
            H_Tra_dPt              [CenIndex] [RapIndex] [PVzIndex] ->Write();
            Sep_A_Num->cd();
            H_Rap_A                [CenIndex] [RapIndex] [PVzIndex] ->Write();
            Sep_B_Num->cd();
            H_Rap_B                [CenIndex] [RapIndex] [PVzIndex] ->Write();
            Sep_ALL_A_Num->cd();
            H_Rap_K_A              [CenIndex] [RapIndex] [PVzIndex] ->Write();
            if(RapIndex==0){
                Sep_ALL_B_Num->cd();
                H_Rap_K_B          [CenIndex]            [PVzIndex] ->Write();
            }
        }
        ALL_kStar->cd();
        H_ALL_Kstar                           [RapIndex] ->Write();
        H_ALL_Mix_Kstar                       [RapIndex] ->Write();
        H_ALL_Tra_Kstar                       [RapIndex] ->Write();
        ALL_dRap->cd();     
        H_ALL_dRap                            [RapIndex] ->Write();
        H_ALL_Mix_dRap                        [RapIndex] ->Write();
        H_ALL_Tra_dRap                        [RapIndex] ->Write();
        ALL_dPt->cd();     
        H_ALL_dPt                             [RapIndex] ->Write();
        H_ALL_Mix_dPt                         [RapIndex] ->Write();
        H_ALL_Tra_dPt                         [RapIndex] ->Write();
        ALL_Mass->cd();     
        H_ALL_Mass                            [RapIndex] ->Write();
        H_ALL_Mix_Mass                        [RapIndex] ->Write();
        H_ALL_Tra_Mass                        [RapIndex] ->Write();
        ALL_A_Num->cd();
        H_ALL_Rap_A                           [RapIndex] ->Write();
        ALL_B_Num->cd();
        H_ALL_Rap_B                           [RapIndex] ->Write();
        ALL_ALL_A_Num->cd();
        H_ALL_Rap_K_A                         [RapIndex] ->Write();
    }
    ALL_ALL_B_Num->cd();
    H_ALL_Rap_K_B->Write();
    fileA->Close();
    return;
}

std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize, TString DataName)
{
    //This is 329
    std::vector<int> Result;Result.clear();
    // int CentralityListSize = sizeof(CentralityList)/sizeof(CentralityList[0]);
    if (DataName == "pAu_39_String") {
        // data from CentralityDefinition.C 
        int NchTable[21] = { 1000 , 78 , 69 , 62 , 57 , 53 , 49 , 45 , 41 , 38 , 35 , 32 , 29 , 26 , 23 , 21 , 18 , 15 , 12 , 8 ,  0};
        int CenTable[21] = {    0 ,  5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 55 , 60 , 65 , 70 , 75 , 80 , 85 , 90 ,95 ,100};
        for (int i=0;i<CentralityListSize;i++) {
            for (int j=0;j<21;j++){
                if (CenTable[j] == CentralityList[i]) {
                    Result.push_back(NchTable[j]);
                    break;
                }
            }
        }
    }
    if (DataName == "pAu_62_String") {
        // data from CentralityDefinition.C 
        int NchTable[21] = { 1000 , 100 , 88 , 79 , 72 , 66 , 60 , 55 , 51 , 46 , 42 , 38 , 34 , 31 , 28 , 25 , 22 , 18 , 15 ,  9 ,  0};
        int CenTable[21] = {    0 ,   5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 55 , 60 , 65 , 70 , 75 , 80 , 85 , 90 , 95 ,100};
        for (int i=0;i<CentralityListSize;i++) {
            for (int j=0;j<21;j++){
                if (CenTable[j] == CentralityList[i]) {
                    Result.push_back(NchTable[j]);
                    break;
                }
            }
        }
    }
    if (DataName == "dAu_200_String") {
        // data from CentralityDefinition.C 
        int NchTable[21] = { 1000 , 306 , 271 , 246 , 225 , 205 , 188 , 171 , 154 , 139 , 123 , 109 , 95 , 83 , 71 , 60 , 50 , 42 , 34 , 24 ,  0};
        int CenTable[21] = {    0 ,   5 ,  10 ,  15 ,  20 ,  25 ,  30 ,  35 ,  40 ,  45 ,  50 ,  55 , 60 , 65 , 70 , 75 , 80 , 85 , 90 , 95 ,100};
        for (int i=0;i<CentralityListSize;i++) {
            for (int j=0;j<21;j++){
                if (CenTable[j] == CentralityList[i]) {
                    Result.push_back(NchTable[j]);
                    break;
                }
            }
        }
    }
    if (DataName == "AuAu_14p6_String") {
        // data from CentralityDefinition.C 
        int NchTable[21] = { 5000 , 1738, 1551, 1387, 1242, 1112,  994,  887,  789,  700,  619,  545, 477, 415, 359, 308, 262, 220, 181, 140,  0};
        int CenTable[21] = {    0 ,   5 ,  10 ,  15 ,  20 ,  25 ,  30 ,  35 ,  40 ,  45 ,  50 ,  55 , 60 , 65 , 70 , 75 , 80 , 85 , 90 , 95 ,100};
        for (int i=0;i<CentralityListSize;i++) {
            for (int j=0;j<21;j++){
                if (CenTable[j] == CentralityList[i]) {
                    Result.push_back(NchTable[j]);
                    break;
                }
            }
        }
    }
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

float* GetPairMassAndKstar(float p1x , float p1y , float p1z , float p2x , float p2y , float p2z , float p3x , float p3y , float p3z , float AMass , float BMass , float CMass) {
    float E1 = pow(p1x*p1x+p1y*p1y+p1z*p1z+AMass*AMass,0.5);
    float E2 = pow(p2x*p2x+p2y*p2y+p2z*p2z+BMass*BMass,0.5);
    float E3 = pow(p3x*p3x+p3y*p3y+p3z*p3z+CMass*CMass,0.5);
    float Tot_E = E1+E2+E3;
    float beta[3] = { -(p1x+p2x+p3x)/Tot_E , -(p1y+p2y+p3y)/Tot_E , -(p1z+p2z+p3z)/Tot_E };
    float beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    float gamma = 1.0 / std::sqrt(1.0 - beta2);
    float gamma2 = (beta2 > 0) ? (gamma - 1.0) / beta2 : 0.0;

    float bp1 = beta[0]*p1x + beta[1]*p1y + beta[2]*p1z;
    float bp2 = beta[0]*p2x + beta[1]*p2y + beta[2]*p2z;
    float bp3 = beta[0]*p3x + beta[1]*p3y + beta[2]*p3z;

    // float New_Px = p1x + gamma2 * bp1 * beta[0] + gamma * beta[0] * E1;
    // float New_Py = p1y + gamma2 * bp1 * beta[1] + gamma * beta[1] * E1;
    // float New_Pz = p1z + gamma2 * bp1 * beta[2] + gamma * beta[2] * E1;
    float New_Px = (p1x - p2x) + gamma2 * (bp1-bp2) * beta[0] + gamma * beta[0] * (E1-E2);
    float New_Py = (p1y - p2y) + gamma2 * (bp1-bp2) * beta[1] + gamma * beta[1] * (E1-E2);
    float New_Pz = (p1z - p2z) + gamma2 * (bp1-bp2) * beta[2] + gamma * beta[2] * (E1-E2);

    float* MassAndKstar = new float[4];
    MassAndKstar[0] = (gamma * (E1 + bp1 + E2 + bp2 + E3 + bp3));
    MassAndKstar[1] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    New_Px = (p1x - p3x) + gamma2 * (bp1-bp3) * beta[0] + gamma * beta[0] * (E1-E3);
    New_Py = (p1y - p3y) + gamma2 * (bp1-bp3) * beta[1] + gamma * beta[1] * (E1-E3);
    New_Pz = (p1z - p3z) + gamma2 * (bp1-bp3) * beta[2] + gamma * beta[2] * (E1-E3);
    MassAndKstar[2] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    New_Px = (p2x - p3x) + gamma2 * (bp2-bp3) * beta[0] + gamma * beta[0] * (E2-E3);
    New_Py = (p2y - p3y) + gamma2 * (bp2-bp3) * beta[1] + gamma * beta[1] * (E2-E3);
    New_Pz = (p2z - p3z) + gamma2 * (bp2-bp3) * beta[2] + gamma * beta[2] * (E2-E3);
    MassAndKstar[3] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    return MassAndKstar;
}

float* GetPairMassAndKstar(float p1x , float p1y , float p1z , float p2x , float p2y , float p2z , float p3x , float p3y , float p3z , float p4x , float p4y , float p4z , float AMass , float BMass , float CMass , float DMass) {
    float E1 = pow(p1x*p1x+p1y*p1y+p1z*p1z+AMass*AMass,0.5);
    float E2 = pow(p2x*p2x+p2y*p2y+p2z*p2z+BMass*BMass,0.5);
    float E3 = pow(p3x*p3x+p3y*p3y+p3z*p3z+CMass*CMass,0.5);
    float E4 = pow(p4x*p4x+p4y*p4y+p4z*p4z+DMass*DMass,0.5);
    float Tot_E = E1+E2+E3+E4;
    float beta[3] = { -(p1x+p2x+p3x+p4x)/Tot_E , -(p1y+p2y+p3y+p4y)/Tot_E , -(p1z+p2z+p3z+p4z)/Tot_E };
    float beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    float gamma = 1.0 / std::sqrt(1.0 - beta2);
    float gamma2 = (beta2 > 0) ? (gamma - 1.0) / beta2 : 0.0;

    float bp1 = beta[0]*p1x + beta[1]*p1y + beta[2]*p1z;
    float bp2 = beta[0]*p2x + beta[1]*p2y + beta[2]*p2z;
    float bp3 = beta[0]*p3x + beta[1]*p3y + beta[2]*p3z;
    float bp4 = beta[0]*p4x + beta[1]*p4y + beta[2]*p4z;

    // float New_Px = p1x + gamma2 * bp1 * beta[0] + gamma * beta[0] * E1;
    // float New_Py = p1y + gamma2 * bp1 * beta[1] + gamma * beta[1] * E1;
    // float New_Pz = p1z + gamma2 * bp1 * beta[2] + gamma * beta[2] * E1;
    float New_Px = (p1x - p2x) + gamma2 * (bp1-bp2) * beta[0] + gamma * beta[0] * (E1-E2);
    float New_Py = (p1y - p2y) + gamma2 * (bp1-bp2) * beta[1] + gamma * beta[1] * (E1-E2);
    float New_Pz = (p1z - p2z) + gamma2 * (bp1-bp2) * beta[2] + gamma * beta[2] * (E1-E2);

    float* MassAndKstar = new float[7];
    MassAndKstar[0] = (gamma * (E1 + bp1 + E2 + bp2 + E3 + bp3 + E4 + bp4));
    MassAndKstar[1] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    New_Px = (p1x - p3x) + gamma2 * (bp1-bp3) * beta[0] + gamma * beta[0] * (E1-E3);
    New_Py = (p1y - p3y) + gamma2 * (bp1-bp3) * beta[1] + gamma * beta[1] * (E1-E3);
    New_Pz = (p1z - p3z) + gamma2 * (bp1-bp3) * beta[2] + gamma * beta[2] * (E1-E3);
    MassAndKstar[2] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    New_Px = (p2x - p3x) + gamma2 * (bp2-bp3) * beta[0] + gamma * beta[0] * (E2-E3);
    New_Py = (p2y - p3y) + gamma2 * (bp2-bp3) * beta[1] + gamma * beta[1] * (E2-E3);
    New_Pz = (p2z - p3z) + gamma2 * (bp2-bp3) * beta[2] + gamma * beta[2] * (E2-E3);
    MassAndKstar[3] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    New_Px = (p2x - p4x) + gamma2 * (bp2-bp4) * beta[0] + gamma * beta[0] * (E2-E4);
    New_Py = (p2y - p4y) + gamma2 * (bp2-bp4) * beta[1] + gamma * beta[1] * (E2-E4);
    New_Pz = (p2z - p4z) + gamma2 * (bp2-bp4) * beta[2] + gamma * beta[2] * (E2-E4);
    MassAndKstar[4] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    New_Px = (p3x - p4x) + gamma2 * (bp3-bp4) * beta[0] + gamma * beta[0] * (E3-E4);
    New_Py = (p3y - p4y) + gamma2 * (bp3-bp4) * beta[1] + gamma * beta[1] * (E3-E4);
    New_Pz = (p3z - p4z) + gamma2 * (bp3-bp4) * beta[2] + gamma * beta[2] * (E3-E4);
    MassAndKstar[5] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    New_Px = (p1x - p4x) + gamma2 * (bp1-bp4) * beta[0] + gamma * beta[0] * (E1-E4);
    New_Py = (p1y - p4y) + gamma2 * (bp1-bp4) * beta[1] + gamma * beta[1] * (E1-E4);
    New_Pz = (p1z - p4z) + gamma2 * (bp1-bp4) * beta[2] + gamma * beta[2] * (E1-E4);
    MassAndKstar[6] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
    return MassAndKstar;
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

bool IfInVector(int Num , std::vector<int> V)
{
    for (int i=0;i<V.size();i++) {
        if (Num == V.at(i)){
            return true;
        }
    }
    return false;
}

Double_t massList(int PID, TString DataName)
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

Double_t massListSigma(int PID, TString DataName)
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

void print(Event Temp)
{
	cout<<"EventID: "<<Temp.eventID<<endl;
    cout<<"    Particle A:"<<Temp.A_particles.size()<<endl;
    cout<<"        px       py       pz       mass"<<endl;
    for(int i=0;i<Temp.A_particles.size();i++){
        cout<<"     "<<(i+1)<<"  "<<Temp.A_particles[i].px<<" "<<Temp.A_particles[i].py<<" "<<Temp.A_particles[i].pz<<" "<<Temp.A_particles[i].mass<<endl;
    }
    cout<<"    Particle B:"<<Temp.B_particles.size()<<endl;
    cout<<"        px       py       pz       mass"<<endl;
    for(int i=0;i<Temp.B_particles.size();i++){
        cout<<"     "<<(i+1)<<"  "<<Temp.B_particles[i].px<<" "<<Temp.B_particles[i].py<<" "<<Temp.B_particles[i].pz<<" "<<Temp.B_particles[i].mass<<endl;
    }
    return ;
}

float CenCorr(float Vz, TString DataName)
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