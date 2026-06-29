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
#include "TH2.h"
#include "TProfile.h"
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
#include "TProfile.h"
// #endif
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <stdio.h>
using namespace std;

// 两体关联
// 使用这个编译：
// singularity exec -e --env DISPLAY=$DISPLAY -B /direct -B /gpfs -B /star -B /cvmfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh
// g++ -O2 -std=c++11 S_Two.cpp -o S_One `root-config --cflags --libs`

// 定义粒子结构体
struct ArmParticle {
    float   px;       // x方向动量
    float   py;       // y方向动量
    float   pz;       // z方向动量
    float   mass;     // 质量
    float   eta;      // 赝快度
    float   y;        // 快度
    double  pt;       // 横向动量
    bool    IsRecord; // 是否被记录
    int     TreeID;   // ID in one event
    double  p;        // 三动量绝对值
    double  E;        // 能量
    std::vector<int>   ParentID; // Parent Particle ID in one event
    
    ArmParticle()
        : px(0), py(0), pz(0), mass(0),
          eta(0), y(0), pt(0),
          IsRecord(false), TreeID(0) {}
    
    // 构造函数
    ArmParticle(double _px, double _py, double _pz, double _mass, int _TreeID) 
        : px(_px), py(_py), pz(_pz), mass(_mass), TreeID(_TreeID) {
        // 计算赝快度、快度和横向动量
        pt = sqrt(px*px + py*py);
        p = sqrt(pt*pt + pz*pz);
        E = sqrt(p*p+mass*mass);
        eta = -1.0*log(tan(0.5*(acos(pz/p))));
        y = 0.5 * log((E + pz) / (E - pz));
        IsRecord = false;
    }
    
    // 计算能量
    float energy() const {
        return sqrt(px*px + py*py + pz*pz + mass*mass);
    }
    
    // 转换为四动量
    TLorentzVector lorentzVector() const {
        return TLorentzVector(px, py, pz, energy());
    }

    // ArmParticle(const ArmParticle& other)
    //     : px(other.px), py(other.py), pz(other.pz),
    //       mass(other.mass), eta(other.eta), y(other.y), pt(other.pt),
    //       IsRecord(other.IsRecord), TreeID(other.TreeID),
    //       ParentID(other.ParentID) {}

    // ArmParticle& operator=(const ArmParticle& other) {
    //     if (this != &other) {
    //         px = other.px;
    //         py = other.py;
    //         pz = other.pz;
    //         mass = other.mass;
    //         eta = other.eta;
    //         y = other.y;
    //         pt = other.pt;
    //         IsRecord = other.IsRecord;
    //         TreeID = other.TreeID;
    //         ParentID = other.ParentID;
    //     }
    //     return *this;
    // }

};

// ROOT 5（特别是 ROOT 5.34.39）的字典机制有些老旧。
// 它会自动为 std::vector<T> 生成迭代器类型的字典（如 random_access_iterator<T,long>），但 C++11 之后这些类型模板已不再定义，因此报错。
#if defined(__CINT__) || defined(__CLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ArmParticle+;
// #pragma link C++ class std::vector<ArmParticle>+;
#endif
// 这告诉 ROOT：
// “只生成 ArmParticle 和 std::vector<ArmParticle> 的字典，
// 不要去尝试生成任何 random_access_iterator 之类的模板。”


#define Pi 3.1415926535898
#define HowMuchEventMixing 10

// int CentralityBin[] = {0 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 55 , 60 , 65 , 70 , 75 , 80 , 85 , 90 , 95 , 100};// %
// int CentralityBin[] = {0 , 10 , 20 , 30 , 40 , 50 , 60 , 70 , 80 , 90 , 100};// %
int CentralityBin[] = {0 , 10 , 20 , 30 , 40 , 50 , 60 , 70 , 80};// %
// const float PVzBin[] = {-45.0 , -35.0 , -25.0 , -15.0 , -5.0 , 5.0 , 15.0 , 25.0 , 35.0 , 45.0 , 55.0}; // Primary Vertex Z (cm) d+Au@200 GeV RUN 21 : -45 ~ 55 cm
const float PVzBin[] = {-80.0 , -70.0 , -60.0 , -50.0 , -40.0 , -30.0 , -20.0 , -10.0 , 0.0 , 10.0 , 20.0 , 30.0 , 40.0 , 50.0 , 60.0}; // Primary Vertex Z (cm) Au+Au@19.6 GeV RUN 19 
const float yBin[]  = {-10000.0 , 10000.0}; // B_y
const float AyCut[] = {-10000.0 , 10000.0}; // A_y
int FeedDown[] = { 0 };
const float EtaCut[] = {-1.5 , 1.5}; // EtaCut for both A and B
const float MassSigmaWidth = 3.0;
const float Sideband_MassSigmaSta   = 5.0;
const float Sideband_MassSigmaEnd   = 8.0;

const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1
const Int_t PVzBinNum = sizeof(PVzBin)/sizeof(PVzBin[0]) - 1; // -1
const Int_t yBinNum = sizeof(yBin)/sizeof(yBin[0]) - 1; // -1
const Int_t FeedDownNum = sizeof(FeedDown)/sizeof(FeedDown[0]);

void print(std::vector<int> Temp);
void print(std::vector<float> Temp);
std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize, TString DataName);
bool IfInVector(int Num , const std::vector<int>& V);
std::vector<int> GetDaughterPDGLit(int ID);
Double_t massList(int PID, TString DataName);
Double_t massListSigma(int PID, TString DataName);
inline bool GetSide(
    const ArmParticle& A,
    const ArmParticle& B,
    float& cosPhiOut,
    float& phiOut,
    bool IfRemoveFeedPair,
    const std::vector<float>& MotherMass,
    const std::vector<float>& MotherMassSigma,
    float MassSigmaWidth);
inline bool GetSide(
        const ArmParticle& A,
        const ArmParticle& B,
        TH1D& H_P_tot,
        TH1D& H_beta ,
        double& cosPhiOut,
        double& phiOut,
        bool IfRemoveFeedPair,
        const std::vector<float>& MotherMass,
        const std::vector<float>& MotherMassSigma,
        float MassSigmaWidth);
inline bool GetAngle(
    const ArmParticle& A,
    const ArmParticle& B,
    TH1D& H_P_tot,
    TH1D& H_beta,
    double& cosPhiOut,
    double& phiOut,
    bool IfRemoveFeedPair,
    const std::vector<float>& MotherMass,
    const std::vector<float>& MotherMassSigma,
    float MassSigmaWidth);
float CenCorr(float Vz, TString DataName);


// 定义事件结构体
struct Event {
    int eventID;                    // 事件ID
    std::vector<ArmParticle> A_particles;  // A类粒子 主粒子
    std::vector<ArmParticle> B_particles;  // B类粒子
    
    Event() : eventID(-1) {}
    // 构造函数
    Event(int _eventID) 
        : eventID(_eventID) {}
};

void print(Event Temp);

void S_Two(
    TString MidName,
    TString DataName,
    int OutputFileIndex,
    TString OutMidName,
    int A_PDG,
    int B_PDG,
    int If_SideBand_A,
    int If_SideBand_B,
    int Mode,
    int SP_ME,
    int RecordingMethod = 0,
    int CutID = 0
) {
    std::cout<<"Start MM.cpp"<<std::endl;

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
    float drap , dpt;
    std::vector<float> Side_Store , drap_Store , dpt_Store , mass_Store;
    std::vector<int>   IfRecorded;
    bool Is2Body = true;
    float NNch , Eta;
    TString TreeName = "hadronTree";

    bool Is_SideBand_A = true;
    bool Is_SideBand_B = true;
    if (If_SideBand_A == 0) Is_SideBand_A = false;
    if (If_SideBand_B == 0) Is_SideBand_B = false;

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
    int Aid , Bid , Cid , Did;
    float APx  , BPx  , CPx ;
    float APy  , BPy  , CPy ;
    float APz  , BPz  , CPz ;
    float APt  , BPt  , CPt ;
    float ARap , BRap , CRap;
    float BMass = massList(B_PDG, DataName)           , AMass = massList(A_PDG, DataName)          ;
    float BMassSigma = massListSigma(B_PDG, DataName) , AMassSigma = massListSigma(A_PDG, DataName);
    std::vector<std::vector<int> > D_ParID;
    bool IsSame;
    float  P_B , kStar;
    float  phi , CosPhi;
    double d_phi , d_CosPhi;
    double dRap;

    // //                                    centrality    A_Rapidity   PrimaryVertex
    // std::vector<Event>    EventPool         [50]           [50]          [50];
    // std::vector<ArmParticle> A_Array                       [50]              , B_Array;
    // std::vector<ArmParticle> A_List                        [50]              , B_List ;
    // TH1D                 *H1D_Side           [50]           [50]          [50];
    // TH1D                 *H1D_ALL_Side                      [50]     ;
    // TH1D                 *H1D_Mix_Side_ABC       [50]           [50]          [50];
    // TH1D                 *H1D_ALL_Mix_Side_ABC                  [50]     ;
    // TH1D                 *H_Tra_Side       [50]           [50]          [50];
    // TH1D                 *H_ALL_Tra_Side                  [50]     ;
    // TH1D                 *H_dRap            [50]           [50]          [50];
    // TH1D                 *H_ALL_dRap                       [50]     ;
    // TH1D                 *H_Mix_dRap        [50]           [50]          [50];
    // TH1D                 *H_ALL_Mix_dRap                   [50]     ;
    // TH1D                 *H_Tra_dRap        [50]           [50]          [50];
    // TH1D                 *H_ALL_Tra_dRap                   [50]     ;
    // TH1D                 *H_dPt             [50]           [50]          [50];
    // TH1D                 *H_ALL_dPt                        [50]     ;
    // TH1D                 *H_Mix_dPt         [50]           [50]          [50];
    // TH1D                 *H_ALL_Mix_dPt                    [50]     ;
    // TH1D                 *H_Tra_dPt         [50]           [50]          [50];
    // TH1D                 *H_ALL_Tra_dPt                    [50]     ;
    // TH1D                 *H_ALL_Mass                       [50]     ;
    // TH1D                 *H_ALL_Mix_Mass                   [50]     ;
    // TH1D                 *H_ALL_Tra_Mass                   [50]     ;
    // TH1D                 *H_Rap_A           [50]           [50]          [50];
    // TH1D                 *H_ALL_Rap_A                      [50]     ;
    // TH1D                 *H_Rap_K_A         [50]           [50]          [50];
    // TH1D                 *H_ALL_Rap_K_A                    [50]     ;
    // TH1D                 *H_Rap_B           [50]           [50]          [50];
    // TH1D                 *H_ALL_Rap_B                      [50]     ;
    // TH1D                 *H_Rap_K_B         [50]                         [50];
    // TH1D                 *H_ALL_Rap_K_B                             ;
    // // Used for test
    // TH2F                 *H_ALL_dRap_ARp                   [50]     ;
    // TH2F                 *H_ALL_Mix_dRap_ARp               [50]     ;
    // TH2F                 *H_Rap_A_B         [50]           [50]          [50];
    // TH2F                 *H_ALL_Rap_A_B                    [50]     ;
    // TH2F                 *H_Mix_Rap_A_B     [50]           [50]          [50];
    // TH2F                 *H_ALL_Mix_Rap_A_B                [50]     ;
    std::vector<std::vector<std::vector<std::vector<Event>>>> EventPool;
    std::vector<std::vector<ArmParticle>> A_Array(CentralityBinNum);
    std::vector<ArmParticle> B_Array;
    std::vector<ArmParticle> C_Array;
    std::vector<std::vector<ArmParticle>> A_List(yBinNum);
    std::vector<ArmParticle> B_List;
    std::vector<ArmParticle> C_List;
    std::vector<TH1D*>                                               H_ALL               ;
    std::vector<TH1D*>                                               H_ALL_Mix           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     H                   ;
    std::vector<std::vector<std::vector<TH1D*>>>                     H_Mix               ;
    std::vector<TH1D*>                                               H_ALL_Cos           ;
    std::vector<TH1D*>                                               H_ALL_Mix_Cos       ;
    std::vector<std::vector<std::vector<TH1D*>>>                     H_Cos               ;
    std::vector<std::vector<std::vector<TH1D*>>>                     H_Mix_Cos           ;
    std::vector<TH1D*>                                               H_dRap_ALL          ;
    std::vector<TH1D*>                                               H_dRap_ALL_Mix      ;
    std::vector<std::vector<std::vector<TH1D*>>>                     H_dRap              ;
    std::vector<std::vector<std::vector<TH1D*>>>                     H_dRap_Mix          ;

    std::vector<std::vector<std::vector<TProfile*>>> B_A_Num_Ratio;
    std::vector<TProfile*> B_A_Num_Ratio_ALL;
    std::vector<TH2D*> B_A_Num_2D;
    
    TH1D* H_P_tot     = new TH1D("H_P_tot","H_P_tot",200,0,10);
    TH1D* H_beta      = new TH1D("H_beta" ,"H_beta" ,500,0,2);
    TH1D* H_Mix_P_tot = new TH1D("H_Mix_P_tot","H_Mix_P_tot",200,0,10);
    TH1D* H_Mix_beta  = new TH1D("H_Mix_beta" ,"H_Mix_beta" ,500,0,2);


    if (RecordingMethod == 0) {
        EventPool.resize(CentralityBinNum);
        H_ALL       .resize(yBinNum, nullptr);
        H_ALL_Mix    .resize(yBinNum, nullptr);
        H      .resize(CentralityBinNum);
        H_Mix   .resize(CentralityBinNum);
        H_ALL_Cos       .resize(yBinNum, nullptr);
        H_ALL_Mix_Cos    .resize(yBinNum, nullptr);
        H_Cos      .resize(CentralityBinNum);
        H_Mix_Cos   .resize(CentralityBinNum);
        H_dRap_ALL       .resize(yBinNum, nullptr);
        H_dRap_ALL_Mix    .resize(yBinNum, nullptr);
        H_dRap      .resize(CentralityBinNum);
        H_dRap_Mix   .resize(CentralityBinNum);
        B_A_Num_Ratio      .resize(CentralityBinNum);
        B_A_Num_Ratio_ALL       .resize(yBinNum, nullptr);
        B_A_Num_2D       .resize(yBinNum, nullptr);
        for (i = 0; i < CentralityBinNum; i++) {
            EventPool[i].resize(yBinNum);
            H              [i].resize(yBinNum);
            H_Mix          [i].resize(yBinNum);
            H_Cos          [i].resize(yBinNum);
            H_Mix_Cos      [i].resize(yBinNum);
            H_dRap         [i].resize(yBinNum);
            H_dRap_Mix     [i].resize(yBinNum);
            B_A_Num_Ratio  [i].resize(yBinNum);
            for (j = 0; j < yBinNum; j++) {
                EventPool[i][j].resize(PVzBinNum);
                H                     [i][j].resize(PVzBinNum, nullptr);
                H_Mix                 [i][j].resize(PVzBinNum, nullptr);
                H_Cos                 [i][j].resize(PVzBinNum, nullptr);
                H_Mix_Cos             [i][j].resize(PVzBinNum, nullptr);
                H_dRap                [i][j].resize(PVzBinNum, nullptr);
                H_dRap_Mix            [i][j].resize(PVzBinNum, nullptr);
                B_A_Num_Ratio         [i][j].resize(PVzBinNum, nullptr);
            }
        }
    }
    for (j = 0; j < CentralityBinNum; j++) {
        A_Array[j].resize(HowMuchEventMixing);
    }
    for (j = 0; j < yBinNum; j++) {
        A_List[j].resize(HowMuchEventMixing);
    }

    ArmParticle           A(0,0,0,0,0), B(0,0,0,0,0), C(0,0,0,0,0), D(0,0,0,0,0);
    Event                 TempEvent(0);

    int SideBinNum = 100;
    float SideSta = 0 , SideEnd = Pi;
    
    int dRapBinNum = 500;
    float dRapSta = -5 , dRapEnd = 5;
    
    int SRapBinNum = 1000;
    float SRapSta = -10 , SRapEnd = 10;
    
    int dPtBinNum = 200;
    float dPtSta = 0 , dPtEnd = 10;
    
    int MBinNum = 1000 , MBinPar = 100;
    float MSta = floor((AMass + BMass)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;
    cout<<"Mass Region: [ "<<MSta<<" , "<<MEnd<<" ], BinNum = "<<MBinNum<<". "<<endl;

    if (RecordingMethod == 0) {
        for (RapIndex=0;RapIndex<yBinNum;RapIndex++) {
            for (CenIndex=0;CenIndex<CentralityBinNum;CenIndex++) {
                for (PVzIndex=0;PVzIndex<PVzBinNum;PVzIndex++) {
                    H               [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("H_%d_%d_%d"         ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    H_Mix           [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("H_Mix_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("Mix, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    H_Cos           [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("H_Cos_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    H_Mix_Cos       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("H_Mix_Cos_%d_%d_%d" ,CenIndex,RapIndex,PVzIndex), Form("Mix, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    H_dRap          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("H_dRap_%d_%d_%d"         ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dRapBinNum,dRapSta,dRapEnd);
                    H_dRap_Mix      [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("H_dRap_Mix_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("Mix, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),dRapBinNum,dRapSta,dRapEnd);
                    B_A_Num_Ratio   [CenIndex] [RapIndex] [PVzIndex] = new TProfile(Form("B_A_Num_Ratio_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("B Num average vs. A Num, [%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),20,0,20,0,20);
                    
                }
            }
            H_ALL                  [RapIndex] = new TH1D(Form("H_ALL_%d"      ,          RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            H_ALL_Mix              [RapIndex] = new TH1D(Form("H_ALL_Mix_%d"  ,          RapIndex), Form("ALL Mix,  %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            H_ALL_Cos              [RapIndex] = new TH1D(Form("H_ALL_Cos_%d"      ,      RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            H_ALL_Mix_Cos          [RapIndex] = new TH1D(Form("H_ALL_Mix_Cos_%d"  ,      RapIndex), Form("ALL Mix,  %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            H_dRap_ALL             [RapIndex] = new TH1D(Form("H_dRap_ALL_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),dRapBinNum,dRapSta,dRapEnd);
            H_dRap_ALL_Mix         [RapIndex] = new TH1D(Form("H_dRap_ALL_Mix_%d"  ,     RapIndex), Form("ALL Mix,  %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),dRapBinNum,dRapSta,dRapEnd);
            B_A_Num_Ratio_ALL      [RapIndex] = new TProfile(Form("B_A_Num_Ratio_ALL_%d"  ,     RapIndex), Form("B Num average vs. A Num,  %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),20,0,20,0,20);
            B_A_Num_2D             [RapIndex] = new TH2D(Form("B_A_Num_2D_%d"  ,     RapIndex), Form("B Num vs. A Num,  %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),20,0,20,20,0,20);
            
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

    // =========================
    // 打开并逐行读取文件
    // =========================
    TChain *hadronTree = new TChain(TreeName);
    std::ifstream infile(MidName.Data());

    if (!infile.is_open()) {
        std::cerr << "Error: cannot open file " << MidName << std::endl;
        return;
    }

    std::string line;
    int lineCount = 0;

    while (std::getline(infile, line)) {
        lineCount++;

        // 跳过空行（可选）
        if (line.empty()) continue;

        // 跳过注释（可选，比如 # 开头）
        if (line[0] == '#') continue;

        // 处理每一行
        std::cout << "Line " << lineCount << ": " << line << std::endl;

        // 如果你需要转成 TString：
        TString tline(line);
        TFile *f = TFile::Open(line.c_str());

        if (!f || f->IsZombie()) {
            std::cerr << "Bad file: " << line << std::endl;
            continue;
        }
        
        // 检查 tree 是否存在
        TTree *t = (TTree*)f->Get(TreeName);
        
        if (!t) {
            std::cerr << "No tree " << TreeName << " in " << line << std::endl;
            f->Close();
            // delete f;
            continue;
        }
        
        // 只有通过检查才加入
        hadronTree->Add(line.c_str());
        
        f->Close();
        delete f;

    }

    infile.close();

    std::cout << "Total lines read: " << lineCount << std::endl;

    // TChain *hadronTree = new TChain(TreeName);
    // for(i=StartFileIndex;i <= EndFileIndex;i++){
    //     TString filename = MidName;
    //     filename+=i;
    //     filename+=".root";
    //     hadronTree->Add(filename);
    // }
    Int_t PDGMult  ;
    Int_t refMult  ;
    Int_t grefMult ;
    Int_t EventID  ;
    Int_t RunID    ;
    Int_t TriggerID;
    Int_t Nch      ;
    float PVz      ;

    hadronTree->SetBranchAddress("PDGMult"  ,&PDGMult  );
    hadronTree->SetBranchAddress("refMult"  ,&refMult  );
    // hadronTree->SetBranchAddress("grefMult" ,&grefMult );
    // hadronTree->SetBranchAddress("EventID"  ,&EventID  );
    // hadronTree->SetBranchAddress("RunID"    ,&RunID    );
    // hadronTree->SetBranchAddress("TriggerID",&TriggerID);
    // hadronTree->SetBranchAddress("Nch"      ,&Nch      );
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


        D_ParID.clear();
        B_List.clear();
        C_List.clear();
        TempEvent.eventID = EntriesID;
        TempEvent.A_particles.clear();
        TempEvent.B_particles.clear();
        for (size_t st=0;st<MatchedRap.size();st++) {
            A_Array[MatchedRap.at(st)].clear();
            A_List [MatchedRap.at(st)].clear();
        }
        MatchedRap.clear();
        // 定Centrality
        CenIndex = -1;
        for (k=0;k<CentralityBinNum;k++){
            // NNch = CenCorr(PVz, DataName) * Nch;
            NNch = CenCorr(PVz, DataName) * refMult;
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
                if ((!Is_SideBand_A && (fabs(InvariantMass->at(i) - AMass) <= MassSigmaWidth*AMassSigma)) || 
                    ( Is_SideBand_A && ((fabs(InvariantMass->at(i) - AMass) >= Sideband_MassSigmaSta*AMassSigma) && (fabs(InvariantMass->at(i) - AMass) <= Sideband_MassSigmaEnd*AMassSigma))))
                {

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

                    A = ArmParticle(mix_px->at(i),mix_py->at(i),mix_pz->at(i),AMass,i);
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
                    continue;
                }
            }
            else if (PDG->at(i) == B_PDG) {
                if ((!Is_SideBand_B && (fabs(InvariantMass->at(i) - BMass) <= MassSigmaWidth*BMassSigma)) || 
                    ( Is_SideBand_B && ((fabs(InvariantMass->at(i) - BMass) >= Sideband_MassSigmaSta*BMassSigma) && (fabs(InvariantMass->at(i) - BMass) <= Sideband_MassSigmaEnd*BMassSigma))))
                {

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

                    B = ArmParticle(mix_px->at(i),mix_py->at(i),mix_pz->at(i),BMass,i);
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
                        D_ParID.push_back(Temp);
                        // IfFoundOmega = true;
                        // cout<<"Found Omega"<<endl;
                    }
                }
            }
        }
        // 筛选A、B、C粒子：优先级：D（母粒子） > A > B
        // 筛选B粒子
        for (Bid=0;Bid<B_List.size();Bid++) {
            IfRecord = true;
            // 如果B、D有血缘关系，不记录B
            for (Did = 0;Did < D_ParID.size();Did++) {
                if (IfInVector(B_List[Bid].TreeID , D_ParID.at(Did))) {IfRecord = false;break;}
            }
            // 如果B、A有血缘关系，不记录B
            if (IfRecord) {
                for (i=0;i<MatchedRap.size();i++) {
                    for (Aid=0;Aid<A_List[MatchedRap[i]].size();Aid++) {
                        if (IfInVector(B_List[Bid].TreeID , A_List[MatchedRap[i]][Aid].ParentID)) {IfRecord = false;break;}
                    }
                }
            }
            if (IfRecord) TempEvent.B_particles.push_back(B_List[Bid]);
        }
        // 筛选A粒子
        for (i=0;i<MatchedRap.size();i++) {
            for (Aid=0;Aid<A_List[MatchedRap[i]].size();Aid++) {
                IfRecord = true;
                // 如果A、D有血缘关系，不记录A
                for (Did = 0;Did < D_ParID.size();Did++) {
                    if (IfInVector(A_List[MatchedRap[i]][Aid].TreeID , D_ParID.at(Did))) {IfRecord = false;break;}
                }
                if (IfRecord) A_Array[MatchedRap[i]].push_back(A_List[MatchedRap[i]][Aid]);
            }
        }
        // if (TempEvent.B_particles.size() >= HowMuchEventMixing+1) continue;
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
                        for (int Aid = 0; Aid < HowMuchEventMixing + 1; ++Aid) {

                            auto& eventA = EventPool[CenIndex][RapIndex][PVzIndex][Aid];
                            const auto& A_particles = eventA.A_particles;
                        
                            for (int Bid = 0; Bid < HowMuchEventMixing + 1; ++Bid) {
                        
                                auto& eventB = EventPool[CenIndex][RapIndex][PVzIndex][Bid];
                                const auto& B_particles = eventB.B_particles;

                                if (Aid == Bid) {
                                    IsSame = true;
                                }else{
                                    IsSame = false;
                                }

                                if (IsSame) {
                                    B_A_Num_Ratio      [CenIndex][RapIndex][PVzIndex]->Fill(A_particles.size(),B_particles.size());
                                    B_A_Num_Ratio_ALL            [RapIndex]          ->Fill(A_particles.size(),B_particles.size());
                                    B_A_Num_2D                   [RapIndex]          ->Fill(A_particles.size(),B_particles.size());
                                }
                        
                                for (const auto& A : A_particles) {
                        
                                    for (const auto& B : B_particles) {
                                        
                                        if (IsSame) {
                                            if (GetSide(A,B , *H_P_tot, *H_beta, d_CosPhi,d_phi, IfRemoveFeedPair, MotherMass, MotherMassSigma, MassSigmaWidth)){
                                            // if (GetAngle(A,B , *H_P_tot, *H_beta, d_CosPhi,d_phi, IfRemoveFeedPair, MotherMass, MotherMassSigma, MassSigmaWidth)){
                                                H                [CenIndex][RapIndex][PVzIndex]->Fill(d_phi);
                                                H_ALL                      [RapIndex]->Fill(d_phi);
                                                H_Cos            [CenIndex][RapIndex][PVzIndex]->Fill(d_CosPhi);
                                                H_ALL_Cos                  [RapIndex]->Fill(d_CosPhi);
                                                dRap = (A.y > 0.0) ? (B.y - A.y) : (A.y - B.y);
                                                H_dRap_ALL                 [RapIndex]->Fill(dRap);
                                                H_dRap           [CenIndex][RapIndex][PVzIndex]->Fill(dRap);
                                            }
                                        }else{
                                            if (GetSide(A,B , *H_P_tot, *H_beta, d_CosPhi,d_phi, IfRemoveFeedPair, MotherMass, MotherMassSigma, MassSigmaWidth)){
                                            // if (GetAngle(A,B , *H_Mix_P_tot, *H_Mix_beta, d_CosPhi,d_phi, IfRemoveFeedPair, MotherMass, MotherMassSigma, MassSigmaWidth)){
                                                H_Mix            [CenIndex][RapIndex][PVzIndex]->Fill(d_phi);
                                                H_ALL_Mix                  [RapIndex]->Fill(d_phi);
                                                H_Mix_Cos        [CenIndex][RapIndex][PVzIndex]->Fill(d_CosPhi);
                                                H_ALL_Mix_Cos              [RapIndex]->Fill(d_CosPhi);
                                                dRap = (A.y > 0.0) ? (B.y - A.y) : (A.y - B.y);
                                                H_dRap_ALL_Mix             [RapIndex]->Fill(dRap);
                                                H_dRap_Mix       [CenIndex][RapIndex][PVzIndex]->Fill(dRap);
                                            }
                                        }
                                        ++AccumSameNum;
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
    OutputFileName += "H_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileA = new TFile(OutputFileName, "RECREATE");
    TDirectory *folder_Side     = fileA->mkdir("Side");
    TDirectory *ALL_Side        = folder_Side->mkdir("ALL");
    TDirectory *Sep_Side        = folder_Side->mkdir("Sep");
    TDirectory *folder_dRap     = fileA->mkdir("dRap");
    TDirectory *ALL_dRap        = folder_dRap->mkdir("ALL");
    TDirectory *Sep_dRap        = folder_dRap->mkdir("Sep");
    TDirectory *folder_BAR     = fileA->mkdir("B_A_Num");
    TDirectory *ALL_BAR        = folder_BAR->mkdir("ALL");
    TDirectory *Sep_BAR        = folder_BAR->mkdir("Sep");
    fileA->cd();
    H_P_tot->Write();
    H_beta ->Write();
    H_Mix_P_tot->Write();
    H_Mix_beta ->Write();
    for (RapIndex=0;RapIndex<yBinNum;RapIndex++) {
        for (CenIndex=0;CenIndex<CentralityBinNum;CenIndex++) {
            for (PVzIndex=0;PVzIndex<PVzBinNum;PVzIndex++) {
                Sep_Side->cd();
                H                    [CenIndex] [RapIndex] [PVzIndex] ->Write();
                H_Mix                [CenIndex] [RapIndex] [PVzIndex] ->Write();
                H_Cos                [CenIndex] [RapIndex] [PVzIndex] ->Write();
                H_Mix_Cos            [CenIndex] [RapIndex] [PVzIndex] ->Write();
                Sep_dRap->cd();
                H_dRap               [CenIndex] [RapIndex] [PVzIndex] ->Write();
                H_dRap_Mix           [CenIndex] [RapIndex] [PVzIndex] ->Write();
                Sep_BAR->cd();
                B_A_Num_Ratio        [CenIndex] [RapIndex] [PVzIndex] ->Write();
            }
        }
        ALL_Side->cd();
        H_ALL                                   [RapIndex] ->Write();
        H_ALL_Mix                               [RapIndex] ->Write();
        H_ALL_Cos                               [RapIndex] ->Write();
        H_ALL_Mix_Cos                           [RapIndex] ->Write();
        ALL_dRap->cd();
        H_dRap_ALL                              [RapIndex] ->Write();
        H_dRap_ALL_Mix                          [RapIndex] ->Write();
        ALL_BAR->cd();
        B_A_Num_Ratio_ALL                       [RapIndex] ->Write();
        B_A_Num_2D                              [RapIndex] ->Write();
    }
    fileA->Close();
    cout<<"FINISH!"<<endl;
    return;
}

int main(int argc, char** argv) {
    // 检查参数数量
    if(argc < 12) {
        std::cerr << "Usage: " << argv[0] 
                  << " MidName DataName OutputFileIndex OutMidName"
                  << " A_PDG B_PDG If_SideBand_A If_SideBand_B Mode SP_ME [CutID]" << std::endl;
        return 1;
    }

    S_Two(
        TString(argv[1]),
        TString(argv[2]),
        atoi(argv[3]),
        TString(argv[4]),
        atoi(argv[5]),
        atoi(argv[6]),
        atoi(argv[7]),
        atoi(argv[8]),
        atoi(argv[9]),
        atoi(argv[10]),
        atoi(argv[11]),
        (argc > 12 ? atoi(argv[12]) : 0)
    );

    return 0;
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
    if (DataName == "AuAu_19_19") {
        // data from https://drupal.star.bnl.gov/STAR/system/files/19p6GeVCentrality_v1.pdf
        int NchTable[21] = { 500 , 296 , 243 , 201 , 165 , 135 , 110 , 88 , 70 , 55 , 43 , 32 , 24 , 18 , 13 , 9 , 6};
        int CenTable[21] = {   0 ,   5 ,  10 ,  15 ,  20 ,  25 ,  30 , 35 , 40 , 45 , 50 , 55 , 60 , 65 , 70 , 75 ,80};
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

inline bool GetSide(
    const ArmParticle& A,
    const ArmParticle& B,
    float& cosPhiOut,
    float& phiOut,
    bool IfRemoveFeedPair,
    const std::vector<float>& MotherMass,
    const std::vector<float>& MotherMassSigma,
    float MassSigmaWidth)
{
    //--------------------------------------------------
    // Total four momentum
    //--------------------------------------------------

    const double Px = double(A.px) + double(B.px);
    const double Py = double(A.py) + double(B.py);
    const double Pz = double(A.pz) + double(B.pz);
    const double E  = double(A.E ) + double(B.E );

    //--------------------------------------------------
    // Invariant mass
    //--------------------------------------------------

    const double M2 =
        E*E
      - Px*Px
      - Py*Py
      - Pz*Pz;

    if (M2 <= 0.0)
        return false;

    //--------------------------------------------------
    // Feed-down rejection
    //--------------------------------------------------

    if (IfRemoveFeedPair) {

        const double M = std::sqrt(M2);

        for (size_t i = 0; i < MotherMass.size(); ++i) {

            if (std::fabs(M - MotherMass[i])
                < MassSigmaWidth * MotherMassSigma[i]) {

                return false;
            }
        }
    }

    //--------------------------------------------------
    // beta
    //--------------------------------------------------

    const double invE = 1.0 / E;

    const double betaX = -Px * invE;
    const double betaY = -Py * invE;
    const double betaZ = -Pz * invE;

    const double beta2 =
        betaX*betaX +
        betaY*betaY +
        betaZ*betaZ;

    if (beta2 < 1e-20 || beta2 >= 1.0)
        return false;

    const double betaAbs =
        std::sqrt(beta2);

    //--------------------------------------------------
    // beta direction
    //--------------------------------------------------

    const double invBeta =
        1.0 / betaAbs;

    const double nx = -betaX * invBeta;
    const double ny = -betaY * invBeta;
    const double nz = -betaZ * invBeta;

    //--------------------------------------------------
    // longitudinal momentum
    //--------------------------------------------------

    const double pPar =
          double(B.px)*nx
        + double(B.py)*ny
        + double(B.pz)*nz;

    //--------------------------------------------------
    // total momentum squared
    //--------------------------------------------------

    const double p2 =
          double(B.px)*double(B.px)
        + double(B.py)*double(B.py)
        + double(B.pz)*double(B.pz);

    //--------------------------------------------------
    // transverse momentum squared
    //--------------------------------------------------

    const double pPerp2 =
        std::max(0.0,
                 p2 - pPar*pPar);

    //--------------------------------------------------
    // gamma
    //--------------------------------------------------

    const double gamma =
        1.0 / std::sqrt(1.0 - beta2);

    //--------------------------------------------------
    // boosted longitudinal momentum
    //--------------------------------------------------

    const double pParStar =
        gamma * (pPar - betaAbs * double(B.E));

    //--------------------------------------------------
    // numerically stable cos(phi)
    //--------------------------------------------------

    const double denom =
        pParStar * pParStar;

    double cosPhi;

    if (denom <= 1e-30) {

        cosPhi = 0.0;

    } else {

        const double ratio =
            pPerp2 / denom;

        cosPhi =
            ((pParStar >= 0.0) ? 1.0 : -1.0)
            /
            std::sqrt(1.0 + ratio);
    }

    //--------------------------------------------------
    // Clamp
    //--------------------------------------------------

    cosPhi =
        std::max(-1.0,
        std::min( 1.0, cosPhi));

    //--------------------------------------------------
    // output
    //--------------------------------------------------

    cosPhiOut = float(cosPhi);
    phiOut    = float(std::acos(cosPhi));

    return true;
}

inline bool GetSide(
    const ArmParticle& A,
    const ArmParticle& B,
    TH1D& H_P_tot,
    TH1D& H_beta,
    double& cosPhiOut,
    double& phiOut,
    bool IfRemoveFeedPair,
    const std::vector<float>& MotherMass,
    const std::vector<float>& MotherMassSigma,
    float MassSigmaWidth)
{
    const double AE = A.E;
    const double BE = B.E;
    const double TotE = AE + BE;
    const double p[3] = {A.px+B.px , A.py+B.py , A.pz+B.pz};
    double beta[4] = { -(p[0])/TotE , -(p[1])/TotE , -(p[2])/TotE , 0.0};
    beta[3] = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];

    double P_tot = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);

    H_P_tot.Fill(P_tot);
    H_beta .Fill(sqrt(beta[3]) );

    // if (sqrt(beta[3]) > 1) {
    //     cout<<"#########################################"<<endl;
    //     cout<<"APx   = "<<A.px<<endl;
    //     cout<<"APy   = "<<A.py<<endl;
    //     cout<<"APz   = "<<A.pz<<endl;
    //     cout<<"AMass = "<<A.mass<<endl;
    //     cout<<"BPx   = "<<B.px<<endl;
    //     cout<<"BPy   = "<<B.py<<endl;
    //     cout<<"BPz   = "<<B.pz<<endl;
    //     cout<<"BMass = "<<B.mass<<endl;
    //     cout<<"AE    = "<<AE<<endl;
    //     cout<<"BE    = "<<BE<<endl;
    //     cout<<"TotE=AE  + BE  = "<<TotE<<endl;
    //     cout<<"Px = APx + BPx = "<<p[0]<<endl;
    //     cout<<"Py = APy + BPy = "<<p[1]<<endl;
    //     cout<<"Pz = APz + BPz = "<<p[2]<<endl;
    //     cout<<"P_tot          = "<<sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])<<endl;
    //     cout<<"beta[0]        = "<<beta[0]<<endl;
    //     cout<<"beta[1]        = "<<beta[1]<<endl;
    //     cout<<"beta[2]        = "<<beta[2]<<endl;
    //     cout<<"beta = P_tot/E = "<<sqrt(beta[3])<<endl;
    //     cout<<"#########################################"<<endl;
    // }

    const double gamma  = 1.0/(sqrt(1-beta[3]));
    const double gamma2 = 1.0/(sqrt(1-beta[3])*(1+sqrt(1-beta[3])));

    const double bpB = beta[0]*B.px + beta[1]*B.py + beta[2]*B.pz;

    const double New_BPx = B.px + gamma2*beta[0]*bpB + gamma*beta[0]*BE;
    const double New_BPy = B.py + gamma2*beta[1]*bpB + gamma*beta[1]*BE;
    const double New_BPz = B.pz + gamma2*beta[2]*bpB + gamma*beta[2]*BE;

    cosPhiOut = (New_BPx*(p[0])+New_BPy*(p[1])+New_BPz*(p[2])) / (sqrt(New_BPx*New_BPx+New_BPy*New_BPy+New_BPz*New_BPz)*P_tot);
    phiOut    = std::acos(cosPhiOut);

    return true;
}
// {
//     //--------------------------------------------------
//     // Total four momentum (promoted to double)
//     //--------------------------------------------------

//     const double Px = double(A.px) + double(B.px);
//     const double Py = double(A.py) + double(B.py);
//     const double Pz = double(A.pz) + double(B.pz);

//     //----------------------------------------------------------------------
//     // Recompute single-particle energies in DOUBLE precision from the raw
//     // px,py,pz,mass — do NOT use the stored float A.E / B.E.  For highly
//     // relativistic particles (p >> m) the float E loses m² in the mantissa
//     // and rounds to E ≈ p, which makes P_tot / E_tot ≥ 1 (unphysical).
//     //----------------------------------------------------------------------

//     const double EA =
//         std::sqrt( double(A.px)*double(A.px)
//                  + double(A.py)*double(A.py)
//                  + double(A.pz)*double(A.pz)
//                  + double(A.mass)*double(A.mass) );

//     const double EB =
//         std::sqrt( double(B.px)*double(B.px)
//                  + double(B.py)*double(B.py)
//                  + double(B.pz)*double(B.pz)
//                  + double(B.mass)*double(B.mass) );

//     const double E = EA + EB;

//     //--------------------------------------------------
//     // Invariant mass (mass-based, avoids E² - P² cancellation)
//     //   M² = m₁² + m₂² + 2(E₁E₂ - p₁·p₂)
//     //--------------------------------------------------

//     const double M2 =
//           double(A.mass) * double(A.mass)
//         + double(B.mass) * double(B.mass)
//         + 2.0 * ( EA * EB
//                 - double(A.px) * double(B.px)
//                 - double(A.py) * double(B.py)
//                 - double(A.pz) * double(B.pz) );

//     if (M2 <= 0.0)
//         return false;

//     const double M = std::sqrt(M2);

//     //--------------------------------------------------
//     // Feed-down rejection (was missing from original)
//     //--------------------------------------------------

//     if (IfRemoveFeedPair) {

//         for (size_t i = 0; i < MotherMass.size(); ++i) {

//             if (std::fabs(M - MotherMass[i])
//                 < MassSigmaWidth * MotherMassSigma[i]) {

//                 return false;
//             }
//         }
//     }

//     //--------------------------------------------------
//     // Total momentum magnitude & beta
//     //--------------------------------------------------

//     const double P_tot = std::sqrt(Px*Px + Py*Py + Pz*Pz);

//     const double beta = P_tot / E;             // |β| = |P_tot| / E_tot

//     H_P_tot.Fill(P_tot);
//     H_beta .Fill(beta );

//     // if (beta > 1) {
//     //     cout<<"#########################################"<<endl;
//     //     cout<<"APx   = "<<A.px<<endl;
//     //     cout<<"APy   = "<<A.py<<endl;
//     //     cout<<"APz   = "<<A.pz<<endl;
//     //     cout<<"AMass = "<<A.mass<<endl;
//     //     cout<<"BPx   = "<<B.px<<endl;
//     //     cout<<"BPy   = "<<B.py<<endl;
//     //     cout<<"BPz   = "<<B.pz<<endl;
//     //     cout<<"BMass = "<<B.mass<<endl;
//     //     cout<<"AE    = "<<EA<<endl;
//     //     cout<<"BE    = "<<EB<<endl;
//     //     cout<<"E  = AE  + BE  = "<<E<<endl;
//     //     cout<<"Px = APx + BPx = "<<Px<<endl;
//     //     cout<<"Py = APy + BPy = "<<Py<<endl;
//     //     cout<<"Pz = APz + BPz = "<<Pz<<endl;
//     //     cout<<"P_tot          = "<<P_tot<<endl;
//     //     cout<<"beta = P_tot/E = "<<beta<<endl;
//     //     cout<<"#########################################"<<endl;
//     // }

//     if (beta < 1e-20 || beta >= 1.0)
//         return false;

//     //--------------------------------------------------
//     // Boost direction unit vector  n̂ = +P̂_tot
//     //--------------------------------------------------

//     const double invP = 1.0 / P_tot;
//     const double nx = Px * invP;
//     const double ny = Py * invP;
//     const double nz = Pz * invP;

//     //--------------------------------------------------
//     // Decompose B momentum: parallel + perpendicular
//     //--------------------------------------------------

//     const double pPar =
//           double(B.px) * nx
//         + double(B.py) * ny
//         + double(B.pz) * nz;

//     const double p2 =
//           double(B.px) * double(B.px)
//         + double(B.py) * double(B.py)
//         + double(B.pz) * double(B.pz);

//     const double pPerp2 =
//         std::max(0.0, p2 - pPar * pPar);

//     //--------------------------------------------------
//     // gamma (via invariant mass: γ = E/M — avoids 1-β² cancellation)
//     //--------------------------------------------------

//     const double gamma = E / M;

//     //--------------------------------------------------
//     // Boosted parallel momentum:  p'∥ = γ (p∥ - β E_B)
//     // Perpendicular component is Lorentz invariant: p'⊥ = p⊥
//     //--------------------------------------------------

//     const double pParStar =
//         gamma * (pPar - beta * EB);

//     //--------------------------------------------------
//     // Numerically stable  cos(θ*) = sign(p'∥) / √(1 + p'⊥²/p'∥²)
//     //--------------------------------------------------

//     double cosPhi;
//     const double denom = pParStar * pParStar;

//     if (denom <= 1e-30) {

//         cosPhi = 0.0;

//     } else {

//         const double ratio = pPerp2 / denom;

//         cosPhi =
//             ((pParStar >= 0.0) ? 1.0 : -1.0)
//             /
//             std::sqrt(1.0 + ratio);
//     }

//     //--------------------------------------------------
//     // Clamp to [-1, 1] (floating-point edge cases)
//     //--------------------------------------------------

//     cosPhi = std::max(-1.0, std::min(1.0, cosPhi));

//     cosPhiOut = cosPhi;
//     phiOut    = std::acos(cosPhi);

//     return true;
// }


inline bool GetAngle(
    const ArmParticle& A,
    const ArmParticle& B,
    TH1D& H_P_tot,
    TH1D& H_beta,
    double& cosPhiOut,
    double& phiOut,
    bool IfRemoveFeedPair,
    const std::vector<float>& MotherMass,
    const std::vector<float>& MotherMassSigma,
    float MassSigmaWidth)
{
    double beta[4] = { -A.px/A.E , -A.py/A.E , -A.pz/A.E , 0.0};
    beta[3] = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];

    H_P_tot.Fill(A.p);
    H_beta .Fill(sqrt(beta[3]) );

    // if (sqrt(beta[3]) > 1) {
    //     cout<<"#########################################"<<endl;
    //     cout<<"APx   = "<<A.px<<endl;
    //     cout<<"APy   = "<<A.py<<endl;
    //     cout<<"APz   = "<<A.pz<<endl;
    //     cout<<"AMass = "<<A.mass<<endl;
    //     cout<<"BPx   = "<<B.px<<endl;
    //     cout<<"BPy   = "<<B.py<<endl;
    //     cout<<"BPz   = "<<B.pz<<endl;
    //     cout<<"BMass = "<<B.mass<<endl;
    //     cout<<"AE    = "<<AE<<endl;
    //     cout<<"BE    = "<<BE<<endl;
    //     cout<<"TotE=AE  + BE  = "<<TotE<<endl;
    //     cout<<"Px = APx + BPx = "<<p[0]<<endl;
    //     cout<<"Py = APy + BPy = "<<p[1]<<endl;
    //     cout<<"Pz = APz + BPz = "<<p[2]<<endl;
    //     cout<<"P_tot          = "<<sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])<<endl;
    //     cout<<"beta[0]        = "<<beta[0]<<endl;
    //     cout<<"beta[1]        = "<<beta[1]<<endl;
    //     cout<<"beta[2]        = "<<beta[2]<<endl;
    //     cout<<"beta = P_tot/E = "<<sqrt(beta[3])<<endl;
    //     cout<<"#########################################"<<endl;
    // }

    const double gamma  = 1.0/(sqrt(1-beta[3]));
    const double gamma2 = 1.0/(sqrt(1-beta[3])*(1+sqrt(1-beta[3])));

    const double bpB = beta[0]*B.px + beta[1]*B.py + beta[2]*B.pz;

    const double New_BPx = B.px + gamma2*beta[0]*bpB + gamma*beta[0]*B.E;
    const double New_BPy = B.py + gamma2*beta[1]*bpB + gamma*beta[1]*B.E;
    const double New_BPz = B.pz + gamma2*beta[2]*bpB + gamma*beta[2]*B.E;

    cosPhiOut = (New_BPx*A.px+New_BPy*A.py+New_BPz*A.pz) / (sqrt(New_BPx*New_BPx+New_BPy*New_BPy+New_BPz*New_BPz)*A.p);
    phiOut    = std::acos(cosPhiOut);
    if (cosPhiOut<-1.0) {cosPhiOut = -1.0;phiOut = Pi ;}
    if (cosPhiOut> 1.0) {cosPhiOut =  1.0;phiOut = 0.0;}

    return true;
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

bool IfInVector(int Num , const std::vector<int>& V)
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
    if (DataName == "AuAu_19_19"){
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
                Result = 1.6726;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 1.6726;
                break;
            case 3312 :// XiFitMass
                Result = 1.3221;
                break;
            case -3312 :// XiBarFitMass
                Result = 1.3221;
                break;
            case 3122 :// LambdaFitMass
                Result = 1.1159;
                break;
            case -3122 :// LambdaBarFitMass
                Result = 1.1158;
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
    if (DataName == "AuAu_19_19"){
        switch (PID)
        {
            case 3334 :// OmegaFitMass
                Result = 0.0018;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 0.0018;
                break;
            case 1003314 :// XiRPdgMass
                Result = 0.0029;
                break;
            case -1003314 :// XiRPdgMass
                Result = 0.0024;
                break;
            case 3312 :// XiFitMass
                Result = 0.0017;
                break;
            case -3312 :// XiBarFitMass
                Result = 0.0017;
                break;
            case 3122 :// LambdaFitMass
                Result = 0.0013;
                break;
            case -3122 :// LambdaBarFitMass
                Result = 0.0013;
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
    if (DataName == "AuAu_19_19") {// data from https://drupal.star.bnl.gov/STAR/system/files/19p6GeVCentrality_v1.pdf
        if      (Vz < -135.0) {
            return 0.998171;
        }
        else if (Vz < -125.0) {
            return 0.99364;
        }
        else if (Vz < -115.0) {
            return 0.991578;
        }
        else if (Vz < -105.0) {
            return 0.990252;
        }
        else if (Vz < -95.0) {
            return 0.990494;
        }
        else if (Vz < -85.0) {
            return 0.990065;
        }
        else if (Vz < -75.0) {
            return 0.990332;
        }
        else if (Vz < -65.0) {
            return 0.996478;
        }
        else if (Vz < -55.0) {
            return 0.999687;
        }
        else if (Vz < -45.0) {
            return 0.998645;
        }
        else if (Vz < -35.0) {
            return 0.993835;
        }
        else if (Vz < -25.0) {
            return 0.996273;
        }
        else if (Vz < -15.0) {
            return 0.998307;
        }
        else if (Vz < -5.0) {
            return 0.999295;
        }
        else if (Vz < 5.0) {
            return 1.0;
        }
        else if (Vz < 15.0) {
            return 1.00056;
        }
        else if (Vz < 25.0) {
            return 1.00019;
        }
        else if (Vz < 35.0) {
            return 0.999894;
        }
        else if (Vz < 45.0) {
            return 0.998907;
        }
        else if (Vz < 55.0) {
            return 1.00468;
        }
        else if (Vz < 65.0) {
            return 1.0055;
        }
        else if (Vz < 75.0) {
            return 1.00142;
        }
        else if (Vz < 85.0) {
            return 0.996247;
        }
        else if (Vz < 95.0) {
            return 0.995789;
        }
        else if (Vz < 105.0) {
            return 0.996513;
        }
        else if (Vz < 115.0) {
            return 0.996928;
        }
        else if (Vz < 125.0) {
            return 0.998196;
        }
        else if (Vz < 135.0) {
            return 1.00097;
        }
        else                {
            return 1.00699;
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
