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
#include <array>
#include <cmath>
using namespace std;

// 四体关联
// 使用这个编译：
// singularity exec -e --env DISPLAY=$DISPLAY -B /direct -B /gpfs -B /star -B /cvmfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh
// g++ -O2 -std=c++11 S_Four.cpp -o S_One `root-config --cflags --libs`

using Vec3 = std::array<double, 3>;
using Vec2 = std::array<double, 2>;

double dot(const Vec3& a, const Vec3& b);

Vec3 cross(const Vec3& a, const Vec3& b);

Vec3 normalize(const Vec3& v);

// 辅助函数：计算最优旋转后的投影坐标
bool computeRotatedProjections(const Vec3& b, const Vec3& c, const Vec3& d, const Vec3& n, double &ThetaB, double &ThetaC, double &ThetaD);

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
    ArmParticle(float _px, float _py, float _pz, float _mass, int _TreeID) 
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
    const ArmParticle& C,
    const ArmParticle& D,
    TH1D& H_P_tot,
    TH1D& H_beta,
    double& BthetaOut,
    double& CthetaOut,
    double& DthetaOut,
    double& BphiOut,
    double& CphiOut,
    double& DphiOut,
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
    std::vector<ArmParticle> C_particles;  // C类粒子
    std::vector<ArmParticle> D_particles;  // D类粒子
    
    Event() : eventID(-1) {}
    // 构造函数
    Event(int _eventID) 
        : eventID(_eventID) {}
};

void print(Event Temp);

void S_Four(
    TString MidName,
    TString DataName,
    int OutputFileIndex,
    TString OutMidName,
    int A_PDG,
    int B_PDG,
    int C_PDG,
    int D_PDG,
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

    bool Check_C_D = false;
    if (D_PDG == C_PDG) {
        Check_C_D = true;
    }

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
    int Aid , Bid , Cid , Did , Mid;
    float APx  , BPx  , CPx  , DPx ;
    float APy  , BPy  , CPy  , DPy ;
    float APz  , BPz  , CPz  , DPz ;
    float APt  , BPt  , CPt  , DPt ;
    float ARap , BRap , CRap , DRap;
    double ATheta , BTheta , CTheta , DTheta;
    float BMass = massList(B_PDG, DataName)           , AMass = massList(A_PDG, DataName)           , CMass = massList(C_PDG, DataName)           , DMass = massList(D_PDG, DataName)          ;
    float BMassSigma = massListSigma(B_PDG, DataName) , AMassSigma = massListSigma(A_PDG, DataName) , CMassSigma = massListSigma(C_PDG, DataName) , DMassSigma = massListSigma(D_PDG, DataName);
    std::vector<std::vector<int> > Mother_ParID;
    bool IsSame;
    float   P_B , kStar;
    double  Bphi , Btheta;
    double  Cphi , Ctheta;
    double  Dphi , Dtheta;
    enum MixType {
        SAME,       // ABCD
    
        ABC_D,      // ABC|D
        ABD_C,      // ABD|C
        ACD_B,      // ACD|B
        BCD_A,      // BCD|A
    
        AB_CD,      // AB|CD
        AC_BD,      // AC|BD
        AD_BC,      // AD|BC
    
        AB_C_D,     // AB|C|D
        AC_B_D,     // AC|B|D
        AD_B_C,     // AD|B|C
        BC_A_D,     // BC|A|D
        BD_A_C,     // BD|A|C
        CD_A_B,     // CD|A|B
    
        A_B_C_D     // A|B|C|D
    };

    std::vector<std::vector<std::vector<std::vector<Event>>>> EventPool;
    std::vector<std::vector<ArmParticle>> A_Array(CentralityBinNum);
    std::vector<ArmParticle> B_Array;
    std::vector<ArmParticle> C_Array;
    std::vector<std::vector<ArmParticle>> A_List(yBinNum);
    std::vector<ArmParticle> B_List;
    std::vector<ArmParticle> C_List;
    std::vector<ArmParticle> D_List;
    std::vector<TH1D*>                                               B_ALL_phi_ABCD             ;
    std::vector<TH1D*>                                               B_ALL_phi_A_B_C_D          ;
    std::vector<TH1D*>                                               B_ALL_phi_AB_CD            ;
    std::vector<TH1D*>                                               B_ALL_phi_AC_BD            ;
    std::vector<TH1D*>                                               B_ALL_phi_AD_BC            ;
    std::vector<TH1D*>                                               B_ALL_phi_ABC_D            ;
    std::vector<TH1D*>                                               B_ALL_phi_ABD_C            ;
    std::vector<TH1D*>                                               B_ALL_phi_ACD_B            ;
    std::vector<TH1D*>                                               B_ALL_phi_BCD_A            ;
    std::vector<TH1D*>                                               B_ALL_phi_AB_C_D           ;
    std::vector<TH1D*>                                               B_ALL_phi_AC_B_D           ;
    std::vector<TH1D*>                                               B_ALL_phi_AD_B_C           ;
    std::vector<TH1D*>                                               B_ALL_phi_BC_A_D           ;
    std::vector<TH1D*>                                               B_ALL_phi_BD_A_C           ;
    std::vector<TH1D*>                                               B_ALL_phi_CD_A_B           ;
    
    std::vector<TH1D*>                                               C_ALL_phi_ABCD             ;
    std::vector<TH1D*>                                               C_ALL_phi_A_B_C_D          ;
    std::vector<TH1D*>                                               C_ALL_phi_AB_CD            ;
    std::vector<TH1D*>                                               C_ALL_phi_AC_BD            ;
    std::vector<TH1D*>                                               C_ALL_phi_AD_BC            ;
    std::vector<TH1D*>                                               C_ALL_phi_ABC_D            ;
    std::vector<TH1D*>                                               C_ALL_phi_ABD_C            ;
    std::vector<TH1D*>                                               C_ALL_phi_ACD_B            ;
    std::vector<TH1D*>                                               C_ALL_phi_BCD_A            ;
    std::vector<TH1D*>                                               C_ALL_phi_AB_C_D           ;
    std::vector<TH1D*>                                               C_ALL_phi_AC_B_D           ;
    std::vector<TH1D*>                                               C_ALL_phi_AD_B_C           ;
    std::vector<TH1D*>                                               C_ALL_phi_BC_A_D           ;
    std::vector<TH1D*>                                               C_ALL_phi_BD_A_C           ;
    std::vector<TH1D*>                                               C_ALL_phi_CD_A_B           ;
    
    std::vector<TH1D*>                                               D_ALL_phi_ABCD             ;
    std::vector<TH1D*>                                               D_ALL_phi_A_B_C_D          ;
    std::vector<TH1D*>                                               D_ALL_phi_AB_CD            ;
    std::vector<TH1D*>                                               D_ALL_phi_AC_BD            ;
    std::vector<TH1D*>                                               D_ALL_phi_AD_BC            ;
    std::vector<TH1D*>                                               D_ALL_phi_ABC_D            ;
    std::vector<TH1D*>                                               D_ALL_phi_ABD_C            ;
    std::vector<TH1D*>                                               D_ALL_phi_ACD_B            ;
    std::vector<TH1D*>                                               D_ALL_phi_BCD_A            ;
    std::vector<TH1D*>                                               D_ALL_phi_AB_C_D           ;
    std::vector<TH1D*>                                               D_ALL_phi_AC_B_D           ;
    std::vector<TH1D*>                                               D_ALL_phi_AD_B_C           ;
    std::vector<TH1D*>                                               D_ALL_phi_BC_A_D           ;
    std::vector<TH1D*>                                               D_ALL_phi_BD_A_C           ;
    std::vector<TH1D*>                                               D_ALL_phi_CD_A_B           ;

    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_ABCD                 ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_AB_CD                ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_AC_BD                ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_AD_BC                ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_A_B_C_D              ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_ABC_D           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_ABD_C           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_ACD_B           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_BCD_A           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_AB_C_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_AC_B_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_AD_B_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_BC_A_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_BD_A_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_phi_CD_A_B          ;

    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_ABCD                 ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_AB_CD                ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_AC_BD                ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_AD_BC                ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_A_B_C_D              ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_ABC_D           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_ABD_C           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_ACD_B           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_BCD_A           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_AB_C_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_AC_B_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_AD_B_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_BC_A_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_BD_A_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_phi_CD_A_B          ;

    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_ABCD                 ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_AB_CD                ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_AC_BD                ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_AD_BC                ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_A_B_C_D              ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_ABC_D           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_ABD_C           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_ACD_B           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_BCD_A           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_AB_C_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_AC_B_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_AD_B_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_BC_A_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_BD_A_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_phi_CD_A_B          ;

    std::vector<TH1D*>                                               B_ALL_theta_ABCD         ;
    std::vector<TH1D*>                                               B_ALL_theta_AB_CD        ;
    std::vector<TH1D*>                                               B_ALL_theta_AC_BD        ;
    std::vector<TH1D*>                                               B_ALL_theta_AD_BC        ;
    std::vector<TH1D*>                                               B_ALL_theta_A_B_C_D      ;
    std::vector<TH1D*>                                               B_ALL_theta_ABC_D           ;
    std::vector<TH1D*>                                               B_ALL_theta_ABD_C           ;
    std::vector<TH1D*>                                               B_ALL_theta_ACD_B           ;
    std::vector<TH1D*>                                               B_ALL_theta_BCD_A           ;
    std::vector<TH1D*>                                               B_ALL_theta_AB_C_D          ;
    std::vector<TH1D*>                                               B_ALL_theta_AC_B_D          ;
    std::vector<TH1D*>                                               B_ALL_theta_AD_B_C          ;
    std::vector<TH1D*>                                               B_ALL_theta_BC_A_D          ;
    std::vector<TH1D*>                                               B_ALL_theta_BD_A_C          ;
    std::vector<TH1D*>                                               B_ALL_theta_CD_A_B          ;

    std::vector<TH1D*>                                               C_ALL_theta_ABCD         ;
    std::vector<TH1D*>                                               C_ALL_theta_AB_CD        ;
    std::vector<TH1D*>                                               C_ALL_theta_AC_BD        ;
    std::vector<TH1D*>                                               C_ALL_theta_AD_BC        ;
    std::vector<TH1D*>                                               C_ALL_theta_A_B_C_D      ;
    std::vector<TH1D*>                                               C_ALL_theta_ABC_D           ;
    std::vector<TH1D*>                                               C_ALL_theta_ABD_C           ;
    std::vector<TH1D*>                                               C_ALL_theta_ACD_B           ;
    std::vector<TH1D*>                                               C_ALL_theta_BCD_A           ;
    std::vector<TH1D*>                                               C_ALL_theta_AB_C_D          ;
    std::vector<TH1D*>                                               C_ALL_theta_AC_B_D          ;
    std::vector<TH1D*>                                               C_ALL_theta_AD_B_C          ;
    std::vector<TH1D*>                                               C_ALL_theta_BC_A_D          ;
    std::vector<TH1D*>                                               C_ALL_theta_BD_A_C          ;
    std::vector<TH1D*>                                               C_ALL_theta_CD_A_B          ;

    std::vector<TH1D*>                                               D_ALL_theta_ABCD         ;
    std::vector<TH1D*>                                               D_ALL_theta_AB_CD        ;
    std::vector<TH1D*>                                               D_ALL_theta_AC_BD        ;
    std::vector<TH1D*>                                               D_ALL_theta_AD_BC        ;
    std::vector<TH1D*>                                               D_ALL_theta_A_B_C_D      ;
    std::vector<TH1D*>                                               D_ALL_theta_ABC_D           ;
    std::vector<TH1D*>                                               D_ALL_theta_ABD_C           ;
    std::vector<TH1D*>                                               D_ALL_theta_ACD_B           ;
    std::vector<TH1D*>                                               D_ALL_theta_BCD_A           ;
    std::vector<TH1D*>                                               D_ALL_theta_AB_C_D          ;
    std::vector<TH1D*>                                               D_ALL_theta_AC_B_D          ;
    std::vector<TH1D*>                                               D_ALL_theta_AD_B_C          ;
    std::vector<TH1D*>                                               D_ALL_theta_BC_A_D          ;
    std::vector<TH1D*>                                               D_ALL_theta_BD_A_C          ;
    std::vector<TH1D*>                                               D_ALL_theta_CD_A_B          ;

    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_ABCD             ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_AB_CD            ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_AC_BD            ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_AD_BC            ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_A_B_C_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_ABC_D           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_ABD_C           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_ACD_B           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_BCD_A           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_AB_C_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_AC_B_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_AD_B_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_BC_A_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_BD_A_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     B_theta_CD_A_B          ;

    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_ABCD             ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_AB_CD            ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_AC_BD            ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_AD_BC            ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_A_B_C_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_ABC_D           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_ABD_C           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_ACD_B           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_BCD_A           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_AB_C_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_AC_B_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_AD_B_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_BC_A_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_BD_A_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     C_theta_CD_A_B          ;

    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_ABCD             ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_AB_CD            ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_AC_BD            ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_AD_BC            ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_A_B_C_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_ABC_D           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_ABD_C           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_ACD_B           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_BCD_A           ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_AB_C_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_AC_B_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_AD_B_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_BC_A_D          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_BD_A_C          ;
    std::vector<std::vector<std::vector<TH1D*>>>                     D_theta_CD_A_B          ;

    TH1D* H_P_tot     = new TH1D("H_P_tot","H_P_tot",200,0,10);
    TH1D* H_beta      = new TH1D("H_beta" ,"H_beta" ,500,0,2);

    if (RecordingMethod == 0) {
        EventPool.resize(CentralityBinNum);
        B_ALL_phi_ABCD       .resize(yBinNum, nullptr);
        B_ALL_phi_A_B_C_D    .resize(yBinNum, nullptr);
        B_ALL_phi_AB_CD     .resize(yBinNum, nullptr);
        B_ALL_phi_AC_BD     .resize(yBinNum, nullptr);
        B_ALL_phi_AD_BC     .resize(yBinNum, nullptr);
        B_ALL_phi_ABC_D     .resize(yBinNum, nullptr);
        B_ALL_phi_ABD_C     .resize(yBinNum, nullptr);
        B_ALL_phi_ACD_B     .resize(yBinNum, nullptr);
        B_ALL_phi_BCD_A     .resize(yBinNum, nullptr);
        B_ALL_phi_AB_C_D    .resize(yBinNum, nullptr);
        B_ALL_phi_AC_B_D    .resize(yBinNum, nullptr);
        B_ALL_phi_AD_B_C    .resize(yBinNum, nullptr);
        B_ALL_phi_BC_A_D    .resize(yBinNum, nullptr);
        B_ALL_phi_BD_A_C    .resize(yBinNum, nullptr);
        B_ALL_phi_CD_A_B    .resize(yBinNum, nullptr);
        B_phi_ABCD      .resize(CentralityBinNum);
        B_phi_AB_CD      .resize(CentralityBinNum);
        B_phi_AC_BD      .resize(CentralityBinNum);
        B_phi_AD_BC      .resize(CentralityBinNum);
        B_phi_A_B_C_D   .resize(CentralityBinNum);
        B_phi_ABC_D     .resize(CentralityBinNum);
        B_phi_ABD_C     .resize(CentralityBinNum);
        B_phi_ACD_B     .resize(CentralityBinNum);
        B_phi_BCD_A     .resize(CentralityBinNum);
        B_phi_AB_C_D    .resize(CentralityBinNum);
        B_phi_AC_B_D    .resize(CentralityBinNum);
        B_phi_AD_B_C    .resize(CentralityBinNum);
        B_phi_BC_A_D    .resize(CentralityBinNum);
        B_phi_BD_A_C    .resize(CentralityBinNum);
        B_phi_CD_A_B    .resize(CentralityBinNum);
        B_ALL_theta_ABCD       .resize(yBinNum, nullptr);
        B_ALL_theta_AB_CD      .resize(yBinNum, nullptr);
        B_ALL_theta_AC_BD      .resize(yBinNum, nullptr);
        B_ALL_theta_AD_BC      .resize(yBinNum, nullptr);
        B_ALL_theta_A_B_C_D    .resize(yBinNum, nullptr);
        B_ALL_theta_ABC_D      .resize(yBinNum, nullptr);
        B_ALL_theta_ABD_C      .resize(yBinNum, nullptr);
        B_ALL_theta_ACD_B      .resize(yBinNum, nullptr);
        B_ALL_theta_BCD_A      .resize(yBinNum, nullptr);
        B_ALL_theta_AB_C_D     .resize(yBinNum, nullptr);
        B_ALL_theta_AC_B_D     .resize(yBinNum, nullptr);
        B_ALL_theta_AD_B_C     .resize(yBinNum, nullptr);
        B_ALL_theta_BC_A_D     .resize(yBinNum, nullptr);
        B_ALL_theta_BD_A_C     .resize(yBinNum, nullptr);
        B_ALL_theta_CD_A_B     .resize(yBinNum, nullptr);
        B_theta_ABCD      .resize(CentralityBinNum);
        B_theta_AB_CD     .resize(CentralityBinNum);
        B_theta_AC_BD     .resize(CentralityBinNum);
        B_theta_AD_BC     .resize(CentralityBinNum);
        B_theta_A_B_C_D   .resize(CentralityBinNum);
        B_theta_ABC_D     .resize(CentralityBinNum);
        B_theta_ABD_C     .resize(CentralityBinNum);
        B_theta_ACD_B     .resize(CentralityBinNum);
        B_theta_BCD_A     .resize(CentralityBinNum);
        B_theta_AB_C_D    .resize(CentralityBinNum);
        B_theta_AC_B_D    .resize(CentralityBinNum);
        B_theta_AD_B_C    .resize(CentralityBinNum);
        B_theta_BC_A_D    .resize(CentralityBinNum);
        B_theta_BD_A_C    .resize(CentralityBinNum);
        B_theta_CD_A_B    .resize(CentralityBinNum);
        
        C_ALL_phi_ABCD       .resize(yBinNum, nullptr);
        C_ALL_phi_A_B_C_D    .resize(yBinNum, nullptr);
        C_ALL_phi_AB_CD     .resize(yBinNum, nullptr);
        C_ALL_phi_AC_BD     .resize(yBinNum, nullptr);
        C_ALL_phi_AD_BC     .resize(yBinNum, nullptr);
        C_ALL_phi_ABC_D     .resize(yBinNum, nullptr);
        C_ALL_phi_ABD_C     .resize(yBinNum, nullptr);
        C_ALL_phi_ACD_B     .resize(yBinNum, nullptr);
        C_ALL_phi_BCD_A     .resize(yBinNum, nullptr);
        C_ALL_phi_AB_C_D    .resize(yBinNum, nullptr);
        C_ALL_phi_AC_B_D    .resize(yBinNum, nullptr);
        C_ALL_phi_AD_B_C    .resize(yBinNum, nullptr);
        C_ALL_phi_BC_A_D    .resize(yBinNum, nullptr);
        C_ALL_phi_BD_A_C    .resize(yBinNum, nullptr);
        C_ALL_phi_CD_A_B    .resize(yBinNum, nullptr);
        C_phi_ABCD      .resize(CentralityBinNum);
        C_phi_AB_CD      .resize(CentralityBinNum);
        C_phi_AC_BD      .resize(CentralityBinNum);
        C_phi_AD_BC      .resize(CentralityBinNum);
        C_phi_A_B_C_D   .resize(CentralityBinNum);
        C_phi_ABC_D     .resize(CentralityBinNum);
        C_phi_ABD_C     .resize(CentralityBinNum);
        C_phi_ACD_B     .resize(CentralityBinNum);
        C_phi_BCD_A     .resize(CentralityBinNum);
        C_phi_AB_C_D    .resize(CentralityBinNum);
        C_phi_AC_B_D    .resize(CentralityBinNum);
        C_phi_AD_B_C    .resize(CentralityBinNum);
        C_phi_BC_A_D    .resize(CentralityBinNum);
        C_phi_BD_A_C    .resize(CentralityBinNum);
        C_phi_CD_A_B    .resize(CentralityBinNum);
        C_ALL_theta_ABCD       .resize(yBinNum, nullptr);
        C_ALL_theta_AB_CD      .resize(yBinNum, nullptr);
        C_ALL_theta_AC_BD      .resize(yBinNum, nullptr);
        C_ALL_theta_AD_BC      .resize(yBinNum, nullptr);
        C_ALL_theta_A_B_C_D    .resize(yBinNum, nullptr);
        C_ALL_theta_ABC_D      .resize(yBinNum, nullptr);
        C_ALL_theta_ABD_C      .resize(yBinNum, nullptr);
        C_ALL_theta_ACD_B      .resize(yBinNum, nullptr);
        C_ALL_theta_BCD_A      .resize(yBinNum, nullptr);
        C_ALL_theta_AB_C_D     .resize(yBinNum, nullptr);
        C_ALL_theta_AC_B_D     .resize(yBinNum, nullptr);
        C_ALL_theta_AD_B_C     .resize(yBinNum, nullptr);
        C_ALL_theta_BC_A_D     .resize(yBinNum, nullptr);
        C_ALL_theta_BD_A_C     .resize(yBinNum, nullptr);
        C_ALL_theta_CD_A_B     .resize(yBinNum, nullptr);
        C_theta_ABCD      .resize(CentralityBinNum);
        C_theta_AB_CD     .resize(CentralityBinNum);
        C_theta_AC_BD     .resize(CentralityBinNum);
        C_theta_AD_BC     .resize(CentralityBinNum);
        C_theta_A_B_C_D   .resize(CentralityBinNum);
        C_theta_ABC_D     .resize(CentralityBinNum);
        C_theta_ABD_C     .resize(CentralityBinNum);
        C_theta_ACD_B     .resize(CentralityBinNum);
        C_theta_BCD_A     .resize(CentralityBinNum);
        C_theta_AB_C_D    .resize(CentralityBinNum);
        C_theta_AC_B_D    .resize(CentralityBinNum);
        C_theta_AD_B_C    .resize(CentralityBinNum);
        C_theta_BC_A_D    .resize(CentralityBinNum);
        C_theta_BD_A_C    .resize(CentralityBinNum);
        C_theta_CD_A_B    .resize(CentralityBinNum);
        
        D_ALL_phi_ABCD       .resize(yBinNum, nullptr);
        D_ALL_phi_A_B_C_D    .resize(yBinNum, nullptr);
        D_ALL_phi_AB_CD     .resize(yBinNum, nullptr);
        D_ALL_phi_AC_BD     .resize(yBinNum, nullptr);
        D_ALL_phi_AD_BC     .resize(yBinNum, nullptr);
        D_ALL_phi_ABC_D     .resize(yBinNum, nullptr);
        D_ALL_phi_ABD_C     .resize(yBinNum, nullptr);
        D_ALL_phi_ACD_B     .resize(yBinNum, nullptr);
        D_ALL_phi_BCD_A     .resize(yBinNum, nullptr);
        D_ALL_phi_AB_C_D    .resize(yBinNum, nullptr);
        D_ALL_phi_AC_B_D    .resize(yBinNum, nullptr);
        D_ALL_phi_AD_B_C    .resize(yBinNum, nullptr);
        D_ALL_phi_BC_A_D    .resize(yBinNum, nullptr);
        D_ALL_phi_BD_A_C    .resize(yBinNum, nullptr);
        D_ALL_phi_CD_A_B    .resize(yBinNum, nullptr);
        D_phi_ABCD      .resize(CentralityBinNum);
        D_phi_AB_CD      .resize(CentralityBinNum);
        D_phi_AC_BD      .resize(CentralityBinNum);
        D_phi_AD_BC      .resize(CentralityBinNum);
        D_phi_A_B_C_D   .resize(CentralityBinNum);
        D_phi_ABC_D     .resize(CentralityBinNum);
        D_phi_ABD_C     .resize(CentralityBinNum);
        D_phi_ACD_B     .resize(CentralityBinNum);
        D_phi_BCD_A     .resize(CentralityBinNum);
        D_phi_AB_C_D    .resize(CentralityBinNum);
        D_phi_AC_B_D    .resize(CentralityBinNum);
        D_phi_AD_B_C    .resize(CentralityBinNum);
        D_phi_BC_A_D    .resize(CentralityBinNum);
        D_phi_BD_A_C    .resize(CentralityBinNum);
        D_phi_CD_A_B    .resize(CentralityBinNum);
        D_ALL_theta_ABCD       .resize(yBinNum, nullptr);
        D_ALL_theta_AB_CD      .resize(yBinNum, nullptr);
        D_ALL_theta_AC_BD      .resize(yBinNum, nullptr);
        D_ALL_theta_AD_BC      .resize(yBinNum, nullptr);
        D_ALL_theta_A_B_C_D    .resize(yBinNum, nullptr);
        D_ALL_theta_ABC_D      .resize(yBinNum, nullptr);
        D_ALL_theta_ABD_C      .resize(yBinNum, nullptr);
        D_ALL_theta_ACD_B      .resize(yBinNum, nullptr);
        D_ALL_theta_BCD_A      .resize(yBinNum, nullptr);
        D_ALL_theta_AB_C_D     .resize(yBinNum, nullptr);
        D_ALL_theta_AC_B_D     .resize(yBinNum, nullptr);
        D_ALL_theta_AD_B_C     .resize(yBinNum, nullptr);
        D_ALL_theta_BC_A_D     .resize(yBinNum, nullptr);
        D_ALL_theta_BD_A_C     .resize(yBinNum, nullptr);
        D_ALL_theta_CD_A_B     .resize(yBinNum, nullptr);
        D_theta_ABCD      .resize(CentralityBinNum);
        D_theta_AB_CD     .resize(CentralityBinNum);
        D_theta_AC_BD     .resize(CentralityBinNum);
        D_theta_AD_BC     .resize(CentralityBinNum);
        D_theta_A_B_C_D   .resize(CentralityBinNum);
        D_theta_ABC_D     .resize(CentralityBinNum);
        D_theta_ABD_C     .resize(CentralityBinNum);
        D_theta_ACD_B     .resize(CentralityBinNum);
        D_theta_BCD_A     .resize(CentralityBinNum);
        D_theta_AB_C_D    .resize(CentralityBinNum);
        D_theta_AC_B_D    .resize(CentralityBinNum);
        D_theta_AD_B_C    .resize(CentralityBinNum);
        D_theta_BC_A_D    .resize(CentralityBinNum);
        D_theta_BD_A_C    .resize(CentralityBinNum);
        D_theta_CD_A_B    .resize(CentralityBinNum);
        for (i = 0; i < CentralityBinNum; i++) {
            EventPool[i].resize(yBinNum);
            B_phi_ABCD           [i].resize(yBinNum);
            B_phi_AB_CD           [i].resize(yBinNum);
            B_phi_AC_BD           [i].resize(yBinNum);
            B_phi_AD_BC           [i].resize(yBinNum);
            B_phi_A_B_C_D       [i].resize(yBinNum);
            B_phi_ABC_D           [i].resize(yBinNum);
            B_phi_ABD_C           [i].resize(yBinNum);
            B_phi_ACD_B           [i].resize(yBinNum);
            B_phi_BCD_A           [i].resize(yBinNum);
            B_phi_AB_C_D          [i].resize(yBinNum);
            B_phi_AC_B_D          [i].resize(yBinNum);
            B_phi_AD_B_C          [i].resize(yBinNum);
            B_phi_BC_A_D          [i].resize(yBinNum);
            B_phi_BD_A_C          [i].resize(yBinNum);
            B_phi_CD_A_B          [i].resize(yBinNum);
            B_theta_ABCD       [i].resize(yBinNum);
            B_theta_AB_CD      [i].resize(yBinNum);
            B_theta_AC_BD      [i].resize(yBinNum);
            B_theta_AD_BC      [i].resize(yBinNum);
            B_theta_A_B_C_D   [i].resize(yBinNum);
            B_theta_ABC_D      [i].resize(yBinNum);
            B_theta_ABD_C      [i].resize(yBinNum);
            B_theta_ACD_B      [i].resize(yBinNum);
            B_theta_BCD_A      [i].resize(yBinNum);
            B_theta_AB_C_D     [i].resize(yBinNum);
            B_theta_AC_B_D     [i].resize(yBinNum);
            B_theta_AD_B_C     [i].resize(yBinNum);
            B_theta_BC_A_D     [i].resize(yBinNum);
            B_theta_BD_A_C     [i].resize(yBinNum);
            B_theta_CD_A_B     [i].resize(yBinNum);
            
            C_phi_ABCD           [i].resize(yBinNum);
            C_phi_AB_CD           [i].resize(yBinNum);
            C_phi_AC_BD           [i].resize(yBinNum);
            C_phi_AD_BC           [i].resize(yBinNum);
            C_phi_A_B_C_D       [i].resize(yBinNum);
            C_phi_ABC_D           [i].resize(yBinNum);
            C_phi_ABD_C           [i].resize(yBinNum);
            C_phi_ACD_B           [i].resize(yBinNum);
            C_phi_BCD_A           [i].resize(yBinNum);
            C_phi_AB_C_D          [i].resize(yBinNum);
            C_phi_AC_B_D          [i].resize(yBinNum);
            C_phi_AD_B_C          [i].resize(yBinNum);
            C_phi_BC_A_D          [i].resize(yBinNum);
            C_phi_BD_A_C          [i].resize(yBinNum);
            C_phi_CD_A_B          [i].resize(yBinNum);
            C_theta_ABCD       [i].resize(yBinNum);
            C_theta_AB_CD      [i].resize(yBinNum);
            C_theta_AC_BD      [i].resize(yBinNum);
            C_theta_AD_BC      [i].resize(yBinNum);
            C_theta_A_B_C_D   [i].resize(yBinNum);
            C_theta_ABC_D      [i].resize(yBinNum);
            C_theta_ABD_C      [i].resize(yBinNum);
            C_theta_ACD_B      [i].resize(yBinNum);
            C_theta_BCD_A      [i].resize(yBinNum);
            C_theta_AB_C_D     [i].resize(yBinNum);
            C_theta_AC_B_D     [i].resize(yBinNum);
            C_theta_AD_B_C     [i].resize(yBinNum);
            C_theta_BC_A_D     [i].resize(yBinNum);
            C_theta_BD_A_C     [i].resize(yBinNum);
            C_theta_CD_A_B     [i].resize(yBinNum);
            
            D_phi_ABCD           [i].resize(yBinNum);
            D_phi_AB_CD           [i].resize(yBinNum);
            D_phi_AC_BD           [i].resize(yBinNum);
            D_phi_AD_BC           [i].resize(yBinNum);
            D_phi_A_B_C_D       [i].resize(yBinNum);
            D_phi_ABC_D           [i].resize(yBinNum);
            D_phi_ABD_C           [i].resize(yBinNum);
            D_phi_ACD_B           [i].resize(yBinNum);
            D_phi_BCD_A           [i].resize(yBinNum);
            D_phi_AB_C_D          [i].resize(yBinNum);
            D_phi_AC_B_D          [i].resize(yBinNum);
            D_phi_AD_B_C          [i].resize(yBinNum);
            D_phi_BC_A_D          [i].resize(yBinNum);
            D_phi_BD_A_C          [i].resize(yBinNum);
            D_phi_CD_A_B          [i].resize(yBinNum);
            D_theta_ABCD       [i].resize(yBinNum);
            D_theta_AB_CD      [i].resize(yBinNum);
            D_theta_AC_BD      [i].resize(yBinNum);
            D_theta_AD_BC      [i].resize(yBinNum);
            D_theta_A_B_C_D   [i].resize(yBinNum);
            D_theta_ABC_D      [i].resize(yBinNum);
            D_theta_ABD_C      [i].resize(yBinNum);
            D_theta_ACD_B      [i].resize(yBinNum);
            D_theta_BCD_A      [i].resize(yBinNum);
            D_theta_AB_C_D     [i].resize(yBinNum);
            D_theta_AC_B_D     [i].resize(yBinNum);
            D_theta_AD_B_C     [i].resize(yBinNum);
            D_theta_BC_A_D     [i].resize(yBinNum);
            D_theta_BD_A_C     [i].resize(yBinNum);
            D_theta_CD_A_B     [i].resize(yBinNum);
            for (j = 0; j < yBinNum; j++) {
                EventPool[i][j].resize(PVzBinNum);
                B_phi_ABCD                     [i][j].resize(PVzBinNum, nullptr);
                B_phi_AB_CD                    [i][j].resize(PVzBinNum, nullptr);
                B_phi_AC_BD                    [i][j].resize(PVzBinNum, nullptr);
                B_phi_AD_BC                    [i][j].resize(PVzBinNum, nullptr);
                B_phi_A_B_C_D                  [i][j].resize(PVzBinNum, nullptr);
                B_phi_ABC_D                    [i][j].resize(PVzBinNum, nullptr);
                B_phi_ABD_C                    [i][j].resize(PVzBinNum, nullptr);
                B_phi_ACD_B                    [i][j].resize(PVzBinNum, nullptr);
                B_phi_BCD_A                    [i][j].resize(PVzBinNum, nullptr);
                B_phi_AB_C_D                   [i][j].resize(PVzBinNum, nullptr);
                B_phi_AC_B_D                   [i][j].resize(PVzBinNum, nullptr);
                B_phi_AD_B_C                   [i][j].resize(PVzBinNum, nullptr);
                B_phi_BC_A_D                   [i][j].resize(PVzBinNum, nullptr);
                B_phi_BD_A_C                   [i][j].resize(PVzBinNum, nullptr);
                B_phi_CD_A_B                   [i][j].resize(PVzBinNum, nullptr);
                B_theta_ABCD                   [i][j].resize(PVzBinNum, nullptr);
                B_theta_AB_CD                  [i][j].resize(PVzBinNum, nullptr);
                B_theta_AC_BD                  [i][j].resize(PVzBinNum, nullptr);
                B_theta_AD_BC                  [i][j].resize(PVzBinNum, nullptr);
                B_theta_A_B_C_D                [i][j].resize(PVzBinNum, nullptr);
                B_theta_ABC_D                  [i][j].resize(PVzBinNum, nullptr);
                B_theta_ABD_C                  [i][j].resize(PVzBinNum, nullptr);
                B_theta_ACD_B                  [i][j].resize(PVzBinNum, nullptr);
                B_theta_BCD_A                  [i][j].resize(PVzBinNum, nullptr);
                B_theta_AB_C_D                 [i][j].resize(PVzBinNum, nullptr);
                B_theta_AC_B_D                 [i][j].resize(PVzBinNum, nullptr);
                B_theta_AD_B_C                 [i][j].resize(PVzBinNum, nullptr);
                B_theta_BC_A_D                 [i][j].resize(PVzBinNum, nullptr);
                B_theta_BD_A_C                 [i][j].resize(PVzBinNum, nullptr);
                B_theta_CD_A_B                 [i][j].resize(PVzBinNum, nullptr);
                
                C_phi_ABCD                     [i][j].resize(PVzBinNum, nullptr);
                C_phi_AB_CD                    [i][j].resize(PVzBinNum, nullptr);
                C_phi_AC_BD                    [i][j].resize(PVzBinNum, nullptr);
                C_phi_AD_BC                    [i][j].resize(PVzBinNum, nullptr);
                C_phi_A_B_C_D                  [i][j].resize(PVzBinNum, nullptr);
                C_phi_ABC_D                    [i][j].resize(PVzBinNum, nullptr);
                C_phi_ABD_C                    [i][j].resize(PVzBinNum, nullptr);
                C_phi_ACD_B                    [i][j].resize(PVzBinNum, nullptr);
                C_phi_BCD_A                    [i][j].resize(PVzBinNum, nullptr);
                C_phi_AB_C_D                   [i][j].resize(PVzBinNum, nullptr);
                C_phi_AC_B_D                   [i][j].resize(PVzBinNum, nullptr);
                C_phi_AD_B_C                   [i][j].resize(PVzBinNum, nullptr);
                C_phi_BC_A_D                   [i][j].resize(PVzBinNum, nullptr);
                C_phi_BD_A_C                   [i][j].resize(PVzBinNum, nullptr);
                C_phi_CD_A_B                   [i][j].resize(PVzBinNum, nullptr);
                C_theta_ABCD                   [i][j].resize(PVzBinNum, nullptr);
                C_theta_AB_CD                  [i][j].resize(PVzBinNum, nullptr);
                C_theta_AC_BD                  [i][j].resize(PVzBinNum, nullptr);
                C_theta_AD_BC                  [i][j].resize(PVzBinNum, nullptr);
                C_theta_A_B_C_D                [i][j].resize(PVzBinNum, nullptr);
                C_theta_ABC_D                  [i][j].resize(PVzBinNum, nullptr);
                C_theta_ABD_C                  [i][j].resize(PVzBinNum, nullptr);
                C_theta_ACD_B                  [i][j].resize(PVzBinNum, nullptr);
                C_theta_BCD_A                  [i][j].resize(PVzBinNum, nullptr);
                C_theta_AB_C_D                 [i][j].resize(PVzBinNum, nullptr);
                C_theta_AC_B_D                 [i][j].resize(PVzBinNum, nullptr);
                C_theta_AD_B_C                 [i][j].resize(PVzBinNum, nullptr);
                C_theta_BC_A_D                 [i][j].resize(PVzBinNum, nullptr);
                C_theta_BD_A_C                 [i][j].resize(PVzBinNum, nullptr);
                C_theta_CD_A_B                 [i][j].resize(PVzBinNum, nullptr);
                
                D_phi_ABCD                     [i][j].resize(PVzBinNum, nullptr);
                D_phi_AB_CD                    [i][j].resize(PVzBinNum, nullptr);
                D_phi_AC_BD                    [i][j].resize(PVzBinNum, nullptr);
                D_phi_AD_BC                    [i][j].resize(PVzBinNum, nullptr);
                D_phi_A_B_C_D                  [i][j].resize(PVzBinNum, nullptr);
                D_phi_ABC_D                    [i][j].resize(PVzBinNum, nullptr);
                D_phi_ABD_C                    [i][j].resize(PVzBinNum, nullptr);
                D_phi_ACD_B                    [i][j].resize(PVzBinNum, nullptr);
                D_phi_BCD_A                    [i][j].resize(PVzBinNum, nullptr);
                D_phi_AB_C_D                   [i][j].resize(PVzBinNum, nullptr);
                D_phi_AC_B_D                   [i][j].resize(PVzBinNum, nullptr);
                D_phi_AD_B_C                   [i][j].resize(PVzBinNum, nullptr);
                D_phi_BC_A_D                   [i][j].resize(PVzBinNum, nullptr);
                D_phi_BD_A_C                   [i][j].resize(PVzBinNum, nullptr);
                D_phi_CD_A_B                   [i][j].resize(PVzBinNum, nullptr);
                D_theta_ABCD                   [i][j].resize(PVzBinNum, nullptr);
                D_theta_AB_CD                  [i][j].resize(PVzBinNum, nullptr);
                D_theta_AC_BD                  [i][j].resize(PVzBinNum, nullptr);
                D_theta_AD_BC                  [i][j].resize(PVzBinNum, nullptr);
                D_theta_A_B_C_D                [i][j].resize(PVzBinNum, nullptr);
                D_theta_ABC_D                  [i][j].resize(PVzBinNum, nullptr);
                D_theta_ABD_C                  [i][j].resize(PVzBinNum, nullptr);
                D_theta_ACD_B                  [i][j].resize(PVzBinNum, nullptr);
                D_theta_BCD_A                  [i][j].resize(PVzBinNum, nullptr);
                D_theta_AB_C_D                 [i][j].resize(PVzBinNum, nullptr);
                D_theta_AC_B_D                 [i][j].resize(PVzBinNum, nullptr);
                D_theta_AD_B_C                 [i][j].resize(PVzBinNum, nullptr);
                D_theta_BC_A_D                 [i][j].resize(PVzBinNum, nullptr);
                D_theta_BD_A_C                 [i][j].resize(PVzBinNum, nullptr);
                D_theta_CD_A_B                 [i][j].resize(PVzBinNum, nullptr);
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
                    B_phi_ABCD           [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_ABCD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_AB_CD          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_AB_CD_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_AC_BD          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_AC_BD_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_AD_BC          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_AD_BC_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_A_B_C_D        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_A_B_C_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_ABC_D          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_ABC_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_ABD_C          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_ABD_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_ACD_B          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_ACD_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_BCD_A          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_BCD_A_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_AB_C_D         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_AB_C_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_AC_B_D         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_AC_B_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_AD_B_C         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_AD_B_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_BC_A_D         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_BC_A_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_BD_A_C         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_BD_A_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_phi_CD_A_B         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_phi_CD_A_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    B_theta_ABCD         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_ABCD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_AB_CD        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_AB_CD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_AC_BD        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_AC_BD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_AD_BC        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_AD_BC_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_A_B_C_D      [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_A_B_C_D_%d_%d_%d" ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_ABC_D        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_ABC_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_ABD_C        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_ABD_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_ACD_B        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_ACD_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_BCD_A        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_BCD_A_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_AB_C_D       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_AB_C_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_AC_B_D       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_AC_B_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_AD_B_C       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_AD_B_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_BC_A_D       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_BC_A_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_BD_A_C       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_BD_A_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_theta_CD_A_B       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("B_theta_CD_A_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    B_phi_ABCD           [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_AB_CD          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_AC_BD          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_AD_BC          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_A_B_C_D        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_ABC_D          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_ABD_C          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_ACD_B          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_BCD_A          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_AB_C_D         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_AC_B_D         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_AD_B_C         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_BC_A_D         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_BD_A_C         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_phi_CD_A_B         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B phi");
                    B_theta_ABCD         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_AB_CD        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_AC_BD        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_AD_BC        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_A_B_C_D      [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_ABC_D        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_ABD_C        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_ACD_B        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_BCD_A        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_AB_C_D       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_AC_B_D       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_AD_B_C       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_BC_A_D       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_BD_A_C       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");
                    B_theta_CD_A_B       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("B theta");

                    C_phi_ABCD           [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_ABCD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_AB_CD          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_AB_CD_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_AC_BD          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_AC_BD_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_AD_BC          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_AD_BC_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_A_B_C_D        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_A_B_C_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_ABC_D          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_ABC_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_ABD_C          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_ABD_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_ACD_B          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_ACD_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_BCD_A          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_BCD_A_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_AB_C_D         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_AB_C_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_AC_B_D         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_AC_B_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_AD_B_C         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_AD_B_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_BC_A_D         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_BC_A_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_BD_A_C         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_BD_A_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_phi_CD_A_B         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_phi_CD_A_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    C_theta_ABCD         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_ABCD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_AB_CD        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_AB_CD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_AC_BD        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_AC_BD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_AD_BC        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_AD_BC_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_A_B_C_D      [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_A_B_C_D_%d_%d_%d" ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_ABC_D        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_ABC_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_ABD_C        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_ABD_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_ACD_B        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_ACD_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_BCD_A        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_BCD_A_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_AB_C_D       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_AB_C_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_AC_B_D       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_AC_B_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_AD_B_C       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_AD_B_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_BC_A_D       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_BC_A_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_BD_A_C       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_BD_A_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_theta_CD_A_B       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("C_theta_CD_A_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    C_phi_ABCD           [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_AB_CD          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_AC_BD          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_AD_BC          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_A_B_C_D        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_ABC_D          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_ABD_C          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_ACD_B          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_BCD_A          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_AB_C_D         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_AC_B_D         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_AD_B_C         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_BC_A_D         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_BD_A_C         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_phi_CD_A_B         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C phi");
                    C_theta_ABCD         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_AB_CD        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_AC_BD        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_AD_BC        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_A_B_C_D      [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_ABC_D        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_ABD_C        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_ACD_B        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_BCD_A        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_AB_C_D       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_AC_B_D       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_AD_B_C       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_BC_A_D       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_BD_A_C       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");
                    C_theta_CD_A_B       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("C theta");

                    D_phi_ABCD           [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_ABCD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_AB_CD          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_AB_CD_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_AC_BD          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_AC_BD_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_AD_BC          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_AD_BC_%d_%d_%d"    ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_A_B_C_D        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_A_B_C_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_ABC_D          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_ABC_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_ABD_C          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_ABD_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_ACD_B          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_ACD_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_BCD_A          [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_BCD_A_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_AB_C_D         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_AB_C_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_AC_B_D         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_AC_B_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_AD_B_C         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_AD_B_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_BC_A_D         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_BC_A_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_BD_A_C         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_BD_A_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_phi_CD_A_B         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_phi_CD_A_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,SideSta,SideEnd);
                    D_theta_ABCD         [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_ABCD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_AB_CD        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_AB_CD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_AC_BD        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_AC_BD_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_AD_BC        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_AD_BC_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_A_B_C_D      [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_A_B_C_D_%d_%d_%d" ,CenIndex,RapIndex,PVzIndex), Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"   ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_ABC_D        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_ABC_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_ABD_C        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_ABD_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_ACD_B        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_ACD_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_BCD_A        [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_BCD_A_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_AB_C_D       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_AB_C_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_AC_B_D       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_AC_B_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_AD_B_C       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_AD_B_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_BC_A_D       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_BC_A_D_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_BD_A_C       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_BD_A_C_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_theta_CD_A_B       [CenIndex] [RapIndex] [PVzIndex] = new TH1D(Form("D_theta_CD_A_B_%d_%d_%d"     ,CenIndex,RapIndex,PVzIndex),Form("[%d,%d]/100, %f<A_y<%f, %f<PV_z<%f"       ,CentralityBin[CenIndex],CentralityBin[CenIndex+1],yBin[RapIndex],yBin[RapIndex+1],PVzBin[PVzIndex],PVzBin[PVzIndex+1]),SideBinNum,-1,1);
                    D_phi_ABCD           [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_AB_CD          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_AC_BD          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_AD_BC          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_A_B_C_D        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_ABC_D          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_ABD_C          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_ACD_B          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_BCD_A          [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_AB_C_D         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_AC_B_D         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_AD_B_C         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_BC_A_D         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_BD_A_C         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_phi_CD_A_B         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D phi");
                    D_theta_ABCD         [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_AB_CD        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_AC_BD        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_AD_BC        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_A_B_C_D      [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_ABC_D        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_ABD_C        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_ACD_B        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_BCD_A        [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_AB_C_D       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_AC_B_D       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_AD_B_C       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_BC_A_D       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_BD_A_C       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                    D_theta_CD_A_B       [CenIndex] [RapIndex] [PVzIndex]->GetXaxis()->SetTitle("D theta");
                }
            }
            B_ALL_phi_ABCD               [RapIndex] = new TH1D(Form("B_ALL_phi_ABCD_%d"      ,          RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_A_B_C_D            [RapIndex] = new TH1D(Form("B_ALL_phi_A_B_C_D_%d"  ,          RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_AB_CD              [RapIndex] = new TH1D(Form("B_ALL_phi_AB_CD_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_AC_BD              [RapIndex] = new TH1D(Form("B_ALL_phi_AC_BD_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_AD_BC              [RapIndex] = new TH1D(Form("B_ALL_phi_AD_BC_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_ABC_D              [RapIndex] = new TH1D(Form("B_ALL_phi_ABC_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_ABD_C              [RapIndex] = new TH1D(Form("B_ALL_phi_ABD_C_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_ACD_B              [RapIndex] = new TH1D(Form("B_ALL_phi_ACD_B_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_BCD_A              [RapIndex] = new TH1D(Form("B_ALL_phi_BCD_A_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_AB_C_D             [RapIndex] = new TH1D(Form("B_ALL_phi_AB_C_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_AC_B_D             [RapIndex] = new TH1D(Form("B_ALL_phi_AC_B_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_AD_B_C             [RapIndex] = new TH1D(Form("B_ALL_phi_AD_B_C_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_BC_A_D             [RapIndex] = new TH1D(Form("B_ALL_phi_BC_A_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_BD_A_C             [RapIndex] = new TH1D(Form("B_ALL_phi_BD_A_C_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_phi_CD_A_B             [RapIndex] = new TH1D(Form("B_ALL_phi_CD_A_B_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            B_ALL_theta_ABCD             [RapIndex] = new TH1D(Form("B_ALL_theta_ABCD_%d"      ,      RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_AB_CD            [RapIndex] = new TH1D(Form("B_ALL_theta_AB_CD_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_AC_BD            [RapIndex] = new TH1D(Form("B_ALL_theta_AC_BD_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_AD_BC            [RapIndex] = new TH1D(Form("B_ALL_theta_AD_BC_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_A_B_C_D          [RapIndex] = new TH1D(Form("B_ALL_theta_A_B_C_D_%d"  ,      RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_ABC_D            [RapIndex] = new TH1D(Form("B_ALL_theta_ABC_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_ABD_C            [RapIndex] = new TH1D(Form("B_ALL_theta_ABD_C_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_ACD_B            [RapIndex] = new TH1D(Form("B_ALL_theta_ACD_B_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_BCD_A            [RapIndex] = new TH1D(Form("B_ALL_theta_BCD_A_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_AB_C_D           [RapIndex] = new TH1D(Form("B_ALL_theta_AB_C_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_AC_B_D           [RapIndex] = new TH1D(Form("B_ALL_theta_AC_B_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_AD_B_C           [RapIndex] = new TH1D(Form("B_ALL_theta_AD_B_C_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_BC_A_D           [RapIndex] = new TH1D(Form("B_ALL_theta_BC_A_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_BD_A_C           [RapIndex] = new TH1D(Form("B_ALL_theta_BD_A_C_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_theta_CD_A_B           [RapIndex] = new TH1D(Form("B_ALL_theta_CD_A_B_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            B_ALL_phi_ABCD               [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_A_B_C_D            [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_AB_CD              [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_AC_BD              [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_AD_BC              [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_ABC_D              [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_ABD_C              [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_ACD_B              [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_BCD_A              [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_AB_C_D             [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_AC_B_D             [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_AD_B_C             [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_BC_A_D             [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_BD_A_C             [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_phi_CD_A_B             [RapIndex]->GetXaxis()->SetTitle("B phi");
            B_ALL_theta_ABCD             [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_AB_CD            [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_AC_BD            [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_AD_BC            [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_A_B_C_D          [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_ABC_D            [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_ABD_C            [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_ACD_B            [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_BCD_A            [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_AB_C_D           [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_AC_B_D           [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_AD_B_C           [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_BC_A_D           [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_BD_A_C           [RapIndex]->GetXaxis()->SetTitle("B theta");
            B_ALL_theta_CD_A_B           [RapIndex]->GetXaxis()->SetTitle("B theta");

            C_ALL_phi_ABCD               [RapIndex] = new TH1D(Form("C_ALL_phi_ABCD_%d"      ,          RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_A_B_C_D            [RapIndex] = new TH1D(Form("C_ALL_phi_A_B_C_D_%d"  ,          RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_AB_CD              [RapIndex] = new TH1D(Form("C_ALL_phi_AB_CD_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_AC_BD              [RapIndex] = new TH1D(Form("C_ALL_phi_AC_BD_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_AD_BC              [RapIndex] = new TH1D(Form("C_ALL_phi_AD_BC_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_ABC_D              [RapIndex] = new TH1D(Form("C_ALL_phi_ABC_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_ABD_C              [RapIndex] = new TH1D(Form("C_ALL_phi_ABD_C_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_ACD_B              [RapIndex] = new TH1D(Form("C_ALL_phi_ACD_B_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_BCD_A              [RapIndex] = new TH1D(Form("C_ALL_phi_BCD_A_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_AB_C_D             [RapIndex] = new TH1D(Form("C_ALL_phi_AB_C_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_AC_B_D             [RapIndex] = new TH1D(Form("C_ALL_phi_AC_B_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_AD_B_C             [RapIndex] = new TH1D(Form("C_ALL_phi_AD_B_C_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_BC_A_D             [RapIndex] = new TH1D(Form("C_ALL_phi_BC_A_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_BD_A_C             [RapIndex] = new TH1D(Form("C_ALL_phi_BD_A_C_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_phi_CD_A_B             [RapIndex] = new TH1D(Form("C_ALL_phi_CD_A_B_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            C_ALL_theta_ABCD             [RapIndex] = new TH1D(Form("C_ALL_theta_ABCD_%d"      ,      RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_AB_CD            [RapIndex] = new TH1D(Form("C_ALL_theta_AB_CD_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_AC_BD            [RapIndex] = new TH1D(Form("C_ALL_theta_AC_BD_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_AD_BC            [RapIndex] = new TH1D(Form("C_ALL_theta_AD_BC_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_A_B_C_D          [RapIndex] = new TH1D(Form("C_ALL_theta_A_B_C_D_%d"  ,      RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_ABC_D            [RapIndex] = new TH1D(Form("C_ALL_theta_ABC_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_ABD_C            [RapIndex] = new TH1D(Form("C_ALL_theta_ABD_C_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_ACD_B            [RapIndex] = new TH1D(Form("C_ALL_theta_ACD_B_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_BCD_A            [RapIndex] = new TH1D(Form("C_ALL_theta_BCD_A_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_AB_C_D           [RapIndex] = new TH1D(Form("C_ALL_theta_AB_C_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_AC_B_D           [RapIndex] = new TH1D(Form("C_ALL_theta_AC_B_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_AD_B_C           [RapIndex] = new TH1D(Form("C_ALL_theta_AD_B_C_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_BC_A_D           [RapIndex] = new TH1D(Form("C_ALL_theta_BC_A_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_BD_A_C           [RapIndex] = new TH1D(Form("C_ALL_theta_BD_A_C_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_theta_CD_A_B           [RapIndex] = new TH1D(Form("C_ALL_theta_CD_A_B_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            C_ALL_phi_ABCD               [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_A_B_C_D            [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_AB_CD              [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_AC_BD              [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_AD_BC              [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_ABC_D              [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_ABD_C              [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_ACD_B              [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_BCD_A              [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_AB_C_D             [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_AC_B_D             [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_AD_B_C             [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_BC_A_D             [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_BD_A_C             [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_phi_CD_A_B             [RapIndex]->GetXaxis()->SetTitle("C phi");
            C_ALL_theta_ABCD             [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_AB_CD            [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_AC_BD            [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_AD_BC            [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_A_B_C_D          [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_ABC_D            [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_ABD_C            [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_ACD_B            [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_BCD_A            [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_AB_C_D           [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_AC_B_D           [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_AD_B_C           [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_BC_A_D           [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_BD_A_C           [RapIndex]->GetXaxis()->SetTitle("C theta");
            C_ALL_theta_CD_A_B           [RapIndex]->GetXaxis()->SetTitle("C theta");
            
            D_ALL_phi_ABCD               [RapIndex] = new TH1D(Form("D_ALL_phi_ABCD_%d"      ,          RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_A_B_C_D            [RapIndex] = new TH1D(Form("D_ALL_phi_A_B_C_D_%d"  ,          RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_AB_CD              [RapIndex] = new TH1D(Form("D_ALL_phi_AB_CD_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_AC_BD              [RapIndex] = new TH1D(Form("D_ALL_phi_AC_BD_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_AD_BC              [RapIndex] = new TH1D(Form("D_ALL_phi_AD_BC_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_ABC_D              [RapIndex] = new TH1D(Form("D_ALL_phi_ABC_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_ABD_C              [RapIndex] = new TH1D(Form("D_ALL_phi_ABD_C_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_ACD_B              [RapIndex] = new TH1D(Form("D_ALL_phi_ACD_B_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_BCD_A              [RapIndex] = new TH1D(Form("D_ALL_phi_BCD_A_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_AB_C_D             [RapIndex] = new TH1D(Form("D_ALL_phi_AB_C_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_AC_B_D             [RapIndex] = new TH1D(Form("D_ALL_phi_AC_B_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_AD_B_C             [RapIndex] = new TH1D(Form("D_ALL_phi_AD_B_C_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_BC_A_D             [RapIndex] = new TH1D(Form("D_ALL_phi_BC_A_D_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_BD_A_C             [RapIndex] = new TH1D(Form("D_ALL_phi_BD_A_C_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_phi_CD_A_B             [RapIndex] = new TH1D(Form("D_ALL_phi_CD_A_B_%d"  ,           RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,SideSta,SideEnd);
            D_ALL_theta_ABCD             [RapIndex] = new TH1D(Form("D_ALL_theta_ABCD_%d"      ,      RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_AB_CD            [RapIndex] = new TH1D(Form("D_ALL_theta_AB_CD_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_AC_BD            [RapIndex] = new TH1D(Form("D_ALL_theta_AC_BD_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_AD_BC            [RapIndex] = new TH1D(Form("D_ALL_theta_AD_BC_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_A_B_C_D          [RapIndex] = new TH1D(Form("D_ALL_theta_A_B_C_D_%d"  ,      RapIndex), Form("ALL %f<A_y<%f"  ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_ABC_D            [RapIndex] = new TH1D(Form("D_ALL_theta_ABC_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_ABD_C            [RapIndex] = new TH1D(Form("D_ALL_theta_ABD_C_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_ACD_B            [RapIndex] = new TH1D(Form("D_ALL_theta_ACD_B_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_BCD_A            [RapIndex] = new TH1D(Form("D_ALL_theta_BCD_A_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_AB_C_D           [RapIndex] = new TH1D(Form("D_ALL_theta_AB_C_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_AC_B_D           [RapIndex] = new TH1D(Form("D_ALL_theta_AC_B_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_AD_B_C           [RapIndex] = new TH1D(Form("D_ALL_theta_AD_B_C_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_BC_A_D           [RapIndex] = new TH1D(Form("D_ALL_theta_BC_A_D_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_BD_A_C           [RapIndex] = new TH1D(Form("D_ALL_theta_BD_A_C_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_theta_CD_A_B           [RapIndex] = new TH1D(Form("D_ALL_theta_CD_A_B_%d"      ,     RapIndex),Form("ALL %f<A_y<%f"      ,yBin[RapIndex],yBin[RapIndex+1]),SideBinNum,-1,1);
            D_ALL_phi_ABCD               [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_A_B_C_D            [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_AB_CD              [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_AC_BD              [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_AD_BC              [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_ABC_D              [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_ABD_C              [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_ACD_B              [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_BCD_A              [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_AB_C_D             [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_AC_B_D             [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_AD_B_C             [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_BC_A_D             [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_BD_A_C             [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_phi_CD_A_B             [RapIndex]->GetXaxis()->SetTitle("D phi");
            D_ALL_theta_ABCD             [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_AB_CD            [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_AC_BD            [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_AD_BC            [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_A_B_C_D          [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_ABC_D            [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_ABD_C            [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_ACD_B            [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_BCD_A            [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_AB_C_D           [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_AC_B_D           [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_AD_B_C           [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_BC_A_D           [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_BD_A_C           [RapIndex]->GetXaxis()->SetTitle("D theta");
            D_ALL_theta_CD_A_B           [RapIndex]->GetXaxis()->SetTitle("D theta");
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


        Mother_ParID.clear();
        B_List.clear();
        C_List.clear();
        TempEvent.eventID = EntriesID;
        TempEvent.A_particles.clear();
        TempEvent.B_particles.clear();
        TempEvent.C_particles.clear();
        TempEvent.D_particles.clear();
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
            else if ((PDG->at(i) == B_PDG)) {
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
            else if ((PDG->at(i) == C_PDG)) {
                if (fabs(InvariantMass->at(i) - CMass) <= MassSigmaWidth*CMassSigma) {

                    if (IfRemoveHighTPCsigma) {
                        if (abs(C_PDG) == 321) {
                            if (fabs(nSigmaKaon->at(i))>1) continue;
                        }
                    }
                    if (IfRemoveLownHits) {
                        if ((abs(C_PDG) == 321) || (abs(C_PDG) == 211) || (abs(C_PDG) == 2212)) {
                            if (nHitsFit->at(i) < 20) continue;
                        }
                    }
                    if (IfCutHighDCA) {
                        if ((abs(C_PDG) == 321) || (abs(C_PDG) == 211) || (abs(C_PDG) == 2212)) {
                            if ( (0 > dcatopv->at(i)) || (dcatopv->at(i) > 0.5)) continue;
                        }
                    }

                    C = ArmParticle(mix_px->at(i),mix_py->at(i),mix_pz->at(i),CMass,i);
                    C.ParentID.clear();C.ParentID.push_back(i);
                    for (k=ParentSta->at(i);k<=ParentEnd->at(i);k++){
                        C.ParentID.push_back(ParentList->at(k));
                    }
                    if (IfRemoveSpliteMerge) {
                        for (k=SE_ParentSta->at(i);k<=SE_ParentEnd->at(i);k++){
                            C.ParentID.push_back(SE_ParentList->at(k));
                        }
                        for (k=ME_ParentSta->at(i);k<=ME_ParentEnd->at(i);k++){
                            C.ParentID.push_back(ME_ParentList->at(k));
                        }
                    }
                    if ((C.eta < EtaCut[0]) || (C.eta > EtaCut[1])) continue;
                    // TempEvent.B_particles.push_back(C);
                    C_List.push_back(C);
                    continue;
                }
            }
            else if ((PDG->at(i) == D_PDG)) {
                if (fabs(InvariantMass->at(i) - DMass) <= MassSigmaWidth*DMassSigma) {

                    if (IfRemoveHighTPCsigma) {
                        if (abs(D_PDG) == 321) {
                            if (fabs(nSigmaKaon->at(i))>1) continue;
                        }
                    }
                    if (IfRemoveLownHits) {
                        if ((abs(D_PDG) == 321) || (abs(D_PDG) == 211) || (abs(D_PDG) == 2212)) {
                            if (nHitsFit->at(i) < 20) continue;
                        }
                    }
                    if (IfCutHighDCA) {
                        if ((abs(D_PDG) == 321) || (abs(D_PDG) == 211) || (abs(D_PDG) == 2212)) {
                            if ( (0 > dcatopv->at(i)) || (dcatopv->at(i) > 0.5)) continue;
                        }
                    }

                    D = ArmParticle(mix_px->at(i),mix_py->at(i),mix_pz->at(i),DMass,i);
                    D.ParentID.clear();D.ParentID.push_back(i);
                    for (k=ParentSta->at(i);k<=ParentEnd->at(i);k++){
                        D.ParentID.push_back(ParentList->at(k));
                    }
                    if (IfRemoveSpliteMerge) {
                        for (k=SE_ParentSta->at(i);k<=SE_ParentEnd->at(i);k++){
                            D.ParentID.push_back(SE_ParentList->at(k));
                        }
                        for (k=ME_ParentSta->at(i);k<=ME_ParentEnd->at(i);k++){
                            D.ParentID.push_back(ME_ParentList->at(k));
                        }
                    }
                    if ((D.eta < EtaCut[0]) || (D.eta > EtaCut[1])) continue;
                    // TempEvent.B_particles.push_back(C);
                    D_List.push_back(D);
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
                        Mother_ParID.push_back(Temp);
                        // IfFoundOmega = true;
                        // cout<<"Found Omega"<<endl;
                    }
                }
            }
        }
        // 筛选A、B、C、D粒子：优先级：M（母粒子） > A > B > C > D
        // 筛选D粒子
        for (Did=0;Did<D_List.size();Did++) {
            IfRecord = true;
            // 如果D、M有血缘关系，不记录D
            for (Mid = 0;Mid < Mother_ParID.size();Mid++) {
                if (IfInVector(D_List[Did].TreeID , Mother_ParID.at(Mid))) {IfRecord = false;break;}
            }
            // 如果D、C有血缘关系，不记录D
            if (IfRecord) {
                for (Cid = 0;Cid < C_List.size();Cid++) {
                    if (IfInVector(D_List[Did].TreeID , C_List[Cid].ParentID)) {IfRecord = false;break;}
                }
            }
            // 如果D、B有血缘关系，不记录D
            if (IfRecord) {
                for (Bid = 0;Bid < B_List.size();Bid++) {
                    if (IfInVector(D_List[Did].TreeID , B_List[Bid].ParentID)) {IfRecord = false;break;}
                }
            }
            // 如果D、A有血缘关系，不记录D
            if (IfRecord) {
                for (i=0;i<MatchedRap.size();i++) {
                    for (Aid=0;Aid<A_List[MatchedRap[i]].size();Aid++) {
                        if (IfInVector(D_List[Did].TreeID , A_List[MatchedRap[i]][Aid].ParentID)) {IfRecord = false;break;}
                    }
                }
            }
            if (IfRecord) TempEvent.D_particles.push_back(D_List[Did]);
        }
        // 筛选C粒子
        for (Cid=0;Cid<C_List.size();Cid++) {
            IfRecord = true;
            // 如果C、D有血缘关系，不记录C
            for (Did = 0;Did < Mother_ParID.size();Did++) {
                if (IfInVector(C_List[Cid].TreeID , Mother_ParID.at(Did))) {IfRecord = false;break;}
            }
            // 如果C、B有血缘关系，不记录C
            if (IfRecord) {
                for (Bid = 0;Bid < B_List.size();Bid++) {
                    if (IfInVector(C_List[Cid].TreeID , B_List[Bid].ParentID)) {IfRecord = false;break;}
                }
            }
            // 如果C、A有血缘关系，不记录C
            if (IfRecord) {
                for (i=0;i<MatchedRap.size();i++) {
                    for (Aid=0;Aid<A_List[MatchedRap[i]].size();Aid++) {
                        if (IfInVector(C_List[Cid].TreeID , A_List[MatchedRap[i]][Aid].ParentID)) {IfRecord = false;break;}
                    }
                }
            }
            if (IfRecord) TempEvent.C_particles.push_back(C_List[Cid]);
        }
        // 筛选B粒子
        for (Bid=0;Bid<B_List.size();Bid++) {
            IfRecord = true;
            // 如果B、D有血缘关系，不记录B
            for (Did = 0;Did < Mother_ParID.size();Did++) {
                if (IfInVector(B_List[Bid].TreeID , Mother_ParID.at(Did))) {IfRecord = false;break;}
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
                for (Did = 0;Did < Mother_ParID.size();Did++) {
                    if (IfInVector(A_List[MatchedRap[i]][Aid].TreeID , Mother_ParID.at(Did))) {IfRecord = false;break;}
                }
                if (IfRecord) A_Array[MatchedRap[i]].push_back(A_List[MatchedRap[i]][Aid]);
            }
        }

        if (Check_C_D) {
            for(Cid=0;Cid<TempEvent.C_particles.size();Cid++) TempEvent.D_particles.push_back(TempEvent.C_particles[Cid]);
        }

        // if (TempEvent.B_particles.size() >= HowMuchEventMixing+1) continue;
        // 确保同时记录到A、B、...粒子
        if (MatchedRap.size() == 0) continue;                                                        // 有A粒子
        if (TempEvent.B_particles.size() == 0) continue;                                             // 有B粒子
        if (TempEvent.C_particles.size() == 0) continue;                                             // 有C粒子
        if (TempEvent.D_particles.size() == 0) continue;                                             // 有C粒子
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
                        
                                for (int Cid = 0; Cid < HowMuchEventMixing + 1; ++Cid) {
                        
                                    auto& eventC = EventPool[CenIndex][RapIndex][PVzIndex][Cid];
                                    const auto& C_particles = eventC.C_particles;

                                    for (int Did = 0; Did < HowMuchEventMixing + 1; ++Did) {
                                        auto& eventD = EventPool[CenIndex][RapIndex][PVzIndex][Did];
                                        const auto& D_particles = eventD.D_particles;
                                        //==================================================
                                        // Determine mixing type
                                        //==================================================
                            
                                        MixType mixType;

                                        if (Aid == Bid) {
                                        
                                            if (Bid == Cid) {
                                        
                                                if (Cid == Did)
                                                    mixType = SAME;          // ABCD
                                                else
                                                    mixType = ABC_D;         // ABC|D
                                        
                                            } else {
                                        
                                                if (Cid == Did)
                                                    mixType = AB_CD;         // AB|CD
                                                else if (Aid == Did)
                                                    mixType = ABD_C;         // ABD|C
                                                else
                                                    mixType = AB_C_D;        // AB|C|D
                                            }
                                        
                                        }
                                        else if (Aid == Cid) {
                                        
                                            if (Cid == Did)
                                                mixType = ACD_B;             // ACD|B
                                            else if (Bid == Did)
                                                mixType = AC_BD;             // AC|BD
                                            else
                                                mixType = AC_B_D;            // AC|B|D
                                        
                                        }
                                        else if (Aid == Did) {
                                        
                                            if (Bid == Cid)
                                                mixType = AD_BC;             // AD|BC
                                            else
                                                mixType = AD_B_C;            // AD|B|C
                                        
                                        }
                                        else if (Bid == Cid) {
                                        
                                            if (Cid == Did)
                                                mixType = BCD_A;             // BCD|A
                                            else
                                                mixType = BC_A_D;            // BC|A|D
                                        
                                        }
                                        else if (Bid == Did) {
                                        
                                            mixType = BD_A_C;                // BD|A|C
                                        
                                        }
                                        else if (Cid == Did) {
                                        
                                            mixType = CD_A_B;                // CD|A|B
                                        
                                        }
                                        else {
                                        
                                            mixType = A_B_C_D;               // A|B|C|D
                                        
                                        }
                                                        
                                        //==================================================
                                        // Select histogram pointers ONCE
                                        //==================================================
                            
                                        TH1* B_phiLocal = nullptr;
                                        TH1* B_phiGlobal = nullptr;
                                        TH1* B_thetaLocal = nullptr;
                                        TH1* B_thetaGlobal = nullptr;
                                        TH1* C_phiLocal = nullptr;
                                        TH1* C_phiGlobal = nullptr;
                                        TH1* C_thetaLocal = nullptr;
                                        TH1* C_thetaGlobal = nullptr;
                                        TH1* D_phiLocal = nullptr;
                                        TH1* D_phiGlobal = nullptr;
                                        TH1* D_thetaLocal = nullptr;
                                        TH1* D_thetaGlobal = nullptr;
                            
                                        switch (mixType) {
                            
                                            case A_B_C_D:
                                                B_phiLocal      = B_phi_A_B_C_D    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_A_B_C_D  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_A_B_C_D          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_A_B_C_D        [RapIndex];
                                                C_phiLocal      = C_phi_A_B_C_D    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_A_B_C_D  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_A_B_C_D          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_A_B_C_D        [RapIndex];
                                                D_phiLocal      = D_phi_A_B_C_D    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_A_B_C_D  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_A_B_C_D          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_A_B_C_D        [RapIndex];
                                                break;
                            
                                            case AB_CD:
                                                B_phiLocal      = B_phi_AB_CD    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_AB_CD  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_AB_CD          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_AB_CD        [RapIndex];
                                                C_phiLocal      = C_phi_AB_CD    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_AB_CD  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_AB_CD          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_AB_CD        [RapIndex];
                                                D_phiLocal      = D_phi_AB_CD    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_AB_CD  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_AB_CD          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_AB_CD        [RapIndex];
                                                break;
                            
                                            case AC_BD:
                                                B_phiLocal      = B_phi_AC_BD    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_AC_BD  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_AC_BD          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_AC_BD        [RapIndex];
                                                C_phiLocal      = C_phi_AC_BD    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_AC_BD  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_AC_BD          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_AC_BD        [RapIndex];
                                                D_phiLocal      = D_phi_AC_BD    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_AC_BD  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_AC_BD          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_AC_BD        [RapIndex];
                                                break;
                            
                                            case AD_BC:
                                                B_phiLocal      = B_phi_AD_BC    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_AD_BC  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_AD_BC          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_AD_BC        [RapIndex];
                                                C_phiLocal      = C_phi_AD_BC    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_AD_BC  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_AD_BC          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_AD_BC        [RapIndex];
                                                D_phiLocal      = D_phi_AD_BC    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_AD_BC  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_AD_BC          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_AD_BC        [RapIndex];
                                                break;
                            
                                            case ABC_D:
                                                B_phiLocal      = B_phi_ABC_D    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_ABC_D  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_ABC_D          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_ABC_D        [RapIndex];
                                                C_phiLocal      = C_phi_ABC_D    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_ABC_D  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_ABC_D          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_ABC_D        [RapIndex];
                                                D_phiLocal      = D_phi_ABC_D    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_ABC_D  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_ABC_D          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_ABC_D        [RapIndex];
                                                break;
                            
                                            case ABD_C:
                                                B_phiLocal      = B_phi_ABD_C    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_ABD_C  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_ABD_C          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_ABD_C        [RapIndex];
                                                C_phiLocal      = C_phi_ABD_C    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_ABD_C  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_ABD_C          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_ABD_C        [RapIndex];
                                                D_phiLocal      = D_phi_ABD_C    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_ABD_C  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_ABD_C          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_ABD_C        [RapIndex];
                                                break;
                            
                                            case ACD_B:
                                                B_phiLocal      = B_phi_ACD_B    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_ACD_B  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_ACD_B          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_ACD_B        [RapIndex];
                                                C_phiLocal      = C_phi_ACD_B    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_ACD_B  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_ACD_B          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_ACD_B        [RapIndex];
                                                D_phiLocal      = D_phi_ACD_B    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_ACD_B  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_ACD_B          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_ACD_B        [RapIndex];
                                                break;
                            
                                            case BCD_A:
                                                B_phiLocal      = B_phi_BCD_A    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_BCD_A  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_BCD_A          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_BCD_A        [RapIndex];
                                                C_phiLocal      = C_phi_BCD_A    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_BCD_A  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_BCD_A          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_BCD_A        [RapIndex];
                                                D_phiLocal      = D_phi_BCD_A    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_BCD_A  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_BCD_A          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_BCD_A        [RapIndex];
                                                break;

                            
                                            case AB_C_D:
                                                B_phiLocal      = B_phi_AB_C_D    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_AB_C_D  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_AB_C_D          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_AB_C_D        [RapIndex];
                                                C_phiLocal      = C_phi_AB_C_D    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_AB_C_D  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_AB_C_D          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_AB_C_D        [RapIndex];
                                                D_phiLocal      = D_phi_AB_C_D    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_AB_C_D  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_AB_C_D          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_AB_C_D        [RapIndex];
                                                break;

                            
                                            case AC_B_D:
                                                B_phiLocal      = B_phi_AC_B_D    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_AC_B_D  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_AC_B_D          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_AC_B_D        [RapIndex];
                                                C_phiLocal      = C_phi_AC_B_D    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_AC_B_D  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_AC_B_D          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_AC_B_D        [RapIndex];
                                                D_phiLocal      = D_phi_AC_B_D    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_AC_B_D  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_AC_B_D          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_AC_B_D        [RapIndex];
                                                break;
                            
                                            case AD_B_C:
                                                B_phiLocal      = B_phi_AD_B_C    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_AD_B_C  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_AD_B_C          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_AD_B_C        [RapIndex];
                                                C_phiLocal      = C_phi_AD_B_C    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_AD_B_C  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_AD_B_C          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_AD_B_C        [RapIndex];
                                                D_phiLocal      = D_phi_AD_B_C    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_AD_B_C  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_AD_B_C          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_AD_B_C        [RapIndex];
                                                break;
                            
                                            case BC_A_D:
                                                B_phiLocal      = B_phi_BC_A_D    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_BC_A_D  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_BC_A_D          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_BC_A_D        [RapIndex];
                                                C_phiLocal      = C_phi_BC_A_D    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_BC_A_D  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_BC_A_D          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_BC_A_D        [RapIndex];
                                                D_phiLocal      = D_phi_BC_A_D    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_BC_A_D  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_BC_A_D          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_BC_A_D        [RapIndex];
                                                break;
                            
                                            case BD_A_C:
                                                B_phiLocal      = B_phi_BD_A_C    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_BD_A_C  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_BD_A_C          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_BD_A_C        [RapIndex];
                                                C_phiLocal      = C_phi_BD_A_C    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_BD_A_C  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_BD_A_C          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_BD_A_C        [RapIndex];
                                                D_phiLocal      = D_phi_BD_A_C    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_BD_A_C  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_BD_A_C          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_BD_A_C        [RapIndex];
                                                break;
                            
                                            case CD_A_B:
                                                B_phiLocal      = B_phi_CD_A_B    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_CD_A_B  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_CD_A_B          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_CD_A_B        [RapIndex];
                                                C_phiLocal      = C_phi_CD_A_B    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_CD_A_B  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_CD_A_B          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_CD_A_B        [RapIndex];
                                                D_phiLocal      = D_phi_CD_A_B    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_CD_A_B  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_CD_A_B          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_CD_A_B        [RapIndex];
                                                break;

                                            case SAME:
                                                B_phiLocal      = B_phi_ABCD    [CenIndex][RapIndex][PVzIndex];
                                                B_thetaLocal    = B_theta_ABCD  [CenIndex][RapIndex][PVzIndex];
                                                B_phiGlobal     = B_ALL_phi_ABCD          [RapIndex];
                                                B_thetaGlobal   = B_ALL_theta_ABCD        [RapIndex];
                                                C_phiLocal      = C_phi_ABCD    [CenIndex][RapIndex][PVzIndex];
                                                C_thetaLocal    = C_theta_ABCD  [CenIndex][RapIndex][PVzIndex];
                                                C_phiGlobal     = C_ALL_phi_ABCD          [RapIndex];
                                                C_thetaGlobal   = C_ALL_theta_ABCD        [RapIndex];
                                                D_phiLocal      = D_phi_ABCD    [CenIndex][RapIndex][PVzIndex];
                                                D_thetaLocal    = D_theta_ABCD  [CenIndex][RapIndex][PVzIndex];
                                                D_phiGlobal     = D_ALL_phi_ABCD          [RapIndex];
                                                D_thetaGlobal   = D_ALL_theta_ABCD        [RapIndex];
                                                break;
                                        }

                                        for (const auto& A : A_particles) {

                                            for (j=0;j<B_particles.size();j++) {
                                                const auto& B = B_particles[j];
                                
                                                for (k=0;k<C_particles.size();k++) {

                                                    const auto& C = C_particles[k];

                                                    for (l=0;l<D_particles.size();l++) {

                                                        if ((Check_C_D) && (Cid == Did)) {
                                                            if (k == l) continue;
                                                            if (IfInVector(C_particles[k].TreeID , D.ParentID)) continue;
                                                        }
                                                        const auto& D = D_particles[k];
                                                
                                                        if (GetSide(
                                                            A,
                                                            B,
                                                            C,
                                                            D,
                                                            *H_P_tot,
                                                            *H_beta,
                                                            Btheta,
                                                            Ctheta,
                                                            Dtheta,
                                                            Bphi,
                                                            Cphi,
                                                            Dphi,
                                                            IfRemoveFeedPair,
                                                            MotherMass,
                                                            MotherMassSigma,
                                                            MassSigmaWidth))
                                                        {
                                                            B_phiLocal     ->Fill(Bphi);
                                                            B_thetaLocal   ->Fill(Btheta);
                                                            B_phiGlobal    ->Fill(Bphi);
                                                            B_thetaGlobal  ->Fill(Btheta);
                                                            C_phiLocal     ->Fill(Cphi);
                                                            C_thetaLocal   ->Fill(Ctheta);
                                                            C_phiGlobal    ->Fill(Cphi);
                                                            C_thetaGlobal  ->Fill(Ctheta);
                                                            D_phiLocal     ->Fill(Dphi);
                                                            D_thetaLocal   ->Fill(Dtheta);
                                                            D_phiGlobal    ->Fill(Dphi);
                                                            D_thetaGlobal  ->Fill(Dtheta);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if (mixType == SAME) {
                                            ++AccumSameNum;
                                        }
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
    fileA->cd();
    H_P_tot->Write();
    H_beta ->Write();
    for (RapIndex=0;RapIndex<yBinNum;RapIndex++) {
        for (CenIndex=0;CenIndex<CentralityBinNum;CenIndex++) {
            for (PVzIndex=0;PVzIndex<PVzBinNum;PVzIndex++) {
                Sep_Side->cd();
                B_phi_ABCD           [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_AB_CD          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_AC_BD          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_AD_BC          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_A_B_C_D        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_ABC_D          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_ABD_C          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_ACD_B          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_BCD_A          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_AB_C_D         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_AC_B_D         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_AD_B_C         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_BC_A_D         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_BD_A_C         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_phi_CD_A_B         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_ABCD         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_AB_CD        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_AC_BD        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_AD_BC        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_A_B_C_D      [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_ABC_D        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_ABD_C        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_ACD_B        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_BCD_A        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_AB_C_D       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_AC_B_D       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_AD_B_C       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_BC_A_D       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_BD_A_C       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                B_theta_CD_A_B       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                
                C_phi_ABCD           [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_AB_CD          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_AC_BD          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_AD_BC          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_A_B_C_D        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_ABC_D          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_ABD_C          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_ACD_B          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_BCD_A          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_AB_C_D         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_AC_B_D         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_AD_B_C         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_BC_A_D         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_BD_A_C         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_phi_CD_A_B         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_ABCD         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_AB_CD        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_AC_BD        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_AD_BC        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_A_B_C_D      [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_ABC_D        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_ABD_C        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_ACD_B        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_BCD_A        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_AB_C_D       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_AC_B_D       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_AD_B_C       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_BC_A_D       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_BD_A_C       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                C_theta_CD_A_B       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                
                D_phi_ABCD           [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_AB_CD          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_AC_BD          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_AD_BC          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_A_B_C_D        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_ABC_D          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_ABD_C          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_ACD_B          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_BCD_A          [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_AB_C_D         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_AC_B_D         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_AD_B_C         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_BC_A_D         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_BD_A_C         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_phi_CD_A_B         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_ABCD         [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_AB_CD        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_AC_BD        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_AD_BC        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_A_B_C_D      [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_ABC_D        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_ABD_C        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_ACD_B        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_BCD_A        [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_AB_C_D       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_AC_B_D       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_AD_B_C       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_BC_A_D       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_BD_A_C       [CenIndex] [RapIndex] [PVzIndex] ->Write();
                D_theta_CD_A_B       [CenIndex] [RapIndex] [PVzIndex] ->Write();
            }
        }
        ALL_Side->cd();
        B_ALL_phi_ABCD               [RapIndex] ->Write();
        B_ALL_phi_A_B_C_D            [RapIndex] ->Write();
        B_ALL_phi_AB_CD              [RapIndex] ->Write();
        B_ALL_phi_AC_BD              [RapIndex] ->Write();
        B_ALL_phi_AD_BC              [RapIndex] ->Write();
        B_ALL_phi_ABC_D              [RapIndex] ->Write();
        B_ALL_phi_ABD_C              [RapIndex] ->Write();
        B_ALL_phi_ACD_B              [RapIndex] ->Write();
        B_ALL_phi_BCD_A              [RapIndex] ->Write();
        B_ALL_phi_AB_C_D             [RapIndex] ->Write();
        B_ALL_phi_AC_B_D             [RapIndex] ->Write();
        B_ALL_phi_AD_B_C             [RapIndex] ->Write();
        B_ALL_phi_BC_A_D             [RapIndex] ->Write();
        B_ALL_phi_BD_A_C             [RapIndex] ->Write();
        B_ALL_phi_CD_A_B             [RapIndex] ->Write();
        B_ALL_theta_ABCD             [RapIndex] ->Write();
        B_ALL_theta_AB_CD            [RapIndex] ->Write();
        B_ALL_theta_AC_BD            [RapIndex] ->Write();
        B_ALL_theta_AD_BC            [RapIndex] ->Write();
        B_ALL_theta_A_B_C_D          [RapIndex] ->Write();
        B_ALL_theta_ABC_D            [RapIndex] ->Write();
        B_ALL_theta_ABD_C            [RapIndex] ->Write();
        B_ALL_theta_ACD_B            [RapIndex] ->Write();
        B_ALL_theta_BCD_A            [RapIndex] ->Write();
        B_ALL_theta_AB_C_D           [RapIndex] ->Write();
        B_ALL_theta_AC_B_D           [RapIndex] ->Write();
        B_ALL_theta_AD_B_C           [RapIndex] ->Write();
        B_ALL_theta_BC_A_D           [RapIndex] ->Write();
        B_ALL_theta_BD_A_C           [RapIndex] ->Write();
        B_ALL_theta_CD_A_B           [RapIndex] ->Write();
        
        C_ALL_phi_ABCD               [RapIndex] ->Write();
        C_ALL_phi_A_B_C_D            [RapIndex] ->Write();
        C_ALL_phi_AB_CD              [RapIndex] ->Write();
        C_ALL_phi_AC_BD              [RapIndex] ->Write();
        C_ALL_phi_AD_BC              [RapIndex] ->Write();
        C_ALL_phi_ABC_D              [RapIndex] ->Write();
        C_ALL_phi_ABD_C              [RapIndex] ->Write();
        C_ALL_phi_ACD_B              [RapIndex] ->Write();
        C_ALL_phi_BCD_A              [RapIndex] ->Write();
        C_ALL_phi_AB_C_D             [RapIndex] ->Write();
        C_ALL_phi_AC_B_D             [RapIndex] ->Write();
        C_ALL_phi_AD_B_C             [RapIndex] ->Write();
        C_ALL_phi_BC_A_D             [RapIndex] ->Write();
        C_ALL_phi_BD_A_C             [RapIndex] ->Write();
        C_ALL_phi_CD_A_B             [RapIndex] ->Write();
        C_ALL_theta_ABCD             [RapIndex] ->Write();
        C_ALL_theta_AB_CD            [RapIndex] ->Write();
        C_ALL_theta_AC_BD            [RapIndex] ->Write();
        C_ALL_theta_AD_BC            [RapIndex] ->Write();
        C_ALL_theta_A_B_C_D          [RapIndex] ->Write();
        C_ALL_theta_ABC_D            [RapIndex] ->Write();
        C_ALL_theta_ABD_C            [RapIndex] ->Write();
        C_ALL_theta_ACD_B            [RapIndex] ->Write();
        C_ALL_theta_BCD_A            [RapIndex] ->Write();
        C_ALL_theta_AB_C_D           [RapIndex] ->Write();
        C_ALL_theta_AC_B_D           [RapIndex] ->Write();
        C_ALL_theta_AD_B_C           [RapIndex] ->Write();
        C_ALL_theta_BC_A_D           [RapIndex] ->Write();
        C_ALL_theta_BD_A_C           [RapIndex] ->Write();
        C_ALL_theta_CD_A_B           [RapIndex] ->Write();
        
        D_ALL_phi_ABCD               [RapIndex] ->Write();
        D_ALL_phi_A_B_C_D            [RapIndex] ->Write();
        D_ALL_phi_AB_CD              [RapIndex] ->Write();
        D_ALL_phi_AC_BD              [RapIndex] ->Write();
        D_ALL_phi_AD_BC              [RapIndex] ->Write();
        D_ALL_phi_ABC_D              [RapIndex] ->Write();
        D_ALL_phi_ABD_C              [RapIndex] ->Write();
        D_ALL_phi_ACD_B              [RapIndex] ->Write();
        D_ALL_phi_BCD_A              [RapIndex] ->Write();
        D_ALL_phi_AB_C_D             [RapIndex] ->Write();
        D_ALL_phi_AC_B_D             [RapIndex] ->Write();
        D_ALL_phi_AD_B_C             [RapIndex] ->Write();
        D_ALL_phi_BC_A_D             [RapIndex] ->Write();
        D_ALL_phi_BD_A_C             [RapIndex] ->Write();
        D_ALL_phi_CD_A_B             [RapIndex] ->Write();
        D_ALL_theta_ABCD             [RapIndex] ->Write();
        D_ALL_theta_AB_CD            [RapIndex] ->Write();
        D_ALL_theta_AC_BD            [RapIndex] ->Write();
        D_ALL_theta_AD_BC            [RapIndex] ->Write();
        D_ALL_theta_A_B_C_D          [RapIndex] ->Write();
        D_ALL_theta_ABC_D            [RapIndex] ->Write();
        D_ALL_theta_ABD_C            [RapIndex] ->Write();
        D_ALL_theta_ACD_B            [RapIndex] ->Write();
        D_ALL_theta_BCD_A            [RapIndex] ->Write();
        D_ALL_theta_AB_C_D           [RapIndex] ->Write();
        D_ALL_theta_AC_B_D           [RapIndex] ->Write();
        D_ALL_theta_AD_B_C           [RapIndex] ->Write();
        D_ALL_theta_BC_A_D           [RapIndex] ->Write();
        D_ALL_theta_BD_A_C           [RapIndex] ->Write();
        D_ALL_theta_CD_A_B           [RapIndex] ->Write();
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
                  << " A_PDG B_PDG C_PDG D_PDG Mode SP_ME [CutID]" << std::endl;
        return 1;
    }

    S_Four(
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
    const ArmParticle& C,
    const ArmParticle& D,
    TH1D& H_P_tot,
    TH1D& H_beta,
    double& BthetaOut,
    double& CthetaOut,
    double& DthetaOut,
    double& BphiOut,
    double& CphiOut,
    double& DphiOut,
    bool IfRemoveFeedPair,
    const std::vector<float>& MotherMass,
    const std::vector<float>& MotherMassSigma,
    float MassSigmaWidth)
{
    const double AE = A.E;
    const double BE = B.E;
    const double CE = C.E;
    const double DE = D.E;
    const double TotE = AE + BE + CE + DE;
    double p[4] = {A.px+B.px+C.px+D.px , A.py+B.py+C.py+D.py , A.pz+B.pz+C.pz+D.pz , 0.0};
    p[3] = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    H_P_tot.Fill(p[3]);
    const double n[3] = {p[0]/p[3] , p[1]/p[3] , p[2]/p[3]};
    double beta[4] = { -(p[0])/TotE , -(p[1])/TotE , -(p[2])/TotE , 0.0};
    beta[3] = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    H_beta.Fill(sqrt(beta[3]));

    const double gamma  = 1.0/(sqrt(1-beta[3]));
    const double gamma2 = 1.0/(sqrt(1-beta[3])*(1+sqrt(1-beta[3])));

    //////////////////////////////////
    // Three body figure
    //////////////////////////////////

    const double bpB = beta[0]*B.px + beta[1]*B.py + beta[2]*B.pz;
    const double bpC = beta[0]*C.px + beta[1]*C.py + beta[2]*C.pz;
    const double bpD = beta[0]*D.px + beta[1]*D.py + beta[2]*D.pz;

    const double New_BPx = B.px + gamma2*beta[0]*bpB + gamma*beta[0]*BE;
    const double New_BPy = B.py + gamma2*beta[1]*bpB + gamma*beta[1]*BE;
    const double New_BPz = B.pz + gamma2*beta[2]*bpB + gamma*beta[2]*BE;
    const double New_CPx = C.px + gamma2*beta[0]*bpC + gamma*beta[0]*CE;
    const double New_CPy = C.py + gamma2*beta[1]*bpC + gamma*beta[1]*CE;
    const double New_CPz = C.pz + gamma2*beta[2]*bpC + gamma*beta[2]*CE;
    const double New_DPx = D.px + gamma2*beta[0]*bpD + gamma*beta[0]*DE;
    const double New_DPy = D.py + gamma2*beta[1]*bpD + gamma*beta[1]*DE;
    const double New_DPz = D.pz + gamma2*beta[2]*bpD + gamma*beta[2]*DE;

    double cosPhiOut = (New_BPx*(n[0])+New_BPy*(n[1])+New_BPz*(n[2])) / (sqrt(New_BPx*New_BPx+New_BPy*New_BPy+New_BPz*New_BPz));
    BphiOut = std::acos(cosPhiOut);
    cosPhiOut = (New_CPx*(n[0])+New_CPy*(n[1])+New_CPz*(n[2])) / (sqrt(New_CPx*New_CPx+New_CPy*New_CPy+New_CPz*New_CPz));
    CphiOut = std::acos(cosPhiOut);
    cosPhiOut = (New_DPx*(n[0])+New_DPy*(n[1])+New_DPz*(n[2])) / (sqrt(New_DPx*New_DPx+New_DPy*New_DPy+New_DPz*New_DPz));
    DphiOut = std::acos(cosPhiOut);

    double v[3] = {-n[1],n[0],0};
    const double Rv = sqrt(n[0]*n[0] + n[1]*n[1]);
    v[0] = v[0]/Rv;v[1] = v[1]/Rv;

    const Vec3 NewB = {New_BPx,New_BPy,New_BPz};
    const Vec3 NewC = {New_CPx,New_CPy,New_CPz};
    const Vec3 NewD = {New_DPx,New_DPy,New_DPz};
    const Vec3 N    = {n[0]   ,n[1]   ,n[2]   };

    if(computeRotatedProjections(NewB, NewC, NewD, N, BthetaOut, CthetaOut, DthetaOut)){
        return true;
    }else{
        return false;
    }

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


double dot(const Vec3& a, const Vec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return { a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] };
}

Vec3 normalize(const Vec3& v) {
    double len = std::sqrt(dot(v, v));
    if (len < 1e-12) return {0,0,0};
    return { v[0]/len, v[1]/len, v[2]/len };
}

// 返回 ThetaB, ThetaC, ThetaD (弧度)
bool computeRotatedProjections(const Vec3& b, const Vec3& c, const Vec3& d, const Vec3& n, double &ThetaB, double &ThetaC, double &ThetaD) {
    // 1. 投影到平面 P
    auto project = [&n](const Vec3& v) -> Vec3 {
        double proj_n = dot(v, n);
        return { v[0] - proj_n * n[0],
                 v[1] - proj_n * n[1],
                 v[2] - proj_n * n[2] };
    };
    Vec3 proj_b = project(b);
    Vec3 proj_c = project(c);
    Vec3 proj_d = project(d);

    // 2. 在平面 P 上建立正交基 (u, v)
    Vec3 a = {1, 0, 0};
    if (std::fabs(dot(n, a)) > 0.9999) {
        a = {0, 1, 0};
    }
    Vec3 u = normalize(cross(n, a));
    Vec3 v = cross(n, u);   // 已在平面内，且与 u 正交

    // 3. 将投影向量转换到二维坐标 (x, y)
    auto to2D = [&](const Vec3& p) -> Vec2 {
        return { dot(p, u), dot(p, v) };
    };
    Vec2 vb0 = to2D(proj_b);
    Vec2 vc0 = to2D(proj_c);
    Vec2 vd0 = to2D(proj_d);

    // 4. 三个固定方向单位向量
    const double sqrt3 = std::sqrt(3.0);
    Vec2 uB = {0.0, 1.0};
    Vec2 uC = {-sqrt3/2.0, -0.5};
    Vec2 uD = { sqrt3/2.0, -0.5};

    // 5. 计算每个点的 Ai, Bi
    auto computeAB = [](const Vec2& v0, const Vec2& u) -> std::pair<double, double> {
        double A = u[0]*v0[0] + u[1]*v0[1];
        double B = -u[0]*v0[1] + u[1]*v0[0];
        return {A, B};
    };
    std::pair<double,double> tmpB = computeAB(vb0, uB);
    double Ab = tmpB.first;
    double Bb = tmpB.second;
    
    std::pair<double,double> tmpC = computeAB(vc0, uC);
    double Ac = tmpC.first;
    double Bc = tmpC.second;
    
    std::pair<double,double> tmpD = computeAB(vd0, uD);
    double Ad = tmpD.first;
    double Bd = tmpD.second;

    // 6. 计算总和
    double sum_A2_minus_B2 = (Ab*Ab - Bb*Bb) + (Ac*Ac - Bc*Bc) + (Ad*Ad - Bd*Bd);
    double sum_AB = Ab*Bb + Ac*Bc + Ad*Bd;

    // 7. 最优旋转角度 θ
    double theta = 0.0;
    double C1 = 0.5 * sum_A2_minus_B2;
    double C2 = sum_AB;
    if (std::fabs(C1) > 1e-12 || std::fabs(C2) > 1e-12) {
        theta = 0.5 * std::atan2(C2, C1);
    }

    // 8. 旋转所有点
    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);
    auto rotate = [cos_t, sin_t](const Vec2& v) -> Vec2 {
        return { v[0]*cos_t - v[1]*sin_t,
                 v[0]*sin_t + v[1]*cos_t };
    };
    Vec2 vb_rot = rotate(vb0);
    Vec2 vc_rot = rotate(vc0);
    Vec2 vd_rot = rotate(vd0);

    // 9. 计算每个向量相对于 uB (0,1) 的夹角
    auto angleFromUB = [](const Vec2& vec) -> double {
        // uB 方向为 (0,1)，夹角 = atan2(x, y)
        return std::atan2(vec[0], vec[1]);
    };
    ThetaB = angleFromUB(vb_rot);
    ThetaC = angleFromUB(vc_rot);
    ThetaD = angleFromUB(vd_rot);

    return true;
}