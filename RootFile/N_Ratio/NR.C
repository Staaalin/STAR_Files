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


#define A_Num_Per_Event 5
#define B_Num_Per_Event 5
#define HowMuchEventMixing 10

struct Particle{
    unsigned short int TreeID;
    float Px;
    float Py;
    float Pz;

    float Rap;
    float Pt;
};

struct ParticlePool{
    unsigned int EvtID;
    bool IfMadeSame = false;
    unsigned short int ListA_Index;
    Particle ListA[A_Num_Per_Event];
    unsigned short int ListB_Index;
    Particle ListB[B_Num_Per_Event];
};

void print(std::vector<int> Temp);
void print(std::vector<float> Temp);
float CenCorr(float Vz,TString Name);
Double_t massList(int PID,TString Name);
Double_t massListSigma(int PID,TString Name);
bool IfInVector(int Num , std::vector<int> V);
bool IfInVector(int Num , std::vector<unsigned short int> V);
bool IfCommonElement(std::vector<int> A , std::vector<int> B);
void DltElement(std::vector<int> &V , int ID);
std::vector<int> GetDaughterPDGLit(int ID);
std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize);
float GetPairMass(float p1x,float p1y,float p1z,float m1,float p2x,float p2y,float p2z,float m2);
void GetPairMassAndKstar(Particle PA , Particle PB , float AMass , float BMass , float (&MassAndKstar)[2]);

// const int CentralityBin[] = {0 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 60 , 70 , 80};// %
const int CentralityBin[] = {0 , 10 , 20 , 40 , 60 , 100};// %
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

    int i , j , k , l , m , n;// used as Index
    int Aid , Bid , Cid;// used as Index
    uint8_t GenI , GenJ , GenK; // used as index generator looping in pool when filling hist
    float tPx , tPy , tPz , tPtSqu , tPt , tRap , tEnergy , PairMass , MassAndKstar[2];
    std::vector<int> Temp;
    std::vector<float> CMass , CMassSigma;
    float C_Mass;
    Particle AnyParticle;
    Particle ParticleA[A_Num_Per_Event] , ParticleB[B_Num_Per_Event];
    Particle ParticleA_Tmp , ParticleB_Tmp;
    float Mass_Store[A_Num_Per_Event*B_Num_Per_Event];
    float Kstar_Store[A_Num_Per_Event*B_Num_Per_Event];
    int   A_Num_Store[A_Num_Per_Event*B_Num_Per_Event];
    int   B_Num_Store[A_Num_Per_Event*B_Num_Per_Event];
    unsigned short int    ParticleASize , ParticleBSize , ParticleCSize; // Particle*Size == Particle*.Size
    unsigned short int    ParticleASizeR, ParticleBSizeR; // Particle*Size after cut
    std::vector<std::vector<unsigned short int> > A_ParID,B_ParID,C_ParID;
    bool A_IfRecord[A_Num_Per_Event],B_IfRecord[B_Num_Per_Event];
    uint8_t A_yIndex[A_Num_Per_Event],B_yIndex[B_Num_Per_Event];
    bool IfMatched[yBinNum];
    bool IfRecord = true , IfRemoveFeedPair = false;
    float BMass = massList(B_PDG)           , AMass = massList(A_PDG);
    float BMassSigma = massListSigma(B_PDG) , AMassSigma = massListSigma(A_PDG);

    for (int i = 0;i < FeedDownNum;i++){
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

    for (int i=0;i<FeedDownNum;i++) {
        if ( IfInVector(A_PDG , GetDaughterPDGLit(FeedDown[i])) && IfInVector(B_PDG , GetDaughterPDGLit(FeedDown[i])) ) IfRemoveFeedPair = true;
    }

    std::vector<int> NchList = GetNchList(CentralityBin , CentralityBinNum+1);     // centrality
    cout<<"NchList = ";
    print(NchList);
    cout<<" "<<endl;

    // ############################################################################################################# //
    // ####                                            Declare Histgram                                         #### //
    // ############################################################################################################# //
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
    uint8_t CenIndex , PVzIndex , yIndex;
    bool RapIndex[15];
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
    int RebinNum[] = {1,2,4,5,10,20,25,50};                    //               RebinNum
    TH1D* H_A_Num_dRap                    [15]                 [15]       [15]    [15];
    TH1D* H_Mix_A_Num_dRap                [15]                 [15]       [15]    [15];
    TH1D* H_Res_A_Num_dRap                [15]                 [15]       [15]    [15];
    TH1D* H_ALL_A_Num_dRap                                     [15]               [15];
    TH1D* H_ALL_Mix_A_Num_dRap                                 [15]               [15];
    TH1D* H_ALL_Res_A_Num_dRap                                 [15]               [15];
    TH1D* H_B_Num_dRap                    [15]                 [15]       [15]    [15];
    TH1D* H_Mix_B_Num_dRap                [15]                 [15]       [15]    [15];
    TH1D* H_Res_B_Num_dRap                [15]                 [15]       [15]    [15];
    TH1D* H_ALL_B_Num_dRap                                     [15]               [15];
    TH1D* H_ALL_Mix_B_Num_dRap                                 [15]               [15];
    TH1D* H_ALL_Res_B_Num_dRap                                 [15]               [15];
    //                                                                          EventPool
    ParticlePool       Tot_Pool           [15]                 [15]       [15]  [HowMuchEventMixing+1];
    uint8_t            Tot_Pool_F_Index   [15]                 [15]       [15] ;                     // The point in pool
    uint8_t            Tot_Pool_Num       [15]                 [15]       [15] ;
    bool               Tot_Pool_IfFilled  [15]                 [15]       [15] ;
    ParticlePool       Tot_Pool_Tmp                            [15]       ; // Temp Store
    ParticlePool       Tot_Pool_TTmp                                      ; // Temp Store
    ParticlePool       Tot_Pool_TTTmp                                     ; // Temp Store

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
    for (i=0;i<CentralityBinNum;i++){
        for (k=0;k<PVzBinNum;k++){
            for (j=0;j<yBinNum;j++){
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
                H_dRap     [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum,dRapSta,dRapEnd);

                HistNameTemp1 = "H_M_";HistNameTemp1+="dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Mixed dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Mix_dRap [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum,dRapSta,dRapEnd);

                HistNameTemp1 = "H_R_";HistNameTemp1+="dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Resed dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Res_dRap [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum,dRapSta,dRapEnd);
                
                HistNameTemp1 = "H_";HistNameTemp1+="dPt_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "dPt, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_dPt      [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,dPtBinNum,dPtSta,dPtEnd);

                HistNameTemp1 = "H_M_";HistNameTemp1+="dPt_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Mixed dPt, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Mix_dPt  [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,dPtBinNum,dPtSta,dPtEnd);

                HistNameTemp1 = "H_R_";HistNameTemp1+="dPt_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Resed dPt, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Res_dPt  [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,dPtBinNum,dPtSta,dPtEnd);
                
                HistNameTemp1 = "H_";HistNameTemp1+="Mass_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Mass, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Mass     [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,MBinNum,MSta,MEnd);

                HistNameTemp1 = "H_M_";HistNameTemp1+="Mass_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Mixed Mass, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Mix_Mass [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,MBinNum,MSta,MEnd);

                HistNameTemp1 = "H_R_";HistNameTemp1+="Mass_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Resed Mass, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Res_Mass [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,MBinNum,MSta,MEnd);
                
                HistNameTemp1 = "H_";HistNameTemp1+="A_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "A_Num, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_A_Num     [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                HistNameTemp1 = "H_M_";HistNameTemp1+="A_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Mixed A_Num, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Mix_A_Num [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                HistNameTemp1 = "H_R_";HistNameTemp1+="A_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Resed A_Num, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Res_A_Num [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);
                
                HistNameTemp1 = "H_";HistNameTemp1+="B_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "B_Num, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_B_Num     [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                HistNameTemp1 = "H_M_";HistNameTemp1+="B_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Mixed B_Num, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Mix_B_Num [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                HistNameTemp1 = "H_R_";HistNameTemp1+="B_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Resed B_Num, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Res_B_Num [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);
                
                HistNameTemp1 = "H_";HistNameTemp1+="Event_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Event_Num, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Event_Num     [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                HistNameTemp1 = "H_R_";HistNameTemp1+="Event_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;
                HistNameTemp2 = "Resed Event_Num, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                H_Res_Event_Num [i][j][k] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);
                
                for (int RebinIndex=0;RebinIndex<(sizeof(RebinNum)/sizeof(RebinNum[0]));RebinIndex++){
                    HistNameTemp1 = "H_";HistNameTemp1+="A_Num_dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                    HistNameTemp2 = "A_Num_dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                    H_A_Num_dRap    [i][j][k][RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);

                    HistNameTemp1 = "H_M_";HistNameTemp1+="A_Num_dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                    HistNameTemp2 = "Mixed A_Num_dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                    H_Mix_A_Num_dRap[i][j][k][RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);

                    HistNameTemp1 = "H_R_";HistNameTemp1+="A_Num_dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                    HistNameTemp2 = "Resed A_Num_dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                    H_Res_A_Num_dRap[i][j][k][RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);
                    
                    HistNameTemp1 = "H_";HistNameTemp1+="B_Num_dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                    HistNameTemp2 = "B_Num_dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                    H_B_Num_dRap    [i][j][k][RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);

                    HistNameTemp1 = "H_M_";HistNameTemp1+="B_Num_dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                    HistNameTemp2 = "Mixed B_Num_dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                    H_Mix_B_Num_dRap[i][j][k][RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);

                    HistNameTemp1 = "H_R_";HistNameTemp1+="B_Num_dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+=k;HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                    HistNameTemp2 = "Resed B_Num_dRap, [";HistNameTemp2+=CentralityBin[i];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[i+1];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[k];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[k+1];HistNameTemp2+=" cm";
                    H_Res_B_Num_dRap[i][j][k][RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);
                }

                if ((i == 0)&&(k == 0)) {
                    HistNameTemp1 = "H_";HistNameTemp1+="Kstar_";HistNameTemp1+="ALL";HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Kstar, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Kstar    [j] = new TH1D(HistNameTemp1,HistNameTemp2,kStarBinNum,kStarSta,kStarEnd);

                    HistNameTemp1 = "H_M_";HistNameTemp1+="Kstar_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Mixed Kstar, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Mix_Kstar[j] = new TH1D(HistNameTemp1,HistNameTemp2,kStarBinNum,kStarSta,kStarEnd);

                    HistNameTemp1 = "H_R_";HistNameTemp1+="Kstar_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Resed Kstar, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Res_Kstar[j] = new TH1D(HistNameTemp1,HistNameTemp2,kStarBinNum,kStarSta,kStarEnd);
                    
                    HistNameTemp1 = "H_";HistNameTemp1+="dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "dRap, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_dRap     [j] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum,dRapSta,dRapEnd);

                    HistNameTemp1 = "H_M_";HistNameTemp1+="dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Mixed dRap, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Mix_dRap [j] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum,dRapSta,dRapEnd);

                    HistNameTemp1 = "H_R_";HistNameTemp1+="dRap_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Resed dRap, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Res_dRap [j] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum,dRapSta,dRapEnd);
                    
                    HistNameTemp1 = "H_";HistNameTemp1+="dPt_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "dPt, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_dPt      [j] = new TH1D(HistNameTemp1,HistNameTemp2,dPtBinNum,dPtSta,dPtEnd);

                    HistNameTemp1 = "H_M_";HistNameTemp1+="dPt_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Mixed dPt, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Mix_dPt  [j] = new TH1D(HistNameTemp1,HistNameTemp2,dPtBinNum,dPtSta,dPtEnd);

                    HistNameTemp1 = "H_R_";HistNameTemp1+="dPt_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Resed dPt, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Res_dPt  [j] = new TH1D(HistNameTemp1,HistNameTemp2,dPtBinNum,dPtSta,dPtEnd);
                    
                    HistNameTemp1 = "H_";HistNameTemp1+="Mass_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Mass, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Mass     [j] = new TH1D(HistNameTemp1,HistNameTemp2,MBinNum,MSta,MEnd);

                    HistNameTemp1 = "H_M_";HistNameTemp1+="Mass_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Mixed Mass, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Mix_Mass [j] = new TH1D(HistNameTemp1,HistNameTemp2,MBinNum,MSta,MEnd);

                    HistNameTemp1 = "H_R_";HistNameTemp1+="Mass_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Resed Mass, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Res_Mass [j] = new TH1D(HistNameTemp1,HistNameTemp2,MBinNum,MSta,MEnd);
                    
                    HistNameTemp1 = "H_";HistNameTemp1+="A_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "A_Num, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_A_Num     [j] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                    HistNameTemp1 = "H_M_";HistNameTemp1+="A_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Mixed A_Num, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Mix_A_Num [j] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                    HistNameTemp1 = "H_R_";HistNameTemp1+="A_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Resed A_Num, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Res_A_Num [j] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);
                    
                    HistNameTemp1 = "H_";HistNameTemp1+="B_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "B_Num, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_B_Num     [j] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                    HistNameTemp1 = "H_M_";HistNameTemp1+="B_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Mixed B_Num, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Mix_B_Num [j] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                    HistNameTemp1 = "H_R_";HistNameTemp1+="B_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Resed B_Num, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Res_B_Num [j] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);
                    
                    HistNameTemp1 = "H_";HistNameTemp1+="Event_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Event_Num, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Event_Num     [j] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);

                    HistNameTemp1 = "H_R_";HistNameTemp1+="Event_Num_";HistNameTemp1+=i;HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";
                    HistNameTemp2 = "Resed Event_Num, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                    HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                    HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                    H_ALL_Res_Event_Num [j] = new TH1D(HistNameTemp1,HistNameTemp2,1,-1,1);
                
                    for (int RebinIndex=0;RebinIndex<(sizeof(RebinNum)/sizeof(RebinNum[0]));RebinIndex++){
                        HistNameTemp1 = "H_";HistNameTemp1+="A_Num_dRap_";HistNameTemp1+="ALL";HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                        HistNameTemp2 = "A_Num_dRap, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                        HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                        HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                        H_A_Num_dRap       [j]   [RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);

                        HistNameTemp1 = "H_M_";HistNameTemp1+="A_Num_dRap_";HistNameTemp1+="ALL";HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                        HistNameTemp2 = "Mixed A_Num_dRap, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                        HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                        HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                        H_Mix_A_Num_dRap   [j]   [RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);

                        HistNameTemp1 = "H_R_";HistNameTemp1+="A_Num_dRap_";HistNameTemp1+="ALL";HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                        HistNameTemp2 = "Resed A_Num_dRap, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                        HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                        HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                        H_Res_A_Num_dRap   [j]   [RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);
                        
                        HistNameTemp1 = "H_";HistNameTemp1+="B_Num_dRap_";HistNameTemp1+="ALL";HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                        HistNameTemp2 = "B_Num_dRap, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                        HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                        HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                        H_B_Num_dRap       [j]   [RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);

                        HistNameTemp1 = "H_M_";HistNameTemp1+="B_Num_dRap_";HistNameTemp1+="ALL";HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                        HistNameTemp2 = "Mixed B_Num_dRap, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                        HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                        HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                        H_Mix_B_Num_dRap   [j]   [RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);

                        HistNameTemp1 = "H_R_";HistNameTemp1+="B_Num_dRap_";HistNameTemp1+="ALL";HistNameTemp1+="_";HistNameTemp1+=j;HistNameTemp1+="_";HistNameTemp1+="ALL";HistNameTemp1+="_R";HistNameTemp1+=RebinNum[RebinIndex];
                        HistNameTemp2 = "Resed B_Num_dRap, [";HistNameTemp2+=CentralityBin[0];HistNameTemp2+="%,";HistNameTemp2+=CentralityBin[CentralityBinNum];HistNameTemp2+="%], ";
                        HistNameTemp2+=yBin[j];HistNameTemp2+="<y_";HistNameTemp2+=B_PDG;HistNameTemp2+="<";HistNameTemp2+=yBin[j+1];HistNameTemp2+=", ";
                        HistNameTemp2+=PVzBin[0];HistNameTemp2+="<PVz<";HistNameTemp2+=PVzBin[PVzBinNum];HistNameTemp2+=" cm";
                        H_Res_B_Num_dRap   [j]   [RebinIndex] = new TH1D(HistNameTemp1,HistNameTemp2,dRapBinNum/RebinNum[RebinIndex],dRapSta,dRapEnd);
                    }
                }
            }
        }
    }
    
    // ############################################################################################################# //
    // ####                                              Events Loop                                            #### //
    // ############################################################################################################# //
    TString TreeName = "hadronTree";
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
    // hadronTree->SetBranchAddress("TriggerID",&TriggerID);
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
    cout << "Events number: " << nentries << endl;
    A_ParID.clear();B_ParID.clear();C_ParID.clear();
    for (i=0;i<PVzBinNum;i++) {
        for (j=0;j<yBinNum;j++) {
            for (k=0;k<PVzBinNum;k++) {
                Tot_Pool_Num       [i]                 [j]       [k] = 0;
                Tot_Pool_F_Index   [i]                 [j]       [k] = -1;
                Tot_Pool_C_Index   [i]                 [j]       [k] = -1;
                Tot_Pool_IfFilled  [i]                 [j]       [k] = false;
            }
        }
    }
    for (int EntriesID = 0 ; EntriesID < nentries ; EntriesID++) {
        hadronTree->GetEntry(EntriesID);

        // Decide Event Index
        CenIndex = -1;
        for (k=0;k<CentralityBinNum;k++){
            NNch = CenCorr(PVz) * Nch;
            // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
            if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                CenIndex = k;
                break;
            }
        }
        if (CenIndex == -1) continue;
        
        PVzIndex = -1;
        for (k=0;k<PVzBinNum;k++){
            // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
            if ((PVzBin[k] <= PVz) && (PVz < PVzBin[k+1])) {
                PVzIndex = k;
                break;
            }
        }
        if (PVzIndex == -1) continue;

        // initialize RapIndex
        for (k=0;k<yBinNum;k++) {
            RapIndex[k] = false;
        }

        ParticleASize = 0;ParticleBSize = 0;ParticleCSize = 0;
        for (j=0;j<PDGMult;j++) {
            if (PDG->at(j) == A_PDG) {
                if (fabs(InvariantMass->at(j) - AMass) <= 3*AMassSigma) {
                    AnyParticle.TreeID = j;
                    ParticleA[ParticleASize] = AnyParticle;
                    ParticleASize++;
                }
                else{continue;}
            }
            else if (PDG->at(j) == B_PDG) {
                if (fabs(InvariantMass->at(j) - BMass) <= 3*AMassSigma) {
                    AnyParticle.TreeID = j;
                    ParticleB[ParticleBSize] = AnyParticle;
                    ParticleBSize++;
                }
                else{continue;}
            }
            else {
                for (int l = 0;l < FeedDownNum;l++) {
                    if ( abs(PDG->at(j)) == FeedDown[l] ) {
                        if ((fabs(InvariantMass->at(j) - CMass.at(l)) > 3*CMassSigma.at(l))) continue;
                        Temp.clear();Temp.push_back(j);
                        for (k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                            Temp.push_back(ParentList->at(k));
                        }
                        C_ParID.push_back(Temp);
                        ParticleCSize++;
                        // IfFoundOmega = true;
                        // cout<<"Found Omega"<<endl;
                    }
                }
            }
        }

        if ((ParticleASize * ParticleBSize) == 0) continue; // if the particle A and B are not found.

        for (i=0;i<ParticleASize;i++) {
            j = ParticleA[i].TreeID;
            Temp.clear();Temp.push_back(j);
            for (k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                Temp.push_back(ParentList->at(k));
            }
            A_ParID.push_back(Temp);
            A_IfRecord[i] = true;
            tPx = mix_px->at(j);
            tPy = mix_py->at(j);
            tPz = mix_pz->at(j);
            ParticleA[i].Px = tPx;
            ParticleA[i].Py = tPy;
            ParticleA[i].Pz = tPz;
            tPtSqu = tPx*tPx + tPy*tPy;
            tEnergy = pow(tPtSqu + tPz*tPz + AMass*AMass,0.5);
            ParticleA[i].Pt = pow(tPtSqu,0.5);
            ParticleA[i].Rap = 0.5*log((tEnergy+tPz)/(tEnergy-tPz));
        }
        for (i=0;i<ParticleBSize;i++) {
            j = ParticleB[i].TreeID;
            Temp.clear();Temp.push_back(j);
            for (k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                Temp.push_back(ParentList->at(k));
            }
            B_ParID.push_back(Temp);
            B_IfRecord[i] = true;
            tPx = mix_px->at(j);
            tPy = mix_py->at(j);
            tPz = mix_pz->at(j);
            ParticleB[i].Px = tPx;
            ParticleB[i].Py = tPy;
            ParticleB[i].Pz = tPz;
            tPtSqu = tPx*tPx + tPy*tPy;
            tEnergy = pow(tPtSqu + tPz*tPz + BMass*BMass,0.5);
            ParticleB[i].Pt = pow(tPtSqu,0.5);
            rap = 0.5*log((tEnergy+tPz)/(tEnergy-tPz));
            ParticleB[i].Rap = rap;
            // Decide B-Rapidity Index
            B_yIndex[i] = -1;
            for (k=0;k<yBinNum;k++){
                if ((yBin[k] <= rap) && (rap < yBin[k+1])) {
                    B_yIndex[i] = k;
                    break;
                }
            }
            if (B_yIndex[i] = -1) B_IfRecord[i] = false;
        }

        // 如果A、B有血缘关系，保留B
        for (Bid = 0;Bid < ParticleBSize;Bid++) {
            for (Aid = 0;Aid < ParticleASize;;Aid++) {
                if (IfInVector(ParticleA[Aid].TreeID , B_ParID.at(Bid))){
                    A_IfRecord[Aid] = false;
                }
                // if (IfCommonElement(A_ParID.at(Aid) , B_ParID.at(Bid))){
                //     A_IfRecord.at(Aid) = 0;
                //     // cout<<"Meet 2!"<<endl;
                //     // cout<<"{ "<<A_PDG<<" } "<<A_TreID.at(Aid)<<" th ";print(A_ParID.at(Aid));
                //     // cout<<"{ "<<B_PDG<<" } "<<B_TreID.at(Bid)<<" th ";print(B_ParID.at(Bid));
                // }
            }
        }
        
        // 如果A、B与C有血缘关系，不记录A和B
        for (Aid = 0;Aid < ParticleASize;Aid++) {
            for (Cid = 0;Cid < ParticleCSize;Cid++) {
                if (IfInVector(ParticleA[Aid].TreeID , C_ParID.at(Cid))) {
                    A_IfRecord[Aid] = false;
                    // cout<<"Meet 3!"<<endl;
                    // cout<<"{ "<<A_PDG<<" } "<<A_TreID.at(Aid)<<" th ";print(A_ParID.at(Aid));
                    // cout<<"{ "<<FeedDown[0]<<" } "<<(C_ParID.at(Cid)).at(0)<<" th ";print(C_ParID.at(Cid));
                }
                // if (IfCommonElement(A_ParID.at(Aid) , C_ParID.at(Cid))){
                //     A_IfRecord.at(Aid) = 0;
                //     // cout<<"Meet 4!"<<endl;
                //     // cout<<"{ "<<A_PDG<<" } "<<A_TreID.at(Aid)<<" th ";print(A_ParID.at(Aid));
                //     // cout<<"{ "<<FeedDown[0]<<" } "<<(C_ParID.at(Cid)).at(0)<<" th ";print(C_ParID.at(Cid));
                // }
            }
        }
        for (Bid = 0;Bid < ParticleBSize;Bid++) {
            for (Cid = 0;Cid < ParticleCSize;Cid++) {
                if (IfInVector(ParticleB[Bid].TreeID , C_ParID.at(Cid))) {
                    B_IfRecord[Bid] = false;
                    // cout<<"Meet 5!"<<endl;
                    // cout<<"{ "<<B_PDG<<" } "<<B_TreID.at(Bid)<<" th ";print(B_ParID.at(Bid));
                    // cout<<"{ "<<FeedDown[0]<<" } "<<(C_ParID.at(Cid)).at(0)<<" th ";print(C_ParID.at(Cid));
                }
                // if (IfCommonElement(B_ParID.at(Bid) , C_ParID.at(Cid))){
                //     B_IfRecord.at(Bid) = 0;
                //     // cout<<"Meet 6!"<<endl;
                //     // cout<<"{ "<<B_PDG<<" } "<<B_TreID.at(Bid)<<" th ";print(B_ParID.at(Bid));
                //     // cout<<"{ "<<FeedDown[0]<<" } "<<(C_ParID.at(Cid)).at(0)<<" th ";print(C_ParID.at(Cid));
                // }
            }
        }
        
        // A rapidity cut
        for (Aid = 0;Aid < ParticleASize;Aid++) {
            tRap = ParticleA[Aid].Rap;
            if ((tRap < AyCut[0]) || (tRap > AyCut[1])){
                A_IfRecord[Aid] = false;
            }
        }

        // A Eta Cut
        for (Aid = 0;Aid < ParticleASize;Aid++) {
            tPt = ParticleA[Aid].Pt;
            tPz = ParticleA[Aid].Pz;
            Eta = -1.0*log(tan(0.5*(acos(tPz/pow(tPt*tPt+tPz*tPz,0.5)))));
            if ((Eta < EtaCut[0]) || (Eta > EtaCut[1])){
                A_IfRecord[Aid] = false;
            }
        }

        // B Eta Cut
        for (Bid = 0;Bid < ParticleBSize;Bid++) {
            tPt = ParticleB[Bid].Pt;
            tPz = ParticleB[Bid].Pz;
            Eta = -1.0*log(tan(0.5*(acos(tPz/pow(tPt*tPt+tPz*tPz,0.5)))));
            if ((Eta < EtaCut[0]) || (Eta > EtaCut[1])){
                B_IfRecord[Bid] = false;
            }
        }

        // C FeedDown Cut
        if (IfRemoveFeedPair) {
            for (Cid = 0;Cid<ParticleCSize;Cid++){
                for (Bid = 0;Bid < ParticleBSize;Bid++) {
                    if (B_IfRecord[Bid]) {
                        C_Mass = CMass[Cid];
                        for (Aid = 0;Aid < ParticleASize;Aid++) {
                            if (A_IfRecord[Aid]) {
                                PairMass = GetPairMass(ParticleA[Aid].Px,ParticleA[Aid].Py,ParticleA[Aid].Pz,AMass,ParticleB[Bid].Px,ParticleB[Bid].Py,ParticleB[Bid].Pz,BMass);
                                if (fabs(PairMass - C_Mass)<3*CMassSigma) {
                                    B_IfRecord[Bid] = false;
                                    A_IfRecord[Aid] = false;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        ParticleASizeR = 0;ParticleBSizeR = 0;
        for (Aid = 0;Aid < ParticleASize;Aid++) {
            if (A_IfRecord[Aid]) ParticleASizeR++;
        }
        for (Bid = 0;Aid < ParticleBSize;Bid++) {
            if (B_IfRecord[Bid]) ParticleBSizeR++;
        if ((ParticleASizeR * ParticleBSizeR) == 0) continue; // if the particle A and B are not found after cut.

        // FIll in the pool
        for (i=0;i<yBinNum;i++) {// initial
            IfMatched   [i]             = false;
            Tot_Pool_Tmp[i].ListA_Index = 0;
            Tot_Pool_Tmp[i].ListB_Index = 0;
        }
        for (Bid = 0;Bid < ParticleBSize;Bid++) {
            if (B_IfRecord[Bid]) {
                yIndex = B_yIndex[Bid];
                IfMatched   [yIndex]                                         = true;
                Tot_Pool_Tmp[yIndex].EvtID                                   = EntriesID;
                Tot_Pool_Tmp[yIndex].ListB[Tot_Pool_Tmp[yIndex].ListB_Index] = ParticleB[Bid];
                Tot_Pool_Tmp[yIndex].ListB_Index++;
            }
        }
        for (i=0;i<yBinNum;i++) {
            if (IfMatched[i]) {
                j = 0;
                for (Aid = 0;Aid < ParticleASize;Aid++) {
                    if (A_IfRecord[Aid]) {
                        Tot_Pool_Tmp[i].ListA[j] = ParticleA[Aid];
                        j++;
                    }
                }
                Tot_Pool_Tmp[i].ListA_Index = j-1;
                Tot_Pool_Num[CenIndex] [i] [PVzIndex]++;
                if (Tot_Pool_Num[CenIndex] [i] [PVzIndex] == (HowMuchEventMixing)) {
                    Tot_Pool_IfFilled[CenIndex] [i] [PVzIndex] = true;
                }
                Tot_Pool [CenIndex] [i] [PVzIndex] [Tot_Pool_Num[CenIndex] [i] [PVzIndex]] = Tot_Pool_Tmp[i];
            }
        }
        // ############################################################################################################# //
        // ####                                          Fill in the Hist                                           #### //
        // ############################################################################################################# //
        for (i=0;i<yBinNum;i++) {
            if (Tot_Pool_IfFilled[CenIndex] [i] [PVzIndex]) {
                for (j=0;j<=HowMuchEventMixing;j++) {
                    Tot_Pool_TTmp = Tot_Pool[CenIndex] [i] [PVzIndex] [j];
                    for (k=0;k<=HowMuchEventMixing;k++) {
                        Tot_Pool_TTTmp = Tot_Pool[CenIndex] [i] [PVzIndex] [k];
                        for (Bid=0;Bid<Tot_Pool_TTmp.ListB_Index;Bid++) {
                            for (Aid=0;Aid<Tot_Pool_TTTmp.ListA_Index;Aid++) {
                                ParticleA_Tmp = Tot_Pool_TTTmp.ListA[Aid];
                                ParticleB_Tmp = Tot_Pool_TTmp.ListB[Bid];
                                GetPairMassAndKstar(ParticleB_Tmp , ParticleA_Tmp , AMass , BMass , MassAndKstar);
                                Mass_Store[m] = MassAndKstar[0];
                                Kstar_Store[m] = MassAndKstar[1];
                                tRap = ParticleA_Tmp.Rap - ParticleB_Tmp.Rap;
                                tPt  = ParticleA_Tmp.Pt  - ParticleB_Tmp.Pt ;
                                if (j == k) {
                                    H_Event_Num      [CenIndex] [i] [PVzIndex]->Fill(0);
                                    H_Res_Event_Num  [CenIndex] [i] [PVzIndex]->Fill(0);
                                    H_ALL_Event_Num             [i]           ->Fill(0);
                                    H_ALL_Res_Event_Num         [i]           ->Fill(0);
                                    H_Mass           [CenIndex] [i] [PVzIndex]->Fill(MassAndKstar[0]);
                                    H_Kstar          [CenIndex] [i] [PVzIndex]->Fill(MassAndKstar[1]);
                                    H_dRap           [CenIndex] [i] [PVzIndex]->Fill(tRap);
                                    H_dPt            [CenIndex] [i] [PVzIndex]->Fill(tPt);
                                    H_Res_Mass       [CenIndex] [i] [PVzIndex]->Fill(MassAndKstar[0]);
                                    H_Res_Kstar      [CenIndex] [i] [PVzIndex]->Fill(MassAndKstar[1]);
                                    H_Res_dRap       [CenIndex] [i] [PVzIndex]->Fill(tRap);
                                    H_Res_dPt        [CenIndex] [i] [PVzIndex]->Fill(tPt);
                                    H_ALL_Mass                  [i]           ->Fill(MassAndKstar[0]);
                                    H_ALL_Kstar                 [i]           ->Fill(MassAndKstar[1]);
                                    H_ALL_dRap                  [i]           ->Fill(tRap);
                                    H_ALL_dPt                   [i]           ->Fill(tPt);
                                    H_ALL_Res_Mass              [i]           ->Fill(MassAndKstar[0]);
                                    H_ALL_Res_Kstar             [i]           ->Fill(MassAndKstar[1]);
                                    H_ALL_Res_dRap              [i]           ->Fill(tRap);
                                    H_ALL_Res_dPt               [i]           ->Fill(tPt);
                                }
                                else {
                                    for (m=0;m<2;m++) {
                                        H_Mix_A_Num_dRap[i] [MassAndKstar[m].RapIndex]->Fill(MassAndKstar[m].Mass);
                                        H_Mix_B_Num_dRap[i] [MassAndKstar[m].RapIndex]->Fill(MassAndKstar[m].Mass);
                                    }
                                }
                            }
                        }
                    }
                }
                Tot_Pool_Num[CenIndex] [i] [PVzIndex] = -1;
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

bool IfInVector(int Num , std::vector<unsigned short int> V)
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

void GetPairMassAndKstar(Particle PA , Particle PB , float AMass , float BMass , float (&MassAndKstar)[2]) {
    float p1x = PA.Px ,p1y = PA.Py , p1z = PA.Pz  , p2x = PB.Px , p2y = PB.Py ,  p2z = PB.Pz;
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

    MassAndKstar[0] = (gamma * (E1 + bp1 + E2 + bp2));
    MassAndKstar[1] = 0.5*pow(New_Px*New_Px+New_Py*New_Py+New_Pz*New_Pz,0.5);
}