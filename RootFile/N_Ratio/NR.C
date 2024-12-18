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
bool IfInVector(int Num , std::vector<unsigned short int> V);
bool IfCommonElement(std::vector<int> A , std::vector<int> B);
void DltElement(std::vector<int> &V , int ID);
std::vector<int> GetDaughterPDGLit(int ID);
std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize);

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
    vector<Particle> ListA;
    vector<Particle> ListB;
};

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

    int i , j , k , l , m , n;// used as Index
    int Aid , Bid , Cid;// used as Index
    std::vector<int> Temp;
    std::vector<float> CMass , CMassSigma;
    Particle AnyParticle;
    std::vector<Particle> ParticleA     , ParticleB;
    unsigned short int    ParticleASize , ParticleBSize , ParticleCSize;
    std::vector<std::vector<unsigned short int> > A_ParID,B_ParID,C_ParID;
    std::vector<uint8_t> A_IfRecord,B_IfRecord;
    bool IfRecord = true , IfRemoveFeedPair = false;
    float BMass = massList(B_PDG)           , AMass = massList(A_PDG);
    float BMassSigma = massListSigma(B_PDG) , AMassSigma = massListSigma(A_PDG);

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
    ParticleA.clear();ParticleB.clear();A_ParID.clear();B_ParID.clear();C_ParID.clear();A_IfRecord.clear();B_IfRecord.clear();
    for (int EntriesID = 0 ; EntriesID < nentries ; EntriesID++) {
        hadronTree->GetEntry(EntriesID);
        for (j=0;j<PDGMult;j++) {
            if (PDG->at(j) == A_PDG) {
                if (fabs(InvariantMass->at(j) - AMass) <= 3*AMassSigma) {
                    AnyParticle.TreeID = j;
                    ParticleA.push_back(AnyParticle);
                }
                else{continue;}
            }
            else if (PDG->at(j) == B_PDG) {
                if (fabs(InvariantMass->at(j) - BMass) <= 3*AMassSigma) {
                    AnyParticle.TreeID = j;
                    ParticleB.push_back(AnyParticle);
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
                        // IfFoundOmega = true;
                        // cout<<"Found Omega"<<endl;
                    }
                }
            }
        }

        ParticleASize = ParticleA.size();ParticleBSize = ParticleB.size();
        if ((ParticleASize * ParticleBSize) == 0) continue; // if the particle A and B are not found.

        for (i=0;i<ParticleASize;i++) {
            j = ParticleA[i].TreeID;
            Temp.clear();Temp.push_back(j);
            for (k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                Temp.push_back(ParentList->at(k));
            }
            A_ParID.push_back(Temp);
            A_IfRecord.push_back(1);
        }
        for (i=0;i<ParticleBSize;i++) {
            j = ParticleB[i].TreeID;
            Temp.clear();Temp.push_back(j);
            for (k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                Temp.push_back(ParentList->at(k));
            }
            B_ParID.push_back(Temp);
            B_IfRecord.push_back(1);
        }
        
        // 如果A、B有血缘关系，保留B
        for (Bid = 0;Bid < ParticleBSize;Bid++) {
            for (Aid = 0;Aid < ParticleASize;;Aid++) {
                if (IfInVector(ParticleA[Aid].TreeID , B_ParID.at(Bid))){
                    A_IfRecord.at(Aid) = 0;
                }
                // if (IfCommonElement(A_ParID.at(Aid) , B_ParID.at(Bid))){
                //     A_IfRecord.at(Aid) = 0;
                //     // cout<<"Meet 2!"<<endl;
                //     // cout<<"{ "<<A_PDG<<" } "<<A_TreID.at(Aid)<<" th ";print(A_ParID.at(Aid));
                //     // cout<<"{ "<<B_PDG<<" } "<<B_TreID.at(Bid)<<" th ";print(B_ParID.at(Bid));
                // }
            }
        }
        
        ParticleCSize = C_ParID.size();
        // 如果A、B与C有血缘关系，不记录A和B
        for (Aid = 0;Aid < ParticleASize;Aid++) {
            for (Cid = 0;Cid < ParticleCSize;Cid++) {
                if (IfInVector(ParticleA[Aid].TreeID , C_ParID.at(Cid))) {
                    A_IfRecord.at(Aid) = 0;
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
                    B_IfRecord.at(Bid) = 0;
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

        for (Aid = 0;Aid < ParticleASize;Aid++) {
            if (A_IfRecord.at(Aid)==0) continue;
            i = ParticleA[Aid].TreeID;
            ParticleA[Aid].Px = mix_px->at(i);
            ParticleA[Aid].Py = mix_py->at(i);
            ParticleA[Aid].Pz = mix_pz->at(i);
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