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

// #define DataName           "pAu_200_15"
// #define DataName           "AuAu_27_18"
// #define DataName           "dAu_200_16"
#define DataName           "dAu_200_21"
// #define DataName           "dAu_62_16"
// #define DataName           "dAu_39_16"
// #define DataName           "dAu_20_16"
// #define DataName           "pp_200_15"
// #define DataName           "OO_200_21"

float CenCorr(float Vz);
std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize);
void print(std::vector<int> Temp);
void print(std::vector<float> Temp);
float GetPairMass(float p1x , float p1y , float p1z , float p2x , float p2y , float p2z , float AMass , float BMass);
Double_t massList(int PID);

#define A_Num_Per_Event 5
#define B_Num_Per_Event 5
#define RotNum 10

void HR(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,TString OutMidName,int Mode = 0) {
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
        std::vector<int>     *DaughtersID        = nullptr;

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
        TBranch *bDaughtersID                    = nullptr;
    
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
            std::vector<int>     *DaughtersID        = NULL;

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
            TBranch *bDaughtersID                    = NULL;

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
            std::vector<int>     *DaughtersID        = 0;
    
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
            TBranch *bDaughtersID                    = 0;

        #endif
    #endif

    const int CentralityBin[] = {0 , 10 , 30 , 50 , 100};// %
    const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1

    int MBinNum = 500 , MBinPar = 50;
    float MSta , MEnd;

    float Mass_Lambda = 1.1161 , Mass_Sigma_Lambda = 0.0020;
    MSta = floor((Mass_Lambda)/0.0005 - 50)*0.0005 , MEnd = floor((Mass_Lambda)/0.0005 + 50)*0.0005;
    TH3D* H_Lambda       = new TH3D("H_Lambda" ,"Lambda_Distribution"        , 40,-2,2 , 200,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH3D* H_Lambdab      = new TH3D("H_Lambdab","LambdaBar_Distribution"     , 40,-2,2 , 200,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH2D* H_ALL_Lambda   = new TH2D("H_ALL_Lambda" ,"Lambda_Distribution"    , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALL_Lambdab  = new TH2D("H_ALL_Lambdab","LambdaBar_Distribution" , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALLr_Lambda  = new TH2D("H_ALLr_Lambda" ,"Lambda_Distribution"   , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALLr_Lambdab = new TH2D("H_ALLr_Lambdab","LambdaBar_Distribution", 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_Lambda_Pt_y  = new TH2D("H_Lambda_Pt_y" ,"Lambda_Pt_y"           , 272,-1.7,1.7, 224,0,2.8);
    TH2D* H_Lambdab_Pt_y = new TH2D("H_Lambdab_Pt_y","LambdaBar_Pt_y"        , 272,-1.7,1.7, 224,0,2.8);

    float Mass_Xi     = 1.3223 , Mass_Sigma_Xi = 0.0024;
    MSta = floor((Mass_Xi)/0.0005 - 50)*0.0005 , MEnd = floor((Mass_Xi)/0.0005 + 50)*0.0005;
    TH3D* H_Xi           = new TH3D("H_Xi"     ,"Xi_Distribution"            , 40,-2,2 , 200,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH3D* H_Xib          = new TH3D("H_Xib"    ,"XiBar_Distribution"         , 40,-2,2 , 200,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH2D* H_ALL_Xi       = new TH2D("H_ALL_Xi"     ,"Xi_Distribution"        , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALL_Xib      = new TH2D("H_ALL_Xib"    ,"XiBar_Distribution"     , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALLr_Xi      = new TH2D("H_ALLr_Xi"    ,"Xi_Distribution"        , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALLr_Xib     = new TH2D("H_ALLr_Xib"   ,"XiBar_Distribution"     , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_Xi_Pt_y      = new TH2D("H_Xi_Pt_y" ,"Xi_Pt_y"                   , 136,-1.7,1.7, 112,0,2.8);
    TH2D* H_Xib_Pt_y     = new TH2D("H_Xib_Pt_y","XiBar_Pt_y"                , 136,-1.7,1.7, 112,0,2.8);

    float Mass_Omega  = 1.6725 , Mass_Sigma_Omega = 0.0029;
    MSta = floor((Mass_Omega)/0.0005 - 50)*0.0005 , MEnd = floor((Mass_Omega)/0.0005 + 50)*0.0005;
    TH3D* H_Omega       = new TH3D("H_Omega"  ,"Omega_Distribution"         , 40,-2,2 , 200,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH3D* H_Omegab      = new TH3D("H_Omegab" ,"OmegaBar_Distribution"      , 40,-2,2 , 200,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH2D* H_ALL_Omega   = new TH2D("H_ALL_Omega"  ,"Omega_Distribution"     , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALL_Omegab  = new TH2D("H_ALL_Omegab" ,"OmegaBar_Distribution"  , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALLr_Omega  = new TH2D("H_ALLr_Omega" ,"Omega_Distribution"     , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALLr_Omegab = new TH2D("H_ALLr_Omegab","OmegaBar_Distribution"  , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_Omega_Pt_y  = new TH2D("H_Omega_Pt_y" ,"Omega_Pt_y"             , 68,-1.7,1.7, 56,0,2.8);
    TH2D* H_Omegab_Pt_y = new TH2D("H_Omegab_Pt_y","OmegaBar_Pt_y"          , 68,-1.7,1.7, 56,0,2.8);

    float Mass_Kaon = 0.493677;
    MSta = floor((Mass_Kaon)/0.0005 - 50)*0.0005 , MEnd = floor((Mass_Kaon)/0.0005 + 50)*0.0005;
    TH3D* H_Kaon        = new TH3D("H_Kaon" ,"Kaon_Distribution"            , 40,-2,2 , 200,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH3D* H_Kaonb       = new TH3D("H_Kaonb","KaonBar_Distribution"         , 40,-2,2 , 200,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH2D* H_ALL_Kaon    = new TH2D("H_ALL_Kaon" ,"Kaon_Distribution"        , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_ALL_Kaonb   = new TH2D("H_ALL_Kaonb","KaonBar_Distribution"     , 40,-2,2 , 200,MSta,MEnd);
    TH2D* H_Kaon_Pt_y   = new TH2D("H_Kaon_Pt_y" ,"Kaon_Pt_y"               , 200,-1,1, 160,0,1.6);
    TH2D* H_Kaonb_Pt_y  = new TH2D("H_Kaonb_Pt_y","KaonBar_Pt_y"            , 200,-1,1, 160,0,1.6);

    int i , j , k , l , m , n;// used as Index
    int CenIndex;
    int NNch , DID , D1id , D2id;
    float tEnergy , rap , Pt , Pz_T , Mass_T;
    float AMass , BMass , CMass[RotNum];
    float APx , APy , APz , BPx , BPy , BPz;
    float APr , BPr , Theta;

    TRandom3 randGen;

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
    hadronTree->SetBranchAddress("nHitsFit"     ,&nHitsFit     ,&bnHitsFit     );
    hadronTree->SetBranchAddress("nHitsMax"     ,&nHitsMax     ,&bnHitsMax     );
    hadronTree->SetBranchAddress("ParentList"   ,&ParentList   ,&bParentList   );
    hadronTree->SetBranchAddress("ParentSta"    ,&ParentSta    ,&bParentSta    );
    hadronTree->SetBranchAddress("ParentEnd"    ,&ParentEnd    ,&bParentEnd    );
    hadronTree->SetBranchAddress("SE_ParentList",&SE_ParentList,&bSE_ParentList   );
    hadronTree->SetBranchAddress("SE_ParentSta" ,&SE_ParentSta ,&bSE_ParentSta    );
    hadronTree->SetBranchAddress("SE_ParentEnd" ,&SE_ParentEnd ,&bSE_ParentEnd    );
    hadronTree->SetBranchAddress("ME_ParentList",&ME_ParentList,&bME_ParentList   );
    hadronTree->SetBranchAddress("ME_ParentSta" ,&ME_ParentSta ,&bME_ParentSta    );
    hadronTree->SetBranchAddress("ME_ParentEnd" ,&ME_ParentEnd ,&bME_ParentEnd    );
    hadronTree->SetBranchAddress("DaughtersID"  ,&DaughtersID  ,&bDaughtersID     );

    std::vector<int> NchList = GetNchList(CentralityBin , CentralityBinNum+1);     // centrality
    cout<<"NchList = ";
    print(NchList);
    cout<<" "<<endl;

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;
    
    time_t time_start;
    time_t time_now;
    time(&time_start);
    clock_t Tstart = clock();

    for (int EntriesID = 0 ; EntriesID < nentries ; EntriesID++) {
        hadronTree->GetEntry(EntriesID);
    
        if ((EntriesID+1)%5000 == 0) {
            time(&time_now);
            int time_diff = (int)difftime(time_now, time_start);
            cout << time_diff/60 << "min " << time_diff%60 << "s: ";
            long long microseconds = (clock() - Tstart)/10000;
            std::cout << "Microseconds: " << microseconds << "  ";
            cout<<"Calculating Event "<<(EntriesID+1)<<"/"<<nentries<<endl;
            Tstart = clock();
        }
    
        for (i=0;i<PDGMult;i++){
            Mass_T = InvariantMass->at(i);
            if ((Mass_T<0)) continue;
            Pz_T = mix_pz->at(i);
            if      (PDG->at(i) == 321) {
                Pt = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2),0.5);
                tEnergy = pow(Pt*Pt + Pz_T*Pz_T + Mass_Kaon*Mass_Kaon,0.5);
                rap     = 0.5*log((tEnergy+Pz_T)/(tEnergy-Pz_T));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Kaon->Fill(rap,Mass_T,CenIndex);
                H_ALL_Kaon->Fill(rap,Mass_T);
                H_Kaon_Pt_y->Fill(rap,Pt);
                continue;
            }
            if (PDG->at(i) == -321) {
                Pt = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2),0.5);
                tEnergy = pow(Pt*Pt + Pz_T*Pz_T + Mass_Kaon*Mass_Kaon,0.5);
                rap     = 0.5*log((tEnergy+Pz_T)/(tEnergy-Pz_T));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Kaonb->Fill(rap,Mass_T,CenIndex);
                H_ALL_Kaonb->Fill(rap,Mass_T);
                H_Kaonb_Pt_y->Fill(rap,Pt);
                continue;
            }
            if (PDG->at(i) == 3122) {
                Pt = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2),0.5);
                tEnergy = pow(Pt*Pt + Pz_T*Pz_T + Mass_Lambda*Mass_Lambda,0.5);
                rap     = 0.5*log((tEnergy+Pz_T)/(tEnergy-Pz_T));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Lambda->Fill(rap,Mass_T,CenIndex);
                H_ALL_Lambda->Fill(rap,Mass_T);
                if (fabs(Mass_T - Mass_Lambda) <= 3*Mass_Sigma_Lambda) {
                    H_Lambda_Pt_y->Fill(rap,Pt);
                }
            }
            else if (PDG->at(i) == -3122) {
                Pt = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2),0.5);
                tEnergy = pow(Pt*Pt + Pz_T*Pz_T + Mass_Lambda*Mass_Lambda,0.5);
                rap     = 0.5*log((tEnergy+Pz_T)/(tEnergy-Pz_T));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Lambdab->Fill(rap,Mass_T,CenIndex);
                H_ALL_Lambdab->Fill(rap,Mass_T);
                if (fabs(Mass_T - Mass_Lambda) <= 3*Mass_Sigma_Lambda) {
                    H_Lambdab_Pt_y->Fill(rap,Pt);
                }
            }
            else if (PDG->at(i) == 3312) {
                Pt = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2),0.5);
                tEnergy = pow(Pt*Pt + Pz_T*Pz_T + Mass_Xi*Mass_Xi,0.5);
                rap     = 0.5*log((tEnergy+Pz_T)/(tEnergy-Pz_T));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Xi->Fill(rap,Mass_T,CenIndex);
                H_ALL_Xi->Fill(rap,Mass_T);
                if (fabs(Mass_T - Mass_Xi) <= 3*Mass_Sigma_Xi) {
                    H_Xi_Pt_y->Fill(rap,Pt);
                }
            }
            else if (PDG->at(i) == -3312) {
                Pt = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2),0.5);
                tEnergy = pow(Pt*Pt + Pz_T*Pz_T + Mass_Xi*Mass_Xi,0.5);
                rap     = 0.5*log((tEnergy+Pz_T)/(tEnergy-Pz_T));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Xib->Fill(rap,Mass_T,CenIndex);
                H_ALL_Xib->Fill(rap,Mass_T);
                if (fabs(Mass_T - Mass_Xi) <= 3*Mass_Sigma_Xi) {
                    H_Xib_Pt_y->Fill(rap,Pt);
                }
            }
            else if (PDG->at(i) == 3334) {
                Pt = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2),0.5);
                tEnergy = pow(Pt*Pt + Pz_T*Pz_T + Mass_Omega*Mass_Omega,0.5);
                rap     = 0.5*log((tEnergy+Pz_T)/(tEnergy-Pz_T));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Omega->Fill(rap,Mass_T,CenIndex);
                H_ALL_Omega->Fill(rap,Mass_T);
                if (fabs(Mass_T - Mass_Omega) <= 3*Mass_Sigma_Omega) {
                    H_Omega_Pt_y->Fill(rap,Pt);
                }
            }
            else if (PDG->at(i) == -3334) {
                Pt = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2),0.5);
                tEnergy = pow(Pt*Pt + Pz_T*Pz_T + Mass_Omega*Mass_Omega,0.5);
                rap     = 0.5*log((tEnergy+Pz_T)/(tEnergy-Pz_T));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Omegab->Fill(rap,Mass_T,CenIndex);
                H_ALL_Omegab->Fill(rap,Mass_T);
                if (fabs(Mass_T - Mass_Omega) <= 3*Mass_Sigma_Omega) {
                    H_Omegab_Pt_y->Fill(rap,Pt);
                }
            }
            else{ continue; }

            DID = DaughtersID->at(i);
            if (DID <= 1000) continue;
            cout<<"DID = "<<DID<<endl;
            cout<<"PDGMult = "<<PDGMult<<endl;
            D1id = DID/1000;D2id = DID%1000;
            cout<<"D1id = "<<D1id<<endl;
            cout<<"D2id = "<<D2id<<endl;
            cout<<"___________"<<endl;
            if ((PDG->at(D1id) == -1) || (PDG->at(D2id) == -1)) continue;
            AMass = massList(PDG->at(D1id));BMass = massList(PDG->at(D2id));
            APx = mix_px->at(D1id);
            APy = mix_py->at(D1id);
            APz = mix_pz->at(D1id);
            BPx = mix_px->at(D2id);
            BPy = mix_py->at(D2id);
            BPz = mix_pz->at(D2id);
            APr = pow(APx*APx+APy*APy,0.5);
            BPr = pow(BPx*BPx+BPy*BPy,0.5);
            for (j=0;j<RotNum;j++) {
                Theta = randGen.Rndm() * 2 * 3.1415926535898;
                APx = APr*sin(Theta);
                APy = APr*cos(Theta);
                Theta = randGen.Rndm() * 2 * 3.1415926535898;
                BPx = BPr*sin(Theta);
                BPy = BPr*cos(Theta);
                CMass[j] = GetPairMass(APx , APy , APz , BPx , BPy , BPz , AMass , BMass);
            }
            if (PDG->at(i) ==  3122) {for (j=0;j<RotNum;j++) {H_ALLr_Lambda ->Fill(rap,CMass[j]);} continue;}
            if (PDG->at(i) == -3122) {for (j=0;j<RotNum;j++) {H_ALLr_Lambdab->Fill(rap,CMass[j]);} continue;}
            if (PDG->at(i) ==  3312) {for (j=0;j<RotNum;j++) {H_ALLr_Xi     ->Fill(rap,CMass[j]);} continue;}
            if (PDG->at(i) == -3312) {for (j=0;j<RotNum;j++) {H_ALLr_Xib    ->Fill(rap,CMass[j]);} continue;}
            if (PDG->at(i) ==  3334) {for (j=0;j<RotNum;j++) {H_ALLr_Omega  ->Fill(rap,CMass[j]);} continue;}
            if (PDG->at(i) == -3334) {for (j=0;j<RotNum;j++) {H_ALLr_Omegab ->Fill(rap,CMass[j]);} continue;}
        }
    }
    TString OutputFileName = OutMidName;
    OutputFileName += "H_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileA = new TFile(OutputFileName, "RECREATE");
    fileA->cd();
    H_Lambda       -> Write();
    H_Lambdab      -> Write();
    H_ALL_Lambda   -> Write();
    H_ALL_Lambdab  -> Write();
    H_ALLr_Lambda  -> Write();
    H_ALLr_Lambdab -> Write();
    H_Lambda_Pt_y  -> Write();
    H_Lambdab_Pt_y -> Write();
    H_Xi           -> Write();
    H_Xib          -> Write();
    H_ALL_Xi       -> Write();
    H_ALL_Xib      -> Write();
    H_ALLr_Xi      -> Write();
    H_ALLr_Xib     -> Write();
    H_Xi_Pt_y      -> Write();
    H_Xib_Pt_y     -> Write();
    H_Omega        -> Write();
    H_Omegab       -> Write();
    H_ALL_Omega    -> Write();
    H_ALL_Omegab   -> Write();
    H_ALLr_Omega   -> Write();
    H_ALLr_Omegab  -> Write();
    H_Omega_Pt_y   -> Write();
    H_Omegab_Pt_y  -> Write();
    H_Kaon         -> Write();
    H_Kaonb        -> Write();
    H_ALL_Kaon     -> Write();
    H_ALL_Kaonb    -> Write();
    H_Kaon_Pt_y    -> Write();
    H_Kaonb_Pt_y   -> Write();
    fileA->Close();
    return;
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

float GetPairMass(float p1x , float p1y , float p1z , float p2x , float p2y , float p2z , float AMass , float BMass) {
    float E1 = pow(p1x*p1x+p1y*p1y+p1z*p1z+AMass*AMass,0.5);
    float E2 = pow(p2x*p2x+p2y*p2y+p2z*p2z+BMass*BMass,0.5);
    float Tot_E = E1+E2;
    float beta[3] = { -(p1x+p2x)/Tot_E , -(p1y+p2y)/Tot_E , -(p1z+p2z)/Tot_E };
    float beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    float gamma = 1.0 / std::sqrt(1.0 - beta2);
    float gamma2 = (beta2 > 0) ? (gamma - 1.0) / beta2 : 0.0;

    float bp1 = beta[0]*p1x + beta[1]*p1y + beta[2]*p1z;
    float bp2 = beta[0]*p2x + beta[1]*p2y + beta[2]*p2z;

    return gamma * (E1 + bp1 + E2 + bp2);
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