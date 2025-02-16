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

#define A_Num_Per_Event 5
#define B_Num_Per_Event 5
#define HowMuchEventMixing 10

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

    const int CentralityBin[] = {0 , 10 , 30 , 50 , 100};// %
    const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1

    int MBinNum = 500 , MBinPar = 50;
    float MSta , MEnd;

    float Mass_Lambda = 1.1161;
    MSta = floor((Mass_Lambda)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;
    TH3D* H_Lambda  = new TH3D("H_Lambda" ,"Lambda_Distribution"    , 40,-2,2 , 100,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH3D* H_Lambdab = new TH3D("H_Lambdab","LambdaBar_Distribution" , 40,-2,2 , 100,MSta,MEnd , CentralityBinNum,0,CentralityBinNum);
    TH2D* H_ALL_Lambda  = new TH2D("H_ALL_Lambda" ,"Lambda_Distribution"    , 40,-2,2 , 100,MSta,MEnd);
    TH2D* H_ALL_Lambdab = new TH2D("H_ALL_Lambdab","LambdaBar_Distribution" , 40,-2,2 , 100,MSta,MEnd);

    float Mass_Xi     = 1.3223;
    MSta = floor((Mass_Xi)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;
    TH3D* H_Xi      = new TH3D("H_Xi"     ,"Xi_Distribution"        , 40,-2,2 , 100,MSta,MEnd         , CentralityBinNum,0,CentralityBinNum);
    TH3D* H_Xib     = new TH3D("H_Xib"    ,"XiBar_Distribution"     , 40,-2,2 , 100,MSta,MEnd         , CentralityBinNum,0,CentralityBinNum);
    TH2D* H_ALL_Xi      = new TH2D("H_ALL_Xi"     ,"Xi_Distribution"        , 40,-2,2 , 100,MSta,MEnd);
    TH2D* H_ALL_Xib     = new TH2D("H_ALL_Xib"    ,"XiBar_Distribution"     , 40,-2,2 , 100,MSta,MEnd);

    float Mass_Omega  = 1.6725;
    MSta = floor((Mass_Omega)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;
    TH3D* H_Omega   = new TH3D("H_Omega"  ,"Omega_Distribution"     , 40,-2,2 , 100,MSta,MEnd   , CentralityBinNum,0,CentralityBinNum);
    TH3D* H_Omegab  = new TH3D("H_Omegab" ,"OmegaBar_Distribution"  , 40,-2,2 , 100,MSta,MEnd   , CentralityBinNum,0,CentralityBinNum);
    TH2D* H_ALL_Omega   = new TH2D("H_ALL_Omega"  ,"Omega_Distribution"     , 40,-2,2 , 100,MSta,MEnd);
    TH2D* H_ALL_Omegab  = new TH2D("H_ALL_Omegab" ,"OmegaBar_Distribution"  , 40,-2,2 , 100,MSta,MEnd);

    int i , j , k , l , m , n;// used as Index
    int CenIndex;
    int NNch;
    float tEnergy , rap;

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
    hadronTree->SetBranchAddress("dcatopv"      ,&dcatopv      ,&bdcatopv      );
    // hadronTree->SetBranchAddress("nSigmaProton" ,&nSigmaProton ,&bnSigmaProton );
    // hadronTree->SetBranchAddress("nSigmaPion"   ,&nSigmaPion   ,&bnSigmaPion   );
    hadronTree->SetBranchAddress("nSigmaKaon"   ,&nSigmaKaon   ,&bnSigmaKaon   );
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
    
        if ((EntriesID+1)%200 == 0) {
            time(&time_now);
            int time_diff = (int)difftime(time_now, time_start);
            cout << time_diff/60 << "min " << time_diff%60 << "s: ";
            long long microseconds = (clock() - Tstart)/10000;
            std::cout << "Microseconds: " << microseconds << "  ";
            cout<<"Calculating Event "<<(EntriesID+1)<<"/"<<nentries<<endl;
            Tstart = clock();
        }
    
        for (i=0;i<PDGMult;i++){
            if       (PDG->at(i) == 3122) {
                tEnergy = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2) + pow(mix_pz->at(i),2) + Mass_Lambda*Mass_Lambda,0.5);
                rap     = 0.5*log((tEnergy+mix_pz->at(i))/(tEnergy-mix_pz->at(i)));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Lambda->Fill(rap,InvariantMass->at(i),CenIndex);
                H_ALL_Lambda->Fill(rap,InvariantMass->at(i));
            }
            else if (PDG->at(i) == -3122) {
                tEnergy = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2) + pow(mix_pz->at(i),2) + Mass_Lambda*Mass_Lambda,0.5);
                rap     = 0.5*log((tEnergy+mix_pz->at(i))/(tEnergy-mix_pz->at(i)));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Lambdab->Fill(rap,InvariantMass->at(i),CenIndex);
                H_ALL_Lambdab->Fill(rap,InvariantMass->at(i));
            }
            else if (PDG->at(i) == 3312) {
                tEnergy = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2) + pow(mix_pz->at(i),2) + Mass_Xi*Mass_Xi,0.5);
                rap     = 0.5*log((tEnergy+mix_pz->at(i))/(tEnergy-mix_pz->at(i)));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Xi->Fill(rap,InvariantMass->at(i),CenIndex);
                H_ALL_Xi->Fill(rap,InvariantMass->at(i));
            }
            else if (PDG->at(i) == -3312) {
                tEnergy = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2) + pow(mix_pz->at(i),2) + Mass_Xi*Mass_Xi,0.5);
                rap     = 0.5*log((tEnergy+mix_pz->at(i))/(tEnergy-mix_pz->at(i)));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Xib->Fill(rap,InvariantMass->at(i),CenIndex);
                H_ALL_Xib->Fill(rap,InvariantMass->at(i));
            }
            else if (PDG->at(i) == 3334) {
                tEnergy = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2) + pow(mix_pz->at(i),2) + Mass_Omega*Mass_Omega,0.5);
                rap     = 0.5*log((tEnergy+mix_pz->at(i))/(tEnergy-mix_pz->at(i)));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Omega->Fill(rap,InvariantMass->at(i),CenIndex);
                H_ALL_Omega->Fill(rap,InvariantMass->at(i));
            }
            else if (PDG->at(i) == -3334) {
                tEnergy = pow(pow(mix_px->at(i),2) + pow(mix_py->at(i),2) + pow(mix_pz->at(i),2) + Mass_Omega*Mass_Omega,0.5);
                rap     = 0.5*log((tEnergy+mix_pz->at(i))/(tEnergy-mix_pz->at(i)));
                CenIndex = 0;
                for (k=0;k<CentralityBinNum;k++){
                    NNch = CenCorr(PVz) * Nch;
                    // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
                    if ((NchList.at(k) >= NNch) && (NNch > NchList.at(k+1))) {
                        CenIndex = k;
                        break;
                    }
                }
                H_Omegab->Fill(rap,InvariantMass->at(i),CenIndex);
                H_ALL_Omegab->Fill(rap,InvariantMass->at(i));
            }
        }
    }
    TString OutputFileName = OutMidName;
    OutputFileName += "H_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileA = new TFile(OutputFileName, "RECREATE");
    fileA->cd();
    TH3D* H_Lambda      -> Write();
    TH3D* H_Lambdab     -> Write();
    TH2D* H_ALL_Lambda  -> Write();
    TH2D* H_ALL_Lambdab -> Write();
    TH3D* H_Xi          -> Write();
    TH3D* H_Xib         -> Write();
    TH2D* H_ALL_Xi      -> Write();
    TH2D* H_ALL_Xib     -> Write();
    TH3D* H_Omega       -> Write();
    TH3D* H_Omegab      -> Write();
    TH2D* H_ALL_Omega   -> Write();
    TH2D* H_ALL_Omegab  -> Write();
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