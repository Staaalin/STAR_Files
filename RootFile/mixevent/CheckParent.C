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

void CheckParent(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,TString OutMidName) {

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

    const int ParName[] = { 3122 , 3312 , 3334 , -3122 , -3312 , -3334 };
    const int ParNameNum = sizeof(ParName)/sizeof(ParName[0]);
    vector<float> ML , MLS;
    for (int i=0;i<ParNameNum;i++) {
        ML .push_back(     massList(ParName[i]));
        MLS.push_back(massListSigma(ParName[i]));
    }

    int ParID = -1 , H1Index , H2Index;

    const int MBinNum = 500 , MBinPar = 50;
    float MSta , MEnd;

    TH1D* H_M [3][20];
    TString HistName1 , HistName2;

    for (int i=0;i<ParNameNum;i++) {
        MSta = floor(ML.at(i)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;
        for (int j = 0;j<3;j++){
            HistName1 = "HM_";
            HistName1 += ParName[i];
            HistName1 += "_";
            if (j == 0) HistName1 += "NoCut";
            if (j == 1) HistName1 += "HaveParent";
            if (j == 2) HistName1 += "StrictParent";
            H_M [j] [i] = new TH1D(HistName1,HistName1,MBinNum,MSta,MEnd);
        }
    }

    TString TreeName = "hadronTree";
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
            cout<<"Calculating Event "<<(EntriesID+1)<<"/"<<nentries<<endl;
            Tstart = clock();
        }

        for (int j=0;j<PDGMult;j++) {
            H2Index = -1;
            for (int k=0;k<ParNameNum;k++) {
                if (PDG->at(j) == ParName[k]) H2Index = k;
            }
            if (H2Index == -1) continue;
            
            H_M[0][H2Index] -> Fill(InvariantMass->at(j));
            if (ParentEnd->at(j) >= ParentSta->at(j)) {
                H_M[1][H2Index] -> Fill(InvariantMass->at(j));
            }
        }
    }

    TString OutputFileName = OutMidName;
    OutputFileName += "H_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileA = new TFile(OutputFileName, "RECREATE");
    fileA->cd();
    
    for (int i=0;i<ParNameNum;i++) {
        for (int j = 0;j<3;j++){
            H_M [j] [i]->Write();
        }
    }
    fileA->Close();


    return;
}