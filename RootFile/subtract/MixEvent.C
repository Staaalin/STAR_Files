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
#include "../../KFParticle4Lambda/StRoot/StKFParticleAnalysisMaker/TriggerList.h"
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

int CentralityBin[] = {0,10,20,40,60,80,100};// %
#define CentralityBinNum 6

float PtBin[] = {0 , 0.5 , 1.0 , 1.5 , 2.0 , 3.0 , 10.0}; // Pt
#define PtBinNum 6

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
            case 211 :
                Result = 0.13957;
                break;
            case -211 :
                Result = 0.13957;
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
                Result = 0;
        }
    }
    return Result;
}

std::vector<int> GetNchList(int CentralityList[])
{
    std::vector<int> Result;Result.resize(0);
    int CentralityListSize = sizeof(CentralityList)/sizeof(CentralityList[0]);
    if (DataName == "dAu_200_21") {
        // data from https://drupal.star.bnl.gov/STAR/system/files/pwg5.pdf
        int NchTable[21] = [ 10000 , 55 , 47 , 42 , 38 , 35 , 32 , 29 , 26 , 24 , 21 , 19 , 17 , 15 , 13 , 11 , 9 , 7 , 6 , 4 ,  0];
        int CenTable[21] = [     0 ,  5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 45 , 50 , 55 , 60 , 65 , 70 , 75 ,80 ,85 ,90 ,95 ,100];
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

void MixEvent(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,
              float Energy)
{

    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) 
        std::vector<Int_t> *PDGMult              = nullptr;
        std::vector<Int_t> *refMult              = nullptr;
        std::vector<Int_t> *grefMult             = nullptr;
        std::vector<Int_t> *EventID              = nullptr;
        std::vector<Int_t> *RunID                = nullptr;
        std::vector<Int_t> *TriggerID            = nullptr;
        
        TBranch *bPDGMult                        = nullptr;
        TBranch *brefMult                        = nullptr;
        TBranch *bgrefMult                       = nullptr;
        TBranch *bEventID                        = nullptr;
        TBranch *bRunID                          = nullptr;
        TBranch *bTriggerID                      = nullptr;

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
    
    #else
        #if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
            std::vector<Int_t> *PDGMult              = NULL;
            std::vector<Int_t> *refMult              = NULL;
            std::vector<Int_t> *grefMult             = NULL;
            std::vector<Int_t> *EventID              = NULL;
            std::vector<Int_t> *RunID                = NULL;
            std::vector<Int_t> *TriggerID            = NULL;

            TBranch *bPDGMult                        = NULL;
            TBranch *brefMult                        = NULL;
            TBranch *bgrefMult                       = NULL;
            TBranch *bEventID                        = NULL;
            TBranch *bRunID                          = NULL;
            TBranch *bTriggerID                      = NULL;

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

        #else
            std::vector<Int_t> *PDGMult              = 0;
            std::vector<Int_t> *refMult              = 0;
            std::vector<Int_t> *grefMult             = 0;
            std::vector<Int_t> *EventID              = 0;
            std::vector<Int_t> *RunID                = 0;
            std::vector<Int_t> *TriggerID            = 0;
            
            TBranch *bPDGMult                        = 0;
            TBranch *brefMult                        = 0;
            TBranch *bgrefMult                       = 0;
            TBranch *bEventID                        = 0;
            TBranch *bRunID                          = 0;
            TBranch *bTriggerID                      = 0;
    
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

        #endif
    #endif

    // TString FileName = "output_";
    // cout<<"Start Running"<<endl;
    // cout<<StartFileIndex<<endl;
    // cout<<EndFileIndex<<endl;



    //load data  
    TChain *hadronTree = new TChain("hadronTree");
    for(int i=StartFileIndex;i <= EndFileIndex;i++){
        TString filename = MidName;
        filename+=i;
        filename+=".root";
        hadronTree->Add(filename);
        // cout<<"Add "<<filename<<" Successfully"<<endl;
    }
    
    Int_t mult,npart;

    hadronTree->SetBranchAddress("PDGMult"  ,&PDGMult  );
    hadronTree->SetBranchAddress("refMult"  ,&refMult  );
    hadronTree->SetBranchAddress("grefMult" ,&grefMult );
    hadronTree->SetBranchAddress("EventID"  ,&EventID  );
    hadronTree->SetBranchAddress("RunID"    ,&RunID    );
    hadronTree->SetBranchAddress("TriggerID",&TriggerID);
    
    hadronTree->SetBranchAddress("PDG"          ,&PDG          ,&bPDG          );
    hadronTree->SetBranchAddress("mix_px"       ,&mix_px       ,&bmix_px       );
    hadronTree->SetBranchAddress("mix_py"       ,&mix_py       ,&bmix_py       );
    hadronTree->SetBranchAddress("mix_pz"       ,&mix_pz       ,&bmix_pz       );
    hadronTree->SetBranchAddress("QA_eta"       ,&QA_eta       ,&bQA_eta       );
    hadronTree->SetBranchAddress("dEdx"         ,&dEdx         ,&bdEdx         );
    hadronTree->SetBranchAddress("m2"           ,&m2           ,&bm2           );
    hadronTree->SetBranchAddress("dcatopv"      ,&dcatopv      ,&bdcatopv      );
    hadronTree->SetBranchAddress("nSigmaProton" ,&nSigmaProton ,&bnSigmaProton );
    hadronTree->SetBranchAddress("nSigmaPion"   ,&nSigmaPion   ,&bnSigmaPion   );
    hadronTree->SetBranchAddress("nSigmaKaon"   ,&nSigmaKaon   ,&bnSigmaKaon   );
    hadronTree->SetBranchAddress("InvariantMass",&InvariantMass,&bInvariantMass);
    hadronTree->SetBranchAddress("Decay_Length" ,&Decay_Length ,&bDecay_Length );
    hadronTree->SetBranchAddress("Chi2"         ,&Chi2         ,&bChi2         );

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;
    
    TLorentzVector p0,p3;
    double kstar, rap;
    TVector3 BetaTemp;

    int KaonpPID = 321,KaonmPID = -321,PionpPID = 211,PionmPID = -211,LambdaPID = 3122,LambdabPID = -3122,XiPID = 3312,XibPID = -3312,OmegaPID = 3334,OmegabPID = -3334;

    std::vector<int> NchList = GetNchList(CentralityBin);     // centrality
    std::vector<int> TriggerBin;                              // trigger
    float PtBin[] = {0 , 0.5 , 1.0 , 1.5 , 2.0 , 3.0 , 10.0}; // Pt
    //                                   centrality       trigger      Pt
    std::vector<Float_t> A_Px        [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> A_Py        [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> A_Pz        [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> B_Px        [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> B_Py        [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> B_Pz        [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> Mix_A_Px    [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> Mix_A_Py    [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> Mix_A_Pz    [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> Mix_B_Px    [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> Mix_B_Py    [CentralityBinNum]    [30]    [PtBinNum];
    std::vector<Float_t> Mix_B_Pz    [CentralityBinNum]    [30]    [PtBinNum];
    TH1D* H_Kstar                    [CentralityBinNum]    [30]    [PtBinNum];
    TH1D* H_Mix_Kstar                [CentralityBinNum]    [30]    [PtBinNum];

    for (int i=0;i<CentralityBinNum;i++){
        TString HistName1 = "H_";
        TString HistName2 = "Cen: [";
        HistName1 += i;HistName1 += "_";
        HistName2 += CentralityBin[i];HistName2 += "% , ";
        HistName2 += CentralityBin[i+1];HistName2 += "%]";
        for (int j=0;i<30;j++){
            HistName1 += j;HistName1 += "_";
            HistName2 += CentralityBin[i];HistName2 += "% , ";
            HistName2 += CentralityBin[i+1];HistName2 += "%]";
            for (int k=0;k<=PtBinNum;k++)
        }
    }

    for (int i=0;i<nentries;i++){
        // if (i > 15) {break;}
        hadronTree->GetEntry(i);
        cout<<"Calculating Event "<<(i+1)<<"/"<<nentries<<endl;
        // cout<<mult<<endl;
        // if(b>7){continue;}
        // cout<<"There OK"<<endl;

    }

    TString OutputFileName = "Cor_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *file = new TFile(OutputFileName, "RECREATE");

    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
        dNdy_Kaonp[Itr]->Write();
        dNdy_Kaonm[Itr]->Write();
        dNdy_Pionp[Itr]->Write();
        dNdy_Pionm[Itr]->Write();
    }
    for (int Itr = 0;Itr < HyperonPhaseNum;Itr++) {
        dNdy_Lambda[Itr]  ->Write();
        dNdy_Lambdab[Itr] ->Write();
        dNdy_Omega[Itr]   ->Write();
        dNdy_Omegab[Itr]  ->Write();
        dNdy_Xi[Itr]      ->Write();
        dNdy_Xib[Itr]     ->Write();
    }
    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
        for (int Jtr = 0;Jtr < HyperonPhaseNum;Jtr++) {

            k_Kp_Omega   [Itr][Jtr]->Write();
            k_Kp_Omegab  [Itr][Jtr]->Write();
            k_Km_Omega   [Itr][Jtr]->Write();
            k_Km_Omegab  [Itr][Jtr]->Write();
            k_Kp_Xi      [Itr][Jtr]->Write();
            k_Kp_Xib     [Itr][Jtr]->Write();
            k_Km_Xi      [Itr][Jtr]->Write();
            k_Km_Xib     [Itr][Jtr]->Write();
            k_Kp_Lambda  [Itr][Jtr]->Write();
            k_Kp_Lambdab [Itr][Jtr]->Write();
            k_Km_Lambda  [Itr][Jtr]->Write();
            k_Km_Lambdab [Itr][Jtr]->Write();
            k_Pip_Omega  [Itr][Jtr]->Write();
            k_Pip_Omegab [Itr][Jtr]->Write();
            k_Pim_Omega  [Itr][Jtr]->Write();
            k_Pim_Omegab [Itr][Jtr]->Write();
            k_Pip_Xi     [Itr][Jtr]->Write();
            k_Pip_Xib    [Itr][Jtr]->Write();
            k_Pim_Xi     [Itr][Jtr]->Write();
            k_Pim_Xib    [Itr][Jtr]->Write();
            k_Pip_Lambda [Itr][Jtr]->Write();
            k_Pip_Lambdab[Itr][Jtr]->Write();
            k_Pim_Lambda [Itr][Jtr]->Write();
            k_Pim_Lambdab[Itr][Jtr]->Write();
        }
    }

    file->Write();

    return;
}
