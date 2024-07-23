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

int CentralityBin[] = {0,25,50,75,100};// %
#define CentralityBinNum 4

float PtBin[] = {0 , 1.0 , 10.0}; // Pt
#define PtBinNum 2

float yBin[] = {-1.0 , 0.0 , 1.0}; // B_y
#define yBinNum 2

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

void MixEvent(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,
              int A_PDG,int B_PDG)
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
    Int_t PDGMult  ;
    Int_t refMult  ;
    Int_t grefMult ;
    Int_t EventID  ;
    Int_t RunID    ;
    Int_t TriggerID;

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
    
    double kstar, rap;
    TVector3 BetaTemp;

    int KaonpPID = 321,KaonmPID = -321,PionpPID = 211,PionmPID = -211,LambdaPID = 3122,LambdabPID = -3122,XiPID = 3312,XibPID = -3312,OmegaPID = 3334,OmegabPID = -3334;

    std::vector<int> NchList = GetNchList(CentralityBin);     // centrality
    //                                        centrality          B_y        B_Pt
    std::vector<Float_t> A_Px        ;
    std::vector<Float_t> A_Py        ;
    std::vector<Float_t> A_Pz        ;
    std::vector<Float_t> A_EvtID     ;
    std::vector<Float_t> A_TreID     ;
    std::vector<Float_t> B_Px        ;
    std::vector<Float_t> B_Py        ;
    std::vector<Float_t> B_Pz        ;
    std::vector<Float_t> B_EvtID     ;
    std::vector<Float_t> B_TreID     ;
    float                Mix_A_Px         [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    float                Mix_A_Py         [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    float                Mix_A_Pz         [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    float                Mix_A_EvtID      [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    float                Mix_A_TreID      [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    float                Mix_B_Px         [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    float                Mix_B_Py         [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    float                Mix_B_Pz         [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    float                Mix_B_EvtID      [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    float                Mix_B_TreID      [CentralityBinNum]   [yBinNum]  [PtBinNum] [500];
    int                  Mix_event_Num    [CentralityBinNum]   [yBinNum]  [PtBinNum];
    int                  Mix_A_Num        [CentralityBinNum]   [yBinNum]  [PtBinNum];
    int                  Mix_B_Num        [CentralityBinNum]   [yBinNum]  [PtBinNum];
    TH1D* H_Kstar                         [CentralityBinNum]   [yBinNum]  [PtBinNum];
    TH1D* H_Mix_Kstar                     [CentralityBinNum]   [yBinNum]  [PtBinNum];

    for (int i=0;i<CentralityBinNum;i++){
        for (int j=0;j<yBinNum;j++){
            for (int k=0;k<PtBinNum;k++){
                TString HistName1 = "H_";
                TString HistName2 = "Cen: [";
                HistName1 += i;HistName1 += "_";
                HistName2 += CentralityBin[i];HistName2 += "% , ";
                HistName2 += CentralityBin[i+1];HistName2 += "%], ";
                HistName1 += j;HistName1 += "_";
                HistName2 += yBin[j];HistName2 += " < y";HistName2 += B_PDG;HistName2 += " <  ";
                HistName2 += yBin[j+1];HistName2 += ", ";
                HistName1 += k;
                HistName2 += PtBin[k];HistName2 += " < Pt";HistName2 += B_PDG;HistName2 += " <  ";
                HistName2 += PtBin[k+1];
                TString HistName1s = HistName1;
                TString HistName2s = HistName2;
                HistName1s += "_S";
                HistName1 += "_M";
                HistName2 += ", Mix";
                H_Kstar[i][j][k] = new TH1D(HistName1s,HistName2s,500,0,10);
                H_Mix_Kstar[i][j][k] = new TH1D(HistName1,HistName2,500,0,10);
                Mix_event_Num[i][j][k] = 0;
                Mix_A_Num[i][j][k] = 0;
                Mix_B_Num[i][j][k] = 0;
            }
        }
    }

    for (int i=0;i<nentries;i++){
        // if (i > 15) {break;}
        hadronTree->GetEntry(i);
        if ((i+1)%200 == 0) cout<<"Calculating Event "<<(i+1)<<"/"<<nentries<<endl;
        // cout<<mult<<endl;
        // if(b>7){continue;}
        // cout<<"There OK"<<endl;
        A_Px.resize(0);B_Px.resize(0);
        A_Py.resize(0);B_Py.resize(0);
        A_Pz.resize(0);B_Pz.resize(0);
        A_EvtID.resize(0);B_EvtID.resize(0);
        A_TreID.resize(0);B_TreID.resize(0);

        for (int j=0;j<PDGMult;j++){
            if (PDG->at(j) == A_PDG) {
                A_Px.push_back(mix_px->at(j));
                A_Py.push_back(mix_py->at(j));
                A_Pz.push_back(mix_pz->at(j));
                A_EvtID.push_back(i);
                A_TreID.push_back(j);
            }
            if (PDG->at(j) == B_PDG) {
                B_Px.push_back(mix_px->at(j));
                B_Py.push_back(mix_py->at(j));
                B_Pz.push_back(mix_pz->at(j));
                B_EvtID.push_back(i);
                B_TreID.push_back(j);
            }
        }

        // cout<<"B_Px.size() = "<<B_Px.size()<<endl;
        if (B_Px.size() == 1){
            // cout<<"FOUND"<<endl;

            int CenIndex = -1;
            for (int k=0;k<CentralityBinNum;k++){
                if ((NchList[k] <= refMult) && (refMult < NchList[k+1])) {
                    CenIndex = k;
                    break;
                }
            }

            float tEnergy = pow(pow(B_Px[0],2) + pow(B_Py[0],2) + pow(B_Pz[0],2) + pow(massList(B_PDG),2),0.5);
            rap = 0.5*log((tEnergy+B_Pz[0])/(tEnergy-B_Pz[0]));
            // cout<<"px = "<<B_Px[0]<<" , py = "<<B_Py[0]<<" , pz = "<<B_Pz[0]<<endl;
            // cout<<"rap = "<<rap<<endl;
            int RapIndex = -1;
            for (int k=0;k<yBinNum;k++){
                if ((yBin[k] <= rap) && (rap < yBin[k+1])) {
                    RapIndex = k;
                    break;
                }
            }

            float B_Pt = pow(pow(B_Px[0],2) + pow(B_Py[0],2),0.5);
            int PtIndex = -1;
            for (int k=0;k<PtBinNum;k++){
                if ((PtBin[k] <= B_Pt) && (B_Pt < PtBin[k+1])) {
                    PtIndex = k;
                    break;
                }
            }

            // cout<<"CenIndex = "<<CenIndex<<" , "<<"RapIndex = "<<RapIndex<<" , "<<"PtIndex = "<<PtIndex<<endl;
            if ((CenIndex == -1) || (RapIndex == -1) || (PtIndex == -1)) {
                continue;
            }

            // Fill Single Event
            for (int j=0;j<B_Px.size();j++){
                TLorentzVector p1;
                p1.SetXYZM(B_Px[j],B_Py[j],B_Pz[j],massList(B_PDG));
                for (int k=0;k<A_Px.size();k++){
                    TLorentzVector p2,p3,p4 = p1;
                    p2.SetXYZM(A_Px[k],A_Py[k],A_Pz[k],massList(A_PDG));
                    p3 = p4 + p2;
                    p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                    kstar = 0.5 * (p4 - p2).Rho();
                    H_Kstar[CenIndex][RapIndex][PtIndex]->Fill(kstar);
                }
            }

            // Fill Mix Event
            // cout<<"Here is OK 1"<<endl;
            for (int j=0;j<B_Px.size();j++){
                Mix_B_Px[CenIndex][RapIndex][PtIndex][Mix_B_Num[CenIndex][RapIndex][PtIndex]] = B_Px[j];
                Mix_B_Py[CenIndex][RapIndex][PtIndex][Mix_B_Num[CenIndex][RapIndex][PtIndex]] = B_Py[j];
                Mix_B_Pz[CenIndex][RapIndex][PtIndex][Mix_B_Num[CenIndex][RapIndex][PtIndex]] = B_Pz[j];
                Mix_B_EvtID[CenIndex][RapIndex][PtIndex][Mix_B_Num[CenIndex][RapIndex][PtIndex]] = B_EvtID[j];
                Mix_B_TreID[CenIndex][RapIndex][PtIndex][Mix_B_Num[CenIndex][RapIndex][PtIndex]] = B_TreID[j];
                Mix_B_Num[CenIndex][RapIndex][PtIndex]++;
            }
            for (int j=0;j<A_Px.size();j++){
                Mix_A_Px[CenIndex][RapIndex][PtIndex][Mix_A_Num[CenIndex][RapIndex][PtIndex]] = A_Px[j];
                Mix_A_Py[CenIndex][RapIndex][PtIndex][Mix_A_Num[CenIndex][RapIndex][PtIndex]] = A_Py[j];
                Mix_A_Pz[CenIndex][RapIndex][PtIndex][Mix_A_Num[CenIndex][RapIndex][PtIndex]] = A_Pz[j];
                Mix_A_EvtID[CenIndex][RapIndex][PtIndex][Mix_A_Num[CenIndex][RapIndex][PtIndex]] = A_EvtID[j];
                Mix_A_TreID[CenIndex][RapIndex][PtIndex][Mix_A_Num[CenIndex][RapIndex][PtIndex]] = A_TreID[j];
                Mix_A_Num[CenIndex][RapIndex][PtIndex]++;
            }
            Mix_event_Num[CenIndex][RapIndex][PtIndex]++;
            // cout<<"Here is OK 2"<<endl;

            if (Mix_event_Num[CenIndex][RapIndex][PtIndex] == 10){
                for (int j=0;j<Mix_B_Num[CenIndex][RapIndex][PtIndex];j++){
                    TLorentzVector p1;
                    p1.SetXYZM(Mix_B_Px[CenIndex][RapIndex][PtIndex][j],Mix_B_Py[CenIndex][RapIndex][PtIndex][j],Mix_B_Pz[CenIndex][RapIndex][PtIndex][j],massList(B_PDG));
                    for (int k=0;k<Mix_A_Num[CenIndex][RapIndex][PtIndex];k++){
                        TLorentzVector p2,p3,p4 = p1;
                        p2.SetXYZM(Mix_A_Px[CenIndex][RapIndex][PtIndex][k],Mix_A_Py[CenIndex][RapIndex][PtIndex][k],Mix_A_Pz[CenIndex][RapIndex][PtIndex][k],massList(A_PDG));
                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        H_Mix_Kstar[CenIndex][RapIndex][PtIndex]->Fill(kstar);
                    }
                }
                Mix_event_Num[CenIndex][RapIndex][PtIndex] = 0;
                Mix_A_Num[CenIndex][RapIndex][PtIndex] = 0;
                Mix_B_Num[CenIndex][RapIndex][PtIndex] = 0;
            }
        }
        
    }

    TString OutputFileName = "Cor_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *file = new TFile(OutputFileName, "RECREATE");


    for (int i=0;i<CentralityBinNum;i++){
        for (int j=0;j<yBinNum;j++){
            for (int k=0;k<PtBinNum;k++){
                if (Mix_event_Num[i][j][k] != 0) cout<<"["<<i<<","<<j<<","<<k<<"] remain "<<Mix_event_Num[i][j][k]<<endl;
                H_Kstar[i][j][k]->Write();
                H_Mix_Kstar[i][j][k]->Write();
            }
        }
    }

    // Writing remaining pool

    int buffer_size = 5000000;
    int BPDGMult  ;
    int BCrefMult ;
    int BCgrefMult;
    int BevtID    ;
    int BrunID    ;
    int BTriggerID;
    std::vector<int> BPDG               ;BPDG            .resize(0);
    std::vector<float> Bpx              ;Bpx             .resize(0);
    std::vector<float> Bpy              ;Bpy             .resize(0);
    std::vector<float> Bpz              ;Bpz             .resize(0);
    std::vector<float> BQA_eta          ;BQA_eta         .resize(0);
    std::vector<float> BQA_dEdx         ;BQA_dEdx        .resize(0);
    std::vector<float> BQA_m2           ;BQA_m2          .resize(0);
    std::vector<float> BQA_DCA_V0_PV    ;BQA_DCA_V0_PV   .resize(0);
    std::vector<float> BQA_nSigmaProton ;BQA_nSigmaProton.resize(0);
    std::vector<float> BQA_nSigmaPion   ;BQA_nSigmaPion  .resize(0);
    std::vector<float> BQA_nSigmaKaon   ;BQA_nSigmaKaon  .resize(0);
    std::vector<float> BInvariantMass   ;BInvariantMass  .resize(0);
    std::vector<float> BQA_Decay_Length ;BQA_Decay_Length.resize(0);
    std::vector<float> BQA_Chi2         ;BQA_Chi2        .resize(0);
    BhadronTree = new TTree("hadronTree", "Tree_STAR");
    BhadronTree->Branch("PDGMult"            ,&BPDGMult             ,"PDGMult/I"                           );
    BhadronTree->Branch("refMult"            ,&BCrefMult            ,"refMult/I"                           );
    BhadronTree->Branch("grefMult"           ,&BCgrefMult           ,"grefMult/I"                          );
    BhadronTree->Branch("EventID"            ,&BevtID               ,"EventID/I"                           );
    BhadronTree->Branch("RunID"              ,&BrunID               ,"RunID/I"                             );
    BhadronTree->Branch("TriggerID"          ,&BTriggerID           ,"TriggerID/I"                         );
    BhadronTree->Branch("PDG"                ,&BPDG                 );
    BhadronTree->Branch("mix_px"             ,&Bpx                  );
    BhadronTree->Branch("mix_py"             ,&Bpy                  );
    BhadronTree->Branch("mix_pz"             ,&Bpz                  );
    BhadronTree->Branch("QA_eta"             ,&BQA_eta                 );

    // Used for PID QA
    BhadronTree->Branch("dEdx"               ,&BQA_dEdx              );
    BhadronTree->Branch("m2"                 ,&BQA_m2                );
    BhadronTree->Branch("dcatopv"            ,&BQA_DCA_V0_PV         );
    BhadronTree->Branch("nSigmaProton"       ,&BQA_nSigmaProton      );
    BhadronTree->Branch("nSigmaPion"         ,&BQA_nSigmaPion        );
    BhadronTree->Branch("nSigmaKaon"         ,&BQA_nSigmaKaon        );
    
    // Used for Reconstruction QA
    BhadronTree->Branch("InvariantMass"      ,&BInvariantMass        );
    BhadronTree->Branch("Decay_Length"       ,&BQA_Decay_Length      );
    BhadronTree->Branch("Chi2"               ,&BQA_Chi2              );

    std::vector<Int_t> Mix_EvtID;
    std::vector<std::vector<Int_t> > Mix_TreID;
    for (int i=0;i<CentralityBinNum;i++){
        for (int j=0;j<yBinNum;j++){
            for (int k=0;k<PtBinNum;k++){

                for (int m=0;m<Mix_B_Num[i][j][k];m++){
                    int nIndex = -1;
                    for (int n=0;n<Mix_EvtID.size();n++){
                        if (Mix_B_EvtID[i][j][k][m] == Mix_EvtID[n]){
                            nIndex = n;
                            break;
                        }
                    }
                    if (nIndex == -1){
                        Mix_EvtID.push_back(Mix_B_EvtID[i][j][k][m]);
                        std::vector<Int_t> Temp;
                        Mix_TreID.push_back(Temp);
                        nIndex = Mix_EvtID.size() - 1;
                    }
                    Mix_TreID[nIndex].push_back(Mix_B_TreID[i][j][k][m]);
                }
                for (int m=0;m<Mix_A_Num[i][j][k];m++){
                    int nIndex = -1;
                    for (int n=0;n<Mix_EvtID.size();n++){
                        if (Mix_A_EvtID[i][j][k][m] == Mix_EvtID[n]){
                            nIndex = n;
                            break;
                        }
                    }
                    if (nIndex == -1){
                        Mix_EvtID.push_back(Mix_A_EvtID[i][j][k][m]);
                        std::vector<Int_t> Temp;
                        Mix_TreID.push_back(Temp);
                        nIndex = Mix_EvtID.size() - 1;
                    }
                    Mix_TreID[nIndex].push_back(Mix_A_TreID[i][j][k][m]);
                }

            }
        }
    }
    for (int i=0;i<Mix_EvtID.size();i++){
        hadronTree->GetEntry(Mix_EvtID[i]);
        BPDGMult   = PDGMult  ;
        BCrefMult  = refMult  ;
        BCgrefMult = grefMult ;
        BevtID     = EventID  ;
        BrunID     = RunID    ;
        BTriggerID = TriggerID;
        for (int j=0;j<Mix_TreID[i].size();j++){
            BPDG            .emplace_back(PDG            [Mix_TreID[i][j]]);
            Bpx             .emplace_back(px             [Mix_TreID[i][j]]);
            Bpy             .emplace_back(py             [Mix_TreID[i][j]]);
            Bpz             .emplace_back(pz             [Mix_TreID[i][j]]);
            BQA_eta         .emplace_back(QA_eta         [Mix_TreID[i][j]]);
            BQA_dEdx        .emplace_back(QA_dEdx        [Mix_TreID[i][j]]);
            BQA_m2          .emplace_back(QA_m2          [Mix_TreID[i][j]]);
            BQA_DCA_V0_PV   .emplace_back(QA_DCA_V0_PV   [Mix_TreID[i][j]]);
            BQA_nSigmaProton.emplace_back(QA_nSigmaProton[Mix_TreID[i][j]]);
            BQA_nSigmaPion  .emplace_back(QA_nSigmaPion  [Mix_TreID[i][j]]);
            BQA_nSigmaKaon  .emplace_back(QA_nSigmaKaon  [Mix_TreID[i][j]]);
            BInvariantMass  .emplace_back(InvariantMass  [Mix_TreID[i][j]]);
            BQA_Decay_Length.emplace_back(QA_Decay_Length[Mix_TreID[i][j]]);
            BQA_Chi2        .emplace_back(QA_Chi2        [Mix_TreID[i][j]]);
        }
        BhadronTree->Fill();
    }

    BhadronTree->Write();
    file->Write();

    return;
}
