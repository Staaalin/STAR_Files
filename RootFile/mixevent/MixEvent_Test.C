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

int CentralityBin[] = {0,10,20,40,80};// %
#define CentralityBinNum 4 // -1

float PtBin[] = {0 , 10.0}; // Pt
#define PtBinNum 1 // -1

float yBin[] = {-5.0 , 0.0 , 0.5 , 5.0}; // B_y
#define yBinNum 3 // -1

TString KindBin[] = {"Mid","Sid"}
#define KindNum 2
TString PatternBin[] = {"AMBM","AMBS","ASBM"};
#define Pattern 3 // 0:A middle B middle , 1:A middle B sideband , 2:A sideband B middle
// Pattern应当大于KindNum

void print(std::vector<int> Temp)
{
	cout<<"{";
    for (int i = 0;i<Temp.size();i++){
		cout<<" "<<Temp[i];
		if (i != (Temp.size() - 1)) cout<<" ,"; 
	}
	cout<<" }"<<endl;
    return ;
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
                Result = 100;
        }
    }
    return Result;
}

bool IfInVector(int Num , std::vector<int> V)
{
    for (int i=0;i<V.size();i++) {
        if (Num == V[i]){
            return true;
        }
    }
    return false;
}

bool IfCommonElement(std::vector<int> A , std::vector<int> B)
{
    for (int i=0;i<A.size();i++){
        for (int j=0;j<B.size();j++){
            if (A[i] == B[j]) return true;
        }
    }
    return false;
}

void DltElement(std::vector<int> &V , int ID)
{
    std::vector<int> V_T;V_T.resize(0);
    for (int i=0;i<V.size();i++){
        V_T.push_back(V[i]);
    }
    V.resize(0);
    for (int i=0;i<V_T.size();i++){
        if (i == ID) continue;
        V.push_back(V_T[i]);
    }
    return;
}

std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize)
{
    std::vector<int> Result;Result.resize(0);
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

void MixEvent_Test(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,TString OutMidName,
              int A_PDG,int B_PDG,int Mode = 0) // Mode = 0: PDGMult 为vector长度
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

    // TString FileName = "output_";
    // cout<<"Start Running"<<endl;
    // cout<<StartFileIndex<<endl;
    // cout<<EndFileIndex<<endl;


    double kstar, rap;
    TVector3 BetaTemp;

    int KaonpPID = 321,KaonmPID = -321,PionpPID = 211,PionmPID = -211,LambdaPID = 3122,LambdabPID = -3122,XiPID = 3312,XibPID = -3312,OmegaPID = 3334,OmegabPID = -3334;

    std::vector<int> NchList = GetNchList(CentralityBin , CentralityBinNum+1);     // centrality
    cout<<"NchList = ";
    print(NchList);
    cout<<" "<<endl;
    //                                        centrality          B_y        B_Pt
    std::vector<Float_t>              A_Px        ;
    std::vector<Float_t>              A_Py        ;
    std::vector<Float_t>              A_Pz        ;
    std::vector<Int_t>                A_EvtID     ;
    std::vector<Int_t>                A_TreID     ;
    std::vector<std::vector<int> >    A_ParID     ;
    std::vector<Float_t>              A_Mass      ;
    std::vector<TString>              A_Kind      ;
    std::vector<Float_t>              B_Px        ;
    std::vector<Float_t>              B_Py        ;
    std::vector<Float_t>              B_Pz        ;
    std::vector<Int_t>                B_EvtID     ;
    std::vector<Int_t>                B_TreID     ;
    std::vector<std::vector<int> >    B_ParID     ;
    std::vector<Float_t>              B_Mass      ;
    std::vector<TString>              B_Kind      ;
    std::vector<float> Mix_A_Px           [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    std::vector<float> Mix_A_Py           [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    std::vector<float> Mix_A_Pz           [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    std::vector<float> Mix_A_EvtID        [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    std::vector<float> Mix_A_TreID        [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    std::vector<float> Mix_B_Px           [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    std::vector<float> Mix_B_Py           [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    std::vector<float> Mix_B_Pz           [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    std::vector<float> Mix_B_EvtID        [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    std::vector<float> Mix_B_TreID        [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    int   Mix_event_Num                   [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    //
    TH1D* H_Kstar                         [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
    TH1D* H_Mix_Kstar                     [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];

    for (int i=0;i<CentralityBinNum;i++){
        for (int l=0;l<Pattern;l++){
            for (int k=0;k<PtBinNum;k++){
                for (int j=0;j<yBinNum;j++){
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
                    HistName1s += "_S_";
                    HistName1 += "_M_";
                    HistName1 += PatternBin[l];
                    HistName1s += PatternBin[l];
                    HistName2 += ", Mix, ";
                    HistName2 += PatternBin[l];
                    HistName2s += PatternBin[l];
                    H_Kstar[i][j][k][l] = new TH1D(HistName1s,HistName2s,500,0,10);
                    H_Mix_Kstar[i][j][k][l] = new TH1D(HistName1,HistName2,500,0,10);
                    Mix_event_Num[i][j][k][l] = 0;
                }
            }
        }
    }

    int ReadTreeIDMax;
    if (Mode == 0) {ReadTreeIDMax = 1;}
    else           {ReadTreeIDMax = Pattern;}
    for (int ReadTreeID = 0;ReadTreeID < ReadTreeIDMax;ReadTreeID++){
        //load data  
        TString TreeName;
        if (Mode == 0) {TreeName = "hadronTree";}
        else{
            if (ReadTreeID == 0) {TreeName = "AMBM/hadronTree";}
            if (ReadTreeID == 1) {TreeName = "AMBS/hadronTree";}
            if (ReadTreeID == 2) {TreeName = "ASBM/hadronTree";}
        }
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

        hadronTree->SetBranchAddress("PDGMult"  ,&PDGMult  );
        hadronTree->SetBranchAddress("refMult"  ,&refMult  );
        hadronTree->SetBranchAddress("grefMult" ,&grefMult );
        hadronTree->SetBranchAddress("EventID"  ,&EventID  );
        hadronTree->SetBranchAddress("RunID"    ,&RunID    );
        hadronTree->SetBranchAddress("TriggerID",&TriggerID);
        hadronTree->SetBranchAddress("Nch"      ,&Nch      );
        
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
        hadronTree->SetBranchAddress("ParentList"   ,&ParentList   ,&bParentList   );
        hadronTree->SetBranchAddress("ParentSta"    ,&ParentSta    ,&bParentSta    );
        hadronTree->SetBranchAddress("ParentEnd"    ,&ParentEnd    ,&bParentEnd    );

        const Int_t nentries=hadronTree->GetEntries();
        cout << "file number: " << nentries << endl;
        



        time_t time_start;
        time_t time_now;
        time(&time_start);
        for (int i=0;i<nentries;i++){
            int FoundAB = 0;
            // if (i > 15) {break;}
            hadronTree->GetEntry(i);
            if ((i+1)%100 == 0) {
                cout<<"Calculating Event "<<(i+1)<<"/"<<nentries<<endl;
			    time(&time_now);
			    int time_diff = (int)difftime(time_now, time_start);
                cout << time_diff/60 << "min " << time_diff%60 << "s: " << endl;
            }
            // cout<<mult<<endl;
            // if(b>7){continue;}
            // cout<<"There OK"<<endl;
            A_Px.resize(0);   B_Px.resize(0);
            A_Py.resize(0);   B_Py.resize(0);
            A_Pz.resize(0);   B_Pz.resize(0);
            A_EvtID.resize(0);B_EvtID.resize(0);
            A_TreID.resize(0);B_TreID.resize(0);
            A_ParID.resize(0);B_ParID.resize(0);
            A_Mass.resize(0); B_Mass.resize(0);
            A_Kind.resize(0); B_Kind.resize(0);

            //                               A              B
            TString KindSample[2][2] /*= {{"Mid","Sid"},{"Mid","Sid"}}*/;
            KindSample[0][0] = "Mid";KindSample[0][1] = "Sid";
            KindSample[1][0] = "Mid";KindSample[1][1] = "Sid";
            if (Mode != 0) {
                if (ReadTreeID == 0) {KindSample[0][1] = "Mid";KindSample[1][1] = "Mid";}
                if (ReadTreeID == 1) {KindSample[0][1] = "Mid";KindSample[1][0] = "Sid";}
                if (ReadTreeID == 2) {KindSample[0][0] = "Sid";KindSample[1][1] = "Mid";}
            }
            int LoopSize;
            LoopSize = PDGMult;
            for (int j=0;j<LoopSize;j++){
                if (PDG->at(j) == A_PDG) {
                    if      (fabs(InvariantMass->at(j) - massList(A_PDG)) <= 3*massListSigma(A_PDG)) {A_Kind.push_back("Mid");FoundAB++;}
                    else if (fabs(InvariantMass->at(j) - massList(A_PDG)) <= 6*massListSigma(A_PDG)) {A_Kind.push_back("Sid");FoundAB++;}
                    else{continue;}
                    bool IfStore = false;
                    for (int k = 0;k<2;k++) {
                        if (KindSample[0][k] == A_Kind[A_Kind.size()-1]) {IfStore = true;break;}
                    }
                    if (!IfStore) {A_Kind.resize(A_Kind.size()-1);continue;}
                    A_Mass.push_back(InvariantMass->at(j));
                    A_Px.push_back(mix_px->at(j));
                    A_Py.push_back(mix_py->at(j));
                    A_Pz.push_back(mix_pz->at(j));
                    A_EvtID.push_back(i);
                    A_TreID.push_back(j);
                    std::vector<Int_t> Temp;Temp.resize(0);Temp.push_back(j);
                    for (int k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                        Temp.push_back(ParentList->at(k));
                    }
                    A_ParID.push_back(Temp);
                }
                if (PDG->at(j) == B_PDG) {
                    if      (fabs(InvariantMass->at(j) - massList(B_PDG)) <= 3*massListSigma(B_PDG)) {B_Kind.push_back("Mid");FoundAB++;}
                    else if (fabs(InvariantMass->at(j) - massList(B_PDG)) <= 6*massListSigma(B_PDG)) {B_Kind.push_back("Sid");FoundAB++;}
                    else{continue;}
                    bool IfStore = false;
                    for (int k = 0;k<2;k++) {
                        if (KindSample[1][k] == B_Kind[B_Kind.size()-1]) {IfStore = true;break;}
                    }
                    if (!IfStore) {B_Kind.resize(B_Kind.size()-1);continue;}
                    B_Mass.push_back(InvariantMass->at(j));
                    B_Px.push_back(mix_px->at(j));
                    B_Py.push_back(mix_py->at(j));
                    B_Pz.push_back(mix_pz->at(j));
                    B_EvtID.push_back(i);
                    B_TreID.push_back(j);
                    std::vector<Int_t> Temp;Temp.resize(0);Temp.push_back(j);
                    for (int k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                        Temp.push_back(ParentList->at(k));
                    }
                    B_ParID.push_back(Temp);
                }
            }
            
            if ((A_Px.size() == 0) || (B_Px.size() == 0)) continue;

            int NumABCal = 0;// 这是实际上每个事件会进行相对动量计算的次数
            for (int j=0;j<B_Px.size();j++) {
                for (int k=0;k<A_Px.size();k++) {
                    if (IfCommonElement(A_ParID[k] , B_ParID[j])) continue;
                    NumABCal++;
                }
            }
            if (NumABCal == 0) continue;

            //                  mid   sid      used as binary bool
            int A_Pattern[] = {  0  ,  0  };
            int B_Pattern[] = {  0  ,  0  };

            // cout<<"B_Px.size() = "<<B_Px.size()<<endl;
            for (int B_index = 0;B_index < B_Px.size();B_index++){
                // cout<<"FOUND"<<endl;

                int CenIndex = -1;
                for (int k=0;k<CentralityBinNum;k++){
                    // if ((NchList[k] <= refMult) && (refMult < NchList[k+1])) {
                    if ((NchList[k] >= Nch) && (Nch > NchList[k+1])) {
                        CenIndex = k;
                        break;
                    }
                }

                float tEnergy = pow(pow(B_Px[B_index],2) + pow(B_Py[B_index],2) + pow(B_Pz[B_index],2) + pow(massList(B_PDG),2),0.5);
                rap = 0.5*log((tEnergy+B_Pz[B_index])/(tEnergy-B_Pz[B_index]));
                // cout<<"px = "<<B_Px[0]<<" , py = "<<B_Py[0]<<" , pz = "<<B_Pz[0]<<endl;
                // cout<<"rap = "<<rap<<endl;
                int RapIndex = -1;
                for (int k=0;k<yBinNum;k++){
                    if ((yBin[k] <= rap) && (rap < yBin[k+1])) {
                        RapIndex = k;
                        break;
                    }
                }

                float B_Pt = pow(pow(B_Px[B_index],2) + pow(B_Py[B_index],2),0.5);
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

                A_Pattern[0] = 0;A_Pattern[1] = 0;
                B_Pattern[0] = 0;B_Pattern[1] = 0;

                // Fill Single Event
                TLorentzVector p1;
                p1.SetXYZM(B_Px[B_index],B_Py[B_index],B_Pz[B_index],massList(B_PDG));
                for (int k=0;k<A_Px.size();k++){
                    if (IfCommonElement(A_ParID[k] , B_ParID[B_index])) continue;
                    TLorentzVector p2,p3,p4 = p1;
                    p2.SetXYZM(A_Px[k],A_Py[k],A_Pz[k],massList(A_PDG));
                    p3 = p4 + p2;
                    p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                    kstar = 0.5 * (p4 - p2).Rho();
                    if ((A_Kind[k] == "Mid") && (B_Kind[B_index] == "Mid")){
                        H_Kstar[CenIndex][RapIndex][PtIndex][0]->Fill(kstar);
                        A_Pattern[0] += 1;B_Pattern[0] += 1;
                    }
                    if ((A_Kind[k] == "Sid") && (B_Kind[B_index] == "Mid")){
                        H_Kstar[CenIndex][RapIndex][PtIndex][2]->Fill(kstar);
                        A_Pattern[1] += 1;B_Pattern[0] += 1;
                    }
                    if ((A_Kind[k] == "Mid") && (B_Kind[B_index] == "Sid")){
                        H_Kstar[CenIndex][RapIndex][PtIndex][1]->Fill(kstar);
                        A_Pattern[0] += 1;B_Pattern[1] += 1;
                    }
                }
                // cout<<"Nch = "<<Nch<<endl;

                // Fill Mix Event
                // cout<<"Here is OK 1"<<endl;
                for (int k=0;k<Pattern;k++){
                    TString A_Sample , B_Sample;
                    if ( (k == 0) ){
                        if (!((A_Pattern[0] != 0) && (B_Pattern[0] != 0))) continue;
                        if ((TreeName == "AMBS/hadronTree") || (TreeName == "ASBM/hadronTree")) continue;
                        A_Sample = "Mid";B_Sample = "Mid";
                    }
                    if ( (k == 1) ){
                        if (!((A_Pattern[0] != 0) && (B_Pattern[1] != 0))) continue;
                        if ((TreeName == "AMBM/hadronTree") || (TreeName == "ASBM/hadronTree")) continue;
                        A_Sample = "Mid";B_Sample = "Sid";
                    }
                    if ( (k == 2) ){
                        if (!((A_Pattern[1] != 0) && (B_Pattern[0] != 0))) continue;
                        if ((TreeName == "AMBS/hadronTree") || (TreeName == "AMBM/hadronTree")) continue;
                        A_Sample = "Sid";B_Sample = "Mid";
                    }
                    for (int j=0;j<B_Px.size();j++){
                        if (B_Kind[j] != B_Sample) continue;
                        Mix_B_Px[CenIndex][RapIndex][PtIndex][k]   .push_back(B_Px[j]   );
                        Mix_B_Py[CenIndex][RapIndex][PtIndex][k]   .push_back(B_Py[j]   );
                        Mix_B_Pz[CenIndex][RapIndex][PtIndex][k]   .push_back(B_Pz[j]   );
                        Mix_B_EvtID[CenIndex][RapIndex][PtIndex][k].push_back(B_EvtID[j]);
                        Mix_B_TreID[CenIndex][RapIndex][PtIndex][k].push_back(B_TreID[j]);
                    }
                    for (int j=0;j<A_Px.size();j++){
                        if (A_Kind[j] != A_Sample) continue;
                        Mix_A_Px[CenIndex][RapIndex][PtIndex][k]   .push_back(A_Px[j]   );
                        Mix_A_Py[CenIndex][RapIndex][PtIndex][k]   .push_back(A_Py[j]   );
                        Mix_A_Pz[CenIndex][RapIndex][PtIndex][k]   .push_back(A_Pz[j]   );
                        Mix_A_EvtID[CenIndex][RapIndex][PtIndex][k].push_back(A_EvtID[j]);
                        Mix_A_TreID[CenIndex][RapIndex][PtIndex][k].push_back(A_TreID[j]);
                    }
                    Mix_event_Num[CenIndex][RapIndex][PtIndex][k]++;
                    // cout<<"Here is OK 2"<<endl;

                    if (Mix_event_Num[CenIndex][RapIndex][PtIndex][k] == 10){
                        int Mix_B_Size = Mix_B_Px[CenIndex][RapIndex][PtIndex][k].size();
                        for (int j=0;j<Mix_B_Size;j++){
                            int Mix_A_Size = Mix_A_Px[CenIndex][RapIndex][PtIndex][k].size();
                            TLorentzVector p1;
                            p1.SetXYZM(Mix_B_Px[CenIndex][RapIndex][PtIndex][k][j],Mix_B_Py[CenIndex][RapIndex][PtIndex][k][j],Mix_B_Pz[CenIndex][RapIndex][PtIndex][k][j],massList(B_PDG));
                            for (int l=0;l<Mix_A_Size;l++){
                                TLorentzVector p2,p3,p4 = p1;
                                p2.SetXYZM(Mix_A_Px[CenIndex][RapIndex][PtIndex][k][l],Mix_A_Py[CenIndex][RapIndex][PtIndex][k][l],Mix_A_Pz[CenIndex][RapIndex][PtIndex][k][l],massList(A_PDG));
                                p3 = p4 + p2;
                                p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                                kstar = 0.5 * (p4 - p2).Rho();
                                H_Mix_Kstar[CenIndex][RapIndex][PtIndex][k]->Fill(kstar);
                            }
                        }
                        Mix_event_Num[CenIndex][RapIndex][PtIndex][k] = 0;
                        Mix_B_Px[CenIndex][RapIndex][PtIndex][k]   .resize(0);
                        Mix_B_Py[CenIndex][RapIndex][PtIndex][k]   .resize(0);
                        Mix_B_Pz[CenIndex][RapIndex][PtIndex][k]   .resize(0);
                        Mix_B_EvtID[CenIndex][RapIndex][PtIndex][k].resize(0);
                        Mix_B_TreID[CenIndex][RapIndex][PtIndex][k].resize(0);
                        Mix_A_Px[CenIndex][RapIndex][PtIndex][k]   .resize(0);
                        Mix_A_Py[CenIndex][RapIndex][PtIndex][k]   .resize(0);
                        Mix_A_Pz[CenIndex][RapIndex][PtIndex][k]   .resize(0);
                        Mix_A_EvtID[CenIndex][RapIndex][PtIndex][k].resize(0);
                        Mix_A_TreID[CenIndex][RapIndex][PtIndex][k].resize(0);
                    }
                }
            }
            // if (FoundAB > 0) {cout<<"________________________________________"<<endl;}
        }

        TString OutputFileName = OutMidName;
        OutputFileName += "H_";
        OutputFileName += OutputFileIndex;
        OutputFileName += ".root";
        TFile *fileA = new TFile(OutputFileName, "RECREATE");

        fileA->cd();
        for (int i=0;i<CentralityBinNum;i++){
            for (int l=0;l<Pattern;l++){
                for (int j=0;j<yBinNum;j++){
                    for (int k=0;k<PtBinNum;k++){
                        if (Mix_event_Num[i][j][k][l] != 0) {
                            cout<<"["<<i<<","<<j<<","<<k<<","<<l<<"] remain "<<Mix_event_Num[i][j][k][l]<<" events, "<<endl; 
                        }
                        // if ((H_Kstar[i][j][k][l]->Integral())>0) H_Kstar[i][j][k][l]->Write();
                        // if ((H_Mix_Kstar[i][j][k][l]->Integral())>0) H_Mix_Kstar[i][j][k][l]->Write();
                        if (Mode == 0) H_Kstar[i][j][k][l]->Write();
                        H_Mix_Kstar[i][j][k][l]->Write();
                    }
                }
                // if ((H_ABphi_Bphi    [i][l]->Integral())>0) H_ABphi_Bphi    [i][l]->Write();
                // if ((H_ABphi_By      [i][l]->Integral())>0) H_ABphi_By      [i][l]->Write();
                // if ((H_Mix_ABphi_Bphi[i][l]->Integral())>0) H_Mix_ABphi_Bphi[i][l]->Write();
                // if ((H_Mix_ABphi_By  [i][l]->Integral())>0) H_Mix_ABphi_By  [i][l]->Write();
                if (Mode == 0) H_ABphi_Bphi[i][l]->Write();
                if (Mode == 0) H_ABphi_By  [i][l]->Write();
                H_Mix_ABphi_Bphi[i][l]->Write();
                H_Mix_ABphi_By  [i][l]->Write();
            }
        }
        fileA.Close();

        // Writing remaining pool

        OutputFileName = OutMidName;
        OutputFileName += "T_";
        OutputFileName += OutputFileIndex;
        OutputFileName += ".root";
        TFile *fileB = new TFile(OutputFileName, "RECREATE");
        folder_AMBM = fileB->mkdir("AMBM");
        folder_AMBS = fileB->mkdir("AMBS");
        folder_ASBM = fileB->mkdir("ASBM");
        fileB->cd();
        if (Mode == 0){
            for (int WriteTreeIndex = 0;WriteTreeIndex < Pattern;WriteTreeIndex++) {
                if (WriteTreeIndex == 0) folder_AMBM->cd();
                if (WriteTreeIndex == 1) folder_AMBS->cd();
                if (WriteTreeIndex == 2) folder_ASBM->cd();
                int buffer_size = 5000000;
                int BPDGMult  ;
                int BCrefMult ;
                int BCgrefMult;
                int BevtID    ;
                int BrunID    ;
                int BTriggerID;
                int BNch      ;
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
                std::vector<int> BParentList        ;BParentList     .resize(0);
                std::vector<int> BParentSta         ;BParentSta      .resize(0);
                std::vector<int> BParentEnd         ;BParentEnd      .resize(0);
                BhadronTree = new TTree("hadronTree", "Tree_STAR");
                BhadronTree->Branch("PDGMult"            ,&BPDGMult             ,"PDGMult/I"                           );
                BhadronTree->Branch("refMult"            ,&BCrefMult            ,"refMult/I"                           );
                BhadronTree->Branch("grefMult"           ,&BCgrefMult           ,"grefMult/I"                          );
                BhadronTree->Branch("EventID"            ,&BevtID               ,"EventID/I"                           );
                BhadronTree->Branch("RunID"              ,&BrunID               ,"RunID/I"                             );
                BhadronTree->Branch("TriggerID"          ,&BTriggerID           ,"TriggerID/I"                         );
                BhadronTree->Branch("Nch"                ,&BNch                 ,"Nch/I"                               );
                BhadronTree->Branch("PDG"                ,&BPDG                 );
                BhadronTree->Branch("mix_px"             ,&Bpx                  );
                BhadronTree->Branch("mix_py"             ,&Bpy                  );
                BhadronTree->Branch("mix_pz"             ,&Bpz                  );
                BhadronTree->Branch("QA_eta"             ,&BQA_eta              );

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
                
                // Used for restore corralated information
                BhadronTree->Branch("ParentList"         ,&BParentList     );
                BhadronTree->Branch("ParentSta"          ,&BParentSta      );
                BhadronTree->Branch("ParentEnd"          ,&BParentEnd      );

                std::vector<Int_t> Mix_EvtID;
                for (int i=0;i<CentralityBinNum;i++){
                    for (int j=0;j<yBinNum;j++){
                        for (int k=0;k<PtBinNum;k++){
                            for (int m=0;m<Mix_B_Px[i][j][k][WriteTreeIndex].size();m++){
                                int nIndex = -1;
                                for (int n=0;n<Mix_EvtID.size();n++){
                                    if (Mix_B_EvtID[i][j][k][WriteTreeIndex][m] == Mix_EvtID[n]){
                                        nIndex = n;
                                        break;
                                    }
                                }
                                if (nIndex == -1){
                                    Mix_EvtID.push_back(Mix_B_EvtID[i][j][k][WriteTreeIndex][m]);
                                    nIndex = Mix_EvtID.size() - 1;
                                }
                            }
                            for (int m=0;m<Mix_A_Px[i][j][k][WriteTreeIndex].size();m++){
                                int nIndex = -1;
                                for (int n=0;n<Mix_EvtID.size();n++){
                                    if (Mix_A_EvtID[i][j][k][WriteTreeIndex][m] == Mix_EvtID[n]){
                                        nIndex = n;
                                        break;
                                    }
                                }
                                if (nIndex == -1){
                                    Mix_EvtID.push_back(Mix_A_EvtID[i][j][k][WriteTreeIndex][m]);
                                    nIndex = Mix_EvtID.size() - 1;
                                }
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
                    BNch       = Nch      ;
                    for (int j=0;j<PDGMult;j++){
                        BPDG            .push_back(PDG          ->at(j));
                        Bpx             .push_back(mix_px       ->at(j));
                        Bpy             .push_back(mix_py       ->at(j));
                        Bpz             .push_back(mix_pz       ->at(j));
                        BQA_eta         .push_back(QA_eta       ->at(j));
                        BQA_dEdx        .push_back(dEdx         ->at(j));
                        BQA_m2          .push_back(m2           ->at(j));
                        BQA_DCA_V0_PV   .push_back(dcatopv      ->at(j));
                        BQA_nSigmaProton.push_back(nSigmaProton ->at(j));
                        BQA_nSigmaPion  .push_back(nSigmaPion   ->at(j));
                        BQA_nSigmaKaon  .push_back(nSigmaKaon   ->at(j));
                        BInvariantMass  .push_back(InvariantMass->at(j));
                        BQA_Decay_Length.push_back(Decay_Length ->at(j));
                        BQA_Chi2        .push_back(Chi2         ->at(j));
                        BParentSta      .push_back(ParentSta    ->at(j));
                        BParentEnd      .push_back(ParentEnd    ->at(j));
                    }
                    for (int j=0;j<ParentList.size();j++){
                        BParentList     .push_back(ParentList   ->at(j));
                    }
                    BhadronTree->Fill();
                    BPDG            .resize(0);
                    Bpx             .resize(0);
                    Bpy             .resize(0);
                    Bpz             .resize(0);
                    BQA_eta         .resize(0);
                    BQA_dEdx        .resize(0);
                    BQA_m2          .resize(0);
                    BQA_DCA_V0_PV   .resize(0);
                    BQA_nSigmaProton.resize(0);
                    BQA_nSigmaPion  .resize(0);
                    BQA_nSigmaKaon  .resize(0);
                    BInvariantMass  .resize(0);
                    BQA_Decay_Length.resize(0);
                    BQA_Chi2        .resize(0);
                    BParentList     .resize(0);
                    BParentSta      .resize(0);
                    BParentEnd      .resize(0);
                }

                BhadronTree->Write();
            }
        }
        fileB->Write();
        fileB->Close();
    }

    return;
}
