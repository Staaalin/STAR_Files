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
    std::vector<int> V_T;V_T.clear();
    for (int i=0;i<V.size();i++){
        V_T.push_back(V[i]);
    }
    V.clear();
    for (int i=0;i<V_T.size();i++){
        if (i == ID) continue;
        V.push_back(V_T[i]);
    }
    return;
}

std::vector<int> GetNchList(int CentralityList[] , int CentralityListSize)
{
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

void MixEvent_Vector(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,TString OutMidName,
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
    // used as array
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
    // used as value, ***[0] must be used
    std::vector<int>   Mix_event_Num      [CentralityBinNum]   [yBinNum]  [PtBinNum] [Pattern];
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
                    Mix_event_Num[i][j][k][l].push_back(0);
                }
            }
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
    for (int EventID = 0 ; EventID < nentries ; EventID++){
        cout<<"EventID = "<<EventID<<endl;
        hadronTree->GetEntry(EventID);
        if ((EventID+1)%50 == 0) {
            cout<<"Calculating Event "<<(EventID+1)<<"/"<<nentries<<endl;
            time(&time_now);
            int time_diff = (int)difftime(time_now, time_start);
            cout << time_diff/60 << "min " << time_diff%60 << "s: " << endl;
        }
        // cout<<mult<<endl;
        // if(b>7){continue;}
        // cout<<"There OK"<<endl;
        A_Px.clear();   B_Px.clear();
        A_Py.clear();   B_Py.clear();
        A_Pz.clear();   B_Pz.clear();
        A_EvtID.clear();B_EvtID.clear();
        A_TreID.clear();B_TreID.clear();
        A_ParID.clear();B_ParID.clear();
        A_Mass.clear(); B_Mass.clear();
        A_Kind.clear(); B_Kind.clear();

        for (int j=0;j<PDGMult;j++){
            if (PDG->at(j) == A_PDG) {
                if      (fabs(InvariantMass->at(j) - massList(A_PDG)) <= 3*massListSigma(A_PDG)) {A_Kind.push_back("Mid");}
                else if (fabs(InvariantMass->at(j) - massList(A_PDG)) <= 6*massListSigma(A_PDG)) {A_Kind.push_back("Sid");}
                else{continue;}
                A_Mass.push_back(InvariantMass->at(j));
                A_Px.push_back(mix_px->at(j));
                A_Py.push_back(mix_py->at(j));
                A_Pz.push_back(mix_pz->at(j));
                A_EvtID.push_back(i);
                A_TreID.push_back(j);
                std::vector<Int_t> Temp;Temp.clear();Temp.push_back(j);
                for (int k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                    Temp.push_back(ParentList->at(k));
                }
                A_ParID.push_back(Temp);
            }
            if (PDG->at(j) == B_PDG) {
                if      (fabs(InvariantMass->at(j) - massList(B_PDG)) <= 3*massListSigma(B_PDG)) {B_Kind.push_back("Mid");}
                else if (fabs(InvariantMass->at(j) - massList(B_PDG)) <= 6*massListSigma(B_PDG)) {B_Kind.push_back("Sid");}
                else{continue;}
                B_Mass.push_back(InvariantMass->at(j));
                B_Px.push_back(mix_px->at(j));
                B_Py.push_back(mix_py->at(j));
                B_Pz.push_back(mix_pz->at(j));
                B_EvtID.push_back(i);
                B_TreID.push_back(j);
                std::vector<Int_t> Temp;Temp.clear();Temp.push_back(j);
                for (int k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                    Temp.push_back(ParentList->at(k));
                }
                B_ParID.push_back(Temp);
            }
        }

        if ((A_Px.size() == 0) || (B_Px.size() == 0)) {cout<<"yes"<<endl;}

        for (int Bid = 0;Bid < B_Px.size();Bid++) {

            int CenIndex = -1;
            for (int k=0;k<CentralityBinNum;k++){
                // if ((NchList[k] <= refMult) && (refMult < NchList[k+1])) {
                if ((NchList[k] >= Nch) && (Nch > NchList[k+1])) {
                    CenIndex = k;
                    break;
                }
            }

            float tEnergy = pow(pow(B_Px[Bid],2) + pow(B_Py[Bid],2) + pow(B_Pz[Bid],2) + pow(massList(B_PDG),2),0.5);
            rap = 0.5*log((tEnergy+B_Pz[Bid])/(tEnergy-B_Pz[Bid]));
            // cout<<"px = "<<B_Px[0]<<" , py = "<<B_Py[0]<<" , pz = "<<B_Pz[0]<<endl;
            // cout<<"rap = "<<rap<<endl;
            int RapIndex = -1;
            for (int k=0;k<yBinNum;k++){
                if ((yBin[k] <= rap) && (rap < yBin[k+1])) {
                    RapIndex = k;
                    break;
                }
            }

            float B_Pt = pow(pow(B_Px[Bid],2) + pow(B_Py[Bid],2),0.5);
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

        }
        cout<<"2"<<endl;
    }

    return;
}
