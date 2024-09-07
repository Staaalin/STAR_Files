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
#include <ctime>
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

const int CentralityBin[] = {0,10,20,40,80};// %
const float PtBin[] = {0 , 10.0}; // Pt
const float yBin[] = {-5.0 , 0.0 , 0.5 , 5.0}; // B_y

const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1
const Int_t PtBinNum = sizeof(PtBin)/sizeof(PtBin[0]) - 1; // -1
const Int_t yBinNum = sizeof(yBin)/sizeof(yBin[0]) - 1; // -1

switch (CentralityBinNum)
{
case 1:
    #define DCentralityBinNum 1
    break;
case 2:
    #define DCentralityBinNum 2
    break;
case 3:
    #define DCentralityBinNum 3
    break;
case 4:
    #define DCentralityBinNum 4
    break;
case 5:
    #define DCentralityBinNum 5
    break;
case 6:
    #define DCentralityBinNum 6
    break;
case 7:
    #define DCentralityBinNum 7
    break;
case 8:
    #define DCentralityBinNum 8
    break;
case 9:
    #define DCentralityBinNum 9
    break;
case 10:
    #define DCentralityBinNum 10
    break;
default:
    break;
}

switch (PtBinNum)
{
case 1:
    #define DPtBinNum 1
    break;
case 2:
    #define DPtBinNum 2
    break;
case 3:
    #define DPtBinNum 3
    break;
case 4:
    #define DPtBinNum 4
    break;
case 5:
    #define DPtBinNum 5
    break;
case 6:
    #define DPtBinNum 6
    break;
case 7:
    #define DPtBinNum 7
    break;
case 8:
    #define DPtBinNum 8
    break;
case 9:
    #define DPtBinNum 9
    break;
case 10:
    #define DPtBinNum 10
    break;
default:
    break;
}

switch (yBinNum)
{
case 1:
    #define DyBinNum 1
    break;
case 2:
    #define DyBinNum 2
    break;
case 3:
    #define DyBinNum 3
    break;
case 4:
    #define DyBinNum 4
    break;
case 5:
    #define DyBinNum 5
    break;
case 6:
    #define DyBinNum 6
    break;
case 7:
    #define DyBinNum 7
    break;
case 8:
    #define DyBinNum 8
    break;
case 9:
    #define DyBinNum 9
    break;
case 10:
    #define DyBinNum 10
    break;
default:
    break;
}

TString KindBin[] = {"Mid","Sid"}
#define KindNum 2
TString PatternBin[] = {"AMBM","AMBS","ASBM"};
#define Pattern 3 // 0:A middle B middle , 1:A middle B sideband , 2:A sideband B middle
// Pattern应当大于KindNum

int HowMuchEventMixing = 10;

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
    TLorentzVector p1 , p2 , p3 , p4 , p5;

    std::vector<int> NchList = GetNchList(CentralityBin , CentralityBinNum+1);     // centrality
    cout<<"NchList = ";
    print(NchList);
    cout<<" "<<endl;
    //                                        centrality          B_y        B_Pt
    std::vector<Float_t>                  A_Px           ;
    std::vector<Float_t>                  A_Py           ;
    std::vector<Float_t>                  A_Pz           ;
    std::vector<Int_t>                    A_TreID        ;
    std::vector<std::vector<int> >        A_ParID        ;
    std::vector<int>                      A_Kind         ;
    std::vector<Float_t>                  B_Px           ;
    std::vector<Float_t>                  B_Py           ;
    std::vector<Float_t>                  B_Pz           ;
    std::vector<Int_t>                    B_TreID        ;
    std::vector<std::vector<int> >        B_ParID        ;
    std::vector<int>                      B_Kind         ;
    // used as array
    std::vector<float> Mix_A_Px           [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<float> Mix_A_Py           [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<float> Mix_A_Pz           [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<int>   Mix_A_TreID        [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<int>   Mix_A_EvtID        [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<int>   Mix_A_ID                                [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<float> Mix_B_Px           [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<float> Mix_B_Py           [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<float> Mix_B_Pz           [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<int>   Mix_B_TreID        [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<int>   Mix_B_EvtID        [CentralityBinNum]   [yBinNum]  [PtBinNum]  [2] [2] ;
    std::vector<int>   Mix_B_ID                                [yBinNum]  [PtBinNum]  [2] [2] ;
    int                Mix_event_Num      [DCentralityBinNum]  [DyBinNum] [DPtBinNum] [2] [2] ;
    //
    TH1D* H_Kstar                         [DCentralityBinNum]  [DyBinNum] [DPtBinNum] [2] [2] ;
    TH1D* H_Mix_Kstar                     [DCentralityBinNum]  [DyBinNum] [DPtBinNum] [2] [2] ;
    int EventPatternMatch                 [DCentralityBinNum]  [DyBinNum] [DPtBinNum] [2] [2] ;
    // Used for testing
    int TestSum = 0;

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
                    if (l == 0) { // AMBM
                        H_Kstar            [i][j][k][0][0] = new TH1D(HistName1s,HistName2s,500,0,10);
                        H_Mix_Kstar        [i][j][k][0][0] = new TH1D(HistName1,HistName2,500,0,10);
                        Mix_event_Num      [i][j][k][0][0] = 0;
                        if(i==0) Mix_A_ID     [j][k][0][0].push_back(0);
                    }
                    if (l == 1) { // AMBS
                        H_Kstar            [i][j][k][0][1] = new TH1D(HistName1s,HistName2s,500,0,10);
                        H_Mix_Kstar        [i][j][k][0][1] = new TH1D(HistName1,HistName2,500,0,10);
                        Mix_event_Num      [i][j][k][0][1] = 0;
                        if(i==0) Mix_A_ID     [j][k][0][1].push_back(0);
                    }
                    if (l == 2) { // ASBM
                        H_Kstar            [i][j][k][1][0] = new TH1D(HistName1s,HistName2s,500,0,10);
                        H_Mix_Kstar        [i][j][k][1][0] = new TH1D(HistName1,HistName2,500,0,10);
                        Mix_event_Num      [i][j][k][1][0] = 0;
                        if(i==0) Mix_A_ID     [j][k][1][0].push_back(0);
                    }
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
    // hadronTree->SetBranchAddress("refMult"  ,&refMult  );
    // hadronTree->SetBranchAddress("grefMult" ,&grefMult );
    hadronTree->SetBranchAddress("EventID"  ,&EventID  );
    // hadronTree->SetBranchAddress("RunID"    ,&RunID    );
    hadronTree->SetBranchAddress("TriggerID",&TriggerID);
    hadronTree->SetBranchAddress("Nch"      ,&Nch      );
    
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
    
    float BMass = massList(B_PDG)           , AMass = massList(A_PDG);
    float BMassSigma = massListSigma(B_PDG) , AMassSigma = massListSigma(A_PDG);

    clock_t Tstart = clock();
    for (int EntriesID = 0 ; EntriesID < nentries ; EntriesID++){
        hadronTree->GetEntry(EntriesID);
        if ((EntriesID+1)%50 == 0) {
            long long microseconds = (clock() - Tstart)/10000;
            std::cout << "Microseconds: " << microseconds << "  ";
            cout << "Test/Events = " << 1.0*TestSum/50 << "  ";
            cout<<"Calculating Event "<<(EntriesID+1)<<"/"<<nentries<<endl;
            Tstart = clock();
            TestSum = 0;
        }

        A_Px   .resize(0);   B_Px.resize(0);
        A_Py   .resize(0);   B_Py.resize(0);
        A_Pz   .resize(0);   B_Pz.resize(0);
        A_TreID.resize(0);B_TreID.resize(0);
        A_ParID.resize(0);B_ParID.resize(0);
        A_Kind .resize(0); B_Kind.resize(0);

        for (int j=0;j<PDGMult;j++){
            if (PDG->at(j) == A_PDG) {
                if      (fabs(InvariantMass->at(j) - AMass) <= 3*AMassSigma) {A_Kind.push_back(0);}
                else if (fabs(InvariantMass->at(j) - AMass) <= 6*AMassSigma) {A_Kind.push_back(1);}
                else{continue;}
                A_Px.push_back(mix_px->at(j));
                A_Py.push_back(mix_py->at(j));
                A_Pz.push_back(mix_pz->at(j));
                A_TreID.push_back(j);
                std::vector<int> Temp;Temp.clear();Temp.push_back(j);
                for (int k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                    Temp.push_back(ParentList->at(k));
                }
                A_ParID.push_back(Temp);
            }
            if (PDG->at(j) == B_PDG) {
                if      (fabs(InvariantMass->at(j) - BMass) <= 3*BMassSigma) {B_Kind.push_back(0);}
                else if (fabs(InvariantMass->at(j) - BMass) <= 6*BMassSigma) {B_Kind.push_back(1);}
                else{continue;}
                B_Px.push_back(mix_px->at(j));
                B_Py.push_back(mix_py->at(j));
                B_Pz.push_back(mix_pz->at(j));
                B_TreID.push_back(j);
                std::vector<int> Temp;Temp.clear();Temp.push_back(j);
                for (int k=ParentSta->at(j);k<=ParentEnd->at(j);k++){
                    Temp.push_back(ParentList->at(k));
                }
                B_ParID.push_back(Temp);
            }
        }

        if ((A_Px.size() == 0) || (B_Px.size() == 0)) {continue;}
        
        // Event Index
        int CenIndex = -1;
        for (int k=0;k<CentralityBinNum;k++){
            // if ((NchList[k] <= refMult) && (refMult < NchList[k+1])) {
            if ((NchList[k] >= Nch) && (Nch > NchList[k+1])) {
                CenIndex = k;
                break;
            }
        }
        if (CenIndex == -1) continue;

        for (int i = 0;i < CentralityBinNum;i++) {
            for (int j = 0;j < yBinNum;j++) {
                for (int k = 0;k < PtBinNum;k++) {
                    for (int A_Kid = 0;A_Kid < 2;A_Kid++) {
                        for (int B_Kid = 0;B_Kid < 2;B_Kid++) {
                            EventPatternMatch[i]  [j] [k][A_Kid][B_Kid] = 0;
                        }
                    }
                }
            }
        }
        
        for (int i = 0;i < yBinNum;i++) {
            for (int j = 0;j < PtBinNum;j++) {
                Mix_A_ID[i][j] [0][0].resize(0);
                Mix_A_ID[i][j] [1][0].resize(0);
                Mix_A_ID[i][j] [0][1].resize(0);
                Mix_A_ID[i][j] [1][1].resize(0);
            }
        }

        for (int Bid = 0;Bid < B_Px.size();Bid++) {

            float BPx = B_Px[Bid] , BPy = B_Py[Bid] , BPz = B_Pz[Bid];

            // B Index
            float tEnergy = pow(pow(BPx,2) + pow(BPy,2) + pow(BPz,2) + pow(BMass,2),0.5);
            rap = 0.5*log((tEnergy+BPz)/(tEnergy-BPz));
            int RapIndex = -1;
            for (int k=0;k<yBinNum;k++){
                if ((yBin[k] <= rap) && (rap < yBin[k+1])) {
                    RapIndex = k;
                    break;
                }
            }

            float B_Pt = pow(pow(BPx,2) + pow(BPy,2),0.5);
            int PtIndex = -1;
            for (int k=0;k<PtBinNum;k++){
                if ((PtBin[k] <= B_Pt) && (B_Pt < PtBin[k+1])) {
                    PtIndex = k;
                    break;
                }
            }

            if ((RapIndex == -1) || (PtIndex == -1)) {
                continue;
            }

            p2.SetXYZM(BPx,BPy,BPz,BMass);
            int B_Kid = B_Kind[Bid];
            bool IfBFilled = false;
            for (int Aid = 0;Aid < A_Px.size();Aid++) {
                if (IfCommonElement(A_ParID[Aid] , B_ParID[Bid])) continue;
                int A_Kid = A_Kind[Aid];
                p3 = p2;
                p1.SetXYZM(A_Px[Aid],A_Py[Aid],A_Pz[Aid],AMass);
                p4 = p1 + p2;
                // p3.Boost(-p4.BoostVector());p1.Boost(-p4.BoostVector());
                // float TTT = 0.5 * (p3 - p1).Rho();
                H_Kstar[CenIndex][RapIndex][PtIndex][A_Kid][B_Kid]->Fill(0.5 * (p3 - p1).Rho());

                TestSum++;
                // bool IfRecord = true;
                // for (int Cid = 0;Cid < Mix_A_ID[RapIndex][PtIndex] [A_Kid][B_Kid].size();Cid++) {
                //     if (Mix_A_ID[RapIndex][PtIndex] [A_Kid][B_Kid][Cid] == Aid) {
                //         IfRecord = false;
                //         break;
                //     }
                // }
                // if (IfRecord) {
                //     Mix_A_Px   [CenIndex][RapIndex][PtIndex] [A_Kid][B_Kid].push_back(A_Px[Aid]);
                //     Mix_A_Py   [CenIndex][RapIndex][PtIndex] [A_Kid][B_Kid].push_back(A_Py[Aid]);
                //     Mix_A_Pz   [CenIndex][RapIndex][PtIndex] [A_Kid][B_Kid].push_back(A_Pz[Aid]);
                //     Mix_A_EvtID[CenIndex][RapIndex][PtIndex] [A_Kid][B_Kid].push_back(EventID);
                //     EventPatternMatch[CenIndex][RapIndex][PtIndex][A_Kid][B_Kid]++;
                //     Mix_A_ID[RapIndex][PtIndex] [A_Kid][B_Kid].push_back(Aid);
                // }
                // if (!IfBFilled) {
                //     Mix_B_Px   [CenIndex][RapIndex][PtIndex] [A_Kid][B_Kid].push_back(B_Px[Bid]);
                //     Mix_B_Py   [CenIndex][RapIndex][PtIndex] [A_Kid][B_Kid].push_back(B_Py[Bid]);
                //     Mix_B_Pz   [CenIndex][RapIndex][PtIndex] [A_Kid][B_Kid].push_back(B_Pz[Bid]);
                //     Mix_B_EvtID[CenIndex][RapIndex][PtIndex] [A_Kid][B_Kid].push_back(EventID);
                //     EventPatternMatch[CenIndex][RapIndex][PtIndex][A_Kid][B_Kid]++;
                //     IfBFilled = true;
                // }
            }
        }

        // for (int i = 0;i < yBinNum;i++) {
        //     for (int j = 0;j < PtBinNum;j++) {
        //         for (int Aid = 0;Aid < 2;Aid++) {
        //             for (int Bid = 0;Bid < 2;Bid++) {
        //                 if (EventPatternMatch[CenIndex][i][j][Aid][Bid] != 0) {
        //                     Mix_event_Num[CenIndex][i][j][Aid][Bid]++;
        //                     if (Mix_event_Num[CenIndex][i][j][Aid][Bid] == HowMuchEventMixing+1) {
        //                         int Mix_A_Size = Mix_A_Px[CenIndex][i][j][Aid][Bid].size();
        //                         int Mix_B_Size = Mix_B_Px[CenIndex][i][j][Aid][Bid].size();
        //                         for (int Aindex = 0;Aindex < Mix_A_Size;Aindex++) {
        //                             int A_EID = Mix_A_EvtID[CenIndex][i][j][Aid][Bid][Aindex];
        //                             p2.SetXYZM(Mix_A_Px[CenIndex][i][j][Aid][Bid][Aindex],Mix_A_Py[CenIndex][i][j][Aid][Bid][Aindex],Mix_A_Pz[CenIndex][i][j][Aid][Bid][Aindex],AMass);
        //                             for (int Bindex = 0;Bindex < Mix_B_Size;Bindex++) {
        //                                 if (A_EID == Mix_B_EvtID[CenIndex][i][j][Aid][Bid][Bindex]) continue;
        //                                 p3 = p2;
        //                                 p1.SetXYZM(Mix_B_Px[CenIndex][i][j][Aid][Bid][Bindex],Mix_B_Py[CenIndex][i][j][Aid][Bid][Bindex],Mix_B_Pz[CenIndex][i][j][Aid][Bid][Bindex],BMass);
        //                                 p4 = p1 + p2;
        //                                 p3.Boost(-p4.BoostVector());p1.Boost(-p4.BoostVector());
        //                                 H_Mix_Kstar[CenIndex][i][j][Aid][Bid]->Fill(0.5 * (p3 - p1).Rho());
        //                             }
        //                         }
        //                         Mix_event_Num[CenIndex][i][j][Aid][Bid] = 0;
        //                         Mix_A_Px[CenIndex][i][j][Aid][Bid].clear();
        //                         Mix_B_Px[CenIndex][i][j][Aid][Bid].clear();
        //                         Mix_A_Py[CenIndex][i][j][Aid][Bid].clear();
        //                         Mix_B_Py[CenIndex][i][j][Aid][Bid].clear();
        //                         Mix_A_Pz[CenIndex][i][j][Aid][Bid].clear();
        //                         Mix_B_Pz[CenIndex][i][j][Aid][Bid].clear();
        //                         Mix_A_EvtID[CenIndex][i][j][Aid][Bid].clear();
        //                         Mix_B_EvtID[CenIndex][i][j][Aid][Bid].clear();
        //                     }
        //                 }
        //             }
        //         }
        //     }

        // }


    }

    return;
}
