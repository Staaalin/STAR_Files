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
// #include "Math/Vector4D.h"
// #include <Math/PtEtaPhiM4D.h>
// #include <Math/Boost.h>
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

const int CentralityBin[] = {0 , 10 , 20 , 40 , 80};// %
const float yBin[] = {-0.5 , -0.4 , -0.3 , -0.2 , -0.1 , 0.0 ,  0.1 , 0.2 , 0.3 , 0.4 , 0.5}; // B_y
int FeedDown[] = { 3334 , 3312 , 1003314 };

const Int_t CentralityBinNum = sizeof(CentralityBin)/sizeof(CentralityBin[0]) - 1; // -1
const Int_t yBinNum = sizeof(yBin)/sizeof(yBin[0]) - 1; // -1
const Int_t FeedDownNum = sizeof(FeedDown)/sizeof(FeedDown[0]);

TString KindBin[] = {"Mid","Sid"}
#define KindNum 2
TString PatternBin[] = {"AMBM","AMBS","ASBM"};
#define Pattern 3 // 0:A middle B middle , 1:A middle B sideband , 2:A sideband B middle
// Pattern应当大于KindNum

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

bool IfInVector(int Num , std::vector<int> V)
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

void dNdy(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,TString OutMidName,int A_PDG,int B_PDG)
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

    std::vector<Float_t>                  A_Px           ;
    std::vector<Float_t>                  A_Py           ;
    std::vector<Float_t>                  A_Pz           ;
    std::vector<Int_t>                    A_TreID        ;
    std::vector<std::vector<int> >        A_ParID        ;
    std::vector<int>                      A_Kind         ;
    std::vector<Float_t>                  A_Rap          ;
    std::vector<Int_t>                    A_IfRecord     ;
    std::vector<Float_t>                  B_Px           ;
    std::vector<Float_t>                  B_Py           ;
    std::vector<Float_t>                  B_Pz           ;
    std::vector<Int_t>                    B_TreID        ;
    std::vector<std::vector<int> >        B_ParID        ;
    std::vector<int>                      B_Kind         ;
    std::vector<Float_t>                  B_Rap          ;
    std::vector<Int_t>                    B_IfRecord     ;
    std::vector<std::vector<int> >        C_ParID        ; // 用于存储Residal Effect

    int StarBinNum = 100 , StarSta = -1.5 , StarEnd = 1.5;
    StarBinNum = (StarEnd - StarSta)/0.05;

    TH1D* H_All           [20];
    TH1D* H_All_WithNetB  [20];
    TH1D* H_Brap          [20][20];

    for (int i=0;i<CentralityBinNum;i++) {
        TString HistName1All = "H_All_";
        TString HistName2All = "Cen: [ ";
        HistName1All += i;
        HistName2All += CentralityBin[i];HistName2All += " , ";HistName2All += CentralityBin[i+1];HistName2All += " ]";
        H_All[i] = new TH1D(HistName1All,HistName2All,StarBinNum,StarSta,StarEnd);
        HistName1All = "H_All_WithNetB_";
        HistName2All = "Cen: [ ";
        HistName1All += i;
        HistName2All += CentralityBin[i];HistName2All += " , ";HistName2All += CentralityBin[i+1];HistName2All += " ]";
        H_All_WithNetB[i] = new TH1D(HistName1All,HistName2All,StarBinNum,StarSta,StarEnd);
        for (int j=0;j<yBinNum;j++){
            TString HistName1Brap = "H_Brap_";
            TString HistName2Brap = "Cen: [ ";
            HistName1Brap += i;HistName1Brap += "_";
            HistName2Brap += CentralityBin[i];HistName2Brap += " , ";HistName2Brap += CentralityBin[i+1];HistName2Brap += " ] , ";
            HistName1Brap += j;
            HistName2Brap += "y";HistName2Brap += B_PDG;" = [ ";HistName2Brap += yBin[j];HistName2Brap += " , ";HistName2Brap += yBin[i+1];HistName2Brap += " ]";
            H_Brap[i][j] = new TH1D(HistName1Brap,HistName2Brap,StarBinNum,StarSta,StarEnd);
        }
    }

    float tEnergy , APx , APy , APz , BPx , BPy , BPz;
    int A_Kid , B_Kid , Mix_A_Size , Mix_B_Size , A_EID , AidN , BidN;
    int B_Sum , B_Sid;
    int CenIndex , RapIndex;
    std::vector<int> Temp;
    std::vector<float> CMass , CMassSigma;
    bool IfRecord = true , IfRecordThisEvent = true;
    float BMass = massList(B_PDG)           , AMass = massList(A_PDG);
    float BMassSigma = massListSigma(B_PDG) , AMassSigma = massListSigma(A_PDG);

    for (int i = 0;i < FeedDownNum;i++){
        CMass.push_back(massList(FeedDown[i]));
        CMassSigma.push_back(massListSigma(FeedDown[i]));
        if (abs(FeedDown[i]) == A_PDG) {
            FeedDown[i] = 0;
        }
        if (abs(FeedDown[i]) == B_PDG) {
            FeedDown[i] = 0;
        }
    }

    std::vector<int> NchList = GetNchList(CentralityBin , CentralityBinNum+1);     // centrality
    cout<<"NchList = ";
    print(NchList);
    cout<<" "<<endl;

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
    IfRecordThisEvent = true;
    for (int EntriesID = 0 ; EntriesID < nentries ; EntriesID++){
        hadronTree->GetEntry(EntriesID);
        if ((EntriesID+1)%200 == 0) {
            time(&time_now);
            int time_diff = (int)difftime(time_now, time_start);
            cout << time_diff/60 << "min " << time_diff%60 << "s: ";
            long long microseconds = (clock() - Tstart)/10000;
            std::cout << "Microseconds: " << microseconds << "  ";
            cout << "Test/Events = " << 1.0*TestSum/50 << "  ";
            cout<<"Calculating Event "<<(EntriesID+1)<<"/"<<nentries<<endl;
            Tstart = clock();
        }

        A_ParID.resize(0);B_ParID.resize(0);
        A_Rap  .resize(0);  B_Rap.resize(0);
        A_IfRecord.resize(0);B_IfRecord.resize(0);
        C_ParID.resize(0);

        for (int Id = 0; Id < PDGMult; Id++)
        {
            if (PDG->at(Id) == A_PDG) {
                if (fabs(InvariantMass->at(Id) - AMass) > 3*AMassSigma) {continue;}
                // if ((fabs(InvariantMass->at(Id) - AMass) <= 3*AMassSigma) || (fabs(InvariantMass->at(Id) - AMass) > 6*AMassSigma)) {continue;}
                A_Px.push_back(mix_px->at(Id));
                A_Py.push_back(mix_py->at(Id));
                A_Pz.push_back(mix_pz->at(Id));
                A_IfRecord.push_back(1);
                Temp.clear();Temp.push_back(Id);
                for (int k=ParentSta->at(Id);k<=ParentEnd->at(Id);k++){
                    Temp.push_back(ParentList->at(k));
                }
                A_ParID.push_back(Temp);
                tEnergy = pow(pow(mix_px->at(Id),2) + pow(mix_py->at(Id),2) + pow(mix_pz->at(Id),2) + pow(AMass,2),0.5);
                A_Rap.push_back(0.5*log((tEnergy+mix_pz->at(Id))/(tEnergy-mix_pz->at(Id))));
            }
            if (PDG->at(Id) == B_PDG) {
                if (fabs(InvariantMass->at(Id) - BMass) > 3*BMassSigma) {continue;}
                // if ((fabs(InvariantMass->at(Id) - BMass) <= 3*BMassSigma) || (fabs(InvariantMass->at(Id) - BMass) > 6*BMassSigma)) {continue;}
                B_Px.push_back(mix_px->at(Id));
                B_Py.push_back(mix_py->at(Id));
                B_Pz.push_back(mix_pz->at(Id));
                B_IfRecord.push_back(1);
                Temp.clear();Temp.push_back(Id);
                for (int k=ParentSta->at(Id);k<=ParentEnd->at(Id);k++){
                    Temp.push_back(ParentList->at(k));
                }
                B_ParID.push_back(Temp);
                tEnergy = pow(pow(mix_px->at(Id),2) + pow(mix_py->at(Id),2) + pow(mix_pz->at(Id),2) + pow(BMass,2),0.5);
                B_Rap.push_back(0.5*log((tEnergy+mix_pz->at(Id))/(tEnergy-mix_pz->at(Id))));
            }
            if (PDG->at(Id) == -1*B_PDG) {
                IfRecordThisEvent = false;
                break;
            }
            for (int j = 0;j < FeedDownNum;j++) {
                if (abs(PDG->at(Id)) == FeedDown[j]) {
                    if ((fabs(InvariantMass->at(Id) - CMass.at(j)) > 3*CMassSigma.at(j))) continue;
                    Temp.clear();Temp.push_back(Id);
                    for (int k=ParentSta->at(Id);k<=ParentEnd->at(Id);k++){
                        Temp.push_back(ParentList->at(k));
                    }
                    C_ParID.push_back(Temp);
                }
            }
        }

        if (!IfRecordThisEvent) continue;

        // 删除 FeedDown 等血缘关系的粒子
        for (int i = 0;i < A_Rap.size();i++) {
            for (int j = 0;j < B_Rap.size();j++) {
                if (IfCommonElement(A_ParID.at(i) , B_ParID.at(j))) {
                    A_IfRecord.at(i) = 0;
                }
            }
        }
        for (int i = 0;i < C_ParID.size();i++) {
            for (int j = 0;j < A_Rap.size();j++) {
                if (IfCommonElement(C_ParID.at(i) , A_ParID.at(j))) {
                    A_IfRecord.at(j) = 0;
                }
            }
            for (int j = 0;j < B_Rap.size();j++) {
                if (IfCommonElement(C_ParID.at(i) , B_ParID.at(j))) {
                    B_IfRecord.at(j) = 0;
                }
            }
        }

        //
        B_Sum = 0;B_Sid = -1;
        for (int i = 0;i < B_Rap.size();i++) {
            if (B_IfRecord.at(i) != 0) {
                B_Sid = i;
                B_Sum++;
            }
        }
        if (!(B_Sum != 1)) continue;

        // Event Index
        CenIndex = -1;
        for (int k=0;k<CentralityBinNum;k++){
            // if ((NchList.at(k) <= refMult) && (refMult < NchList.at(k+1))) {
            if ((NchList.at(k) >= Nch) && (Nch > NchList.at(k+1))) {
                CenIndex = k;
                break;
            }
        }
        if (CenIndex == -1) continue;

        for (int i=0;i<A_Rap.size();i++) {
            H_All[CenIndex]->Fill(A_Rap.at(j));
        }

        if (!IfRecordThisEvent) continue;

        for (int i=0;i<A_Rap.size();i++) {
            H_All_WithNetB[CenIndex]->Fill(A_Rap.at(j));
        }

        // B_Rap index
        RapIndex = -1
        for (int k=0;k<yBinNum;k++){
            if ((yBin[k] <= B_Rap.at(B_Sid)) && (B_Rap.at(B_Sid) < yBin[k+1])) {
                RapIndex = k;
                break;
            }
        }

        if ((RapIndex == -1)) {continue;}

        for (int i=0;i<A_Rap.size();i++) {
            H_Brap[CenIndex][RapIndex]->Fill(A_Rap.at(j));
        }
        
    }

    TString OutputFileName = OutMidName;
    OutputFileName += "H_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileA = new TFile(OutputFileName, "RECREATE");

    fileA->cd();
    for (int i=0;i<CentralityBinNum;i++){
        H_All[i]->Write();
        H_All_WithNetB[i]->Write();
        for (int j=0;j<yBinNum;j++){
            H_Brap[i][j]->Write();
        }
    }
    fileA.Close();

    cout << "#######################" << endl;
    cout << "# Finish storing Hist #" << endl;
    cout << "#######################" << endl;

    return;

}