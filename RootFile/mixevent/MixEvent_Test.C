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
const float omegaMass = 1.67245;
const float lambdaMass = 1.11568;
const float xiMass = 1.32171;
const float kaonMass = 0.493677;
const int maxTrack = 30000;

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

// void MixEvent_Test(Int_t Aid,Int_t Bid,float Energy,float B_RapMin,float B_RapMax)

void MixEvent_Test(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,TString OutMidName,int A_PDG,int B_PDG,int Mode = 0)
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

    // Int_t Aid = -3122, Bid = 3334;
    const Int_t MultBin39[]   = { 10000 , 47 , 38 , 26 , 0};
    const Int_t MultBinSize39 = sizeof(MultBin39)/sizeof(MultBin39[0]);
    Int_t MultBinSize,FileMax;
    std::vector<int> MultBin;
    MultBinSize = MultBinSize39;
    for(int i=0;i<MultBinSize39;i++){
        MultBin.push_back(MultBin39[i]);
    }
    const Int_t A_NumBin[] = {0,1,2,3,4,5,6,200};
    const Int_t B_NumBin[] = {0,1,2,3,4,200};
    // const Int_t MultBin[] = {200, 225, 250, 275, 300, 350, 400};
    // const Int_t A_NumBin[] = {0,200};
    // const Int_t B_NumBin[] = {0,200};
    const Int_t A_NumBinSize = sizeof(A_NumBin)/sizeof(A_NumBin[0]);
    const Int_t B_NumBinSize = sizeof(B_NumBin)/sizeof(B_NumBin[0]);

    TString midname = "";
    midname = MidName;
    Int_t kBinNum = 1000;
    Float_t kmin = 0;
    Float_t kmax = 10;


    TH1::SetDefaultSumw2("kTRUE");
    TH1D *Ak = new TH1D("Ak", "Ak number", kBinNum, kmin, kmax);
    TH1D *Bk = new TH1D("Bk", "Bk number", kBinNum, kmin, kmax);
    TH1D *Ck = (TH1D *) Ak->Clone();
    Ck->SetName("Ck");
    Ck->SetTitle("Ck");
    TH1D *Dk = (TH1D *) Ak->Clone();
    Dk->SetName("Dk");
    Dk->SetTitle("Dk");
    TH2D *MIX_Con = new TH2D("MIX_Con", "MIX_Con", MultBinSize, 0, MultBinSize, A_NumBinSize, 0, A_NumBinSize);
    
    TH1D *Ay = new TH1D("Ay", "Ay number", kBinNum, -kmax*0.8, kmax*0.8);
    TH1D *By = new TH1D("By", "By number", kBinNum, -kmax*0.8, kmax*0.8);
    TH1D *Cy = (TH1D *) Ay->Clone();
    Cy->SetName("Cy");
    Cy->SetTitle("Cy");

    cout<<"FileMax = "<<FileMax<<endl;
    //load data  
    TChain *hadronTree = new TChain("hadronTree");
    for(int i=StartFileIndex;i <= EndFileIndex;i++){//62 6280; 39 9027
        //TString filename = "data39_Default/afterART7.7a_";//path of data 
        //TString filename = "data39_String/afterART7.7a_";//path of data 
        //TString filename = "data39_Hardon/afterART7.7a_";//path of data 
        //TString filename = "data39_PTurnAu/afterART7.7a_";//path of data
        TString filename = midname;
        filename+=i;
        filename+=".root";
        hadronTree->Add(filename);
        cout<<filename<<endl;
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
    
    TLorentzVector p0,p3;
    float kstar, rap;
    TVector3 BetaTemp;
    std::vector<float> A_px,B_px;
    std::vector<float> A_py,B_py;
    std::vector<float> A_pz,B_pz;
    std::vector<float> A_x,B_x;
    std::vector<float> A_y,B_y;
    std::vector<float> A_z,B_z;
    std::vector<float> A_rap,B_rap;
    std::vector<float> A_mass,B_mass;
    std::vector<TString> A_Kind,B_Kind;
    std::vector<float> Sum_A_px  [MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1],   Sum_B_px[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1];
    std::vector<float> Sum_A_py  [MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1],   Sum_B_py[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1];
    std::vector<float> Sum_A_pz  [MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1],   Sum_B_pz[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1];
    std::vector<float> Sum_A_x   [MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1],    Sum_B_x[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1];
    std::vector<float> Sum_A_y   [MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1],    Sum_B_y[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1];
    std::vector<float> Sum_A_z   [MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1],    Sum_B_z[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1];
    std::vector<float> Sum_A_rap [MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1],  Sum_B_rap[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1];
    std::vector<float> Sum_A_mass[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1], Sum_B_mass[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1];



    Int_t NumA = 0, NumB = 0, SumNumA = 0, SumNumB = 0, BP[5];
    Int_t SumNum[MultBinSize-1][A_NumBinSize-1][B_NumBinSize-1];
    for(int i=0;i<MultBinSize-1;i++){
        for(int k=0;k<A_NumBinSize-1;k++){
            for(int m=0;m<B_NumBinSize-1;m++){
                SumNum[i][k][m] = 0;
            }
        }
    }

    //read data; if the nst event contains particle A, then record n in A_Loc, px in A_px ……; so do B particle
    for (int i=0;i<nentries;i++){
        hadronTree->GetEntry(i);
        if ((i+1)%50 == 0) {
            cout<<"Calculating Event "<<(i+1)<<"/"<<nentries<<endl;
        }
        // if (b >= 1.7){
        //     continue;
        // }
        A_Kind.clear();B_Kind.clear();
        A_px.clear();B_px.clear();
        A_py.clear();B_py.clear();
        A_pz.clear();B_pz.clear();
        A_rap.clear();B_rap.clear();
        A_mass.clear();B_mass.clear();
        NumA = 0;NumB = 0;
        for (int j=0;j < PDGMult;j++){
            
            p0.SetXYZM(mix_px->at(j),mix_py->at(j),mix_pz->at(j),massList(PDG->at(j)));
            rap = p0.Rapidity();

            if(PDG->at(j)==A_PDG || PDG->at(j)==B_PDG){
                //&& rap > 0
                if(PDG->at(j)==A_PDG){
                    if      (fabs(InvariantMass->at(j) - massList(A_PDG)) <= 3*massListSigma(A_PDG)) {A_Kind.push_back("Mid");}
                    else if (fabs(InvariantMass->at(j) - massList(A_PDG)) <= 6*massListSigma(A_PDG)) {continue;}
                    A_px.push_back(mix_px->at(j));
                    A_py.push_back(mix_py->at(j));
                    A_pz.push_back(mix_pz->at(j));
                    A_rap.push_back(rap);
                    A_mass.push_back(InvariantMass->at(j));
                    NumA++;
                }
                else if(PDG->at(j)==B_PDG && rap > B_RapMin && rap < B_RapMax){
                    if      (fabs(InvariantMass->at(j) - massList(B_PDG)) <= 3*massListSigma(B_PDG)) {B_Kind.push_back("Mid");}
                    else if (fabs(InvariantMass->at(j) - massList(B_PDG)) <= 6*massListSigma(B_PDG)) {continue;}
                    B_px.push_back(mix_px->at(j));
                    B_py.push_back(mix_py->at(j));
                    B_pz.push_back(mix_pz->at(j));
                    B_rap.push_back(rap);
                    B_mass.push_back(InvariantMass->at(j));
                    NumB++;
                }
                // cout<<NumA<<","<<NumB<<endl;
            }
        }
        //same
        if(NumB != 1){
            continue;
        }
        if(NumA != 0 && NumB != 0){
            for(int j=0;j<NumA;j++){
                for(int k=0;k<NumB;k++){
                    TLorentzVector p1;
                    p1.SetXYZM(A_px[j],A_py[j],A_pz[j],A_mass[j]);
                    TLorentzVector p2;
                    p2.SetXYZM(B_px[k],B_py[k],B_pz[k],B_mass[k]);
                    // testPx->Fill(B_px[k]);
                    // testPy->Fill(B_py[k]);
                    // testPz->Fill(B_pz[k]);
                    // testPx->Fill(A_px[j]);
                    // testPy->Fill(A_py[j]);
                    // testPz->Fill(A_pz[j]);
                    // testPx->Fill((p1).X());
                    // testPy->Fill((p1).Y());
                    // testPz->Fill((p1).Z());
                    // testPx->Fill((p2).X());
                    // testPy->Fill((p2).Y());
                    // testPz->Fill((p2).Z());
                    // testPx->Fill((p1-p2).X());
                    // testPy->Fill((p1-p2).Y());
                    // testPz->Fill((p1-p2).Z());

                    TLorentzVector p3 = p1 + p2;
                    p1.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                    // testPx->Fill((p1-p2).X());
                    // testPy->Fill((p1-p2).Y());
                    // testPz->Fill((p1-p2).Z());
                    kstar = 0.5 * (p1 - p2).Rho();
                    Ak->Fill(kstar);

                    Ay->Fill(A_rap[j]-B_rap[k]);
                }
            }
        }else{
            //cout<<"Failed"<<endl;
            continue;
        }
        //mix
        for (int j=1;j<MultBinSize;j++){
            if (MultBin[j-1] <= PDGMult && PDGMult < MultBin[j]){
                BP[0] = j-1;
                break;
            }
        }
        for (int j=1;j<A_NumBinSize;j++){
            if (A_NumBin[j-1] <= NumA && NumA < A_NumBin[j]){
                BP[1] = j-1;
                break;
            }
        }
        for (int j=1;j<B_NumBinSize;j++){
            if (B_NumBin[j-1] <= NumB && NumB < B_NumBin[j]){
                BP[2] = j-1;
                break;
            }
        }
        
        for(int j=0;j<NumA;j++){
            Sum_A_x   [BP[0]][BP[1]][BP[2]].push_back(A_x[j]);
            Sum_A_y   [BP[0]][BP[1]][BP[2]].push_back(A_y[j]);
            Sum_A_z   [BP[0]][BP[1]][BP[2]].push_back(A_z[j]);
            Sum_A_px  [BP[0]][BP[1]][BP[2]].push_back(A_px[j]);
            Sum_A_py  [BP[0]][BP[1]][BP[2]].push_back(A_py[j]);
            Sum_A_pz  [BP[0]][BP[1]][BP[2]].push_back(A_pz[j]);
            Sum_A_rap [BP[0]][BP[1]][BP[2]].push_back(A_rap[j]);
            Sum_A_mass[BP[0]][BP[1]][BP[2]].push_back(A_mass[j]);
        }
        for(int j=0;j<NumB;j++){
            Sum_B_x   [BP[0]][BP[1]][BP[2]].push_back(B_x[j]);
            Sum_B_y   [BP[0]][BP[1]][BP[2]].push_back(B_y[j]);
            Sum_B_z   [BP[0]][BP[1]][BP[2]].push_back(B_z[j]);
            Sum_B_px  [BP[0]][BP[1]][BP[2]].push_back(B_px[j]);
            Sum_B_py  [BP[0]][BP[1]][BP[2]].push_back(B_py[j]);
            Sum_B_pz  [BP[0]][BP[1]][BP[2]].push_back(B_pz[j]);
            Sum_B_rap [BP[0]][BP[1]][BP[2]].push_back(B_rap[j]);
            Sum_B_mass[BP[0]][BP[1]][BP[2]].push_back(B_mass[j]);
        }
        SumNum[BP[0]][BP[1]][BP[2]]++;

        if(SumNum[BP[0]][BP[1]][BP[2]] == 10){
            SumNumA = Sum_A_px[BP[0]][BP[1]][BP[2]].size(); SumNumB = Sum_B_px[BP[0]][BP[1]][BP[2]].size();
            for(int j=0;j<SumNumA;j++){
                for(int k=0;k<SumNumB;k++){
                    TLorentzVector p1;
                    p1.SetXYZM(Sum_A_px[BP[0]][BP[1]][BP[2]][j],Sum_A_py[BP[0]][BP[1]][BP[2]][j],Sum_A_pz[BP[0]][BP[1]][BP[2]][j],Sum_A_mass[BP[0]][BP[1]][BP[2]][j]);
                    TLorentzVector p2;
                    p2.SetXYZM(Sum_B_px[BP[0]][BP[1]][BP[2]][k],Sum_B_py[BP[0]][BP[1]][BP[2]][k],Sum_B_pz[BP[0]][BP[1]][BP[2]][k],Sum_B_mass[BP[0]][BP[1]][BP[2]][k]);
                    // testPx->Fill(Sum_B_px[BP][k]);
                    // testPy->Fill(Sum_B_py[BP][k]);
                    // testPz->Fill(Sum_B_pz[BP][k]);
                    // testPx->Fill(Sum_A_px[BP][j]);
                    // testPy->Fill(Sum_A_py[BP][j]);
                    // testPz->Fill(Sum_A_pz[BP][j]);
                    p3 = p1 + p2;
                    p1.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                    // testPx->Fill((p1).X());
                    // testPy->Fill((p1).Y());
                    // testPz->Fill((p1).Z());
                    // testPx->Fill((p2).X());
                    // testPy->Fill((p2).Y());
                    // testPz->Fill((p2).Z());
                    kstar = 0.5 * (p1 - p2).Rho();
                    Bk->Fill(kstar);

                    By->Fill(Sum_A_rap[BP[0]][BP[1]][BP[2]][j]-Sum_B_rap[BP[0]][BP[1]][BP[2]][k]);
                }
            }
            MIX_Con->Fill(BP[0],BP[1]);
            Sum_A_x   [BP[0]][BP[1]][BP[2]].clear();    Sum_B_x[BP[0]][BP[1]][BP[2]].clear();
            Sum_A_y   [BP[0]][BP[1]][BP[2]].clear();    Sum_B_y[BP[0]][BP[1]][BP[2]].clear();
            Sum_A_z   [BP[0]][BP[1]][BP[2]].clear();    Sum_B_z[BP[0]][BP[1]][BP[2]].clear();
            Sum_A_px  [BP[0]][BP[1]][BP[2]].clear();   Sum_B_px[BP[0]][BP[1]][BP[2]].clear();
            Sum_A_py  [BP[0]][BP[1]][BP[2]].clear();   Sum_B_py[BP[0]][BP[1]][BP[2]].clear();
            Sum_A_pz  [BP[0]][BP[1]][BP[2]].clear();   Sum_B_pz[BP[0]][BP[1]][BP[2]].clear();
            Sum_A_rap [BP[0]][BP[1]][BP[2]].clear();  Sum_B_rap[BP[0]][BP[1]][BP[2]].clear();
            Sum_A_mass[BP[0]][BP[1]][BP[2]].clear(); Sum_B_mass[BP[0]][BP[1]][BP[2]].clear();
            SumNum[BP[0]][BP[1]][BP[2]] = 0;
        }
    }

    TString OutputName = OutMidName;
    OutputName += OutputFileIndex;
    OutputName += ".root";

    TFile *file = new TFile(OutputName, "RECREATE");

    Int_t binNorm[2];
    binNorm[0] = Ak->FindBin(3);
    binNorm[1] = Ak->FindBin(4); 
    Double_t factorN1 = 1.0*Bk->Integral(binNorm[0], binNorm[1]) / Ak->Integral(binNorm[0], binNorm[1]);
    Ck->Divide(Ak,Bk,1.0,1.0/factorN1);
    Ck->Sumw2();

    binNorm[0] = Ay->FindBin(3);
    binNorm[1] = Ay->FindBin(4); 
    Double_t factorN2 = 1.0*By->Integral(binNorm[0], binNorm[1]) / Ay->Integral(binNorm[0], binNorm[1]);
    Cy->Divide(Ay,By,1.0,1.0/factorN2);
    Cy->Sumw2();

    Ak->GetXaxis()->SetTitle("k^{*} [GeV]");
    Ak->GetYaxis()->SetTitle("Counts");
    Bk->GetXaxis()->SetTitle("k^{*} [GeV]");
    Bk->GetYaxis()->SetTitle("Counts");
    Ck->GetXaxis()->SetTitle("k^{*} [GeV]");
    Ck->GetYaxis()->SetTitle("C(k^{*})");
    Ck->SetStats(0);

    Ay->GetXaxis()->SetTitle("#Delta y");
    Ay->GetYaxis()->SetTitle("Counts");
    By->GetXaxis()->SetTitle("#Delta y");
    By->GetYaxis()->SetTitle("Counts");
    Cy->GetXaxis()->SetTitle("#Delta y");
    Cy->GetYaxis()->SetTitle("C(#Delta y)");
    Cy->SetStats(0);

    Ak->Write();
    Bk->Write();
    Ck->Write();
    
    Ay->Write();
    By->Write();
    Cy->Write();
    MIX_Con->Write();
    


    file->Write();

    return 0;
}
