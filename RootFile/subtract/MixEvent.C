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

Double_t massList(int PID , TString DataName)
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

Double_t massListSigma(int PID , TString DataName)
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

void MixEvent(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,
              float Energy)
{

    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) 
        std::vector<Int_t>   *id              = nullptr;
        std::vector<Int_t>   *evtID           = nullptr;
        std::vector<Int_t>   *runID           = nullptr;
        std::vector<Float_t> *px              = nullptr;
        std::vector<Float_t> *py              = nullptr;
        std::vector<Float_t> *pz              = nullptr;
        std::vector<Float_t> *mass            = nullptr;
        std::vector<Float_t> *energy          = nullptr;
        
        TBranch *bid              = nullptr;
        TBranch *bevtID           = nullptr;
        TBranch *brunID           = nullptr;
        TBranch *bpx              = nullptr;
        TBranch *bpy              = nullptr;
        TBranch *bpz              = nullptr;
        TBranch *bmass            = nullptr;
        TBranch *benergy          = nullptr;

        std::vector<Float_t> *dEdx               = nullptr;
        std::vector<Float_t> *m2                 = nullptr;
        std::vector<Float_t> *dcatopv            = nullptr;
        std::vector<int>     *hasTOF             = nullptr;
        std::vector<Float_t> *nSigmaProton       = nullptr;
        std::vector<Float_t> *nSigmaPion         = nullptr;
        std::vector<Float_t> *nSigmaKaon         = nullptr;
        std::vector<double>  *zTOF_proton        = nullptr;
        std::vector<double>  *zTOF_pion          = nullptr;
        std::vector<double>  *zTOF_kaon          = nullptr;
	    std::vector<int>     *IfConfuse          = nullptr;
        std::vector<double>  *Decay_Length       = nullptr;
	    std::vector<int>     *IfBadReconstructed = nullptr;

        TBranch *bdEdx               = nullptr;
        TBranch *bm2                 = nullptr;
        TBranch *bdcatopv            = nullptr;
        TBranch *bhasTOF             = nullptr;
        TBranch *bnSigmaProton       = nullptr;
        TBranch *bnSigmaPion         = nullptr;
        TBranch *bnSigmaKaon         = nullptr;
        TBranch *bzTOF_proton        = nullptr;
        TBranch *bzTOF_pion          = nullptr;
        TBranch *bzTOF_kaon          = nullptr;
	    TBranch *bIfConfuse          = nullptr;
	    TBranch *bDecay_Length       = nullptr;
	    TBranch *bIfBadReconstructed = nullptr;
    
    #else
        #if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
            std::vector<Int_t>   *id              = NULL;
            std::vector<Int_t>   *evtID           = NULL;
            std::vector<Int_t>   *runID           = NULL;
            std::vector<Float_t> *px              = NULL;
            std::vector<Float_t> *py              = NULL;
            std::vector<Float_t> *pz              = NULL;
            std::vector<Float_t> *mass            = NULL;
            std::vector<Float_t> *energy          = NULL;
            
            TBranch *bid              = NULL;
            TBranch *bevtID           = NULL;
            TBranch *brunID           = NULL;
            TBranch *bpx              = NULL;
            TBranch *bpy              = NULL;
            TBranch *bpz              = NULL;
            TBranch *bmass            = NULL;
            TBranch *benergy          = NULL;
        
            std::vector<Float_t> *dEdx               = NULL;
            std::vector<Float_t> *m2                 = NULL;
            std::vector<Float_t> *dcatopv            = NULL;
            std::vector<int>     *hasTOF             = NULL;
            std::vector<Float_t> *nSigmaProton       = NULL;
            std::vector<Float_t> *nSigmaPion         = NULL;
            std::vector<Float_t> *nSigmaKaon         = NULL;
            std::vector<double>  *zTOF_proton        = NULL;
            std::vector<double>  *zTOF_pion          = NULL;
            std::vector<double>  *zTOF_kaon          = NULL;
            std::vector<int>     *IfConfuse          = NULL;
            std::vector<double>  *Decay_Length       = NULL;
	        std::vector<int>     *IfBadReconstructed = NULL;

            TBranch *bdEdx               = NULL;
            TBranch *bm2                 = NULL;
            TBranch *bdcatopv            = NULL;
            TBranch *bhasTOF             = NULL;
            TBranch *bnSigmaProton       = NULL;
            TBranch *bnSigmaPion         = NULL;
            TBranch *bnSigmaKaon         = NULL;
            TBranch *bzTOF_proton        = NULL;
            TBranch *bzTOF_pion          = NULL;
            TBranch *bzTOF_kaon          = NULL;
            TBranch *bIfConfuse          = NULL;
            TBranch *bDecay_Length       = NULL;
	        TBranch *bIfBadReconstructed = NULL;

        #else
            std::vector<Int_t>   *id              = 0;
            std::vector<Int_t>   *evtID           = 0;
            std::vector<Int_t>   *runID           = 0;
            std::vector<Float_t> *px              = 0;
            std::vector<Float_t> *py              = 0;
            std::vector<Float_t> *pz              = 0;
            std::vector<Float_t> *mass            = 0;
            std::vector<Float_t> *energy          = 0;
            
            TBranch *bid              = 0;
            TBranch *bevtID           = 0;
            TBranch *brunID           = 0;
            TBranch *bpx              = 0;
            TBranch *bpy              = 0;
            TBranch *bpz              = 0;
            TBranch *bmass            = 0;
            TBranch *benergy          = 0;

            std::vector<Float_t> *dEdx               = 0;
            std::vector<Float_t> *m2                 = 0;
            std::vector<Float_t> *dcatopv            = 0;
            std::vector<int>     *hasTOF             = 0;
            std::vector<Float_t> *nSigmaProton       = 0;
            std::vector<Float_t> *nSigmaPion         = 0;
            std::vector<Float_t> *nSigmaKaon         = 0;
            std::vector<double>  *zTOF_proton        = 0;
            std::vector<double>  *zTOF_pion          = 0;
            std::vector<double>  *zTOF_kaon          = 0;
            std::vector<int>     *IfConfuse          = 0;
            std::vector<double>  *Decay_Length       = 0;
	        std::vector<int>     *IfBadReconstructed = 0;

            TBranch *bdEdx               = 0;
            TBranch *bm2                 = 0;
            TBranch *bdcatopv            = 0;
            TBranch *bhasTOF             = 0;
            TBranch *bnSigmaProton       = 0;
            TBranch *bnSigmaPion         = 0;
            TBranch *bnSigmaKaon         = 0;
            TBranch *bzTOF_proton        = 0;
            TBranch *bzTOF_pion          = 0;
            TBranch *bzTOF_kaon          = 0;
            TBranch *bIfConfuse          = 0;
            TBranch *bDecay_Length       = 0;
	        TBranch *bIfBadReconstructed = 0;

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

    hadronTree->SetBranchAddress("PDGMult",&mult);
    hadronTree->SetBranchAddress("PDG",&id,&bid);
    hadronTree->SetBranchAddress("mix_px",&px,&bpx);
    hadronTree->SetBranchAddress("mix_py",&py,&bpy);
    hadronTree->SetBranchAddress("mix_pz",&pz,&bpz);
    hadronTree->SetBranchAddress("InvariantMass",&mass,&bmass);
    hadronTree->SetBranchAddress("nSigmaProton",&nSigmaProton,&bnSigmaProton);
    hadronTree->SetBranchAddress("nSigmaPion"  ,&nSigmaPion  ,&bnSigmaPion  );
    hadronTree->SetBranchAddress("nSigmaKaon"  ,&nSigmaKaon  ,&bnSigmaKaon  );

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;
    
    TLorentzVector p0,p3;
    double kstar, rap;
    TVector3 BetaTemp;

    int KaonpPID = 321,KaonmPID = -321,PionpPID = 211,PionmPID = -211,LambdaPID = 3122,LambdabPID = -3122,XiPID = 3312,XibPID = -3312,OmegaPID = 3334,OmegabPID = -3334;

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
