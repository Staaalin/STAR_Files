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

void vMass()
{
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) 
        std::vector<Int_t>   *PDG             = nullptr;
        std::vector<Int_t>   *evtID           = nullptr;
        std::vector<Int_t>   *runID           = nullptr;
        std::vector<Float_t> *px              = nullptr;
        std::vector<Float_t> *py              = nullptr;
        std::vector<Float_t> *pz              = nullptr;
        std::vector<Float_t> *InvarentMass    = nullptr;
        std::vector<Float_t> *energy          = nullptr;
        
        TBranch *bPDG             = nullptr;
        TBranch *bevtID           = nullptr;
        TBranch *brunID           = nullptr;
        TBranch *bpx              = nullptr;
        TBranch *bpy              = nullptr;
        TBranch *bpz              = nullptr;
        TBranch *bInvarentMass    = nullptr;
        TBranch *benergy          = nullptr;
    
    #else
        #if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
            std::vector<Int_t>   *PDG             = NULL;
            std::vector<Int_t>   *evtID           = NULL;
            std::vector<Int_t>   *runID           = NULL;
            std::vector<Float_t> *px              = NULL;
            std::vector<Float_t> *py              = NULL;
            std::vector<Float_t> *pz              = NULL;
            std::vector<Float_t> *InvarentMass    = NULL;
            std::vector<Float_t> *energy          = NULL;
            
            TBranch *bPDG             = NULL;
            TBranch *bevtID           = NULL;
            TBranch *brunID           = NULL;
            TBranch *bpx              = NULL;
            TBranch *bpy              = NULL;
            TBranch *bpz              = NULL;
            TBranch *bInvarentMass    = NULL;
            TBranch *benergy          = NULL;
        
        #else
            std::vector<Int_t>   *PDG             = 0;
            std::vector<Int_t>   *evtID           = 0;
            std::vector<Int_t>   *runID           = 0;
            std::vector<Float_t> *px              = 0;
            std::vector<Float_t> *py              = 0;
            std::vector<Float_t> *pz              = 0;
            std::vector<Float_t> *InvarentMass    = 0;
            std::vector<Float_t> *energy          = 0;
            
            TBranch *bPDG             = 0;
            TBranch *bevtID           = 0;
            TBranch *brunID           = 0;
            TBranch *bpx              = 0;
            TBranch *bpy              = 0;
            TBranch *bpz              = 0;
            TBranch *bInvarentMass    = 0;
            TBranch *benergy          = 0;
        #endif
    #endif
    Int_t refMult,grefMult,PDGMult;

    Int_t kBinNum = 2500;
    Float_t kmin = 0;
    Float_t kmax = 5;
	const TString ParticleName[] = { "Lambda" , "Lambdab" , "Omega" };
	const int ParticlePDG[]      = {   3122   ,   -3122   ,   3334  };
    const float FitRange[3][2]   = { {1.0,1.2}, {1.0,1.2} ,{1.6,1.8}};
	int HSize = sizeof(ParticleName)/sizeof(ParticleName[0]);
	// const int HSize = 3;
	// TH1D *HMass[HSize];
    std::vector<TH1D *>    HMass[HSize];
    std::vector<TCanvas *> CMass[HSize];
    std::vector<TF1 *> Vfit[HSize];
	for (int i=0;i<HSize;i++){
		TString HistName1 = "HM";
		TString HistName2 = "The Mass of ";
        TString FitName1 = "FN";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
        FitName1 += ParticleName[i];
		HMass.emplace_back(new TH1D(HistName1, HistName2, kBinNum, kmin, kmax));
		Vfit.emplace_back(new TF1(FitName1, "gaus", FitRange[i][0], FitRange[i][1]));
		
		HMass[i]->GetXaxis()->SetTitle("Mass [GeV]");
		HMass[i]->GetYaxis()->SetTitle("Counts");
	}

    //load data  
    TString midname = "/star/data01/pwg/svianping/output/output_";

    for(int i=1000;i <= 1050;i++){
        TString filename = midname;
        filename+="00";
        filename+=i;
        filename+=".root";
        hadronTree->Add(filename);
        // cout<<filename<<endl;
    }
    
    hadronTree->SetBranchAddress("Mult"          ,&PDGMult);
    hadronTree->SetBranchAddress("refMult"       ,&refMult);
    hadronTree->SetBranchAddress("grefMult"      ,&grefMult);
    hadronTree->SetBranchAddress("PDG"           ,&PDG          ,&bPDG         );
    hadronTree->SetBranchAddress("mix_px"        ,&px           ,&bpx          );
    hadronTree->SetBranchAddress("mix_py"        ,&py           ,&bpy          );
    hadronTree->SetBranchAddress("mix_pz"        ,&pz           ,&bpz          );
    hadronTree->SetBranchAddress("InvarentMass"  ,&InvarentMass ,&bInvarentMass);
    // hadronTree->SetBranchAddress("energy"        ,&energy       ,&benergy      );

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;
    
    //read data
    for (int i=0;i<nentries;i++){
        hadronTree->GetEntry(i);

        if (PDGMult != PDG[0].size()){
            cout<<"Warning! PDGMult = "<<PDGMult<<", but PDG[0].size() = "<<PDG[0].size()<<endl;
        }
		for (int j=0;j<PDGMult;j++){
			for (int k=0;k<HSize;k++){
				if(PDG->at(j) == ParticlePDG[k]){
					HMass[k]->Fill(InvarentMass->at(j));
					break;
				}
			}
		}

    }

    //

    TFile *file = new TFile("/star/u/svianping/STAR_Files/KFParticle4Lambda/output/CheckOutput.root", "RECREATE");


	for (int i=0;i<HSize;i++){
		HMass[i]->Write();
	}

    file->Write();

    return 0;
}
