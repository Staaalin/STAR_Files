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

void readTree()
{
    const Int_t maxMultiplicity = 20000;
    Int_t refMult,grefMult,Mult;
    Int_t PDG[maxMultiplicity],evtID[maxMultiplicity],runID[maxMultiplicity];
    Float_t px[maxMultiplicity],py[maxMultiplicity],pz[maxMultiplicity];
    Float_t InvarentMass[maxMultiplicity],energy[maxMultiplicity];

    Int_t kBinNum = 1000;
    Float_t kmin = 0;
    Float_t kmax = 10;
	TString ParticleName[] = { "Lambda" , "Lambdab" , "Omega" };
	int ParticlePDG[]      = {   3122   ,   -3122   ,   3334  };
	int HSize = sizeof(ParticleName)/sizeof(ParticleName[0]);
	TH1D *HMass[HSize];
	for (int i=0;i<HSize;i++){
		TString HistName1 = "HM";
		TString HistName2 = "The Mass of ";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
		HMass[i] = new TH1D(HistName1, HistName2, kBinNum, kmin, kmax);
		
		HMass[i]->GetXaxis()->SetTitle("Mass [GeV]");
		HMass[i]->GetYaxis()->SetTitle("Counts");
	}

    //load data  
    TString midname = "/star/data01/pwg/svianping/output/output_";

    TChain *hadronTree = new TChain("hadronTree");
    for(int i=2229;i <= 2230;i++){
        TString filename = midname;
        filename+="00";
        filename+=i;
        filename+=".root";
        hadronTree->Add(filename);
        // cout<<filename<<endl;
    }
    
    hadronTree->SetBranchAddress("Mult",&Mult);
    hadronTree->SetBranchAddress("refMult",&refMult);
    hadronTree->SetBranchAddress("grefMult",&grefMult);
    hadronTree->SetBranchAddress("PDG",&PDG);
    hadronTree->SetBranchAddress("px",&px);
    hadronTree->SetBranchAddress("py",&py);
    hadronTree->SetBranchAddress("pz",&pz);
    hadronTree->SetBranchAddress("mass",&InvarentMass);

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;
    
    //read data
    for (int i=0;i<nentries;i++){
        hadronTree->GetEntry(i);

		for (int j=0;j<Mult;j++){
			for (int k=0;k<HSize;k++){
				if(PDG[j] == ParticlePDG[k]){
					HMass[k]->Fill(InvarentMass[j]);
					break;
				}
			}
		}

    }

    TFile *file = new TFile("/star/u/svianping/STAR_Files/KFParticle4Lambda/output/CheckOutput.root", "RECREATE");


	for (int i=0;i<HSize;i++){
		HMass[i]->Write();
	}

    file->Write();

    return 0;
}
