#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <random>
#include "math.h"
#include "string.h"
#include <vector>
#ifndef __CINT__
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
#endif
#include <iostream>
#include <map>
#include <stdio.h>
using namespace std;

void readTree()
{
    const Int_t maxMultiplicity = 20000;
    Int_t refMult,grefMult;
    Int_t PDG[maxMultiplicity],evtID[maxMultiplicity],runID[maxMultiplicity];
    Float_t px[maxMultiplicity],py[maxMultiplicity],pz[maxMultiplicity];
    Float_t mass[maxMultiplicity],energy[maxMultiplicity];

    Int_t kBinNum = 1000;
    Float_t kmin = 0;
    Float_t kmax = 10;
	TString ParticleName[] = { "Lambda" , "Lambdab" , "Omega" };
	int ParticlePDG[]      = {   3122   ,   -3122   ,   3334  };
	int HSize = sizeof(ParticleName)/sizeof(ParticleName[0]);
	TH1D *HMass[HSize];
	for (int i=0;i<ParticleName;i++){
		TString HistName1 = "HM";
		TString HistName2 = "The Mass of ";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
		HMass[i] = new TH1D(HistName1, HistName2, kBinNum, kmin, kmax);
	}

    //load data  
    TString midname = "/star/data01/pwg/svianping/output/output_";

    TChain *hadronTree = new TChain("hadronTree");
    for(int i=2892;i <= 2991;i++){
        TString filename = midname;
        filename+="00";
        filename+=i;
        filename+=".root";
        hadronTree->Add(filename);
        // cout<<filename<<endl;
    }
    
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

		for (int j=0;j<)

    }

    TFile *file = new TFile("CheckOutput.root", "RECREATE");

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
