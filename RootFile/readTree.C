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
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) 
        std::vector<Int_t>   *PDG             = nullptr;
        std::vector<Int_t>   *evtID           = nullptr;
        std::vector<Int_t>   *runID           = nullptr;
        std::vector<Float_t> *px              = nullptr;
        std::vector<Float_t> *py              = nullptr;
        std::vector<Float_t> *pz              = nullptr;
        std::vector<Float_t> *InvariantMass   = nullptr;
        std::vector<Float_t> *energy          = nullptr;
        
        TBranch *bPDG             = nullptr;
        TBranch *bevtID           = nullptr;
        TBranch *brunID           = nullptr;
        TBranch *bpx              = nullptr;
        TBranch *bpy              = nullptr;
        TBranch *bpz              = nullptr;
        TBranch *bInvariantMass   = nullptr;
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
            std::vector<Int_t>   *PDG             = NULL;
            std::vector<Int_t>   *evtID           = NULL;
            std::vector<Int_t>   *runID           = NULL;
            std::vector<Float_t> *px              = NULL;
            std::vector<Float_t> *py              = NULL;
            std::vector<Float_t> *pz              = NULL;
            std::vector<Float_t> *InvariantMass   = NULL;
            std::vector<Float_t> *energy          = NULL;
            
            TBranch *bPDG             = NULL;
            TBranch *bevtID           = NULL;
            TBranch *brunID           = NULL;
            TBranch *bpx              = NULL;
            TBranch *bpy              = NULL;
            TBranch *bpz              = NULL;
            TBranch *bInvariantMass   = NULL;
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
            std::vector<Int_t>   *PDG             = 0;
            std::vector<Int_t>   *evtID           = 0;
            std::vector<Int_t>   *runID           = 0;
            std::vector<Float_t> *px              = 0;
            std::vector<Float_t> *py              = 0;
            std::vector<Float_t> *pz              = 0;
            std::vector<Float_t> *InvariantMass   = 0;
            std::vector<Float_t> *energy          = 0;
            
            TBranch *bPDG             = 0;
            TBranch *bevtID           = 0;
            TBranch *brunID           = 0;
            TBranch *bpx              = 0;
            TBranch *bpy              = 0;
            TBranch *bpz              = 0;
            TBranch *bInvariantMass   = 0;
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
    Int_t refMult,grefMult,PDGMult;

    Int_t kBinNum = 1000;
    Float_t kmin = 0;
    Float_t kmax = 10;
	const TString ParticleName[] = { "Lambda" , "Lambdab" , "Omega"  , "Proton"  , "Protonb" , "Pion" , "Pionb"  , "Kaon" , "Kaonb" , "Un-Identified"};
	const int ParticlePDG[]      = {   3122   ,   -3122   ,   3334   ,   2212    ,  -2212    ,  211   ,  -211    ,  321   ,  -321   ,       0        };
	// const int HSize = sizeof(ParticleName)/sizeof(ParticleName[0]);
	const int HSize = 10;
	TH1D *HMass[HSize];
	TH1D *HBadMass[HSize];
	TH1D *HP[HSize];
	TH1D *HRapidity[HSize];
	TH1D *Hdcatopv[HSize];
	TH1D *HnSigma[HSize];
	TH2D *H_dEdx_p[HSize];
	TH2D *H_sigma_TOF[HSize];
	TH2D *H_sigma_p[HSize];
	TH2D *H_sigma_pT[HSize];
	TH2D *H_Mass_DL[HSize];
	TH2D *H_Mass_dcatopv[HSize];
    // for confused
	TH2D *HC_sigma_TOF[3];
	for (int i=0;i<HSize;i++){
		TString HistName1 = "HM";
		TString HistName2 = "The Mass of ";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
		HMass[i] = new TH1D(HistName1, HistName2, kBinNum, kmin, kmax);
		HMass[i]->GetXaxis()->SetTitle("Mass [GeV]");
		HMass[i]->GetYaxis()->SetTitle("Counts");
        
		TString HistName1 = "HBadM";
		TString HistName2 = "The Mass of Bad";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
		HBadMass[i] = new TH1D(HistName1, HistName2, kBinNum, kmin, kmax);
		HBadMass[i]->GetXaxis()->SetTitle("Mass [GeV]");
		HBadMass[i]->GetYaxis()->SetTitle("Counts");

		TString HistName1 = "HP";
		TString HistName2 = "The Momuntum of ";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
		HP[i] = new TH1D(HistName1, HistName2, kBinNum/10, kmin, 4);
		HP[i]->GetXaxis()->SetTitle("Momuntum [GeV]");
		HP[i]->GetYaxis()->SetTitle("Counts");
        
		TString HistName1 = "HRapidity";
		TString HistName2 = "The Rapidity of ";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
		HRapidity[i] = new TH1D(HistName1, HistName2, 30, -1.5, 1.5);
		HRapidity[i]->GetXaxis()->SetTitle("y");
		HRapidity[i]->GetYaxis()->SetTitle("Counts");
        
		TString HistName1 = "Hdcatopv";
		TString HistName2 = "The dcatopv of ";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
		Hdcatopv[i] = new TH1D(HistName1, HistName2, 60, 0, 3);
		Hdcatopv[i]->GetXaxis()->SetTitle("DCA to PV [cm]");
		Hdcatopv[i]->GetYaxis()->SetTitle("Counts");

		TString HistName1 = "HnSigma";
		TString HistName2 = "The nSigma";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
        HistName2 += " of ";
		HistName2 += ParticleName[i];
		HnSigma[i] = new TH1D(HistName1, HistName2, 600, -3, 3);
		HnSigma[i]->GetXaxis()->SetTitle("nSigma");
		HnSigma[i]->GetYaxis()->SetTitle("Counts");

		TString HistName1 = "H_dEdx_p";
		TString HistName2 = "The dE/dx vs. p of ";
		HistName1 += ParticleName[i];
		HistName2 += ParticleName[i];
		H_dEdx_p[i] = new TH2D(HistName1, HistName2, kBinNum, kmin, kmax, 100, 0, 15);
		H_dEdx_p[i]->GetXaxis()->SetTitle("p [GeV]");
		H_dEdx_p[i]->GetYaxis()->SetTitle("dE/dx [keV/cm]");
        
		TString HistName1 = "H_sigma_TOF";
		HistName1 += ParticleName[i];
		TString HistName2 = "The nSigma";
		HistName2 += ParticleName[i];
        HistName2 += " vs. zTOF";
		HistName2 += ParticleName[i];
        HistName2 += " of ";
		HistName2 += ParticleName[i];
		H_sigma_TOF[i] = new TH2D(HistName1, HistName2, 400, -3, 3, 400, -5, 5);
		H_sigma_TOF[i]->GetXaxis()->SetTitle("nSigma");
		H_sigma_TOF[i]->GetYaxis()->SetTitle("zTOF");
        
		TString HistName1 = "H_sigma_p";
		HistName1 += ParticleName[i];
		TString HistName2 = "The nSigma";
		HistName2 += ParticleName[i];
        HistName2 += " vs. momentum";
		H_sigma_p[i] = new TH2D(HistName1, HistName2, 400, 0, 5, 400, -2.5, 2.5);
		H_sigma_p[i]->GetXaxis()->SetTitle("p [GeV]");
		H_sigma_p[i]->GetYaxis()->SetTitle("nSigma");
        
		TString HistName1 = "H_sigma_pT";
		HistName1 += ParticleName[i];
		TString HistName2 = "The nSigma";
		HistName2 += ParticleName[i];
        HistName2 += " vs. pT";
		H_sigma_pT[i] = new TH2D(HistName1, HistName2, 400, 0, 3, 400, -2.5, 2.5);
		H_sigma_pT[i]->GetXaxis()->SetTitle("pT [GeV]");
		H_sigma_pT[i]->GetYaxis()->SetTitle("nSigma");
        
		TString HistName1 = "H_Mass_DL";
		HistName1 += ParticleName[i];
		TString HistName2 = "The DecayLength";
		HistName2 += ParticleName[i];
        HistName2 += " vs. Mass";
		H_Mass_DL[i] = new TH2D(HistName1, HistName2, 5*kBinNum, kmin, kmax, 400, 0, 50);
		H_Mass_DL[i]->GetXaxis()->SetTitle("Mass [GeV]");
		H_Mass_DL[i]->GetYaxis()->SetTitle("DecayLength [cm]");
        
		TString HistName1 = "H_Mass_dcatopv";
		HistName1 += ParticleName[i];
		TString HistName2 = "The DCA to PV";
		HistName2 += ParticleName[i];
        HistName2 += " vs. Mass";
		H_Mass_dcatopv[i] = new TH2D(HistName1, HistName2, 5*kBinNum, kmin, kmax, 400, 0, 50);
		H_Mass_dcatopv[i]->GetXaxis()->SetTitle("Mass [GeV]");
		H_Mass_dcatopv[i]->GetYaxis()->SetTitle("DCA2PV [cm]");
	}
	const TString ConfusedParticleName[] = {"Proton" , "Pion"  , "Kaon"};
    for (int i=0;i<3;i++) {
		TString HistName1 = "HC_sigma_TOF";
		HistName1 += ConfusedParticleName[i];
		TString HistName2 = "The nSigma";
		HistName2 += ConfusedParticleName[i];
        HistName2 += " vs. zTOF";
		HistName2 += ConfusedParticleName[i];
        HistName2 += " of ";
		HistName2 += ConfusedParticleName[i];
		HC_sigma_TOF[i] = new TH2D(HistName1, HistName2, 200, -35, 35, 200, -4, 4);
		HC_sigma_TOF[i]->GetXaxis()->SetTitle("nSigma");
		HC_sigma_TOF[i]->GetYaxis()->SetTitle("zTOF");
    }

    //load data  
    TChain *hadronTree = new TChain("hadronTree");
    // TString midname = "/star/data01/pwg/svianping/output/output_";
    TString midname = "/star/u/svianping/STAR_Files/KFParticle4Lambda/output_999998.root";

    if (midname != "/star/u/svianping/STAR_Files/KFParticle4Lambda/output_999998.root"){
        for(int i=66090;i <= 66389;i++){
            TString filename = midname;
            filename+=i;
            filename+=".root";
            hadronTree->Add(filename);
            // cout<<filename<<endl;
        }
    }
    else{
        hadronTree->Add(midname);
    }
    // hadronTree->Add("/star/u/svianping/STAR_Files/KFParticle4Lambda/output_999998.root");
    
    hadronTree->SetBranchAddress("PDGMult"       ,&PDGMult);
    hadronTree->SetBranchAddress("refMult"       ,&refMult);
    hadronTree->SetBranchAddress("grefMult"      ,&grefMult);
    hadronTree->SetBranchAddress("PDG"           ,&PDG          ,&bPDG         );
    hadronTree->SetBranchAddress("mix_px"        ,&px           ,&bpx          );
    hadronTree->SetBranchAddress("mix_py"        ,&py           ,&bpy          );
    hadronTree->SetBranchAddress("mix_pz"        ,&pz           ,&bpz          );
    hadronTree->SetBranchAddress("InvariantMass"  ,&InvariantMass ,&bInvariantMass);
    // hadronTree->SetBranchAddress("energy"        ,&energy       ,&benergy      );

	// Used for QA
    hadronTree->SetBranchAddress("dEdx"                 ,&dEdx                   ,&bdEdx                   );
    hadronTree->SetBranchAddress("m2"                   ,&m2                     ,&bm2                     );
    hadronTree->SetBranchAddress("dcatopv"              ,&dcatopv                ,&bdcatopv                );
    hadronTree->SetBranchAddress("hasTOF"               ,&hasTOF                 ,&bhasTOF                 );
    hadronTree->SetBranchAddress("nSigmaProton"         ,&nSigmaProton           ,&bnSigmaProton           );
    hadronTree->SetBranchAddress("nSigmaPion"           ,&nSigmaPion             ,&bnSigmaPion             );
    hadronTree->SetBranchAddress("nSigmaKaon"           ,&nSigmaKaon             ,&bnSigmaKaon             );
    hadronTree->SetBranchAddress("zTOF_proton"          ,&zTOF_proton            ,&bzTOF_proton            );
    hadronTree->SetBranchAddress("zTOF_pion"            ,&zTOF_pion              ,&bzTOF_pion              );
    hadronTree->SetBranchAddress("zTOF_kaon"            ,&zTOF_kaon              ,&bzTOF_kaon              );
	hadronTree->SetBranchAddress("IfConfuse"            ,&IfConfuse              ,&bIfConfuse              );
	hadronTree->SetBranchAddress("Decay_Length"         ,&Decay_Length           ,&bDecay_Length           );
	hadronTree->SetBranchAddress("IfBadReconstructed"   ,&IfBadReconstructed     ,&bIfBadReconstructed     );

    const Int_t nentries=hadronTree->GetEntries();
    cout << "entries number: " << nentries << endl;
    
    //read data
    int ReadedEntries = 0 , ReadedPercent = 0;
    for (int i=0;i<nentries;i++){
        hadronTree->GetEntry(i);

        if (PDGMult != PDG[0].size()){
            cout<<"Warning! PDGMult = "<<PDGMult<<", but PDG[0].size() = "<<PDG[0].size()<<endl;
        }
        if (PDG[0].size() != px[0].size()){
            cout<<"Warning! PDG[0].size() = "<<PDG[0].size()<<", but px[0].size() = "<<px[0].size()<<endl;
        }
		for (int j=0;j<PDGMult;j++){
			for (int k=0;k<HSize;k++){
				if (PDG->at(j) == ParticlePDG[k] && IfConfuse->at(j) != 1){
                    TLorentzVector p0;
                    Float_t Mag = pow(pow(px->at(j),2) + pow(py->at(j),2) + pow(pz->at(j),2),0.5);
                    p0.SetPxPyPzE(px->at(j),py->at(j),pz->at(j),pow(pow(px->at(j),2) + pow(py->at(j),2) + pow(pz->at(j),2) + pow(InvariantMass->at(j),2),0.5));
                    float rap = p0.Rapidity();

                    cout<<IfBadReconstructed->at(j)<<endl;
                    if (IfBadReconstructed->at(j) == 1){
					    HMass[k]->Fill(InvariantMass->at(j));
                    }
                    if (IfBadReconstructed->at(j) == 0){
                        HBadMass[k]->Fill(InvariantMass->at(j));
                    }
					HP[k]->Fill(Mag);
					HRapidity[k]->Fill(rap);
                    Hdcatopv[k]->Fill(dcatopv->at(j));
                    H_Mass_DL[k]->Fill(InvariantMass->at(j),Decay_Length->at(j));// cout<<PDG->at(j)<<"  "<<"InvariantMass = "<<InvariantMass->at(j)<<", Decay_Length = "<<Decay_Length->at(j)<<endl;
                    H_Mass_dcatopv[k]->Fill(InvariantMass->at(j),dcatopv->at(j));
                    if (fabs(PDG->at(j)) == 2212){
                        Hdcatopv[k]->Fill(dcatopv->at(j));
                        H_dEdx_p[k]->Fill(Mag,dEdx->at(j));
                        H_sigma_TOF[k]->Fill(nSigmaProton->at(j),zTOF_proton->at(j));
                        HnSigma[k]->Fill(nSigmaProton->at(j));
                        H_sigma_p[k]->Fill(Mag,nSigmaProton->at(j));
                        H_sigma_pT[k]->Fill(pow(pow(px->at(j),2) + pow(py->at(j),2),0.5),nSigmaProton->at(j));
                    }
                    if (fabs(PDG->at(j)) == 211){
                        Hdcatopv[k]->Fill(dcatopv->at(j));
                        H_dEdx_p[k]->Fill(Mag,dEdx->at(j));
                        H_sigma_TOF[k]->Fill(nSigmaPion->at(j),zTOF_pion->at(j));
                        HnSigma[k]->Fill(nSigmaPion->at(j));
                        H_sigma_p[k]->Fill(Mag,nSigmaPion->at(j));
                        H_sigma_pT[k]->Fill(pow(pow(px->at(j),2) + pow(py->at(j),2),0.5),nSigmaPion->at(j));
                    }
                    if (fabs(PDG->at(j)) == 321){
                        Hdcatopv[k]->Fill(dcatopv->at(j));
                        H_dEdx_p[k]->Fill(Mag,dEdx->at(j));
                        H_sigma_TOF[k]->Fill(nSigmaKaon->at(j),zTOF_kaon->at(j));
                        HnSigma[k]->Fill(nSigmaKaon->at(j));
                        H_sigma_p[k]->Fill(Mag,nSigmaKaon->at(j));
                        H_sigma_pT[k]->Fill(pow(pow(px->at(j),2) + pow(py->at(j),2),0.5),nSigmaKaon->at(j));
                    }
					break;
				}

                if (IfConfuse->at(j) == 1) {
                    TLorentzVector p0;
                    Float_t Mag = pow(pow(px->at(j),2) + pow(py->at(j),2) + pow(pz->at(j),2),0.5);
                    p0.SetPxPyPzE(px->at(j),py->at(j),pz->at(j),pow(pow(px->at(j),2) + pow(py->at(j),2) + pow(pz->at(j),2) + pow(InvariantMass->at(j),2),0.5));
                    float rap = p0.Rapidity();

					HMass[HSize-1]->Fill(InvariantMass->at(j));
					HP[HSize-1]->Fill(Mag);
					HRapidity[HSize-1]->Fill(rap);
                    Hdcatopv[HSize-1]->Fill(dcatopv->at(j));
                    H_dEdx_p[HSize-1]->Fill(Mag,dEdx->at(j));
                    if (fabs(PDG->at(j)) == 2212){
                        HC_sigma_TOF[0]->Fill(nSigmaProton->at(j),zTOF_proton->at(j));
                    }
                    if (fabs(PDG->at(j)) == 211){
                        HC_sigma_TOF[1]->Fill(nSigmaPion->at(j),zTOF_pion->at(j));
                    }
                    if (fabs(PDG->at(j)) == 321){
                        HC_sigma_TOF[2]->Fill(nSigmaKaon->at(j),zTOF_kaon->at(j));
                    }
					break;
                }
			}
		}

        // if((i-ReadedEntries)>nentries/20){
        //     ReadedEntries+=nentries/20;
        //     ReadedPercent+=5;
        //     cout<<ReadedPercent<<"%"<<" events finish"<<endl;
        // }

        // if (i % 5 == 0) // Create a status bar that will auto refresh
        // {
        //     std::string statusBar;
        //     std::cout << i << "th event is read";
        //     statusBar.assign((int)((float)i / nentries * 50.), '#');
        //     statusBar.resize(50, '.');
        //     std::cout << statusBar << "\r ";
        // }

    }
    std::cout<<std::endl;

    TFile *file = new TFile("readTree.root", "RECREATE");


	for (int i=0;i<HSize;i++){
		HMass[i]->Write();
		HBadMass[i]->Write();
		HP[i]->Write();
		HRapidity[i]->Write();
	    Hdcatopv[i]->Write();
        HnSigma[i]->Write();
	    H_dEdx_p[i]->Write();
	    H_sigma_TOF[i]->Write();
	    H_sigma_p[i]->Write();
	    H_sigma_pT[i]->Write();
	    H_Mass_DL[i]->Write();
	    H_Mass_dcatopv[i]->Write();
	}
	for (int i=0;i<3;i++){
        HC_sigma_TOF[i]->Write();
	}

    file->Write();

    return;
}
