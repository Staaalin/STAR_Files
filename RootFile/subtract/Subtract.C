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

void Subtract(TString MidName,int StartFileIndex,int EndFileIndex,int OutputFileIndex,
                   float LambdaFitMass,float LambdaFitMassSigma,
                   float XiFitMass,float XiFitMassSigma,
                   float OmegaFitMass,float OmegaFitMassSigma,
                   float LambdaBarFitMass,float LambdaBarFitMassSigma,
                   float XiBarFitMass,float XiBarFitMassSigma,
                   float OmegaBarFitMass,float OmegaBarFitMassSigma)
{
    // TString FileName = "output_";
    cout<<"Start Running"<<endl;
    cout<<StartFileIndex<<endl;
    cout<<EndFileIndex<<endl;

    map<int, float> massList;
    map<int, float> massListSigma;
    massList.insert(pair<int, float>(321, 0.493677));//Kaon
    massList.insert(pair<int, float>(-321, 0.493677));
    massList.insert(pair<int, float>(311, 0.49765));//Kaon ^0
    massList.insert(pair<int, float>(-311, 0.49765));
    massList.insert(pair<int, float>(211, 0.13957));//Pion
    massList.insert(pair<int, float>(-211, 0.13957));

    // massList.insert(pair<int, float>(3334, 1.67245));//Omega
    // massList.insert(pair<int, float>(-3334, 1.67245));
    massList.insert(pair<int, float>(3334, OmegaFitMass));//Omega
    massList.insert(pair<int, float>(-3334, OmegaBarFitMass));
    massListSigma.insert(pair<int, float>(3334, OmegaFitMassSigma));//Omega
    massListSigma.insert(pair<int, float>(-3334, OmegaBarFitMassSigma));

    // massList.insert(pair<int, float>(3312, 1.32171));//Xi^-
    // massList.insert(pair<int, float>(-3312, 1.32171));
    massList.insert(pair<int, float>(3312, XiFitMass));//Xi^-
    massList.insert(pair<int, float>(-3312, XiBarFitMass));
    massListSigma.insert(pair<int, float>(3312, XiFitMassSigma));//Xi^-
    massListSigma.insert(pair<int, float>(-3312, XiBarFitMassSigma));

    massList.insert(pair<int, float>(3322, 1.3148));//Xi^0
    massList.insert(pair<int, float>(-3322, 1.3148));

    // massList.insert(pair<int, float>(3122, 1.11568));//Lambda
    // massList.insert(pair<int, float>(-3122, 1.11568));
    massList.insert(pair<int, float>(3122, LambdaFitMass));//Lambda
    massList.insert(pair<int, float>(-3122, LambdaBarFitMass));
    massListSigma.insert(pair<int, float>(3122, LambdaFitMassSigma));//Lambda
    massListSigma.insert(pair<int, float>(-3122, LambdaBarFitMassSigma));

    massList.insert(pair<int, float>(3222, 1.18937));//Sigma^+
    massList.insert(pair<int, float>(-3222, 1.18937));

    massList.insert(pair<int, float>(3112, 1.19745));//Sigma^-
    massList.insert(pair<int, float>(-3112, 1.19745));

    // FileMax = 500;
    // Int_t kBinNum = 200;
    // Float_t kmin = 0;
    // Float_t kmax = 20;

    const Int_t maxMultiplicity = 20000;
    Int_t mult,npart;
    Int_t id[maxMultiplicity];
    Float_t px[maxMultiplicity],py[maxMultiplicity],pz[maxMultiplicity];
    Float_t nSigmaPion[maxMultiplicity],nSigmaKaon[maxMultiplicity],z[maxMultiplicity],nSigmaProton[maxMultiplicity],mass[maxMultiplicity],energy[maxMultiplicity],b;

    #define HyperonPhaseNum  2
    #define MesonPhaseNum    2
    TString HyperonPhase[HyperonPhaseNum] = {"m1p1Rap","0p1Rap"}
    TString MesonPhase[MesonPhaseNum] = {"Allp","0p5t2p"}
    TH1D *dNdy_Kaonp  [MesonPhaseNum];
    TH1D *dNdy_Kaonm  [MesonPhaseNum];
    TH1D *dNdy_Pionp  [MesonPhaseNum];
    TH1D *dNdy_Pionm  [MesonPhaseNum];
    TH1D *dNdy_Lambda [HyperonPhaseNum];
    TH1D *dNdy_Lambdab[HyperonPhaseNum];
    TH1D *dNdy_Omega  [HyperonPhaseNum];
    TH1D *dNdy_Omegab [HyperonPhaseNum];
    TH1D *dNdy_Xi     [HyperonPhaseNum];
    TH1D *dNdy_Xib    [HyperonPhaseNum];
    TH1D *k_Kp_Omega  [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Omegab [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Omega  [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Omegab [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Xi     [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Xib    [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Xi     [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Xib    [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Lambda [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Kp_Lambdab[MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Lambda [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Km_Lambdab[MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Omega  [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Omegab [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Omega  [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Omegab [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Xi     [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Xib    [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Xi     [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Xib    [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Lambda [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pip_Lambdab[MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Lambda [MesonPhaseNum][HyperonPhaseNum];
    TH1D *k_Pim_Lambdab[MesonPhaseNum][HyperonPhaseNum];

    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
        TString HistName = "dNdy_Kaonp_";HistName += MesonPhase[Itr];
        dNdy_Kaonp[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Kaonp[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Kaonp[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Kaonm_";HistName += MesonPhase[Itr];
        dNdy_Kaonm[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Kaonm[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Kaonm[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Pionp_";HistName += MesonPhase[Itr];
        dNdy_Pionp[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Pionp[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Pionp[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Pionm_";HistName += MesonPhase[Itr];
        dNdy_Pionm[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Pionm[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Pionm[Itr]->GetYaxis()->SetTitle("Counts");
    }
    for (int Itr = 0;Itr < HyperonPhaseNum;Itr++) {
        TString HistName = "dNdy_Lambda_";HistName += HyperonPhase[Itr];
        dNdy_Lambda[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Lambda[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Lambda[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Lambdab_";HistName += HyperonPhase[Itr];
        dNdy_Lambdab[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Lambdab[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Lambdab[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Omega_";HistName += HyperonPhase[Itr];
        dNdy_Omega[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Omega[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Omega[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Omegab_";HistName += HyperonPhase[Itr];
        dNdy_Omegab[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Omegab[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Omegab[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Xi_";HistName += HyperonPhase[Itr];
        dNdy_Xi[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Xi[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Xi[Itr]->GetYaxis()->SetTitle("Counts");
        HistName = "dNdy_Xib_";HistName += HyperonPhase[Itr];
        dNdy_Xib[Itr] = new TH1D(HistName,HistName,40,-2,2);
        dNdy_Xib[Itr]->GetXaxis()->SetTitle("y");
        dNdy_Xib[Itr]->GetYaxis()->SetTitle("Counts");
    }
    for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
        for (int Jtr = 0;Jtr < HyperonPhaseNum;Jtr++) {
            TString HistName = "k_Kp_Omega_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Omega  [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Omegab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Omegab [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Omega_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Omega  [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Omegab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Omegab [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Xi_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Xi     [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Xib_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Xib    [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Xi_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Xi     [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Xib_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Xib    [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Lambda_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Lambda [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Kp_Lambdab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Kp_Lambdab[Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Lambda_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Lambda [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Km_Lambdab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Km_Lambdab[Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);

            HistName = "k_Pip_Omega_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Omega  [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Omegab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Omegab [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Omega_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Omega  [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Omegab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Omegab [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Xi_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Xi     [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Xib_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Xib    [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Xi_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Xi     [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Xib_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Xib    [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Lambda_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Lambda [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pip_Lambdab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pip_Lambdab[Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Lambda_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Lambda [Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);
            HistName = "k_Pim_Lambdab_";HistName += MesonPhase[Itr];HistName += "_";HistName += HyperonPhase[Jtr];
            k_Pim_Lambdab[Itr][Jtr] = new TH1D(HistName,HistName,100,0,10);

            k_Kp_Omega   [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Omegab  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Omega   [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Omegab  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Xi      [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Xib     [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Xi      [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Xib     [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Lambda  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Lambdab [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Lambda  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Km_Lambdab [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Omega  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Omegab [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Omega  [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Omegab [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Xi     [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Xib    [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Xi     [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Xib    [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Lambda [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pip_Lambdab[Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Lambda [Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Pim_Lambdab[Itr][Jtr]->GetXaxis()->SetTitle("k* [GeV]");
            k_Kp_Omega   [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Omegab  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Omega   [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Omegab  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Xi      [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Xib     [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Xi      [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Xib     [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Lambda  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Kp_Lambdab [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Lambda  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Km_Lambdab [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Omega  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Omegab [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Omega  [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Omegab [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Xi     [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Xib    [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Xi     [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Xib    [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Lambda [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pip_Lambdab[Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Lambda [Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
            k_Pim_Lambdab[Itr][Jtr]->GetYaxis()->SetTitle("dN/dk*");
        }
    }


    //load data  
    TChain *hadronTree = new TChain("hadronTree");
    for(int i=StartFileIndex;i <= EndFileIndex;i++){
        TString filename = MidName;
        filename+=i;
        filename+=".root";
        hadronTree->Add(filename);
        // cout<<"Add "<<filename<<" Successfully"<<endl;
    }
    
    hadronTree->SetBranchAddress("PDGMult",&mult);
    hadronTree->SetBranchAddress("PDG",&id);
    hadronTree->SetBranchAddress("mix_px",&px);
    hadronTree->SetBranchAddress("mix_py",&py);
    hadronTree->SetBranchAddress("mix_pz",&pz);
    hadronTree->SetBranchAddress("InvariantMass",&mass);
    hadronTree->SetBranchAddress("nSigmaProton",&nSigmaProton);
    hadronTree->SetBranchAddress("nSigmaPion",&nSigmaPion);
    hadronTree->SetBranchAddress("nSigmaKaon",&nSigmaKaon);

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;
    
    TLorentzVector p0,p3;
    float kstar, rap;
    TVector3 BetaTemp;

    std::vector<float> Kp_px[MesonPhaseNum],Km_px[MesonPhaseNum],Pip_px[MesonPhaseNum],Pim_px[MesonPhaseNum];
    std::vector<float> Kp_py[MesonPhaseNum],Km_py[MesonPhaseNum],Pip_py[MesonPhaseNum],Pim_py[MesonPhaseNum];
    std::vector<float> Kp_pz[MesonPhaseNum],Km_pz[MesonPhaseNum],Pip_pz[MesonPhaseNum],Pim_pz[MesonPhaseNum];
    std::vector<float> Kp_y [MesonPhaseNum],Km_y [MesonPhaseNum],Pip_y [MesonPhaseNum],Pim_y [MesonPhaseNum];
    std::vector<float> Lambda_px[HyperonPhaseNum],Lambdab_px[HyperonPhaseNum],Omega_px[HyperonPhaseNum],Omegab_px[HyperonPhaseNum],Xi_px[HyperonPhaseNum],Xib_px[HyperonPhaseNum];
    std::vector<float> Lambda_py[HyperonPhaseNum],Lambdab_py[HyperonPhaseNum],Omega_py[HyperonPhaseNum],Omegab_py[HyperonPhaseNum],Xi_py[HyperonPhaseNum],Xib_py[HyperonPhaseNum];
    std::vector<float> Lambda_pz[HyperonPhaseNum],Lambdab_pz[HyperonPhaseNum],Omega_pz[HyperonPhaseNum],Omegab_pz[HyperonPhaseNum],Xi_pz[HyperonPhaseNum],Xib_pz[HyperonPhaseNum];
    std::vector<float> Lambda_y [HyperonPhaseNum],Lambdab_y [HyperonPhaseNum],Omega_y [HyperonPhaseNum],Omegab_y [HyperonPhaseNum],Xi_y [HyperonPhaseNum],Xib_y [HyperonPhaseNum];

    int KaonpPID = 321,KaonmPID = -321,PionpPID = 211,PionmPID = -211,LambdaPID = 3122,LambdabPID = -3122,XiPID = 3312,XibPID = -3312,OmegaPID = 3334,OmegabPID = -3334;

    TString States;

    for (int i=0;i<nentries;i++){
        hadronTree->GetEntry(i);
        // if(b>7){continue;}
        Kp_px.resize(0);Km_px.resize(0);Pip_px.resize(0);Pim_px.resize(0);
        Kp_py.resize(0);Km_py.resize(0);Pip_py.resize(0);Pim_py.resize(0);
        Kp_pz.resize(0);Km_pz.resize(0);Pip_pz.resize(0);Pim_pz.resize(0);
        Kp_y .resize(0);Km_y .resize(0);Pip_y .resize(0);Pim_y .resize(0);
        Lambda_px.resize(0);Lambdab_px.resize(0);Omega_px.resize(0);Omegab_px.resize(0);Xi_px.resize(0);Xib_px.resize(0);
        Lambda_py.resize(0);Lambdab_py.resize(0);Omega_py.resize(0);Omegab_py.resize(0);Xi_py.resize(0);Xib_py.resize(0);
        Lambda_pz.resize(0);Lambdab_pz.resize(0);Omega_pz.resize(0);Omegab_pz.resize(0);Xi_pz.resize(0);Xib_pz.resize(0);
        Lambda_y .resize(0);Lambdab_y .resize(0);Omega_y .resize(0);Omegab_y .resize(0);Xi_y .resize(0);Xib_y .resize(0);
        for (int j=0;j < mult;j++){
            if ((fabs(id[j]) == LambdaPID) || (fabs(id[j]) == XiPID || (fabs(id[j]) == OmegaPID))){
                if ((fabs(mass[j] - massList[id[j]])) > 3*massListSigma[id[j]]) {
                    continue;
                }
                p0.SetXYZM(px[j],py[j],pz[j],massList[id[j]]);
                rap = p0.Rapidity();
                if (-1 > rap || rap > 1) {
                    continue;
                }

                if (rap > 0) {
                    States = "0p1Rap";
                }else{
                    States = "m1p1Rap";
                }
                int HyperonPhaseIndex;
                for (int Itr = 0;Itr < HyperonPhaseNum;Itr++) {
                    if (HyperonPhase[Itr] == States){
                        HyperonPhaseIndex = Itr;
                        break;
                    }
                }

                if (id[j] = LambdaPID) {
                    Lambda_px[HyperonPhaseIndex].push_back(px[j]);
                    Lambda_py[HyperonPhaseIndex].push_back(py[j]);
                    Lambda_pz[HyperonPhaseIndex].push_back(pz[j]);
                    Lambda_y [HyperonPhaseIndex].push_back(rap);
                    dNdy_Lambda[HyperonPhaseIndex]->Fill(rap);
                }
                else if (id[j] = LambdabPID) {
                    Lambdab_px[HyperonPhaseIndex].push_back(px[j]);
                    Lambdab_py[HyperonPhaseIndex].push_back(py[j]);
                    Lambdab_pz[HyperonPhaseIndex].push_back(pz[j]);
                    Lambdab_y [HyperonPhaseIndex].push_back(rap);
                    dNdy_Lambdab[HyperonPhaseIndex]->Fill(rap);
                }
                else if (id[j] = XiPID) {
                    Xi_px[HyperonPhaseIndex].push_back(px[j]);
                    Xi_py[HyperonPhaseIndex].push_back(py[j]);
                    Xi_pz[HyperonPhaseIndex].push_back(pz[j]);
                    Xi_y [HyperonPhaseIndex].push_back(rap);
                    dNdy_Xi[HyperonPhaseIndex]->Fill(rap);
                }
                else if (id[j] = XibPID) {
                    Xib_px[HyperonPhaseIndex].push_back(px[j]);
                    Xib_py[HyperonPhaseIndex].push_back(py[j]);
                    Xib_pz[HyperonPhaseIndex].push_back(pz[j]);
                    Xib_y [HyperonPhaseIndex].push_back(rap);
                    dNdy_Xib[HyperonPhaseIndex]->Fill(rap);
                }
                else if (id[j] = OmegaPID) {
                    Omega_px[HyperonPhaseIndex].push_back(px[j]);
                    Omega_py[HyperonPhaseIndex].push_back(py[j]);
                    Omega_pz[HyperonPhaseIndex].push_back(pz[j]);
                    Omega_y [HyperonPhaseIndex].push_back(rap);
                    dNdy_Omega[HyperonPhaseIndex]->Fill(rap);
                }
                else if (id[j] = OmegabPID) {
                    Omegab_px[HyperonPhaseIndex].push_back(px[j]);
                    Omegab_py[HyperonPhaseIndex].push_back(py[j]);
                    Omegab_pz[HyperonPhaseIndex].push_back(pz[j]);
                    Omegab_y [HyperonPhaseIndex].push_back(rap);
                    dNdy_Omegab[HyperonPhaseIndex]->Fill(rap);
                }
            }

            if ((fabs(id[j]) == KaonpPID) || (fabs(id[j]) == PionpPID)) {
                if ((fabs(id[j]) == KaonpPID) && ((nSigmaKaon[j] < -2) || (nSigmaKaon[j] > 2))) {
                    continue;
                }
                if ((fabs(id[j]) == PionpPID) && ((nSigmaPion[j] < -3) || (nSigmaPion[j] > 3))) {
                    continue;
                }
                p0.SetXYZM(px[j],py[j],pz[j],massList[id[j]]);
                rap = p0.Rapidity();
                float Meg = pow(pow(px[j],2)+pow(py[j],2)+pow(pz[j],2),0.5);
                if ((0.5 < Meg) && (Meg < 2)) {
                    States = "0p5t2p";
                }else{
                    States = "Allp";
                }
                int MesonPhaseIndex;
                for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
                    if (MesonPhase[Itr] == States){
                        MesonPhaseIndex = Itr;
                        break;
                    }
                }
                
                if (id[j] == PionpPID) {
                    Pip_px[MesonPhaseIndex].push_back(px[j]);
                    Pip_py[MesonPhaseIndex].push_back(py[j]);
                    Pip_pz[MesonPhaseIndex].push_back(pz[j]);
                    Pip_y [MesonPhaseIndex].push_back(rap);
                    dNdy_Pionp[MesonPhaseIndex]->Fill(rap);
                }
                else if (id[j] == PionmPID) {
                    Pim_px[MesonPhaseIndex].push_back(px[j]);
                    Pim_py[MesonPhaseIndex].push_back(py[j]);
                    Pim_pz[MesonPhaseIndex].push_back(pz[j]);
                    Pim_y [MesonPhaseIndex].push_back(rap);
                    dNdy_Pionm[MesonPhaseIndex]->Fill(rap);
                }
                else if (id[j] == KaonpPID) {
                    Kp_px[MesonPhaseIndex].push_back(px[j]);
                    Kp_py[MesonPhaseIndex].push_back(py[j]);
                    Kp_pz[MesonPhaseIndex].push_back(pz[j]);
                    Kp_y [MesonPhaseIndex].push_back(rap);
                    dNdy_Kaonp[MesonPhaseIndex]->Fill(rap);
                }
                else if (id[j] == KaonmPID) {
                    Km_px[MesonPhaseIndex].push_back(px[j]);
                    Km_py[MesonPhaseIndex].push_back(py[j]);
                    Km_pz[MesonPhaseIndex].push_back(pz[j]);
                    Km_y [MesonPhaseIndex].push_back(rap);
                    dNdy_Kaomm[MesonPhaseIndex]->Fill(rap);
                }
            }
        }

        for (int Itr = 0;Itr < MesonPhaseNum;Itr++) {
            for (int Jtr = 0;Jtr < HyperonPhaseNum;Jtr++) {
                for (int Ktr = 0;Ktr < Kp_px.size();Ktr++) {
                    TLorentzVector p1;
                    p1.SetXYZM(Kp_px[Ktr],Kp_py[Ktr],Kp_pz[Ktr],massList[KaonpPID]);
                    // float ALength = pow(pow(A_px[j],2)+pow(A_py[j],2)+pow(A_pz[j],2) ,0.5);
                    for (int Ntr = 0;Ntr < Lambda_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambda_px[Ntr],Lambda_py[Ntr],Lambda_pz[Ntr],massList[LambdaPID]);

                        p3 = p4 + p2;
                        // float dpt = fabs(p1.Perp() - p2.Perp());
                        // Apt->Fill(dpt);
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Lambda[Itr][Jtr]->Fill(kstar);
                        // Ay->Fill(fabs(A_y[j] - B_y[k]));

                        // float BLength = pow(pow(B_px[k],2)+pow(B_py[k],2)+pow(B_pz[k],2) ,0.5);
                        // float CosAB = (A_px[j]*B_px[k] + A_py[j]*B_py[k] + A_pz[j]*B_pz[k])/(ALength*BLength);
                        // Aphi->Fill(acos(CosAB));
                    }
                    for (int Ntr = 0;Ntr < Lambdab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambdab_px[Ntr],Lambdab_py[Ntr],Lambdab_pz[Ntr],massList[LambdabPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Lambdab[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xi_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xi_px[Ntr],Xi_py[Ntr],Xi_pz[Ntr],massList[XiPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Xi[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xib_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xib_px[Ntr],Xib_py[Ntr],Xib_pz[Ntr],massList[XibPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Xib[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omega_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omega_px[Ntr],Omega_py[Ntr],Omega_pz[Ntr],massList[OmegaPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Omega[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omegab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omegab_px[Ntr],Omegab_py[Ntr],Omegab_pz[Ntr],massList[OmegabPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Kp_Omegab[Itr][Jtr]->Fill(kstar);
                    }
                }
                for (int Ktr = 0;Ktr < Km_px.size();Ktr++) {
                    TLorentzVector p1;
                    p1.SetXYZM(Km_px[Ktr],Km_py[Ktr],Km_pz[Ktr],massList[KaonmPID]);
                    for (int Ntr = 0;Ntr < Lambda_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambda_px[Ntr],Lambda_py[Ntr],Lambda_pz[Ntr],massList[LambdaPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Lambda[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Lambdab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambdab_px[Ntr],Lambdab_py[Ntr],Lambdab_pz[Ntr],massList[LambdabPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Lambdab[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xi_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xi_px[Ntr],Xi_py[Ntr],Xi_pz[Ntr],massList[XiPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Xi[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xib_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xib_px[Ntr],Xib_py[Ntr],Xib_pz[Ntr],massList[XibPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Xib[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omega_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omega_px[Ntr],Omega_py[Ntr],Omega_pz[Ntr],massList[OmegaPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Omega[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omegab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omegab_px[Ntr],Omegab_py[Ntr],Omegab_pz[Ntr],massList[OmegabPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Km_Omegab[Itr][Jtr]->Fill(kstar);
                    }
                }
                for (int Ktr = 0;Ktr < Pip_px.size();Ktr++) {
                    TLorentzVector p1;
                    p1.SetXYZM(Pip_px[Ktr],Pip_py[Ktr],Pip_pz[Ktr],massList[PionpPID]);
                    for (int Ntr = 0;Ntr < Lambda_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambda_px[Ntr],Lambda_py[Ntr],Lambda_pz[Ntr],massList[LambdaPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Lambda[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Lambdab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambdab_px[Ntr],Lambdab_py[Ntr],Lambdab_pz[Ntr],massList[LambdabPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Lambdab[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xi_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xi_px[Ntr],Xi_py[Ntr],Xi_pz[Ntr],massList[XiPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Xi[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xib_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xib_px[Ntr],Xib_py[Ntr],Xib_pz[Ntr],massList[XibPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Xib[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omega_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omega_px[Ntr],Omega_py[Ntr],Omega_pz[Ntr],massList[OmegaPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Omega[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omegab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omegab_px[Ntr],Omegab_py[Ntr],Omegab_pz[Ntr],massList[OmegabPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pip_Omegab[Itr][Jtr]->Fill(kstar);
                    }
                }
                for (int Ktr = 0;Ktr < Pim_px.size();Ktr++) {
                    TLorentzVector p1;
                    p1.SetXYZM(Pim_px[Ktr],Pim_py[Ktr],Pim_pz[Ktr],massList[PionmPID]);
                    for (int Ntr = 0;Ntr < Lambda_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambda_px[Ntr],Lambda_py[Ntr],Lambda_pz[Ntr],massList[LambdaPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Lambda[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Lambdab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Lambdab_px[Ntr],Lambdab_py[Ntr],Lambdab_pz[Ntr],massList[LambdabPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Lambdab[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xi_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xi_px[Ntr],Xi_py[Ntr],Xi_pz[Ntr],massList[XiPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Xi[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Xib_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Xib_px[Ntr],Xib_py[Ntr],Xib_pz[Ntr],massList[XibPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Xib[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omega_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omega_px[Ntr],Omega_py[Ntr],Omega_pz[Ntr],massList[OmegaPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Omega[Itr][Jtr]->Fill(kstar);
                    }
                    for (int Ntr = 0;Ntr < Omegab_px.size();Ntr++) {
                        TLorentzVector p2,p4 = p1;
                        p2.SetXYZM(Omegab_px[Ntr],Omegab_py[Ntr],Omegab_pz[Ntr],massList[OmegabPID]);

                        p3 = p4 + p2;
                        p4.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                        kstar = 0.5 * (p4 - p2).Rho();
                        k_Pim_Omegab[Itr][Jtr]->Fill(kstar);
                    }
                }
            }
        }

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
