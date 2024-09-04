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
const float omegaMass = 1.67245;
const float lambdaMass = 1.11568;
const float xiMass = 1.32171;
const float kaonMass = 0.493677;
const int maxTrack = 30000;
using namespace std;

void MixEvent_Test(Int_t Aid,Int_t Bid,float Energy,float B_RapMin,float B_RapMax)
{
    map<int, float> massList;
    massList.insert(pair<int, float>(321, 0.493677));//Kaon
    massList.insert(pair<int, float>(-321, 0.493677));

    massList.insert(pair<int, float>(3334, 1.67245));//Omega
    massList.insert(pair<int, float>(-3334, 1.67245));

    massList.insert(pair<int, float>(3312, 1.32171));//Xi^-
    massList.insert(pair<int, float>(-3312, 1.32171));

    massList.insert(pair<int, float>(3322, 1.3148));//Xi^0
    massList.insert(pair<int, float>(-3322, 1.3148));

    massList.insert(pair<int, float>(3122, 1.11568));//Lambda
    massList.insert(pair<int, float>(-3122, 1.11568));

    massList.insert(pair<int, float>(3222, 1.18937));//Sigma^+
    massList.insert(pair<int, float>(-3222, 1.18937));

    massList.insert(pair<int, float>(3112, 1.19745));//Sigma^-
    massList.insert(pair<int, float>(-3112, 1.19745));

    // Int_t Aid = -3122, Bid = 3334;
    const Int_t MultBin39[]   = {0, 8, 12, 16, 20, 25, 31, 42, 300};// 39
    const Int_t MultBin62[]   = {0, 9, 13, 17, 22, 27, 33, 43, 300};// 62
    const Int_t MultBin14p6[] = {398, 448, 498, 548, 598, 1000};// AuAu14.6
    const Int_t MultBinSize39 = sizeof(MultBin39)/sizeof(MultBin39[0]);
    const Int_t MultBinSize62 = sizeof(MultBin62)/sizeof(MultBin62[0]);
    const Int_t MultBinSize14p6 = sizeof(MultBin14p6)/sizeof(MultBin14p6[0]);
    Int_t MultBinSize,FileMax;
    std::vector<int> MultBin;
    if (fabs(Energy - 39.0)<=0.01) {
        MultBinSize = MultBinSize39;
        for(int i=0;i<MultBinSize39;i++){
            MultBin.push_back(MultBin39[i]);
        }
    }
    if (fabs(Energy - 62.0)<=0.01) {
        MultBinSize = MultBinSize62;
        for(int i=0;i<MultBinSize62;i++){
            MultBin.push_back(MultBin62[i]);
        }
    }
    if (fabs(Energy - 14.6)<=0.01) {
        MultBinSize = MultBinSize14p6;
        for(int i=0;i<MultBinSize14p6;i++){
            MultBin.push_back(MultBin14p6[i]);
        }
    }
    const Int_t A_NumBin[] = {0,1,2,3,4,5,6,200};
    const Int_t B_NumBin[] = {0,1,2,3,4,200};
    // const Int_t MultBin[] = {200, 225, 250, 275, 300, 350, 400};
    // const Int_t A_NumBin[] = {0,200};
    // const Int_t B_NumBin[] = {0,200};
    const Int_t A_NumBinSize = sizeof(A_NumBin)/sizeof(A_NumBin[0]);
    const Int_t B_NumBinSize = sizeof(B_NumBin)/sizeof(B_NumBin[0]);

    cout<<"Energy = "<<Energy<<endl;
    TString midname = "";
    if (fabs(Energy - 39.0)<=0.01) {
        midname = "/home/siyuanping/ampt_public/data39_String/Split/SplitFile_";
        // midname = "/home/siyuanping/ampt_public/data39_String/afterART7.7a_";
        FileMax = 9027;
    }
    if (fabs(Energy - 62.0)<=0.01) {
        midname = "/home/siyuanping/ampt_public/data62_String/Split/SplitFile_";
        // midname = "/home/siyuanping/ampt_public/data62_String/afterART7.7a_";
        FileMax = 6280;
    }
    if (fabs(Energy - 14.6)<=0.01) {
        midname = "/home/siyuanping/ampt_public/data14p6_AuAu_String/afterART7.7a_";
        FileMax = 202;
    }
    Int_t kBinNum = 1000;
    Float_t kmin = 0;
    Float_t kmax = 10;

    const Int_t maxMultiplicity = 20000;
    Int_t mult,npart,OriginMult,MultL;
    Int_t id[maxMultiplicity];
    Float_t px[maxMultiplicity],py[maxMultiplicity],pz[maxMultiplicity];
    Float_t x[maxMultiplicity],y[maxMultiplicity],z[maxMultiplicity],t[maxMultiplicity],mass[maxMultiplicity],energy[maxMultiplicity],b;

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
    for(int i=1;i <= FileMax;i++){//62 6280; 39 9027
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
    
    hadronTree->SetBranchAddress("nMultiplicityTree",&mult);
    hadronTree->SetBranchAddress("b",&b);
    hadronTree->SetBranchAddress("id",&id);
    hadronTree->SetBranchAddress("px",&px);
    hadronTree->SetBranchAddress("py",&py);
    hadronTree->SetBranchAddress("pz",&pz);
    hadronTree->SetBranchAddress("mass",&mass);
    hadronTree->SetBranchAddress("energy",&energy); 
    hadronTree->SetBranchAddress("t",t); 
    if ((fabs(Energy - 39.0)<=0.01) || (fabs(Energy - 62.0)<=0.01)) {
        hadronTree->SetBranchAddress("OriginMult",&OriginMult); 
        hadronTree->SetBranchAddress("MultL",&MultL); 
    }

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
        if ((i+1)%500 == 0) {
            cout<<"Calculating Event "<<(i+1)<<"/"<<nentries<<endl;
        }
        // if (b >= 1.7){
        //     continue;
        // }
        A_x.clear();B_x.clear();
        A_y.clear();B_y.clear();
        A_z.clear();B_z.clear();
        A_px.clear();B_px.clear();
        A_py.clear();B_py.clear();
        A_pz.clear();B_pz.clear();
        A_rap.clear();B_rap.clear();
        A_mass.clear();B_mass.clear();
        NumA = 0;NumB = 0;
        for (int j=0;j < mult;j++){
            
            p0.SetXYZM(px[j],py[j],pz[j],mass[j]);
            rap = p0.Rapidity();

            if(id[j]==Aid || id[j]==Bid){
                //&& rap > 0
                if(id[j]==Aid){
                    A_x.push_back(x[j]);
                    A_y.push_back(y[j]);
                    A_z.push_back(z[j]);
                    A_px.push_back(px[j]);
                    A_py.push_back(py[j]);
                    A_pz.push_back(pz[j]);
                    A_rap.push_back(rap);
                    A_mass.push_back(mass[j]);
                    NumA++;
                }
                else if(id[j]==Bid && rap > B_RapMin && rap < B_RapMax){
                    B_x.push_back(x[j]);
                    B_y.push_back(y[j]);
                    B_z.push_back(z[j]);
                    B_px.push_back(px[j]);
                    B_py.push_back(py[j]);
                    B_pz.push_back(pz[j]);
                    B_rap.push_back(rap);
                    B_mass.push_back(mass[j]);
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
            if ((fabs(Energy - 39.0)<=0.01) || (fabs(Energy - 62.0)<=0.01)) {
                if (MultBin[j-1] <= MultL && MultL < MultBin[j]){
                    BP[0] = j-1;
                    break;
                }
            }
            else{
                if (MultBin[j-1] <= mult && mult < MultBin[j]){
                    BP[0] = j-1;
                    break;
                }
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

    TFile *file = new TFile("Cor.root", "RECREATE");

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
