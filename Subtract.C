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

void pAcorSubtract(float Energy, 
                   Int_t Aid, float A_rap_Min, float A_rap_Max, 
                   Int_t Bid, float B_rap_Min, float B_rap_Max, 
                   Int_t Cid, float C_rap_Min, float C_rap_Max, 
                   Int_t Did, float D_rap_Min, float D_rap_Max, 
                   Int_t kBinNum, Float_t kmin, Float_t kmax)
{
    map<int, float> massList;
    massList.insert(pair<int, float>(321, 0.493677));//Kaon
    massList.insert(pair<int, float>(-321, 0.493677));
    massList.insert(pair<int, float>(311, 0.49765));//Kaon ^0
    massList.insert(pair<int, float>(-311, 0.49765));

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

    // Int_t Aid = 321, Bid = -3334, Cid = 321, Did = -3334;
    const float AMass = massList[Aid], BMass = massList[Bid], CMass = massList[Cid], DMass = massList[Did];
    int FileMax = 0;
    TString midname = "";
    if (Energy == 39.0){//62 6280; 39 9027
        midname = "/home/siyuanping/ampt_public/data39_String/Split/SplitFile_";
        // TString midname = "/home/siyuanping/ampt_public/data39_String/afterART7.7a_";
        FileMax = 9027;
    }
    if (Energy == 62.0){//62 6280; 39 9027
        midname = "/home/siyuanping/ampt_public/data62_String/Split/SplitFile_";
        // TString midname = "/home/siyuanping/ampt_public/data62_String/afterART7.7a_";
        FileMax = 6280;
    }
    // FileMax = 500;
    // Int_t kBinNum = 200;
    // Float_t kmin = 0;
    // Float_t kmax = 20;

    const Int_t maxMultiplicity = 20000;
    Int_t mult,npart;
    Int_t id[maxMultiplicity];
    Float_t px[maxMultiplicity],py[maxMultiplicity],pz[maxMultiplicity];
    Float_t x[maxMultiplicity],y[maxMultiplicity],z[maxMultiplicity],t[maxMultiplicity],mass[maxMultiplicity],energy[maxMultiplicity],b;

    TH1::SetDefaultSumw2("kTRUE");
    TH1D *A_dNdy = new TH1D("A_dNdy", "A_dNdy number", 2*kBinNum, -kmax, kmax);
    TH1D *B_dNdy = new TH1D("B_dNdy", "B_dNdy number", 2*kBinNum, -kmax, kmax);
    TH1D *C_dNdy = new TH1D("C_dNdy", "C_dNdy number", 2*kBinNum, -kmax, kmax);
    TH1D *D_dNdy = new TH1D("D_dNdy", "D_dNdy number", 2*kBinNum, -kmax, kmax);
    TH1D *Ak = new TH1D("Ak", "Ak number", kBinNum, kmin, kmax);
    TH1D *Bk = new TH1D("Bk", "Bk number", kBinNum, kmin, kmax);
    TH1D *Apt = new TH1D("Apt", "Apt number", kBinNum, kmin, kmax);
    TH1D *Bpt = new TH1D("Bpt", "Bpt number", kBinNum, kmin, kmax);
    TH1D *Ay = new TH1D("Ay", "Ay number", kBinNum, kmin, kmax);
    TH1D *By = new TH1D("By", "By number", kBinNum, kmin, kmax);
    TH1D *Aphi = new TH1D("Aphi", "Aphi number", 50, 0.0, 3.1415926535897);
    TH1D *Bphi = new TH1D("Bphi", "Bphi number", 50, 0.0, 3.1415926535897);
    TH1D *testPx = new TH1D("testPx", "testPx number", kBinNum, kmin, kmax);
    TH1D *testPy = new TH1D("testPy", "testPy number", kBinNum, kmin, kmax);
    TH1D *testPz = new TH1D("testPz", "testPz number", kBinNum, kmin, kmax);

    //load data  
    TChain *hadronTree = new TChain("hadronTree");
    for(int i=1;i <= FileMax;i++){
        TString filename = midname;
        filename+=i;
        filename+=".root";
        hadronTree->Add(filename);
    }
    
    hadronTree->SetBranchAddress("nMultiplicityTree",&mult);
    hadronTree->SetBranchAddress("b",&b);
    hadronTree->SetBranchAddress("id",&id);
    hadronTree->SetBranchAddress("px",&px);
    hadronTree->SetBranchAddress("py",&py);
    hadronTree->SetBranchAddress("pz",&pz);
    hadronTree->SetBranchAddress("energy",&energy); 
    hadronTree->SetBranchAddress("t",t); 

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;
    
    TLorentzVector p0,p3;
    float kstar, rap;
    TVector3 BetaTemp;
    std::vector<float> A_px, B_px, C_px, D_px;
    std::vector<float> A_py, B_py, C_py, D_py;
    std::vector<float> A_pz, B_pz, C_pz, D_pz;
    std::vector<float> A_y,  B_y,  C_y,  D_y; 

    Int_t NumA = 0, NumB = 0, NumC = 0, NumD = 0, SumNum = 0, SumNumA = 0, SumNumB = 0;//AB-CD

    Int_t num_of_counted_Bid = 0, num_of_counted_Did = 0;

    Int_t num_net_3122 = 0;//Lambda
    Int_t num_net_3222 = 0;//Sigma+
    Int_t num_net_3112 = 0;//Sigma-
    Int_t num_net_3322 = 0;//Xi0
    Int_t num_net_3312 = 0;//Xi-

    Int_t num_3122 = 0;//Lambda
    Int_t num_3222 = 0;//Sigma+
    Int_t num_3112 = 0;//Sigma-
    Int_t num_3322 = 0;//Xi0
    Int_t num_3312 = 0;//Xi-
    Int_t num_m3122 = 0;//Lambda_bar
    Int_t num_m3222 = 0;//Sigma+_bar
    Int_t num_m3112 = 0;//Sigma-_bar
    Int_t num_m3322 = 0;//Xi0_bar
    Int_t num_m3312 = 0;//Xi-_bar



    for (int i=0;i<nentries;i++){
        Int_t num_of_counted_Omega_in_one_event = 0, num_of_counted_Omega_bar_in_one_event = 0;
        hadronTree->GetEntry(i);
        // if(b>7){continue;}
        num_net_3122 = 0;num_net_3222 = 0;num_net_3112 = 0;num_net_3322 = 0;num_net_3312 = 0;
        num_3122 = 0;num_3222 = 0;num_3112 = 0;num_3322 = 0;num_3312 = 0;
        num_m3122 = 0;num_m3222 = 0;num_m3112 = 0;num_m3322 = 0;num_m3312 = 0;
        A_px.clear();B_px.clear();C_px.clear();D_px.clear();
        A_py.clear();B_py.clear();C_py.clear();D_py.clear();
        A_pz.clear();B_pz.clear();C_pz.clear();D_pz.clear();
        A_y.clear(); B_y.clear(); C_y.clear(); D_y.clear(); 
        NumA = 0;NumB = 0;NumC = 0;NumD = 0;
        for (int j=0;j < mult;j++){
            if(id[j]==3122){
                num_net_3122++;
                num_3122++;
            }
            if(id[j]==-3122){
                num_net_3122--;
                num_m3122++;
            }
            if(id[j]==3222){
                num_net_3222++;
                num_3222++;
            }
            if(id[j]==-3222){
                num_net_3222--;
                num_m3222++;
            }
            if(id[j]==3112){
                num_net_3112++;
                num_3112++;
            }
            if(id[j]==-3112){
                num_net_3112--;
                num_m3112++;
            }
            if(id[j]==3322){
                num_net_3322++;
                num_3322++;
            }
            if(id[j]==-3322){
                num_net_3322--;
                num_m3322++;
            }
            if(id[j]==3312){
                num_net_3312++;
                num_3312++;
            }
            if(id[j]==-3312){
                num_net_3312--;
                num_m3312++;
            }
            /////////////////////////////////////////////////
            if(id[j] ==  3334){
                num_of_counted_Omega_in_one_event++;
            }
            if(id[j] == -3334){
                num_of_counted_Omega_bar_in_one_event++;
                //break;//**
            }
            if(num_of_counted_Omega_in_one_event > 1 || num_of_counted_Omega_bar_in_one_event > 1){
                NumA = 0;NumB = 0;NumC = 0;NumD = 0;
                break;
            }
            if(1 <= num_of_counted_Omega_in_one_event && 1 <= num_of_counted_Omega_bar_in_one_event){
                NumA = 0;NumB = 0;NumC = 0;NumD = 0;
                break;
            }

            if(id[j]==Aid || id[j]==Bid || id[j]==Cid || id[j]==Did){
                p0.SetXYZM(px[j],py[j],pz[j],mass[j]);
                rap = p0.Rapidity();

                //&& rap > 0
                if(id[j]==Aid && A_rap_Min <= rap && rap <= A_rap_Max){
                    A_px.push_back(px[j]);
                    A_py.push_back(py[j]);
                    A_pz.push_back(pz[j]);
                    A_y.push_back(rap);
                    NumA++;
                }
                else if(id[j]==Bid && B_rap_Min <= rap && rap <= B_rap_Max){
                    B_px.push_back(px[j]);
                    B_py.push_back(py[j]);
                    B_pz.push_back(pz[j]);
                    B_y.push_back(rap);
                    NumB++;
                }
                if(id[j]==Cid && C_rap_Min <= rap && rap <= C_rap_Max){
                    C_px.push_back(px[j]);
                    C_py.push_back(py[j]);
                    C_pz.push_back(pz[j]);
                    C_y.push_back(rap);
                    NumC++;
                }
                else if(id[j]==Did && D_rap_Min <= rap && rap <= D_rap_Max){
                    D_px.push_back(px[j]);
                    D_py.push_back(py[j]);
                    D_pz.push_back(pz[j]);
                    D_y.push_back(rap);
                    NumD++;
                }
            }
        }
        ////////////////////////////
        // Int_t num_net_s = (num_net_3122 + num_net_3222 + num_net_3112) + 2*(num_net_3322 + num_net_3312);
        // bool Check_Point = false;
        // if(Bid == 3334 && (-1 == num_net_s || -2 == num_net_s)){
        //     Check_Point = true;
        // }
        // if(Bid == -3334 && (1 == num_net_s || 2 == num_net_s)){
        //     Check_Point = true;
        // }
        // if(Check_Point == false){
        //     continue;
        // }
        ////////////////////////////
        // bool Check_Point = false;
        // if(Bid == 3334 && num_m3322 > 0){
        //     Check_Point = true;
        // }
        // if(Bid == -3334 && num_3322 > 0){
        //     Check_Point = true;
        // }
        // if(Bid == 3334 && num_m3312 > 0){
        //     Check_Point = true;
        // }
        // if(Bid == -3334 && num_3312 > 0){
        //     Check_Point = true;
        // }
        // if(Check_Point == false){
        //     continue;
        // }
        ////////////////////////////
        // bool Check_Point = false;
        // if(Bid == 3334 && num_m3122 > 0){
        //     Check_Point = true;
        // }
        // if(Bid == -3334 && num_3122 > 0){
        //     Check_Point = true;
        // }
        // if(Check_Point == false){
        //     continue;
        // }
        ////////////////////////////
        // bool Check_Point = false;
        // if(Bid == 3334 && num_m3322 > 0){
        //     Check_Point = true;
        // }
        // if(Bid == -3334 && num_3322 > 0){
        //     Check_Point = true;
        // }
        // if(Bid == 3334 && num_m3312 > 0){
        //     Check_Point = true;
        // }
        // if(Bid == -3334 && num_3312 > 0){
        //     Check_Point = true;
        // }
        // if(Check_Point == false){
        //     continue;
        // }
        ////////////////////////////

        // Int_t num_net_s = (num_net_3122 + num_net_3222 + num_net_3112) + 2*(num_net_3322 + num_net_3312);
        // bool Check_Point = false;
        // if(Bid == 3334 && (0 == num_net_s)){
        //     Check_Point = true;
        // }
        // if(Bid == -3334 && (0 == num_net_s)){
        //     Check_Point = true;
        // }
        // if(Check_Point == false){
        //     continue;
        // }
        ////////////////////////////

        // Int_t num_net_s = (num_net_3122 + num_net_3222 + num_net_3112) + 2*(num_net_3322 + num_net_3312);
        // bool Check_Point = false;
        // if(Bid == 3334 && num_3122 > 0 && num_m3122 == 0){
        //     Check_Point = true;
        // }
        // if(Bid == -3334){
        //     Check_Point = true;
        // }
        // if(Check_Point == false){
        //     continue;
        // }

        ////////////////////////////

        if(0 == num_of_counted_Omega_in_one_event && 0 == num_of_counted_Omega_bar_in_one_event){
            continue;
        }
        if(NumA != 0 && NumB != 0){
            num_of_counted_Bid++;
            for(int j=0;j<NumA;j++){
                TLorentzVector p1;
                p1.SetXYZM(A_px[j],A_py[j],A_pz[j],AMass);
                float ALength = pow(pow(A_px[j],2)+pow(A_py[j],2)+pow(A_pz[j],2) ,0.5);
                A_dNdy->Fill(A_y[j]);
                for(int k=0;k<NumB;k++){
                    TLorentzVector p2;
                    p2.SetXYZM(B_px[k],B_py[k],B_pz[k],BMass);

                    p3 = p1 + p2;
                    float dpt = fabs(p1.Perp() - p2.Perp());
                    Apt->Fill(dpt);
                    p1.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                    kstar = 0.5 * (p1 - p2).Rho();
                    Ak->Fill(kstar);
                    Ay->Fill(fabs(A_y[j] - B_y[k]));

                    float BLength = pow(pow(B_px[k],2)+pow(B_py[k],2)+pow(B_pz[k],2) ,0.5);
                    float CosAB = (A_px[j]*B_px[k] + A_py[j]*B_py[k] + A_pz[j]*B_pz[k])/(ALength*BLength);
                    Aphi->Fill(acos(CosAB));
                    if (j == 0){
                        B_dNdy->Fill(B_y[k]);
                    }
                }
            }
        }
        if(NumC != 0 && NumD != 0){
            num_of_counted_Did++;
            for(int j=0;j<NumC;j++){
                TLorentzVector p1;
                p1.SetXYZM(C_px[j],C_py[j],C_pz[j],CMass);
                float CLength = pow(pow(C_px[j],2)+pow(C_py[j],2)+pow(C_pz[j],2) ,0.5);
                C_dNdy->Fill(C_y[j]);
                for(int k=0;k<NumD;k++){
                    TLorentzVector p2;
                    p2.SetXYZM(D_px[k],D_py[k],D_pz[k],DMass);

                    p3 = p1 + p2;
                    float dpt = fabs(p1.Perp() - p2.Perp());
                    Bpt->Fill(dpt);
                    p1.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());
                    kstar = 0.5 * (p1 - p2).Rho();
                    Bk->Fill(kstar);
                    By->Fill(fabs(C_y[j] - D_y[k]));

                    float DLength = pow(pow(D_px[k],2)+pow(D_py[k],2)+pow(D_pz[k],2) ,0.5);
                    float CosCD = (C_px[j]*D_px[k] + C_py[j]*D_py[k] + C_pz[j]*D_pz[k])/(CLength*DLength);
                    Bphi->Fill(acos(CosCD));
                    if (j == 0){
                        D_dNdy->Fill(D_y[k]);
                    }
                }
            }
        }
    }


    TFile *file = new TFile("Cor.root", "RECREATE");

    
    TH1D *Ck = (TH1D *) Ak->Clone();
    Ck->SetName("Ck");
    Ck->SetTitle("Ck");
    TH1D *Dk = (TH1D *) Bk->Clone();
    Dk->SetName("Dk");
    Dk->SetTitle("Dk");
    
    TH1D *Cy = (TH1D *) Ay->Clone();
    Cy->SetName("Cy");
    Cy->SetTitle("Cy");
    TH1D *Dy = (TH1D *) By->Clone();
    Dy->SetName("Dy");
    Dy->SetTitle("Dy");
    
    TH1D *Cphi = (TH1D *) Aphi->Clone();
    Cphi->SetName("Cphi");
    Cphi->SetTitle("Cphi");
    TH1D *Dphi = (TH1D *) Bphi->Clone();
    Dphi->SetName("Dphi");
    Dphi->SetTitle("Dphi");

    TH1D *Cpt = (TH1D *) Apt->Clone();
    Cpt->SetName("Cpt");
    Cpt->SetTitle("Cpt");
    TH1D *Dpt = (TH1D *) Bpt->Clone();
    Dpt->SetName("Dpt");
    Dpt->SetTitle("Dpt");
    
    //通过pair归一化
    // Int_t binNorm[2];
    // binNorm[0] = Ak->FindBin(0);
    // binNorm[1] = Ak->FindBin(19.9999); 
    // Double_t factorN1 = Ak->Integral(binNorm[0], binNorm[1]);
    // Ck->Scale(1.0/factorN1);// Ck->Scale(1.0/factorN1)
    // factorN1 = Bk->Integral(binNorm[0], binNorm[1]);
    // Dk->Scale(1.0/factorN1);
    // factorN1 = Apt->Integral(binNorm[0], binNorm[1]);
    // Cpt->Scale(1.0/factorN1);
    // factorN1 = Bpt->Integral(binNorm[0], binNorm[1]);
    // Dpt->Scale(1.0/factorN1);

    //通过Omega和Omega_bar数量归一化
    Ck->Scale(1.0/num_of_counted_Bid);// Ck->Scale(1.0/factorN1)
    Dk->Scale(1.0/num_of_counted_Did);
    Cpt->Scale(1.0/num_of_counted_Bid);
    Dpt->Scale(1.0/num_of_counted_Did);
    Cy->Scale(1.0/num_of_counted_Bid);
    Dy->Scale(1.0/num_of_counted_Did);
    Cphi->Scale(1.0/num_of_counted_Bid);
    Dphi->Scale(1.0/num_of_counted_Did);

    TH1D *Ek = (TH1D *) Ck->Clone();
    Ek->SetName("Ek");
    Ek->Add(Dk, -1.0);
    TH1D *Ey = (TH1D *) Cy->Clone();
    Ey->SetName("Ey");
    Ey->Add(Dy, -1.0);
    TH1D *Ephi = (TH1D *) Cphi->Clone();
    Ephi->SetName("Ephi");
    Ephi->Add(Dphi, -1.0);
    TH1D *Ept = (TH1D *) Cpt->Clone();
    Ept->SetName("Ept");
    Ept->Add(Dpt, -1.0);
    Ak->GetXaxis()->SetTitle("k^{*} [GeV]");
    Ak->GetYaxis()->SetTitle("Counts");
    Bk->GetXaxis()->SetTitle("k^{*} [GeV]");
    Bk->GetYaxis()->SetTitle("Counts");
    Ck->GetXaxis()->SetTitle("k^{*} [GeV]");
    Ck->GetYaxis()->SetTitle("Difference");
    Cpt->GetXaxis()->SetTitle("#Delta P_{T} [GeV]");
    Cpt->GetYaxis()->SetTitle("Difference");

    Ak->Write();
    Bk->Write();
    Ck->Write();
    Dk->Write();
    Ek->Write();
    Ay->Write();
    By->Write();
    Cy->Write();
    Dy->Write();
    Ey->Write();
    Aphi->Write();
    Bphi->Write();
    Cphi->Write();
    Dphi->Write();
    Ephi->Write();
    Apt->Write();
    Bpt->Write();
    Cpt->Write();
    Dpt->Write();
    Ept->Write();
    testPx->Write();
    testPy->Write();
    testPz->Write();
    A_dNdy->Write();
    B_dNdy->Write();
    C_dNdy->Write();
    D_dNdy->Write();
    


    file->Write();

    return 0;
}
