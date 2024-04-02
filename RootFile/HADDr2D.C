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

// void HADDr(int from_FileID,int to_FIleID)
void HADDr2D()
{
    TString midname = "/star/u/svianping/STAR_Files/RootFile/HADDr";

    TList *list1 = new TList;
    TList *list2 = new TList;
    TList *list3 = new TList;
    TList *list4 = new TList;
    TList *list5 = new TList;
    TList *list6 = new TList;
    TList *list7 = new TList;
    TList *list8 = new TList;
    TFile *fileR;
    int Itr = 0;
    for(int i=0;i <= 21;i++){
        TString filename = midname;
        filename+=i;
        filename+=".root";

        // cout<<"FUCK"<<endl;
        TFile *fileR = TFile::Open(filename);
        // TFile *fileR = new TFile(filename,"read");
        TH1F* h1 = (TH1F*)fileR->Get("merge_Lambda");
        TH1F* h2 = (TH1F*)fileR->Get("merge_Lambdab");
        TH1F* h3 = (TH1F*)fileR->Get("merge_Omega");
        TH1F* h4 = (TH1F*)fileR->Get("merge_Omegab");
        TH1F* h5 = (TH1F*)fileR->Get("merge_DaughtersDCA_Lambda");
        TH1F* h6 = (TH1F*)fileR->Get("merge_DaughtersDCA_Lambdab");
        TH1F* h7 = (TH1F*)fileR->Get("merge_DaughtersDCA_Omega");
        TH1F* h8 = (TH1F*)fileR->Get("merge_DaughtersDCA_Omegab");
        list1->Add(h1);
        list2->Add(h2);
        list3->Add(h3);
        list4->Add(h4);
        list5->Add(h5);
        list6->Add(h6);
        list7->Add(h7);
        list8->Add(h8);
        Itr++;
    }
    TH1F *Result[9];
    TH1F *Result[1] = (TH1F*)h1->Clone("Merge_merge_Lambda");
    TH1F *Result[2] = (TH1F*)h2->Clone("Merge_merge_Lambdab");
    TH1F *Result[3] = (TH1F*)h3->Clone("Merge_merge_Omega");
    TH1F *Result[4] = (TH1F*)h4->Clone("Merge_merge_Omegab");
    TH1F *Result[5] = (TH1F*)h5->Clone("Merge_merge_DaughtersDCA_Lambda");
    TH1F *Result[6] = (TH1F*)h6->Clone("Merge_merge_DaughtersDCA_Lambdab");
    TH1F *Result[7] = (TH1F*)h7->Clone("Merge_merge_DaughtersDCA_Omega");
    TH1F *Result[8] = (TH1F*)h8->Clone("Merge_merge_DaughtersDCA_Omegab");
    Result[1]->Reset();
    Result[1]->Merge(list1);
    Result[2]->Reset();
    Result[2]->Merge(list2);
    Result[3]->Reset();
    Result[3]->Merge(list3);
    Result[4]->Reset();
    Result[4]->Merge(list4);
    Result[5]->Reset();
    Result[5]->Merge(list5);
    Result[6]->Reset();
    Result[6]->Merge(list6);
    Result[7]->Reset();
    Result[7]->Merge(list7);
    Result[8]->Reset();
    Result[8]->Merge(list8);
    TFile *file = new TFile("HADDrN.root", "RECREATE");
    Result[1]->Write("Merge_merge_Lambda");
    Result[2]->Write("Merge_merge_Lambdab");
    Result[3]->Write("Merge_merge_Omega");
    Result[4]->Write("Merge_merge_Omegab");
    Result[5]->Write("Merge_merge_DaughtersDCA_Lambda");
    Result[6]->Write("Merge_merge_DaughtersDCA_Lambdab");
    Result[7]->Write("Merge_merge_DaughtersDCA_Omega");
    Result[8]->Write("Merge_merge_DaughtersDCA_Omegab");

    for(int i=1;i<=4;i++){
        TString Name = "c";
        Name += i;
        auto c1 = new TCanvas(Name,Name,600,500);
        gStyle->SetOptStat(0);
        Result[i]->SetFillColor(2);
        Result[i]->SetLineColor(2);
        Result[i+4]->SetFillColor(3);
        Result[i+4]->SetLineColor(3);

        if (i <= 2){
            Result[i]  ->GetXaxis()->SetRangeUser(1.05, 1.8);
            Result[i]  ->GetXaxis()->SetRangeUser(1.05, 1.8);
            Result[i+4]->GetXaxis()->SetRangeUser(1.05, 1.8);
            Result[i+4]->GetXaxis()->SetRangeUser(1.05, 1.8);
        }
        else
        {
            Result[i]  ->GetXaxis()->SetRangeUser(1.6, 2.4);
            Result[i]  ->GetXaxis()->SetRangeUser(1.6, 2.4);
            Result[i+4]->GetXaxis()->SetRangeUser(1.6, 2.4);
            Result[i+4]->GetXaxis()->SetRangeUser(1.6, 2.4);
        }

        Result[i]->Draw();
        Result[i+4]->Draw("sames");

        auto legend = new TLegend(0.1,0.7,0.48,0.9);
        legend->SetHeader("The DCA between Daughts vs. mass"); // option "C" allows to center the header
        legend->AddEntry(Result[i],  "NO CUT","f");
        legend->AddEntry(Result[i+4],"DCA < 0.6 cm","f");
        legend->Draw();
        c1->Write();
    }
    file->Close();

    return;
}
