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

void HADDr_xml(const TString FileList,const TString JobID)
{
    int Filr_Itr = 22; // IMPORTANT
    int Start_File = 57590;
    TString midname  = "/star/data01/pwg/svianping/output/output_";
    TString RecordFileName = "HADDr";

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
     
    RecordFileName += Filr_Itr;
    int File_Start = Start_File + Filr_Itr*400;
    int File_End   = Start_File + (Filr_Itr + 1)*400 - 1;
    if (File_Start >= 66390){return;}
    if (File_End   >= 66390){File_End = 66389;}
    for(int i=File_Start;i <= File_End;i++){
        if (i == 66390){cout<<"END"<<endl;break;}
        TString filename = midname;
        filename+=i;
        filename+=".root";

        // cout<<"FUCK"<<endl;
        TFile *fileR = TFile::Open(filename,"read");
        if (!fileR){continue;}
        // TFile *fileR = new TFile(filename,"read");
        TH1F* h1 = (TH1F*)fileR->Get("HM_Lambda");
        TH1F* h2 = (TH1F*)fileR->Get("HM_Lambdab");
        TH1F* h3 = (TH1F*)fileR->Get("HM_Omega");
        TH1F* h4 = (TH1F*)fileR->Get("HM_Omegab");
        TH1F* h5 = (TH1F*)fileR->Get("HM_DaughtersDCA_Lambda");
        TH1F* h6 = (TH1F*)fileR->Get("HM_DaughtersDCA_Lambdab");
        TH1F* h7 = (TH1F*)fileR->Get("HM_DaughtersDCA_Omega");
        TH1F* h8 = (TH1F*)fileR->Get("HM_DaughtersDCA_Omegab");
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
    TH1F *Result1 = (TH1F*)h1->Clone("merge_Lambda");
    TH1F *Result2 = (TH1F*)h2->Clone("merge_Lambdab");
    TH1F *Result3 = (TH1F*)h3->Clone("merge_Omega");
    TH1F *Result4 = (TH1F*)h4->Clone("merge_Omegab");
    TH1F *Result5 = (TH1F*)h5->Clone("merge_DaughtersDCA_Lambda");
    TH1F *Result6 = (TH1F*)h6->Clone("merge_DaughtersDCA_Lambdab");
    TH1F *Result7 = (TH1F*)h7->Clone("merge_DaughtersDCA_Omega");
    TH1F *Result8 = (TH1F*)h8->Clone("merge_DaughtersDCA_Omegab");
    Result1->Reset();
    Result1->Merge(list1);
    Result2->Reset();
    Result2->Merge(list2);
    Result3->Reset();
    Result3->Merge(list3);
    Result4->Reset();
    Result4->Merge(list4);
    Result5->Reset();
    Result5->Merge(list5);
    Result6->Reset();
    Result6->Merge(list6);
    Result7->Reset();
    Result7->Merge(list7);
    Result8->Reset();
    Result8->Merge(list8);
    RecordFileName += ".root";
    TFile *file = new TFile(RecordFileName, "RECREATE");
    cout<<Result1->Integral(1,1000)<<endl;
    Result1->Write("merge_Lambda");
    Result2->Write("merge_Lambdab");
    Result3->Write("merge_Omega");
    Result4->Write("merge_Omegab");
    Result5->Write("merge_DaughtersDCA_Lambda");
    Result6->Write("merge_DaughtersDCA_Lambdab");
    Result7->Write("merge_DaughtersDCA_Omega");
    Result8->Write("merge_DaughtersDCA_Omegab");
    file->Close();

//////////////////////////////////////////////////////////////////////////////

    // TString midname = "/star/data01/pwg/svianping/output/output_";

    // TH1F* h1,h2,h3,h4,h5,h6,h7,h8;
    // int Itr = 0;
    // for(int i=62690;i <= 66389;i++){ // 66389
    //     TString filename = midname;
    //     filename+=i;
    //     filename+=".root";
    //     if (Itr == 0){
    //         TFile *fileR = TFile::Open(filename);
    //         h1 = (TH1F*)fileR->Get("H_DaughterDCA_LitP1_Mass");
    //         h2 = (TH1F*)fileR->Get("H_DaughterDCA_LitP2_Mass");
    //         h3 = (TH1F*)fileR->Get("H_DaughterDCA_LitP3_Mass");
    //         h4 = (TH1F*)fileR->Get("H_DaughterDCA_LitP4_Mass");
    //         h5 = (TH1F*)fileR->Get("H_DaughterDCA_LitP5_Mass");
    //         h6 = (TH1F*)fileR->Get("H_DaughterDCA_LitP6_Mass");
    //         h7 = (TH1F*)fileR->Get("H_DaughterDCA_LitP8_Mass");
    //         h8 = (TH1F*)fileR->Get("H_DaughterDCA_NOLIM_Mass");
    //         continue;
    //     }

    //     // cout<<"FUCK"<<endl;
    //     TFile *fileR = TFile::Open(filename);
    //     TH1F* t1 = (TH1F*)fileR->Get("H_DaughterDCA_LitP1_Mass");
    //     TH1F* t2 = (TH1F*)fileR->Get("H_DaughterDCA_LitP2_Mass");
    //     TH1F* t3 = (TH1F*)fileR->Get("H_DaughterDCA_LitP3_Mass");
    //     TH1F* t4 = (TH1F*)fileR->Get("H_DaughterDCA_LitP4_Mass");
    //     TH1F* t5 = (TH1F*)fileR->Get("H_DaughterDCA_LitP5_Mass");
    //     TH1F* t6 = (TH1F*)fileR->Get("H_DaughterDCA_LitP6_Mass");
    //     TH1F* t7 = (TH1F*)fileR->Get("H_DaughterDCA_LitP8_Mass");
    //     TH1F* t8 = (TH1F*)fileR->Get("H_DaughterDCA_NOLIM_Mass");
    //     h1->Add(t1,1);
    //     h2->Add(t2,1);
    //     h3->Add(t3,1);
    //     h4->Add(t4,1);
    //     h5->Add(t5,1);
    //     h6->Add(t6,1);
    //     h7->Add(t7,1);
    //     h8->Add(t8,1);
    //     Itr++;
    // }
    // TH1F *Result1 = (TH1F*)h1->Clone("H_DaughterDCA_LitP1_Mass_merge");
    // TH1F *Result2 = (TH1F*)h2->Clone("H_DaughterDCA_LitP2_Mass_merge");
    // TH1F *Result3 = (TH1F*)h3->Clone("H_DaughterDCA_LitP3_Mass_merge");
    // TH1F *Result4 = (TH1F*)h4->Clone("H_DaughterDCA_LitP4_Mass_merge");
    // TH1F *Result5 = (TH1F*)h5->Clone("H_DaughterDCA_LitP5_Mass_merge");
    // TH1F *Result6 = (TH1F*)h6->Clone("H_DaughterDCA_LitP6_Mass_merge");
    // TH1F *Result7 = (TH1F*)h8->Clone("H_DaughterDCA_LitP8_Mass_merge");
    // TH1F *Result8 = (TH1F*)h9->Clone("H_DaughterDCA_NOLIM_Mass_merge");
    // TFile *file = new TFile("HADDr.root", "RECREATE");
    // Result1->Write("H_DaughterDCA_LitP1_Mass_merge");
    // Result2->Write("H_DaughterDCA_LitP2_Mass_merge");
    // Result3->Write("H_DaughterDCA_LitP3_Mass_merge");
    // Result4->Write("H_DaughterDCA_LitP4_Mass_merge");
    // Result5->Write("H_DaughterDCA_LitP5_Mass_merge");
    // Result6->Write("H_DaughterDCA_LitP6_Mass_merge");
    // Result7->Write("H_DaughterDCA_LitP8_Mass_merge");
    // Result8->Write("H_DaughterDCA_NOLIM_Mass_merge");
    // file->Close();

    return;
}
