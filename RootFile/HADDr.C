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

void HADDr()
{

    // TFile *file = new TFile("HADDr.root", "RECREATE");
    
    // TH1F* Result1;
    // TString midname = "/star/data01/pwg/svianping/output/output_";

    // int Itr = 0;
    // for(int i=65690;i <= 66389;i++){
    //     TString filename = midname;
    //     filename+=i;
    //     filename+=".root";

    //         cout<<"FUCK"<<endl;
    //     TFile *fileR = TFile::Open(filename);
    //     TH1F* h1 = (TH1F*)fileR->Get("h_Good_Mass");
    //     if (Itr == 0) {
    //         Result1 = (TH1F*)h1->Clone();
    //     }
    //     else{
    //         Result1->Add((TH1F*)h1->Clone());
    //     }
    //     fileR->Close();
    //     cout<<Itr<<endl;
    //     Itr++;
    // }
    
    // cout<<Result1->Integral(0,5)<<endl;
    // Result1->Write("h_Good_Mass");
    // // file->Write();
    // file->Close();

//////////////////////////////////////////////////////////////////////////////
    
    // TH1F* Result1,Temp;
    TString midname = "/star/data01/pwg/svianping/output/output_";

    TList *list1 = new TList;
    TList *list2 = new TList;
    TList *list3 = new TList;
    TList *list4 = new TList;
    TList *list5 = new TList;
    TList *list6 = new TList;
    TList *list7 = new TList;
    TList *list8 = new TList;
    int Itr = 0;
    for(int i=62690;i <= 66389;i++){ // 66389
        // if (i == 66381){continue;}
        TString filename = midname;
        filename+=i;
        filename+=".root";

        // cout<<"FUCK"<<endl;
        TFile *fileR = TFile::Open(filename);
        TH1F* h1 = (TH1F*)fileR->Get("H_DaughterDCA_LitP1_Mass");
        TH1F* h2 = (TH1F*)fileR->Get("H_DaughterDCA_LitP2_Mass");
        TH1F* h3 = (TH1F*)fileR->Get("H_DaughterDCA_LitP3_Mass");
        TH1F* h4 = (TH1F*)fileR->Get("H_DaughterDCA_LitP4_Mass");
        TH1F* h5 = (TH1F*)fileR->Get("H_DaughterDCA_LitP5_Mass");
        TH1F* h6 = (TH1F*)fileR->Get("H_DaughterDCA_LitP6_Mass");
        TH1F* h7 = (TH1F*)fileR->Get("H_DaughterDCA_LitP8_Mass");
        TH1F* h8 = (TH1F*)fileR->Get("H_DaughterDCA_NOLIM_Mass");
        // TH2D* h3 = (TH2D*)fileR->Get("hH_M_ParentDCA");
        list1->Add(h1);
        list2->Add(h2);
        list3->Add(h3);
        list4->Add(h4);
        list5->Add(h5);
        list6->Add(h6);
        list7->Add(h7);
        list8->Add(h8);
        // if (Itr == 0) {
        //     Result1 = (TH1F*)h1->Clone();
        // }
        // else{
        //     Temp = (TH1F*)h1->Clone();
        //     Result1->Add(Result1,Temp,1.0,1.0);
        // }
        // cout<<Itr<<" : "<<Result1->Integral(1,1000)<<endl;
        // fileR->Close();
        Itr++;
    }
    TH1F *Result1 = (TH1F*)h1->Clone("H_DaughterDCA_LitP1_Mass_merge");
    TH1F *Result2 = (TH1F*)h2->Clone("H_DaughterDCA_LitP2_Mass_merge");
    TH1F *Result3 = (TH1F*)h3->Clone("H_DaughterDCA_LitP3_Mass_merge");
    TH1F *Result4 = (TH1F*)h4->Clone("H_DaughterDCA_LitP4_Mass_merge");
    TH1F *Result5 = (TH1F*)h5->Clone("H_DaughterDCA_LitP5_Mass_merge");
    TH1F *Result6 = (TH1F*)h6->Clone("H_DaughterDCA_LitP6_Mass_merge");
    TH1F *Result7 = (TH1F*)h8->Clone("H_DaughterDCA_LitP8_Mass_merge");
    TH1F *Result8 = (TH1F*)h9->Clone("H_DaughterDCA_NOLIM_Mass_merge");
    // TH2D *Result3 = (TH2D*)h3->Clone("hHM_ParentDCA_Merged");
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
    TFile *file = new TFile("HADDr.root", "RECREATE");
    cout<<Result1->Integral(1,1000)<<endl;
    Result1->Write("H_DaughterDCA_LitP1_Mass_merge");
    Result2->Write("H_DaughterDCA_LitP2_Mass_merge");
    Result3->Write("H_DaughterDCA_LitP3_Mass_merge");
    Result4->Write("H_DaughterDCA_LitP4_Mass_merge");
    Result5->Write("H_DaughterDCA_LitP5_Mass_merge");
    Result6->Write("H_DaughterDCA_LitP6_Mass_merge");
    Result7->Write("H_DaughterDCA_LitP8_Mass_merge");
    Result8->Write("H_DaughterDCA_NOLIM_Mass_merge");
    // file->Write();
    file->Close();

//////////////////////////////////////////////////////////////////////////////

    // TString midname = "/star/data01/pwg/svianping/output/output_";
    // TH1F* Result1,Temp;
    // for(int i=65690;i <= 66389;i++){
    //     TString filename = midname;
    //     filename+=i;
    //     filename+=".root";

    //     TFile *fileR = TFile::Open(filename);
    //     TH1F* h1 = (TH1F*)fileR->Get("h_Good_Mass");
    //     for (int j=0;j<h1->GetNbins();j++){}
    //     fileR->Close();
    // }

    return;
}
