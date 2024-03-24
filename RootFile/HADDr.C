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
    int Itr = 0;
    for(int i=66090;i <= 66389;i++){ // 66389
        // if (i == 66381){continue;}
        TString filename = midname;
        filename+=i;
        filename+=".root";

        // cout<<"FUCK"<<endl;
        TFile *fileR = TFile::Open(filename);
        TH1F* h1 = (TH1F*)fileR->Get("h_Good_Mass");
        TH1F* h2 = (TH1F*)fileR->Get("h_Bad_Mass");
        TH2D* h3 = (TH2D*)fileR->GetPrimitive("hHM_ParentDCA");
        list1->Add(h1);
        list2->Add(h2);
        list3->Add(h3);
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
    TH1F *Result1 = (TH1F*)h1->Clone("h_Good_Mass_Merged");
    TH1F *Result2 = (TH1F*)h2->Clone("h_Bad_Mass_Merged");
    TH2D *Result3 = (TH2D*)h3->Clone("hHM_ParentDCA_Merged");
    Result1->Reset();
    Result1->Merge(list1);
    Result2->Reset();
    Result2->Merge(list2);
    Result3->Reset();
    Result3->Merge(list3);
    TFile *file = new TFile("HADDr.root", "RECREATE");
    cout<<Result1->Integral(1,1000)<<endl;
    Result1->Write("h_Good_Mass_Merged");
    Result2->Write("h_Bad_Mass_Merged");
    Result3->Write("hHM_ParentDCA_Merged");
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
