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
    
    TFile *file = new TFile("HADDr.root", "RECREATE");
    
    TH1F* Result;
    TString midname = "/star/data01/pwg/svianping/output/output_";

    int Itr = 0;
    for(int i=65690;i <= 66389;i++){
        TString filename = midname;
        filename+=i;
        filename+=".root";

        TFile *fileR = new TFile(filename, "READ");
        TH1F* h1 = (TH1F*)fileR->Get("h_Good_Mass");
        if (Itr == 0) {
            Result = (TH1F*)h1->Clone();
        }
        else{
            Result += (TH1F*)h1->Clone();
        }
        fileR->Close();
        Itr++;
        if (Itr/10 == 0){cout<<Itr<<endl;}
    }
    
    Result->Write();
    file->Write();
    file->Close();

    return;
}
