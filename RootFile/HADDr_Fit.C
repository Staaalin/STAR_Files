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

Double_t FitFuction(Double_t *x, Double_t *params){
    // params[0]  多项式常数项
    // params[i]  多项式级数项
    // params[10] 高斯分布幅值
    // params[11] 高斯分布\mu

    Double_t value = params[0];
    for (Int_t i = 1; i < 10; ++i) {
       value += params[i] * TMath::Power(x[0], i); // 多项式的求和
    }
    value += params[10]*TMath::Exp((x[0] - params[11]));
    return value;
}

void HADDr(const TString InputName,const TString OutputName)
{
    TString List_Name[4] = { "HM_Lambda" , "HM_Lambdab" , "HM_Omega" , "HM_Omegab"};
    TString TCan_Name[4] = { "CM_Lambda" , "CM_Lambdab" , "CM_Omega" , "CM_Omegab"};
    TH1F* h[4];
    TCanvas *canvas[4];
    TFile *fileR = TFile::Open(InputName,"read");
    for (int i=0;i<4;i++) {
        h[i] = (TH1F*)fileR->Get(List_Name[i]);
        canvas[i] = new TCanvas(TCan_Name[i] , TCan_Name[i]);
        h[i]->Draw();
    }
    TFile *file = new TFile(OutputName, "RECREATE");

    file->Close();


    return;
}
