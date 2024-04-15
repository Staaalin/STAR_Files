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

const int PolyI = 7;

Double_t FitFuction(Double_t *x, Double_t *params){
    // params[0]  多项式常数项
    // params[i]  多项式级数项
    // params[10] 高斯分布幅值
    // params[11] 高斯分布\mu
    // params[12] 高斯分布\sigma

    Double_t value = params[0];
    for (Int_t i = 1; i < PolyI; ++i) {
       value += params[i] * TMath::Power(x[0], i); // 多项式的求和
    }
    value += params[PolyI]*TMath::Exp(-1.0*TMath::Power((x[0] - params[PolyI+1]),2)/(2*params[PolyI+2]*params[PolyI+2]));
    return value;
}

int XValue2BinID(TH1F* histogram, Double_t x){
    int binCount = histogram->GetNbinsX();
    for (int i = 0;i<binCount;i++){
        if (x <= histogram->GetBinContent(i)){
            return i;
        }
    }
}

Double_t calculateAverage(TH1F* histogram, int startBin, int endBin) {
   // 求和变量和计数器
   Double_t sum = 0;
   int count = 0;

   // 遍历区间内的bin，并求和
   for (int bin = startBin; bin <= endBin; bin++) {
      Double_t binContent = 1.0*histogram->GetBinContent(bin);
      sum += binContent;
      count++;
   }

   // 计算平均值
   Double_t average = sum / (1.0*count);

   return average;
}

void HADDr_Fit(const TString InputName,const TString OutputName)
{
    TString List_Name[4] = { "HM_Lambda" , "HM_Lambdab" , "HM_Omega" , "HM_Omegab"};
    TString TCan_Name[4] = { "CM_Lambda" , "CM_Lambdab" , "CM_Omega" , "CM_Omegab"};
    Double_t  FIT_X_Min[4] = {    1.09     ,    1.09      ,   1.63     ,    1.63    };
    Double_t  FIT_X_Max[4] = {    1.14     ,    1.14      ,   1.71     ,    1.71    };
    Double_t  FIT_X_Wid[4] = {    0.005    ,    0.005     ,   0.009    ,    0.009   };
    Double_t  FIT_X_Mid[4] = {    1.1165   ,    1.1165    ,   1.6725   ,    1.6725  };
    Double_t  FIT_A_Min[4] = {    1.09     ,    1.09      ,   1.635    ,    1.635   };
    Double_t  FIT_A_Max[4] = {    1.10     ,    1.10      ,   1.663    ,    1.663   };
    TH1F* h[4];
    TCanvas *canvas[4];
    TFile *fileR = TFile::Open(InputName,"read");
    cout<<"#####################################################################################"<<endl;
    TFile *file = new TFile(OutputName, "RECREATE");
    for (int i=0;i<4;i++) {
        cout<<"Start in "<<List_Name[i]<<endl;
        h[i] = (TH1F*)fileR->Get(List_Name[i]);
        canvas[i] = new TCanvas(TCan_Name[i] , TCan_Name[i]);
        int maxBinIndex = h[i]->GetMaximumBin();
        Double_t maxBinValue = h[i]->GetBinContent(maxBinIndex);
        cout<<"MaxBin at "<<maxBinIndex<<", and value is "<<maxBinValue<<endl;
        TF1 *customFunction = new TF1("SFunction", FitFuction, FIT_X_Min[i], FIT_X_Max[i], PolyI+3);
        Double_t OrderFCT1 = calculateAverage(h[i] , XValue2BinID(h[i], FIT_A_Min[i]) , XValue2BinID(h[i], FIT_A_Max[i]));
        cout<<"The 1st order poly factor is "<<OrderFCT1<<endl;
        if       (PolyI == 7) {
            customFunction->SetParameters(OrderFCT1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, maxBinValue - OrderFCT1, FIT_X_Mid[i] , 0.5*FIT_X_Wid[i]);
        }else if (PolyI == 5)  {
            customFunction->SetParameters(OrderFCT1, 0.0, 0.1, 0.2, 0.3, maxBinValue - OrderFCT1, FIT_X_Mid[i] , 0.5*FIT_X_Wid[i]);
        }
        h[i]->Draw();
        h[i]->Fit(customFunction);
        customFunction->Draw("same");
        canvas[i]->Draw();
        cout<<"#####################################################################################"<<endl;
    }
    

    file->Close();
    fileR->Close();


    return;
}
