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

Double_t PI = 3.14159265358979383246;

const int PolyI = 7;
// const int PolyI = 4;

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

Double_t PolyFuction(Double_t *x, Double_t *params){
    Double_t value = params[0];
    for (Int_t i = 1; i < PolyI; ++i) {
       value += params[i] * TMath::Power(x[0], i); // 多项式的求和
    }
    return value;
}

Double_t IntPolyFuction(Double_t *x, Double_t *params){
    Double_t value = params[0]*x[0];
    for (Int_t i = 1; i < PolyI; ++i) {
       value += params[i] * TMath::Power(x[0], i+1) / (i+1); // 多项式的求和
    }
    return value;
}

Double_t GausFuction(Double_t *x, Double_t *params){
    Double_t value = params[PolyI]*TMath::Exp(-1.0*TMath::Power((x[0] - params[PolyI+1]),2)/(2*params[PolyI+2]*params[PolyI+2]));
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
    TString   List_Name[4] = { "HM_Lambda" , "HM_Lambdab" , "HM_Omega" , "HM_Omegab"};
    TString   TCan_Name[4] = { "CM_Lambda" , "CM_Lambdab" , "CM_Omega" , "CM_Omegab"};
    Double_t  FIT_X_Min[4] = {    1.090    ,    1.090     ,   1.63     ,    1.63    };
    Double_t  FIT_X_Max[4] = {    1.136    ,    1.136     ,   1.71     ,    1.71    };
    Double_t  FIT_X_Wid[4] = {    0.005    ,    0.005     ,   0.009    ,    0.009   };
    Double_t  FIT_X_Mid[4] = {    1.1165   ,    1.1165    ,   1.6725   ,    1.6725  };
    Double_t  FIT_A_Min[4] = {    1.09     ,    1.09      ,   1.635    ,    1.635   };
    Double_t  FIT_A_Max[4] = {    1.10     ,    1.10      ,   1.663    ,    1.663   };
    Double_t Hist_X_Min[4] = {    1.078    ,    1.078     ,   1.635    ,    1.635   };
    Double_t Hist_X_Max[4] = {    1.178    ,    1.178     ,   1.663    ,    1.663   };
    TH1F* h[4];
    TCanvas *canvas[4];
    TFile *fileR = TFile::Open(InputName,"read");
    cout<<"#####################################################################################"<<endl;
    TFile *file = new TFile(OutputName, "RECREATE");
    for (int i=0;i<4;i++) {
        // if (i > 0) {break;}
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
        }else if (PolyI == 4)  {
            customFunction->SetParameters(OrderFCT1, 0.0, 0.1, 0.2, maxBinValue - OrderFCT1, FIT_X_Mid[i] , 0.5*FIT_X_Wid[i]);
        }
        h[i]->GetXaxis()->SetRangeUser(Hist_X_Min[i], Hist_X_Max[i]);
        h[i]->Draw();
        h[i]->Fit(customFunction, "R");
        customFunction->SetLineColor(kGreen);
        customFunction->Draw("same");
        if       (PolyI == 7) {
            Double_t params[] = { customFunction->GetParameter(0) , customFunction->GetParameter(1) , customFunction->GetParameter(2) ,
                                  customFunction->GetParameter(3) , customFunction->GetParameter(4) , customFunction->GetParameter(5) ,
                                  customFunction->GetParameter(6) , customFunction->GetParameter(7) , customFunction->GetParameter(8) ,
                                  customFunction->GetParameter(9)};
            TF1* PolyF = new TF1 ("PolyFuction",PolyFuction,FIT_X_Min[i], FIT_X_Max[i], PolyI+3);
            PolyF->SetParameters(params);
            PolyF->SetLineColor(kRed);
            PolyF->Draw("same");
            TF1* GausF = new TF1 ("GausFuction",GausFuction,FIT_X_Min[i], FIT_X_Max[i], PolyI+3);
            GausF->SetParameters(params);
            GausF->SetLineColor(kBlue);
            GausF->Draw("same");

            Double_t Signal_Integral = TMath::Erf(3.0/TMath::Power(2.0,0.5))*(customFunction->GetParameter(9))*TMath::Power(2.0*PI,0.5)*customFunction->GetParameter(7);
            Double_t xl[] = {customFunction->GetParameter(8) - 3*customFunction->GetParameter(9)};
            Double_t xr[] = {customFunction->GetParameter(8) + 3*customFunction->GetParameter(9)};
            Double_t BackGr_Integral = IntPolyFuction(xr,params) - IntPolyFuction(xl,params);
            Double_t S_B = 1.0*Signal_Integral/BackGr_Integral;
            TString SBText = "In 3 sigma, S/B = ";
            SBText += S_B;
            cout<<SBText<<endl;
            Signal_Integral = Signal_Integral / (h[i]->GetBinWidth(1));
            TString SText = "In 3 sigma, S = ";
            SText += Signal_Integral;
            cout<<SText<<endl;
            TText *SB_text = new TText(0.7, 0.3, SBText);
            TText *S_text = new TText(0.7, 0.4, SText);
            SB_text->SetTextColor(kRed);
            SB_text->SetTextAlign(22);
            SB_text->SetTextSize(0.04);
            S_text->SetTextColor(kRed);
            S_text->SetTextAlign(22);
            S_text->SetTextSize(0.04);
            SB_text->Draw();
            S_text->Draw();
        }else if (PolyI == 4)  {
        }
        canvas[i]->Draw();
        canvas[i]->Write();
        cout<<"#####################################################################################"<<endl;
    }
    

    file->Close();
    fileR->Close();


    return;
}
