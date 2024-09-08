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

// 计算质心系速度
std::vector<double> calculateBeta(std::vector<double>& p1, std::vector<double>& p2) {
    std::vector<double> Result;
    double totalPx = p1[0] + p2[0];
    double totalPy = p1[1] + p2[1];
    double totalPz = p1[2] + p2[2];
    double totalE = p1[3] + p2[3];
    Result.push_back(totalPx / totalE);
    Result.push_back(totalPy / totalE);
    Result.push_back(totalPz / totalE);
    cout<<"Result = ("<<Result[0]<<" , "<<Result[1]<<" , "<<Result[2]<<" )"<<endl;
    return Result;
}

// 计算 Lorentz boost
std::vector<double> boost(std::vector<double>& p, std::vector<double>& beta) {
    double beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    double gamma = 1.0 / std::sqrt(1.0 - beta2);

    double bp = beta[0]*p[0] + beta[1]*p[1] + beta[2]*p[2];
    double gamma2 = (beta2 > 0) ? (gamma - 1.0) / beta2 : 0.0;

    std::vector<double> boosted;
    boosted.push_back(p[0] + gamma2 * bp * beta[0] + gamma * beta[0] * p[3]);
    boosted.push_back(p[1] + gamma2 * bp * beta[1] + gamma * beta[1] * p[3]);
    boosted.push_back(p[2] + gamma2 * bp * beta[2] + gamma * beta[2] * p[3]);
    boosted.push_back(gamma * (p[3] + bp));

    cout<<"boosted = ("<<boosted[0]<<" , "<<boosted[1]<<" , "<<boosted[2]<<" , "<<boosted[3]<<" )"<<endl;
    return boosted;
}

int Test() {

    int Num = 6;

    float A_Px[] = {  1.0 , -3.4 ,  2.4 ,  9.4 , -4.3 ,  6.6 };
    float A_Py[] = { -1.7 ,  3.4 , -0.4 ,  5.2 ,  8.6 ,  2.2 };
    float A_Pz[] = { -6.1 ,  1.9 ,  7.2 , -2.2 , -1.9 , -8.2 };
    float A_M[]  = {  1.1 ,  9.4 ,  2.9 ,  6.2 ,  8.5 ,  7.7 };
    
    float B_Px[] = {  8.6 ,  2.7 , -1.9 ,  2.1 ,  6.6 , -5.5 };
    float B_Py[] = { -6.7 , -7.3 ,  4.4 , -5.9 , -1.1 ,  0.2 };
    float B_Pz[] = {  4.1 , -0.9 ,  5.4 ,  7.8 ,  5.9 , -4.9 };
    float B_M[]  = {  3.1 ,  2.6 ,  3.4 ,  1.2 ,  4.7 ,  3.6 };

    for (int j = 0;j < Num;j++) {
        for (int k = 0;k < Num;k++) {
            cout<<"################################################################"<<endl;
            if (true) {
                // 定义两个四维矢量
                std::vector<double> p1; // (px, py, pz, E)
                p1.push_back(A_Px[j]);
                p1.push_back(A_Py[j]);
                p1.push_back(A_Pz[j]);
                p1.push_back(A_M[j]);
                std::vector<double> p2;
                p2.push_back(B_Px[k]);
                p2.push_back(B_Py[k]);
                p2.push_back(B_Pz[k]);
                p2.push_back(B_M[k]);

                // 计算质心系速度矢量
                auto beta = calculateBeta(p1, p2);

                // Boost 两个四维矢量到质心系
                std::vector<double> boostedP1 = boost(p1, beta);
                std::vector<double> boostedP2 = boost(p2, beta);

                // 输出结果
                std::cout << "Boosted p1: (" << boostedP1[0] << ", " << boostedP1[1] << ", " << boostedP1[2] << ", " << boostedP1[3] << ")\n";
                std::cout << "Boosted p2: (" << boostedP2[0] << ", " << boostedP2[1] << ", " << boostedP2[2] << ", " << boostedP2[3] << ")\n";
            }
            if (true) {
                TLorentzVector Pp1;
                Pp1.SetXYZM(A_Px[j],A_Py[j],A_Pz[j],A_M[j]);
                TLorentzVector Pp2;
                Pp2.SetXYZM(B_Px[k],B_Py[k],B_Pz[k],B_M[k]);

                TLorentzVector p3 = Pp1 + Pp2;
                Pp1.Boost(-p3.BoostVector());Pp2.Boost(-p3.BoostVector());

                std::cout << "Boosted p1: (" << Pp1.Px() << ", " << Pp1.Py() << ", " << Pp1.Pz() << ", " << Pp1.E() << ")\n";
                std::cout << "Boosted p2: (" << Pp2.Px() << ", " << Pp2.Py() << ", " << Pp2.Pz() << ", " << Pp2.E() << ")\n";
            }
        }
    }

    return 0;
}