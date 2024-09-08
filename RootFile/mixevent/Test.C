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

// 定义一个四维矢量类型
using FourVector = std::array<double, 4>;

// 计算质心系速度
std::array<double, 3> calculateBeta(const FourVector& p1, const FourVector& p2) {
    double totalPx = p1[0] + p2[0];
    double totalPy = p1[1] + p2[1];
    double totalPz = p1[2] + p2[2];
    double totalE = p1[3] + p2[3];

    return {totalPx / totalE, totalPy / totalE, totalPz / totalE};
}

// 计算 Lorentz boost
FourVector boost(const FourVector& p, const std::array<double, 3>& beta) {
    double beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
    double gamma = 1.0 / std::sqrt(1.0 - beta2);

    double bp = beta[0]*p[0] + beta[1]*p[1] + beta[2]*p[2];
    double gamma2 = (beta2 > 0) ? (gamma - 1.0) / beta2 : 0.0;

    FourVector boosted;
    boosted[0] = p[0] + gamma2 * bp * beta[0] + gamma * beta[0] * p[3];
    boosted[1] = p[1] + gamma2 * bp * beta[1] + gamma * beta[1] * p[3];
    boosted[2] = p[2] + gamma2 * bp * beta[2] + gamma * beta[2] * p[3];
    boosted[3] = gamma * (p[3] + bp);

    return boosted;
}

int test() {

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
                FourVector p1 = {A_px[j],A_py[j],A_pz[j],A_mass[j]}; // (px, py, pz, E)
                FourVector p2 = {B_px[k],B_py[k],B_pz[k],B_mass[k]};

                // 计算质心系速度矢量
                auto beta = calculateBeta(p1, p2);

                // Boost 两个四维矢量到质心系
                FourVector boostedP1 = boost(p1, beta);
                FourVector boostedP2 = boost(p2, beta);

                // 输出结果
                std::cout << "Boosted p1: (" << boostedP1[0] << ", " << boostedP1[1] << ", " << boostedP1[2] << ", " << boostedP1[3] << ")\n";
                std::cout << "Boosted p2: (" << boostedP2[0] << ", " << boostedP2[1] << ", " << boostedP2[2] << ", " << boostedP2[3] << ")\n";
            }
            if (true) {
                TLorentzVector p1;
                p1.SetXYZM(A_px[j],A_py[j],A_pz[j],A_mass[j]);
                TLorentzVector p2;
                p2.SetXYZM(B_px[k],B_py[k],B_pz[k],B_mass[k]);

                TLorentzVector p3 = p1 + p2;
                p1.Boost(-p3.BoostVector());p2.Boost(-p3.BoostVector());

                std::cout << "Boosted p1: (" << p1.Px() << ", " << p1.Py() << ", " << p1.Pz() << ", " << p1.E() << ")\n";
                std::cout << "Boosted p2: (" << p2.Px() << ", " << p2.Py() << ", " << p2.Pz() << ", " << p2.E() << ")\n";
            }
        }
    }

    return 0;
}