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


#define A_Num_Per_Event 5
#define B_Num_Per_Event 5
#define Max_Event_Per_Pool 100 // no larger than 255

struct Particle{
    unsigned short int TreeID;
    float Px;
    float Py;
    float Pz;

    float Rap;
    float Pt;
};

struct ParticlePool{
    unsigned int EvtID;
    bool IfMadeSame = false;
    unsigned short int ListA_Index;
    Particle ListA[A_Num_Per_Event];
    unsigned short int ListB_Index;
    Particle ListB[B_Num_Per_Event];
};

float* f() {
    float* MassAndKstar = new float[2];
    MassAndKstar[0] = 5.0;
    MassAndKstar[1] = 10.0;
    return MassAndKstar;
}

void NC()
{
    for (int i=0;i<2;i++) {
        float* R = f();
        cout<<"R = ["<<R[0]<<","<<R[1]<<"]"<<endl;
        delete[] R;
    }

}
