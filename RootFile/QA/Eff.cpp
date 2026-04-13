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
#include "TH2.h"
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
#include "TString.h"
// #endif
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <stdio.h>
using namespace std;

// 使用这个编译：
// singularity exec -e --env DISPLAY=$DISPLAY -B /direct -B /gpfs -B /star -B /cvmfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh
// g++ -O2 -std=c++11 Eff.cpp -o Eff `root-config --cflags --libs`

// 定义粒子结构体
struct ArmParticle {
    float px;       // x方向动量
    float py;       // y方向动量
    float pz;       // z方向动量
    float mass;     // 质量
    float eta;      // 赝快度
    float y;        // 快度
    float pt;       // 横向动量
    bool  IsRecord; // 是否被记录
    int   TreeID;   // ID in one event
    std::vector<int>   ParentID; // Parent Particle ID in one event
    
    ArmParticle()
        : px(0), py(0), pz(0), mass(0),
          eta(0), y(0), pt(0),
          IsRecord(false), TreeID(0) {}
    
    // 构造函数
    ArmParticle(float _px, float _py, float _pz, float _mass, int _TreeID) 
        : px(_px), py(_py), pz(_pz), mass(_mass), TreeID(_TreeID) {
        // 计算赝快度、快度和横向动量
        pt = sqrt(px*px + py*py);
        float p = sqrt(pt*pt + pz*pz);
        float E = sqrt(p*p+mass*mass);
        eta = -1.0*log(tan(0.5*(acos(pz/p))));
        y = 0.5 * log((E + pz) / (E - pz));
        IsRecord = false;
    }
    
    // 计算能量
    float energy() const {
        return sqrt(px*px + py*py + pz*pz + mass*mass);
    }
    
    // 转换为四动量
    TLorentzVector lorentzVector() const {
        return TLorentzVector(px, py, pz, energy());
    }

    ArmParticle(const ArmParticle& other)
        : px(other.px), py(other.py), pz(other.pz),
          mass(other.mass), eta(other.eta), y(other.y), pt(other.pt),
          IsRecord(other.IsRecord), TreeID(other.TreeID),
          ParentID(other.ParentID) {}

    ArmParticle& operator=(const ArmParticle& other) {
        if (this != &other) {
            px = other.px;
            py = other.py;
            pz = other.pz;
            mass = other.mass;
            eta = other.eta;
            y = other.y;
            pt = other.pt;
            IsRecord = other.IsRecord;
            TreeID = other.TreeID;
            ParentID = other.ParentID;
        }
        return *this;
    }

};

// ROOT 5（特别是 ROOT 5.34.39）的字典机制有些老旧。
// 它会自动为 std::vector<T> 生成迭代器类型的字典（如 random_access_iterator<T,long>），但 C++11 之后这些类型模板已不再定义，因此报错。
#if defined(__CINT__) || defined(__CLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ArmParticle+;
// #pragma link C++ class std::vector<ArmParticle>+;
#endif
// 这告诉 ROOT：
// “只生成 ArmParticle 和 std::vector<ArmParticle> 的字典，
// 不要去尝试生成任何 random_access_iterator 之类的模板。”


#define Pi 3.1415926535898
#define HowMuchEventMixing 10

void Eff(
    TString MidName,
    TString DataName,
    int OutputFileIndex,
    TString OutMidName,
    int CutID = 0
) {
    std::cout<<"Start Eff.cpp"<<std::endl;

    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) 

        std::vector<int>     *PDG                = nullptr;
        std::vector<Float_t> *mix_px             = nullptr;
        std::vector<Float_t> *mix_py             = nullptr;
        std::vector<Float_t> *mix_pz             = nullptr;
        std::vector<Float_t> *QA_eta             = nullptr;
        std::vector<Float_t> *dEdx               = nullptr;
        std::vector<Float_t> *m2                 = nullptr;
        std::vector<Float_t> *dcatopv            = nullptr;
        std::vector<Float_t> *nSigmaProton       = nullptr;
        std::vector<Float_t> *nSigmaPion         = nullptr;
        std::vector<Float_t> *nSigmaKaon         = nullptr;
        std::vector<Float_t> *InvariantMass      = nullptr;
        std::vector<Float_t> *Decay_Length       = nullptr;
        std::vector<Float_t> *Chi2               = nullptr;
        std::vector<Float_t> *nHitsFit           = nullptr;
        std::vector<Float_t> *nHitsMax           = nullptr;
        std::vector<int>     *ParentList         = nullptr;
        std::vector<int>     *ParentSta          = nullptr;
        std::vector<int>     *ParentEnd          = nullptr;
        std::vector<int>     *SE_ParentList      = nullptr;
        std::vector<int>     *SE_ParentSta       = nullptr;
        std::vector<int>     *SE_ParentEnd       = nullptr;
        std::vector<int>     *ME_ParentList      = nullptr;
        std::vector<int>     *ME_ParentSta       = nullptr;
        std::vector<int>     *ME_ParentEnd       = nullptr;

        TBranch *bPDG                            = nullptr;
        TBranch *bmix_px                         = nullptr;
        TBranch *bmix_py                         = nullptr;
        TBranch *bmix_pz                         = nullptr;
        TBranch *bQA_eta                         = nullptr;
        TBranch *bdEdx                           = nullptr;
        TBranch *bm2                             = nullptr;
        TBranch *bdcatopv                        = nullptr;
        TBranch *bnSigmaProton                   = nullptr;
        TBranch *bnSigmaPion                     = nullptr;
        TBranch *bnSigmaKaon                     = nullptr;
        TBranch *bInvariantMass                  = nullptr;
        TBranch *bDecay_Length                   = nullptr;
        TBranch *bChi2                           = nullptr;
        TBranch *bnHitsFit                       = nullptr;
        TBranch *bnHitsMax                       = nullptr;
        TBranch *bParentList                     = nullptr;
        TBranch *bParentSta                      = nullptr;
        TBranch *bParentEnd                      = nullptr;
        TBranch *bSE_ParentList                  = nullptr;
        TBranch *bSE_ParentSta                   = nullptr;
        TBranch *bSE_ParentEnd                   = nullptr;
        TBranch *bME_ParentList                  = nullptr;
        TBranch *bME_ParentSta                   = nullptr;
        TBranch *bME_ParentEnd                   = nullptr;
    
    #else
        #if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)

            std::vector<int>     *PDG                = NULL;
            std::vector<Float_t> *mix_px             = NULL;
            std::vector<Float_t> *mix_py             = NULL;
            std::vector<Float_t> *mix_pz             = NULL;
            std::vector<Float_t> *QA_eta             = NULL;
            std::vector<Float_t> *dEdx               = NULL;
            std::vector<Float_t> *m2                 = NULL;
            std::vector<Float_t> *dcatopv            = NULL;
            std::vector<Float_t> *nSigmaProton       = NULL;
            std::vector<Float_t> *nSigmaPion         = NULL;
            std::vector<Float_t> *nSigmaKaon         = NULL;
            std::vector<Float_t> *InvariantMass      = NULL;
            std::vector<Float_t> *Decay_Length       = NULL;
            std::vector<Float_t> *Chi2               = NULL;
            std::vector<Float_t> *nHitsFit           = NULL;
            std::vector<Float_t> *nHitsMax           = NULL;
            std::vector<int>     *ParentList         = NULL;
            std::vector<int>     *ParentSta          = NULL;
            std::vector<int>     *ParentEnd          = NULL;
            std::vector<int>     *SE_ParentList      = NULL;
            std::vector<int>     *SE_ParentSta       = NULL;
            std::vector<int>     *SE_ParentEnd       = NULL;
            std::vector<int>     *ME_ParentList      = NULL;
            std::vector<int>     *ME_ParentSta       = NULL;
            std::vector<int>     *ME_ParentEnd       = NULL;

            TBranch *bPDG                            = NULL;
            TBranch *bmix_px                         = NULL;
            TBranch *bmix_py                         = NULL;
            TBranch *bmix_pz                         = NULL;
            TBranch *bQA_eta                         = NULL;
            TBranch *bdEdx                           = NULL;
            TBranch *bm2                             = NULL;
            TBranch *bdcatopv                        = NULL;
            TBranch *bnSigmaProton                   = NULL;
            TBranch *bnSigmaPion                     = NULL;
            TBranch *bnSigmaKaon                     = NULL;
            TBranch *bInvariantMass                  = NULL;
            TBranch *bDecay_Length                   = NULL;
            TBranch *bChi2                           = NULL;
            TBranch *bnHitsFit                       = NULL;
            TBranch *bnHitsMax                       = NULL;
            TBranch *bParentList                     = NULL;
            TBranch *bParentSta                      = NULL;
            TBranch *bParentEnd                      = NULL;
            TBranch *bSE_ParentList                  = NULL;
            TBranch *bSE_ParentSta                   = NULL;
            TBranch *bSE_ParentEnd                   = NULL;
            TBranch *bME_ParentList                  = NULL;
            TBranch *bME_ParentSta                   = NULL;
            TBranch *bME_ParentEnd                   = NULL;

        #else
    
            std::vector<int>     *PDG                = 0;
            std::vector<Float_t> *mix_px             = 0;
            std::vector<Float_t> *mix_py             = 0;
            std::vector<Float_t> *mix_pz             = 0;
            std::vector<Float_t> *QA_eta             = 0;
            std::vector<Float_t> *dEdx               = 0;
            std::vector<Float_t> *m2                 = 0;
            std::vector<Float_t> *dcatopv            = 0;
            std::vector<Float_t> *nSigmaProton       = 0;
            std::vector<Float_t> *nSigmaPion         = 0;
            std::vector<Float_t> *nSigmaKaon         = 0;
            std::vector<Float_t> *InvariantMass      = 0;
            std::vector<Float_t> *Decay_Length       = 0;
            std::vector<Float_t> *Chi2               = 0;
            std::vector<Float_t> *nHitsFit           = 0;
            std::vector<Float_t> *nHitsMax           = 0;
            std::vector<int>     *ParentList         = 0;
            std::vector<int>     *ParentSta          = 0;
            std::vector<int>     *ParentEnd          = 0;
            std::vector<int>     *SE_ParentList      = 0;
            std::vector<int>     *SE_ParentSta       = 0;
            std::vector<int>     *SE_ParentEnd       = 0;
            std::vector<int>     *ME_ParentList      = 0;
            std::vector<int>     *ME_ParentSta       = 0;
            std::vector<int>     *ME_ParentEnd       = 0;
    
            TBranch *bPDG                            = 0;
            TBranch *bmix_px                         = 0;
            TBranch *bmix_py                         = 0;
            TBranch *bmix_pz                         = 0;
            TBranch *bQA_eta                         = 0;
            TBranch *bdEdx                           = 0;
            TBranch *bm2                             = 0;
            TBranch *bdcatopv                        = 0;
            TBranch *bnSigmaProton                   = 0;
            TBranch *bnSigmaPion                     = 0;
            TBranch *bnSigmaKaon                     = 0;
            TBranch *bInvariantMass                  = 0;
            TBranch *bDecay_Length                   = 0;
            TBranch *bChi2                           = 0;
            TBranch *bnHitsFit                       = 0;
            TBranch *bnHitsMax                       = 0;
            TBranch *bParentList                     = 0;
            TBranch *bParentSta                      = 0;
            TBranch *bParentEnd                      = 0;
            TBranch *bSE_ParentList                  = 0;
            TBranch *bSE_ParentSta                   = 0;
            TBranch *bSE_ParentEnd                   = 0;
            TBranch *bME_ParentList                  = 0;
            TBranch *bME_ParentSta                   = 0;
            TBranch *bME_ParentEnd                   = 0;

        #endif
    #endif

    const int Recording_Particle[] = {321 , -321 , 3122 , -3122 , 3312 , -3312 , 3334 , -3334};
    const Int_t ParticleNum = sizeof(Recording_Particle)/sizeof(Recording_Particle[0]);
    std::vector<int> Particle_Mass , Particle_MassSigma;
    for (int i=0;i<ParticleNum;i++) {
        Particle_Mass.push_back(massList(Recording_Particle[i], DataName));
        Particle_MassSigma.push_back(massListSigma(Recording_Particle[i], DataName));
    }

    ArmParticle           A(0,0,0,0,0), B(0,0,0,0,0), C(0,0,0,0,0), D(0,0,0,0,0);

    bool IfRecord = true , IfRemoveFeedPair = false , IfRemoveSpliteMerge = false , IfRemoveLownHits = false , IfRemoveHighPVz = false , IfRemoveHighTPCsigma = false , IfCutHighDCA = false;

    int kStarBinNum = 400;
    float kStarSta = 0 , kStarEnd = 8;
    
    int dRapBinNum = 500;
    float dRapSta = -5 , dRapEnd = 5;
    
    int SRapBinNum = 1000;
    float SRapSta = -10 , SRapEnd = 10;
    
    int dPtBinNum = 200;
    float dPtSta = 0 , dPtEnd = 4;
    
    int detaBinNum = 200;
    float detaSta = -2 , detaEnd = 2;
    
    int MBinNum = 1000 , MBinPar = 100;
    float MSta = floor((AMass + BMass)/0.0005-MBinPar)*0.0005 , MEnd = MSta + (MBinNum - MBinPar)*0.0005;
    cout<<"Mass Region: [ "<<MSta<<" , "<<MEnd<<" ], BinNum = "<<MBinNum<<". "<<endl;

    std::vector<TH2F*> H_pT_rap;
    std::vector<TH2F*> H_pT_eta;
    H_pT_rap.resize(ParticleNum, nullptr);
    H_pT_eta.resize(ParticleNum, nullptr);
    for(int i=0;i<ParticleNum;i++){
        H_pT_rap[i] = new TH2F(Form("H_pT_rap_%d"   , Recording_Particle[i]),Form("pT vs. rap, %d"  ,Recording_Particle[i]),dPtBinNum,dPtSta,dPtEnd,dRapBinNum,dRapSta,dRapEnd);
        H_pT_eta[i] = new TH2F(Form("H_pT_eta_%d"   , Recording_Particle[i]),Form("pT vs. eta, %d"  ,Recording_Particle[i]),dPtBinNum,dPtSta,dPtEnd,detaBinNum,detaSta,detaEnd);
    }

    cout<<"Histogram initialized!"<<endl;

    // =========================
    // 打开并逐行读取文件
    // =========================
    TChain *hadronTree = new TChain(TreeName);
    std::ifstream infile(MidName.Data());

    if (!infile.is_open()) {
        std::cerr << "Error: cannot open file " << MidName << std::endl;
        return;
    }

    std::string line;
    int lineCount = 0;

    while (std::getline(infile, line)) {
        lineCount++;

        // 跳过空行（可选）
        if (line.empty()) continue;

        // 跳过注释（可选，比如 # 开头）
        if (line[0] == '#') continue;

        // 处理每一行
        std::cout << "Line " << lineCount << ": " << line << std::endl;

        // 如果你需要转成 TString：
        TString tline(line);
        TFile *f = TFile::Open(line.c_str());

        if (!f || f->IsZombie()) {
            std::cerr << "Bad file: " << line << std::endl;
            continue;
        }
        
        // 检查 tree 是否存在
        TTree *t = (TTree*)f->Get(TreeName);
        
        if (!t) {
            std::cerr << "No tree " << TreeName << " in " << line << std::endl;
            f->Close();
            // delete f;
            continue;
        }
        
        // 只有通过检查才加入
        hadronTree->Add(line.c_str());
        
        f->Close();
        delete f;

    }

    infile.close();

    std::cout << "Total lines read: " << lineCount << std::endl;

    // TChain *hadronTree = new TChain(TreeName);
    // for(i=StartFileIndex;i <= EndFileIndex;i++){
    //     TString filename = MidName;
    //     filename+=i;
    //     filename+=".root";
    //     hadronTree->Add(filename);
    // }
    Int_t PDGMult  ;
    Int_t refMult  ;
    Int_t grefMult ;
    Int_t EventID  ;
    Int_t RunID    ;
    Int_t TriggerID;
    Int_t Nch      ;
    float PVz      ;

    hadronTree->SetBranchAddress("PDGMult"  ,&PDGMult  );
    hadronTree->SetBranchAddress("refMult"  ,&refMult  );
    // hadronTree->SetBranchAddress("grefMult" ,&grefMult );
    // hadronTree->SetBranchAddress("EventID"  ,&EventID  );
    // hadronTree->SetBranchAddress("RunID"    ,&RunID    );
    // hadronTree->SetBranchAddress("TriggerID",&TriggerID);
    // hadronTree->SetBranchAddress("Nch"      ,&Nch      );
    hadronTree->SetBranchAddress("PVz"      ,&PVz      );
    
    hadronTree->SetBranchAddress("PDG"          ,&PDG          ,&bPDG          );
    hadronTree->SetBranchAddress("mix_px"       ,&mix_px       ,&bmix_px       );
    hadronTree->SetBranchAddress("mix_py"       ,&mix_py       ,&bmix_py       );
    hadronTree->SetBranchAddress("mix_pz"       ,&mix_pz       ,&bmix_pz       );
    // hadronTree->SetBranchAddress("QA_eta"       ,&QA_eta       ,&bQA_eta       );
    // hadronTree->SetBranchAddress("dEdx"         ,&dEdx         ,&bdEdx         );
    // hadronTree->SetBranchAddress("m2"           ,&m2           ,&bm2           );
    if(IfCutHighDCA) hadronTree->SetBranchAddress("dcatopv"      ,&dcatopv      ,&bdcatopv      );
    // hadronTree->SetBranchAddress("nSigmaProton" ,&nSigmaProton ,&bnSigmaProton );
    // hadronTree->SetBranchAddress("nSigmaPion"   ,&nSigmaPion   ,&bnSigmaPion   );
    if (IfRemoveHighTPCsigma && ((abs(A_PDG) == 321) || (abs(B_PDG) == 321))){
        hadronTree->SetBranchAddress("nSigmaKaon"   ,&nSigmaKaon   ,&bnSigmaKaon   );
    }
    hadronTree->SetBranchAddress("InvariantMass",&InvariantMass,&bInvariantMass);
    // hadronTree->SetBranchAddress("Decay_Length" ,&Decay_Length ,&bDecay_Length );
    // hadronTree->SetBranchAddress("Chi2"         ,&Chi2         ,&bChi2         );
    if (IfRemoveLownHits) {
        hadronTree->SetBranchAddress("nHitsFit"     ,&nHitsFit     ,&bnHitsFit     );
        hadronTree->SetBranchAddress("nHitsMax"     ,&nHitsMax     ,&bnHitsMax     );
    }
    hadronTree->SetBranchAddress("ParentList"   ,&ParentList   ,&bParentList   );
    hadronTree->SetBranchAddress("ParentSta"    ,&ParentSta    ,&bParentSta    );
    hadronTree->SetBranchAddress("ParentEnd"    ,&ParentEnd    ,&bParentEnd    );

    if (CutID == 1) IfRemoveLownHits = true;
    if (CutID == 2) IfRemoveHighPVz = true;
    if (CutID == 3) IfRemoveHighTPCsigma = true;
    if (CutID == 4) IfCutHighDCA = true;

    const Int_t nentries=hadronTree->GetEntries();
    cout << "file number: " << nentries << endl;

    time_t time_start;
    time_t time_now;
    time(&time_start);
    clock_t Tstart = clock();
    for (int EntriesID = 0 ; EntriesID < nentries ; EntriesID++){
        hadronTree->GetEntry(EntriesID);
        if ((EntriesID+1)%20000 == 0) {
            time(&time_now);
            int time_diff = (int)difftime(time_now, time_start);
            cout << time_diff/60 << "min " << time_diff%60 << "s: ";
            long long microseconds = (clock() - Tstart)/10000;
            std::cout << "Microseconds: " << microseconds << "  ";
            cout<<"Calculating Event "<<(EntriesID+1)<<"/"<<nentries<<endl;
            Tstart = clock();
        }

        if (IfRemoveHighPVz) {
            if (!((-25.0 <= PVz) && (PVz < 25.0))) continue;
        }

        // 遍历粒子，筛选A、B、C、D
        for (int i=0;i<PDGMult;i++){
            if (IfRemoveHighTPCsigma) {
                if (abs(PDG->at(i)) == 321) {
                    if (fabs(nSigmaKaon->at(i))>1) continue;
                }
            }
            if (IfRemoveLownHits) {
                if ((abs(PDG->at(i)) == 321) || (abs(PDG->at(i)) == 211) || (abs(PDG->at(i)) == 2212)) {
                    if (nHitsFit->at(i) < 20) continue;
                }
            }
            if (IfCutHighDCA) {
                if ((abs(PDG->at(i)) == 321) || (abs(PDG->at(i)) == 211) || (abs(PDG->at(i)) == 2212)) {
                    if ( (0 > dcatopv->at(i)) || (dcatopv->at(i) > 0.5)) continue;
                }
            }
            for(int j=0;j<ParticleNum;j++) {
                if (PDG->at(i) == Recording_Particle[j]) {
                    if (fabs(InvariantMass->at(i) - Particle_Mass[j]) <= MassSigmaWidth*Particle_MassSigma) {
                        A = ArmParticle(mix_px->at(i),mix_py->at(i),mix_pz->at(i),Particle_Mass[j],i);
                        H_pT_rap[i].Fill(A.pt,A.y);
                        H_pT_eta[i].Fill(A.pt,A.eta);
                    }
                    break;
                }
            }
        }
    }
    // 保存.root文件
    
    TString OutputFileName = OutMidName;
    OutputFileName += "H_";
    OutputFileName += OutputFileIndex;
    OutputFileName += ".root";
    TFile *fileA = new TFile(OutputFileName, "RECREATE");
    std::vector<TDirectory*> folder_Particle;
    folder_Particle.resize(ParticleNum,nullptr);
    for(int i=0;i<ParticleNum;i++) {
        folder_Particle[i]     = fileA->mkdir(TString(Recording_Particle[i]));
        folder_Particle[i]->cd();
        H_pT_rap[i]->Write();
        H_pT_eta[i]->Write();
    }
    fileA->Close();
    return;
}

int main(int argc, char** argv) {
    // 检查参数数量
    if(argc < 5) {
        std::cerr << "Usage: " << argv[0] 
                  << " MidName DataName OutputFileIndex OutMidName CutID" << std::endl;
        return 1;
    }

    Eff(
        TString(argv[1]),
        TString(argv[2]),
        atoi(argv[3]),
        TString(argv[4]),
        atoi(argv[5])
    );

    return 0;
}

Double_t massList(int PID, TString DataName)
{
    Double_t Result;
    if (DataName == "dAu_200_21"){
        switch (PID)
        {
            case 321 :
                Result = 0.493677;
                break;
            case -321 :
                Result = 0.493677;
                break;
            case 310 :
                Result = 0.49794;
                break;
            case 211 :
                Result = 0.13957;
                break;
            case -211 :
                Result = 0.13957;
                break;
            case 1003314 :// XiRPdgMass
                Result = 1.6725;
                break;
            case -1003314 :// XiRPdgMass
                Result = 1.6727;
                break;
            case 3334 :// OmegaFitMass
                Result = 1.6725;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 1.6727;
                break;
            case 3312 :// XiFitMass
                Result = 1.3223;
                break;
            case -3312 :// XiBarFitMass
                Result = 1.3223;
                break;
            case 3122 :// LambdaFitMass
                Result = 1.1161;
                break;
            case -3122 :// LambdaBarFitMass
                Result = 1.1161;
                break;
            default :
                Result = 0;
        }
    }
    if (DataName == "AuAu_19_19"){
        switch (PID)
        {
            case 321 :
                Result = 0.493677;
                break;
            case -321 :
                Result = 0.493677;
                break;
            case 310 :
                Result = 0.49794;
                break;
            case 211 :
                Result = 0.13957;
                break;
            case -211 :
                Result = 0.13957;
                break;
            case 1003314 :// XiRPdgMass
                Result = 1.6725;
                break;
            case -1003314 :// XiRPdgMass
                Result = 1.6727;
                break;
            case 3334 :// OmegaFitMass
                Result = 1.6725;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 1.6727;
                break;
            case 3312 :// XiFitMass
                Result = 1.3223;
                break;
            case -3312 :// XiBarFitMass
                Result = 1.3223;
                break;
            case 3122 :// LambdaFitMass
                Result = 1.1161;
                break;
            case -3122 :// LambdaBarFitMass
                Result = 1.1161;
                break;
            default :
                Result = 0;
        }
    }
    return Result;
}

Double_t massListSigma(int PID, TString DataName)
{
    Double_t Result;
    if (DataName == "dAu_200_21"){
        switch (PID)
        {
            case 3334 :// OmegaFitMass
                Result = 0.0029;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 0.0024;
                break;
            case 1003314 :// XiRPdgMass
                Result = 0.0029;
                break;
            case -1003314 :// XiRPdgMass
                Result = 0.0024;
                break;
            case 3312 :// XiFitMass
                Result = 0.0024;
                break;
            case -3312 :// XiBarFitMass
                Result = 0.0024;
                break;
            case 3122 :// LambdaFitMass
                Result = 0.0020;
                break;
            case -3122 :// LambdaBarFitMass
                Result = 0.0020;
                break;
            default :
                Result = 100;
        }
    }
    if (DataName == "dAu_62_16"){// tbd, used as dAu@200R21
        switch (PID)
        {
            case 3334 :// OmegaFitMass
                Result = 0.0029;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 0.0024;
                break;
            case 1003314 :// XiRPdgMass
                Result = 0.0029;
                break;
            case -1003314 :// XiRPdgMass
                Result = 0.0024;
                break;
            case 3312 :// XiFitMass
                Result = 0.0024;
                break;
            case -3312 :// XiBarFitMass
                Result = 0.0024;
                break;
            case 3122 :// LambdaFitMass
                Result = 0.0020;
                break;
            case -3122 :// LambdaBarFitMass
                Result = 0.0020;
                break;
            default :
                Result = 100;
        }
    }
    if (DataName == "AuAu_19_19"){
        switch (PID)
        {
            case 3334 :// OmegaFitMass
                Result = 0.0029;
                break;
            case -3334 :// OmegaBarFitMass
                Result = 0.0024;
                break;
            case 1003314 :// XiRPdgMass
                Result = 0.0029;
                break;
            case -1003314 :// XiRPdgMass
                Result = 0.0024;
                break;
            case 3312 :// XiFitMass
                Result = 0.0024;
                break;
            case -3312 :// XiBarFitMass
                Result = 0.0024;
                break;
            case 3122 :// LambdaFitMass
                Result = 0.0020;
                break;
            case -3122 :// LambdaBarFitMass
                Result = 0.0020;
                break;
            default :
                Result = 100;
        }
    }
    return Result;
}
