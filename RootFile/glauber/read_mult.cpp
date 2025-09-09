#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

int main(int argc, char** argv) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " <FilePath> <IndexStart> <IndexEnd>" << endl;
        return 1;
    }

    TString FilePath = argv[1];
    int IndexStart = atoi(argv[2]);
    int IndexEnd = atoi(argv[3]);

    // 创建直方图
    TH1D *h = new TH1D("mult", "mult", 500, 0, 500);

    // 循环读取文件
    for (int i = IndexStart; i <= IndexEnd; i++) {
        // 拼接文件名：四位补零
        std::ostringstream oss;
        oss << FilePath << std::setw(4) << std::setfill('0') << i << ".root";
        TString filename = oss.str();

        cout << "Opening file: " << filename << endl;

        // 打开ROOT文件
        TFile *file = TFile::Open(filename);
        if (!file || file->IsZombie()) {
            cout << "Cannot open file: " << filename << endl;
            continue;
        }

        // 获取树
        TTree *tree = (TTree*)file->Get("tree");
        if (!tree) {
            cout << "Tree not found in file: " << filename << endl;
            file->Close();
            continue;
        }

        // 设置分支地址
        UInt_t mult = 0;
        tree->SetBranchAddress("mult", &mult);

        // 读取条目并填充直方图
        Long64_t nentries = tree->GetEntries();
        for (Long64_t j = 0; j < nentries; j++) {
            tree->GetEntry(j);
            h->Fill(mult);
        }

        file->Close();
    }

    // 保存直方图
    TFile *outFile = new TFile("mult.root", "RECREATE");
    h->Write();
    outFile->Close();

    cout << "Histogram saved to mult.root" << endl;
    return 0;
}
