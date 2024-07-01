using namespace std;

#include "stdio.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

void Maker() {
    ifstream inrun;

    StPicoDstMaker* maker = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
	if (! maker) return;
	maker->SetStatus("*",1);
	TChain *tree = maker->chain();
	Long64_t nentries = tree->GetEntries();
	if (nentries <= 0) return;
	Long64_t EventsToRun = TMath::Min(nEvents,nentries);
	cout << nentries << " events in chain " << EventsToRun << " will be read." << endl;
}
