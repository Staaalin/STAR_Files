#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StEpdUtil//StEpdGeom.h"
// ... 你需要的头文件继续加

void CD(int File_Index , const Char_t *inFile = "test.list"); // 原始 CD 函数声明

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: ./runCD [File_Index] [input.list]" << std::endl;
        return 1;
    }

    int fileIndex = atoi(argv[1]);
    const char* inFile = argv[2];

    // 加载 STAR 所需的库
    gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    gSystem->Load("StPicoEvent");
    gSystem->Load("StEpdUtil");
    gSystem->Load("StRefMultCorr");

    // 调用分析函数
    CD(fileIndex, inFile);
    return 0;
}

void CD(int File_Index , const Char_t *inFile = "test.list") {

  const float PI = TMath::Pi();
  Bool_t readEvent;
  StPicoDst *dst;
  StPicoEvent *event;
  int Run	       ;
  TVector3 pV    ;
  float pVz      ;
  float pVx      ;
  float pVy      ;
  float VPDvz    ;
  int RefMult    ;
  int FxtMult    ;
  int TOFMult    ;
  int Ntofmatch  ;
  int NPTracks   ;
  float BBCco    ;
  float ZDCcoin  ;
  int NumCharge  ;
  
  TH2F *hTofMatch_vs_RefMult = new TH2F("hTofMatch_vs_RefMult","nbTofMatch_vs_RefMult",500,0,500,500,0,500);
  hTofMatch_vs_RefMult->GetXaxis()->SetTitle("nBTOFMatch");
  hTofMatch_vs_RefMult->GetYaxis()->SetTitle("RefMult");

  cout<<"Start"<<endl;
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  gSystem->AddIncludePath("-I$STAR/StRoot/StarClassLibrary");
gSystem->Load("StUtilities");
  gSystem->Load("StEpdUtil");
gSystem->Load("StRefMultCorr");
gSystem->Load("StPicoEvent");
gSystem->Load("StPicoDstMaker");
  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();
  cout<<"Finish initting"<<endl;
  if( !picoReader->chain() ) { std::cout << "No chain has been found." << std::endl; }
  Int_t nentries = picoReader->chain()->GetEntries();

  for(int i = 0; i < nentries; i++) {

if((i+1)%1000==0) cout<<"Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
readEvent = picoReader->readPicoEvent(i);
  if( !readEvent ) {
      cout << "Something went wrong, Master! Nothing to analyze..." << endl;
      break;
  }

          // Retrieve picoDst
dst = picoReader->picoDst();
// Retrieve event information
event = dst->event();
  if( !event ) {
      cout << "Something went wrong, Master! Event is hiding from me..." << endl;
      break;
  }

          if( !event->isTrigger(630052) ) continue; // AuAu26p5

Run	  = event->runId();
pV 	  = event->primaryVertex();
pVz	  = pV.Z();
pVx	  = pV.X();
pVy	  = pV.Y();
VPDvz     = event->vzVpd();
RefMult   = event->refMult();
FxtMult   = event->fxtMult();
TOFMult   = event->btofTrayMultiplicity();
Ntofmatch = event->nBTOFMatch();
NPTracks  = dst->numberOfTracks();
BBCco     = event->BBCx();
          ZDCcoin   = event->ZDCx();
          
          
          // Siyuan Ping: Reject Bad Run
          if (Run == 19168001) continue;
          if (Run == 19168002) continue;
          if (Run == 19168003) continue;
          if (Run == 19168016) continue;
          if (Run == 19167054) continue;
          if (Run == 19166003) continue;
          if (Run == 19164002) continue;
          if (Run == 19164021) continue;
          if (Run == 19161031) continue;
          if (Run == 19161032) continue;
          if (Run == 19161033) continue;
          if (Run == 19159042) continue;
          if (Run == 19158053) continue;
          if (Run == 19158054) continue;
          if (Run == 19158055) continue;
          if (Run == 19158056) continue;
          if (Run == 19157033) continue;
          if (Run == 19157034) continue;
          if (Run == 19157035) continue;
          if (Run == 19157036) continue;
          if (Run == 19157037) continue;
          if (Run == 19157038) continue;
          if (Run == 19157039) continue;
          if (Run == 19157040) continue;
          if (Run == 19157041) continue;
          if (Run == 19157042) continue;
          if (Run == 19157043) continue;
          if (Run == 19156034) continue;
          if (Run == 19156035) continue;
          if (Run == 19156036) continue;
          if (Run == 19156038) continue;
          if (Run == 19156039) continue;
          if (Run == 19156069) continue;

          NumCharge = 0;
          for (Int_t iTrack = 0; iTrack < NPTracks; iTrack++) {
                  StPicoTrack *track = dst->track(iTrack);
                  if (! track)            continue;
                  if (! track->charge())  continue;
                  if (! track->isPrimary()) continue;
                  NumCharge++;
          }
          RefMult = NumCharge;

          hTofMatch_vs_RefMult->Fill(Ntofmatch,RefMult);
  }

  TFile *outFile = new TFile("CenDef.root", "RECREATE");
  hTofMatch_vs_RefMult->Write();
  outFile->Close();

}