
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"


class StRefMultCorr;
class CentralityMaker;

//_________________
void RunMaker(const Char_t *inFileName = "test.list") {
  // Next line is not needed if you are not running in a standalone mode
//  gROOT->ProcessLine("#define _VANILLA_ROOT_");
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StEpdUtil");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StPicoDst");
  TString str;
  str = ".x Maker.C++(";
  str += "\")";
//  cout<<str.Data()<<endl;
  gROOT->ProcessLine( str.Data() );
  // Next line should be commented if you run in a batch mode
  gROOT->ProcessLine(".!rm -f Maker_C* ");
  gROOT->ProcessLine(".!rm -f RunMaker_C* ");
}
