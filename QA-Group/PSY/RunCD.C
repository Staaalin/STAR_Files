
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoDst.h"

class StRefMultCorr;
class CentralityMaker;

//_________________
void RunCD(const Char_t *inFile = "test.list") {
  // Next line is not needed if you are not running in a standalone mode
//  gROOT->ProcessLine("#define _VANILLA_ROOT_");
  // gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  // gSystem->Load("StPicoEvent");
  // gSystem->Load("StEpdUtil");
  // gSystem->Load("StRefMultCorr");
  TString str;
  str = ".x CD.C++(";
  str += "\"";
  str += inFile;
  str += "\")";
//  cout<<str.Data()<<endl;
  gROOT->ProcessLine( str.Data() );
}
