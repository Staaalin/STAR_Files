/**
 *  Author: Grigory Nigmatkulov
 *  Date:   July 5, 2018
 *
 *  Description:
 *  This macros takes inFileName argument with a picoDst.root file
 *  or with a list of files (name.lis or name.list). It sets _VANILLA_ROOT_
 *  (necessary for standalone and can be skipped on RACF), loads pre compiled
 *  libStPicoDst.so (from StPicoEvent), compiles and executes a text
 *  PicoDstAnalyzer.C macro with passing inFileName to it, and
 *  cleans up the directory from the compilation products at the end.
 *
 *  Some details:
 *    inFileName - is a name of name.picoDst.root file or a name
 *                 of a name.lis(t) files that contains a list of
 *                 name1.picoDst.root files.
 *    NOTE: inFileName should contain either /absolutePath/inFileName
 *          or /relative2currentDir/inFileName
 *  It is assumed that PicoDstAnalyzer.C is placed in the same
 *  directory where the RunAnalyzer.C is stored.
 **/

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

// #include "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent/StPicoDst.h"

class StRefMultCorr;
class CentralityMaker;

//_________________
void RunAnalyzer_QA(const int cen = 0, const int opt_weight = 0, const Char_t *inFileName = "test.list") {
  // Next line is not needed if you are not running in a standalone mode
//  gROOT->ProcessLine("#define _VANILLA_ROOT_");
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StEpdUtil");
  gSystem->Load("StRefMultCorr");
  TString str;
  str = ".x Gamma_QA.C++(";
  str += cen;
  str += ",";
  str += opt_weight;
  str += ",\"";
  str += inFileName;
  str += "\")";
//  cout<<str.Data()<<endl;
  gROOT->ProcessLine( str.Data() );
  // Next line should be commented if you run in a batch mode
  gROOT->ProcessLine(".!rm -f Gamma_QA_C* ");
  gROOT->ProcessLine(".!rm -f RunAnalyzer_QA_C* ");
}
