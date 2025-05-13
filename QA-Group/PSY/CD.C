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

/// PicoDst headers
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
// #include "StRoot/StPicoEvent/StPicoBTofHit.h"
// #include "StRoot/StPicoEvent/StPicoBTowHit.h"
// #include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
// #include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
// #include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
// #include "StRoot/StEpdUtil/StEpdEpFinder.h"
// #include "StRoot/StRefMultCorr/StRefMultCorr.h"
// #include "StRoot/StRefMultCorr/CentralityMaker.h"
// #include "StRoot/StPicoEvent/StPicoEpdHit.h"
// #include "StRoot/StEpdUtil//StEpdGeom.h"

//class StRefMultCorr;
//class CentralityMaker;
class StPicoDstReader;

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


void CD(const Char_t *inFile = "test.list") {

        TH2F *hTofMatch_vs_RefMult_Roop = new TH2F("hTofMatch_vs_RefMult_Roop","nbTofMatch_vs_RefMult(Calculated from loop)",250,0,250,500,0,500);
        hTofMatch_vs_RefMult_Roop->GetXaxis()->SetTitle("nBTOFMatch");
        hTofMatch_vs_RefMult_Roop->GetYaxis()->SetTitle("RefMult");

        TH2F *hTofMatch_vs_RefMult_Ref = new TH2F("hTofMatch_vs_RefMult_Ref","nbTofMatch_vs_RefMult(Calculated from ->RefMult())",250,0,250,500,0,500);
        hTofMatch_vs_RefMult_Ref->GetXaxis()->SetTitle("nBTOFMatch");
        hTofMatch_vs_RefMult_Ref->GetYaxis()->SetTitle("RefMult");

        TH2F *hTofMatch_vs_RefMult_Fxt = new TH2F("hTofMatch_vs_RefMult_Fxt","nbTofMatch_vs_RefMult(Calculated from ->FxtMult())",250,0,250,500,0,500);
        hTofMatch_vs_RefMult_Fxt->GetXaxis()->SetTitle("nBTOFMatch");
        hTofMatch_vs_RefMult_Fxt->GetYaxis()->SetTitle("RefMult");

        cout<<"Start"<<endl;
        gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
        gSystem->AddIncludePath("-I$STAR/StRoot/StarClassLibrary");
	// gSystem->Load("StUtilities");
        // gSystem->Load("StEpdUtil");
	// gSystem->Load("StRefMultCorr");
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

                if (fabs(pVz-200.0)>2.0) continue;
                if (PVx*PVx+PVy*PVy>4.0) continue;
                
                
                // Siyuan Ping: Reject Bad Run
                // from Run Log assessment
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
                // from QA Rejected
                if (Run == 19159043) continue;
                if (Run == 19159044) continue;
                if (Run == 19159046) continue;
                if (Run == 19160032) continue;
                if (Run == 19160033) continue;
                if (Run == 19160034) continue;
                if (Run == 19160035) continue;
                if (Run == 19160036) continue;
                if (Run == 19160037) continue;
                if (Run == 19160038) continue;
                if (Run == 19160039) continue;
                if (Run == 19160040) continue;
                if (Run == 19160041) continue;
                if (Run == 19160042) continue;
                if (Run == 19160043) continue;
                if (Run == 19160044) continue;
                if (Run == 19161001) continue;
                if (Run == 19161020) continue;
                if (Run == 19161021) continue;
                if (Run == 19161022) continue;
                if (Run == 19161023) continue;
                if (Run == 19161024) continue;
                if (Run == 19161025) continue;
                if (Run == 19161026) continue;
                if (Run == 19161027) continue;
                if (Run == 19161028) continue;
                if (Run == 19161029) continue;
                if (Run == 19161030) continue;
                if (Run == 19161034) continue;
                if (Run == 19161035) continue;
                if (Run == 19161036) continue;
                if (Run == 19161037) continue;
                if (Run == 19161038) continue;
                if (Run == 19161042) continue;
                if (Run == 19162033) continue;
                if (Run == 19162034) continue;
                if (Run == 19164001) continue;
                if (Run == 19164022) continue;
                if (Run == 19164023) continue;
                if (Run == 19164024) continue;
                if (Run == 19164025) continue;
                if (Run == 19167050) continue;
                if (Run == 19167051) continue;
                if (Run == 19167052) continue;
                if (Run == 19167053) continue;
                if (Run == 19168041) continue;
                if (Run == 19168042) continue;

                NumCharge = 0;
                for (Int_t iTrack = 0; iTrack < NPTracks; iTrack++) {
                        StPicoTrack *track = dst->track(iTrack);
                        if (! track)            continue;
                        if (! track->charge())  continue;
                        if (! track->isPrimary()) continue;
                        if (track->gMom().Perp() < 0.06 || track->gMom().Perp() > 2.0) continue;
                        if (fabs(track->gMom().Eta()) > 1.5) continue;
                        if (fabs(track->gMom().Mag()) < 0.1) continue;
                        NumCharge++;
                }

                hTofMatch_vs_RefMult_Roop->Fill(Ntofmatch,NumCharge);
                hTofMatch_vs_RefMult_Ref ->Fill(Ntofmatch,RefMult);
                hTofMatch_vs_RefMult_Fxt ->Fill(Ntofmatch,FxtMult);
        }

        TFile *outFile = new TFile("cen1.v2.root", "RECREATE");
        hTofMatch_vs_RefMult_Roop->Write();
        hTofMatch_vs_RefMult_Ref ->Write();
        hTofMatch_vs_RefMult_Fxt ->Write();
        outFile->Close();

}