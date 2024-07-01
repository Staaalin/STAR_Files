#ifndef StRecenter_h
#define StRecenter_h

#include "StMaker.h"
#include "TString.h"
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StCut;
class StProManger;
class StCorrection;
class StRefMultCorr;

class TProfile;
class TH1F;
class TH2F;

class StRecenter : public StMaker {
    public:
        StRecenter(const char *name, StPicoDstMaker *picoMaker, char* jobid);
        virtual ~StRecenter();

        virtual Int_t Init();
        virtual Int_t Make();
        virtual void  Clear(Option_t *opt="");
        virtual Int_t Finish();

        int Centrality(int gRefMult);

    private:
        StRefMultCorr *mRefMultCorrUtil;
        StPicoDstMaker *mPicoDstMaker;
        StPicoDst      *mPicoDst;
        StPicoEvent *event;
        StCut *mCut;
        StProManger *mProManger;
        StCorrection *mCorrection;

        TString mInPut_Corr_ReCenter;
        TString mOutPut_ReCenter;
        TString mOutPut_Shift;

        TFile *mFile_Corr_ReCenter;
        TFile *mFile_Corr_Shift;

        ClassDef(StRecenter, 1)
};

#endif
