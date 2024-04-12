#ifndef StCut_h
#define StCut_h

#include "TObject.h"
#include "TString.h"
//#include <array>

class StPicoEvent;
class StPicoTrack;
class StPicoBTofPidTraits;
class StRefMultCorr;

using namespace std;

class StCut : public TObject
{
  public:
    StCut();
    ~StCut();

    bool passEventCut(StPicoEvent*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackEP(StPicoTrack*, float);

  private:
    ClassDef(StCut,1)
};
#endif
