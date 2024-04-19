#include "TriggerList.h"
#include <cmath>

std::map<TString, vector<int> > TriggerList;
TriggerList["dAu_200_16"] = {1,
                             2,
                             3,
                             4,
                             5,
                             }

	if ((!mEvent->isTrigger(530003))
	  &&(!mEvent->isTrigger(530002))
	  &&(!mEvent->isTrigger(530806))
	  &&(!mEvent->isTrigger(530101))
	  &&(!mEvent->isTrigger(530102))
	  &&(!mEvent->isTrigger(530201))
	  &&(!mEvent->isTrigger(530202))
	  &&(!mEvent->isTrigger(530213))
	  &&(!mEvent->isTrigger(530851))// bbc
	  &&(!mEvent->isTrigger(530852))// zdce
	  &&(!mEvent->isTrigger(530853))// vpd-10
	  &&(!mEvent->isTrigger(530854))// vpd
	  &&(!mEvent->isTrigger(530855))// zdc
	  &&(!mEvent->isTrigger(530856))// bbcnotac
	  &&(!mEvent->isTrigger(530857))// zdc-notac
	  &&(!mEvent->isTrigger(530861))// bbc
	  &&(!mEvent->isTrigger(530866))// bbcnotac
	  &&(!mEvent->isTrigger(2)) //17132063
	  &&(!mEvent->isTrigger(3))
	  &&(!mEvent->isTrigger(4)) //VPD-5
	  &&(!mEvent->isTrigger(6)) //highMult-VPD-5
	  &&(!mEvent->isTrigger(15))
	  &&(!mEvent->isTrigger(16)) //BHT2-VPD-30
	  &&(!mEvent->isTrigger(17))
	  &&(!mEvent->isTrigger(54))
	  &&(!mEvent->isTrigger(55))
	  &&(!mEvent->isTrigger(56))
	  &&(!mEvent->isTrigger(57))
	  &&(!mEvent->isTrigger(58))
	  &&(!mEvent->isTrigger(59))
	  &&(!mEvent->isTrigger(61)) //vpd-30
	  )return kStOK;