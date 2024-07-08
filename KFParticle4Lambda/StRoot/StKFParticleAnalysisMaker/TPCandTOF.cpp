#include "TriggerList.h"
#include <TString.h>
#include <map>
#include <fstream>
#include <iostream>

static const std::vector<TString> TriggerList::DataNameList = {
	"dAu_200_16",
	"pAu_200_15",
	"dAu_62_16",
	"dAu_39_16",
	"dAu_20_16",
	"AuAu_27_18",
	"dAu_200_21",
	"pp_200_15",
	"OO_200_21"
};

float TriggerList::GetTPCCenter(flout Pt , TString DataName)
{
	if (DataName == "dAu_200_21") {
		if      ( 0.2 <= Pt && Pt < 0.3) {return 0.8}
		else if ( 0.3 <= Pt && Pt < 0.4) {return 0.4}
		else if ( 0.4 <= Pt && Pt < 0.5) {return 0.1}
		else if ( 0.5 <= Pt && Pt < 0.6) {return -0.06}
		else if ( 0.6 <= Pt && Pt < 0.7) {return -0.08}
		else                             {return 0.0}
	}
}

float TriggerList::GetTOF(flout Pt , TString DataName)
{
	if (DataName == "dAu_200_21") {
		if      ( 0.2 <= Pt && Pt < 0.3) {return 0.8}
		else if ( 0.3 <= Pt && Pt < 0.4) {return 0.4}
		else if ( 0.4 <= Pt && Pt < 0.5) {return 0.1}
		else if ( 0.5 <= Pt && Pt < 0.6) {return -0.06}
		else if ( 0.6 <= Pt && Pt < 0.7) {return -0.08}
		else                             {return 0.0}
	}
}