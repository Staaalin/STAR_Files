#include "TriggerList.h"
#include <TString.h>
#include <map>
#include <fstream>
#include <iostream>

	// "dAu_200_16",
	// "pAu_200_15",
	// "dAu_62_16",
	// "dAu_39_16",
	// "dAu_20_16",
	// "AuAu_27_18",
	// "dAu_200_21",
	// "pp_200_15",
	// "OO_200_21"

float KaonTPCCenter(flaot Pt , TString DataName)
{
	if (DataName == "dAu_200_21") {
		if      ( 0.2 <= Pt && Pt < 0.3) {return 0.8  ;}
		else if ( 0.3 <= Pt && Pt < 0.4) {return 0.4  ;}
		else if ( 0.4 <= Pt && Pt < 0.5) {return 0.1  ;}
		else if ( 0.5 <= Pt && Pt < 0.6) {return -0.06;}
		else if ( 0.6 <= Pt && Pt < 0.7) {return -0.08;}
		else                             {return 0.0  ;}
	}
}

std::vector<float> KaonTOFm2(float Pt, TString DataName)
{
    std::vector<float> result;
    if (DataName == "dAu_200_21") {
        if      (0.2 <= Pt && Pt < 0.3) { result = {0.228, 0.261}; }
        else if (0.3 <= Pt && Pt < 0.4) { result = {0.227, 0.259}; }
        else if (0.4 <= Pt && Pt < 0.5) { result = {0.224, 0.261}; }
        else if (0.5 <= Pt && Pt < 0.6) { result = {0.219, 0.264}; }
        else if (0.6 <= Pt && Pt < 0.7) { result = {0.213, 0.268}; }
        else if (0.7 <= Pt && Pt < 0.8) { result = {0.206, 0.273}; }
        else if (0.8 <= Pt && Pt < 0.9) { result = {0.196, 0.282}; }
        else if (0.9 <= Pt && Pt < 1.0) { result = {0.185, 0.290}; }
        else if (1.0 <= Pt && Pt < 1.1) { result = {0.172, 0.299}; }
        else if (1.1 <= Pt && Pt < 1.2) { result = {0.152, 0.313}; }
        else if (1.2 <= Pt && Pt < 1.3) { result = {0.131, 0.326}; }
        else if (1.3 <= Pt && Pt < 1.4) { result = {0.115, 0.335}; }
        else                            { result = {0.115, 0.335}; }
    }
    return result;
}
