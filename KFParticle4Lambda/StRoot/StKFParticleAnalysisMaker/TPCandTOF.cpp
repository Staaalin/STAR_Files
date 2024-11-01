#include "TPCandTOF.h"
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

float TPCandTOF::KaonTPCCenter(float Pt , TString DataName)
{
	if (DataName == "dAu_200_21") {
		if      ( 0.2 <= Pt && Pt < 0.3) {return 0.8  ;}
		else if ( 0.3 <= Pt && Pt < 0.4) {return 0.4  ;}
		else if ( 0.4 <= Pt && Pt < 0.5) {return 0.1  ;}
		else if ( 0.5 <= Pt && Pt < 0.6) {return -0.06;}
		else if ( 0.6 <= Pt && Pt < 0.7) {return -0.08;}
		else                             {return 0.0  ;}
	}
    else{
        return 0.0;
    }
}

std::vector<float> TPCandTOF::KaonTOFm2(float Pt, TString DataName)
{
    std::vector<float> result;
    if (DataName == "dAu_200_21") {
        if      (0.2 <= Pt && Pt < 0.3) { result = {0.229 , 0.267}; }
        else if (0.3 <= Pt && Pt < 0.4) { result = {0.228 , 0.264}; }
        else if (0.4 <= Pt && Pt < 0.5) { result = {0.225 , 0.265}; }
        else if (0.5 <= Pt && Pt < 0.6) { result = {0.221 , 0.268}; }
        else if (0.6 <= Pt && Pt < 0.7) { result = {0.216 , 0.272}; }
        else if (0.7 <= Pt && Pt < 0.8) { result = {0.208 , 0.278}; }
        else if (0.8 <= Pt && Pt < 0.9) { result = {0.202 , 0.284}; }
        else if (0.9 <= Pt && Pt < 1.0) { result = {0.192 , 0.291}; }
        else if (1.0 <= Pt && Pt < 1.1) { result = {0.181 , 0.300}; }
        else if (1.1 <= Pt && Pt < 1.2) { result = {0.170 , 0.309}; }
        else if (1.2 <= Pt && Pt < 1.3) { result = {0.152 , 0.322}; }
        else if (1.3 <= Pt && Pt < 1.4) { result = {0.130 , 0.337}; }
        else                            { result = {0.130 , 0.337}; }
        return result;
    }
    else if (DataName == "dAu_200_16") { // tbd
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
        return result;
    }
    else{
        result = {0.115, 0.335};
        return result;
    }
}

std::vector<float> TPCandTOF::ProtonTOFm2(float Pt, TString DataName)
{
    std::vector<float> result;
    if (DataName == "dAu_200_21") {
        if      (0.2 <= Pt && Pt < 0.3) { result = {0.668 , 1.024}; }
        else if (0.3 <= Pt && Pt < 0.4) { result = {0.797 , 0.959}; }
        else if (0.4 <= Pt && Pt < 0.5) { result = {0.809 , 0.952}; }
        else if (0.5 <= Pt && Pt < 0.6) { result = {0.813 , 0.947}; }
        else if (0.6 <= Pt && Pt < 0.7) { result = {0.813 , 0.947}; }
        else if (0.7 <= Pt && Pt < 0.8) { result = {0.811 , 0.949}; }
        else if (0.8 <= Pt && Pt < 0.9) { result = {0.806 , 0.954}; }
        else if (0.9 <= Pt && Pt < 1.0) { result = {0.799 , 0.959}; }
        else if (1.0 <= Pt && Pt < 1.1) { result = {0.790 , 0.966}; }
        else if (1.1 <= Pt && Pt < 1.2) { result = {0.781 , 0.975}; }
        else if (1.2 <= Pt && Pt < 1.3) { result = {0.770 , 0.984}; }
        else if (1.3 <= Pt && Pt < 1.4) { result = {0.758 , 0.994}; }
        else                            { result = {0.745 , 1.005}; }
        return result;
    }
    else if (DataName == "dAu_200_16") { // tbd
        if      (0.2 <= Pt && Pt < 0.3) { result = {0.668 , 1.024}; }
        else if (0.3 <= Pt && Pt < 0.4) { result = {0.797 , 0.959}; }
        else if (0.4 <= Pt && Pt < 0.5) { result = {0.809 , 0.952}; }
        else if (0.5 <= Pt && Pt < 0.6) { result = {0.813 , 0.947}; }
        else if (0.6 <= Pt && Pt < 0.7) { result = {0.813 , 0.947}; }
        else if (0.7 <= Pt && Pt < 0.8) { result = {0.811 , 0.949}; }
        else if (0.8 <= Pt && Pt < 0.9) { result = {0.806 , 0.954}; }
        else if (0.9 <= Pt && Pt < 1.0) { result = {0.799 , 0.959}; }
        else if (1.0 <= Pt && Pt < 1.1) { result = {0.790 , 0.966}; }
        else if (1.1 <= Pt && Pt < 1.2) { result = {0.781 , 0.975}; }
        else if (1.2 <= Pt && Pt < 1.3) { result = {0.770 , 0.984}; }
        else if (1.3 <= Pt && Pt < 1.4) { result = {0.758 , 0.994}; }
        else                            { result = {0.745 , 1.005}; }
        return result;
    }
    else{
        result = {0.745 , 1.005};
        return result;
    }
}


std::vector<float> TPCandTOF::PionTOFm2(float Pt, TString DataName)
{
    std::vector<float> result;
    if (DataName == "dAu_200_21") {
        if      (0.2 <= Pt && Pt < 0.3) { result = { 0.0150 , 0.0263}; }
        else if (0.3 <= Pt && Pt < 0.4) { result = { 0.0122 , 0.0273}; }
        else if (0.4 <= Pt && Pt < 0.5) { result = { 0.0083 , 0.0301}; }
        else if (0.5 <= Pt && Pt < 0.6) { result = { 0.0034 , 0.0340}; }
        else if (0.6 <= Pt && Pt < 0.7) { result = {-0.0016 , 0.0379}; }
        else if (0.7 <= Pt && Pt < 0.8) { result = {-0.0093 , 0.0439}; }
        else if (0.8 <= Pt && Pt < 0.9) { result = {-0.0169 , 0.0499}; }
        else if (0.9 <= Pt && Pt < 1.0) { result = {-0.0254 , 0.0566}; }
        else if (1.0 <= Pt && Pt < 1.1) { result = {-0.0349 , 0.0640}; }
        else if (1.1 <= Pt && Pt < 1.2) { result = {-0.0453 , 0.0720}; }
        else if (1.2 <= Pt && Pt < 1.3) { result = {-0.0566 , 0.0807}; }
        else if (1.3 <= Pt && Pt < 1.4) { result = {-0.0688 , 0.0897}; }
        else                            { result = {-0.0813 , 0.0985}; }
        return result;
    }
    else if (DataName == "dAu_200_16") { // tbd
        if      (0.2 <= Pt && Pt < 0.3) { result = { 0.0150 , 0.0263}; }
        else if (0.3 <= Pt && Pt < 0.4) { result = { 0.0122 , 0.0273}; }
        else if (0.4 <= Pt && Pt < 0.5) { result = { 0.0083 , 0.0301}; }
        else if (0.5 <= Pt && Pt < 0.6) { result = { 0.0034 , 0.0340}; }
        else if (0.6 <= Pt && Pt < 0.7) { result = {-0.0016 , 0.0379}; }
        else if (0.7 <= Pt && Pt < 0.8) { result = {-0.0093 , 0.0439}; }
        else if (0.8 <= Pt && Pt < 0.9) { result = {-0.0169 , 0.0499}; }
        else if (0.9 <= Pt && Pt < 1.0) { result = {-0.0254 , 0.0566}; }
        else if (1.0 <= Pt && Pt < 1.1) { result = {-0.0349 , 0.0640}; }
        else if (1.1 <= Pt && Pt < 1.2) { result = {-0.0453 , 0.0720}; }
        else if (1.2 <= Pt && Pt < 1.3) { result = {-0.0566 , 0.0807}; }
        else if (1.3 <= Pt && Pt < 1.4) { result = {-0.0688 , 0.0897}; }
        else                            { result = {-0.0813 , 0.0985}; }
        return result;
    }
    else{
        result = {-0.0813 , 0.0985};
        return result;
    }
}
