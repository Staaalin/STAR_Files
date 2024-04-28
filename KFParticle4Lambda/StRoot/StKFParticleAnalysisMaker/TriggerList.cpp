#include "TriggerList.h"
#include <TString.h>
#include <map>
#include <fstream>
#include <iostream>

static const std::vector<TString> TriggerList::DataNameList = {
	"dAu_200_16",
	"pAu_200_15",
	"dAu_62_16"
};

static const std::vector<std::vector<int>> TriggerList::TriggerNameList = {

	{// dAu_200_16
		2,
		3,
		4,
		6,
		15,
		16,
		17,
		54,
		55,
		56,
		57,
		58,
		59,
		61,
		530003,
		530002,
		530806,
		530101,
		530102,
		530201,
		530202,
		530213,
		530851,
		530852,
		530853,
		530854,
		530855,
		530856,
		530857,
		530861,
		530866
	},

	{// pAu_200_15
		20,
		21,
		22,
		27,
		29,
		30,
		31,
		32,
		35,
		48,
		57,
		59,
		510001,
		510003,
		510004,
		510008,
		510009,
		510011,
		510201,
		510202,
		510205,
		510206,
		510301,
		510406,
		510501,
		510602,
		510607,
		510806,
		510807,
		510808,
		510904
	},

	{// dAu_62_16
		3,
		540003,
		540101,
		540112,
		540203,
		540204,
		540603,
		540611,
		540701,
		540806,
		540852,
		540863
	}

};

std::vector<int> TriggerList::GetTriggerList()
{
	int Length = DataNameList.size();
	for (int i=0;i<Length;i++){
		if (DataNameList[i] == TriggerList::DataName){
			std::vector<int> result = TriggerNameList[i];
			return result;
		}
	}
	std::cout<<"NOT FOUND TRIGGER LIST!!"<<std::endl;
}