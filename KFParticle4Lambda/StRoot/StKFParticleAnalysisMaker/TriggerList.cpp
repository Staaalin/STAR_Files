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
	"AuAu_27_18"
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
		540000,
		540002,
		540003,
		540101,
		540102,
		540112,
		540201,
		540203,
		540204,
		540601,
		540602,
		540603,
		540604,
		540611,
		540614,
		540701,
		540702,
		540703,
		540851,
		540852,
		540854,
		540855,
		540856,
		540857,
		540863
	},

	{// dAu_39_16
		560000,
		560001,
		560002,
		560007,
		560008,
		560011,
		560017,
		560018,
		560853
	},

	{// dAu_20_16
		8,
		550001,
		550003,
		550007,
		550009,
		550011,
		550018
	},

	{// AuAu_27_18
		610001,
		610011,
		610021,
		610031,
		610041,
		610051
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