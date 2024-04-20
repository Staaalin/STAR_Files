#include "TriggerList.h"
#include <TString.h>
#include <map>

static const std::vector<TString> TriggerList::DataNameList = {
	"dAu_200_16",
	"pAu_200_15"
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
		-2,
		-3
	}

};

std::vector<int> TriggerList::GetTriggerList()
{
	for (int i=0;i<DataNameList.size();i++){
		if (DataNameList[i] == TriggerList::DataName){
			std::vector<int> result = TriggerNameList[i];
			return result;
		}
	}
	std::cout<<"NOT FOUND TRIGGER LIST!!"<<std::endl;
}