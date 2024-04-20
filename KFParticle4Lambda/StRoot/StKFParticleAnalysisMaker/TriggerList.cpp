#include "TriggerList.h"
#include <TString.h>
#include <map>

TriggerList::DataNameList.clear();
TriggerList::TriggerNameList.clear();

TriggerList::DataNameList.emplace_back("dAu_200_16");
TriggerList::TriggerNameList.emplace_back({
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
								});

TriggerList::DataNameList.emplace_back("pAu_200_15");
TriggerList::TriggerNameList.emplace_back({
								-1,
								-10,
								});

std::vector<int> TriggerList::GetTriggerList()
{
	for (int i=0;i<DataNameList.size();i++){
		if (DataNameList[i] == DataName){
			std::vector<int> result = TriggerNameList[DataName];
			break;
		}
	}
	return result;
}