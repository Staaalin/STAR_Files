#include <TString.h>
#include <map>

class TriggerList
{
    private:
    TString DataName;
    public:
    std::vector<TString> DataNameList;
    std::vector<std::vector<int> > TriggerNameList;
    TriggerList(TString _DataName): DataName(_DataName) {}
    std::vector<int> GetTriggerList();
};
