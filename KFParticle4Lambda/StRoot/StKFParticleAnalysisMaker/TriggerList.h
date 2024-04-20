#include <TString.h>
#include <map>

class TriggerList
{
    private:
    TString DataName;
    std::vector<TString> DataNameList;
    std::vector<std::vector<int> > TriggerNameList;
    public:
    TriggerList(TString _DataName): DataName(_DataName) {}
    std::vector<int> GetTriggerList();
};
