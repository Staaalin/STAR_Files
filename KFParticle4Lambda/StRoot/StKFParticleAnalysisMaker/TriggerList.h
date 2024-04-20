#include <TString.h>
#include <map>

class TriggerList
{
    private:
    TString DataName;
    static const std::vector<TString> DataNameList;
    static const std::vector<std::vector<int> > TriggerNameList;
    public:
    TriggerList(TString _DataName): DataName(_DataName) {}
    std::vector<int> GetTriggerList();
};
