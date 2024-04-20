#include <TString.h>
#include <map>

class TriggerList
{
    private:
    TString DataName;
    std::map<TString,std::vector<int>> DataTriggerList;
    public:
    TriggerList(TString _DataName): DataName(_DataName) {}
    std::vector<int> GetTriggerList();
};
