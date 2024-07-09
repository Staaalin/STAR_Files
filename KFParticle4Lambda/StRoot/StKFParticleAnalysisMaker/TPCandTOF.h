#include <TString.h>
#include <map>

class TPCandTOF
{
    private:
    TString DataName;
    public:
    TPCandTOF(TString _DataName): DataName(_DataName) {}
    float KaonTPCCenter(float Pt , TString DataName);
    std::vector<float> KaonTOFm2(float Pt, TString DataName);
};
