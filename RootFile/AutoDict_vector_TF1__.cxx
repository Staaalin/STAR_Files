#include "vector"
class TF1;
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<TF1*>+;
#pragma link C++ class vector<TF1*>::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<TF1*>::iterator;
#pragma link C++ operators vector<TF1*>::const_iterator;
#pragma link C++ operators vector<TF1*>::reverse_iterator;
#endif
#endif
