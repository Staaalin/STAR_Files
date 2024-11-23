#include "vector"
#include "vector"
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<vector<string,allocator<string> > >+;
#pragma link C++ class vector<vector<string,allocator<string> > >::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<vector<string,allocator<string> > >::iterator;
#pragma link C++ operators vector<vector<string,allocator<string> > >::const_iterator;
#pragma link C++ operators vector<vector<string,allocator<string> > >::reverse_iterator;
#endif
#endif
