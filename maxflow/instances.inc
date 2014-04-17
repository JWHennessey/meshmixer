#include "graph.h"
#include "../MyMesh.hh"

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

// Instantiations: <captype, tcaptype, flowtype>
// IMPORTANT: 
//    flowtype should be 'larger' than tcaptype 
//    tcaptype should be 'larger' than captype


template class Graph<MyMesh, int,int,int>;
template class Graph<MyMesh, short,int,int>;
template class Graph<MyMesh, float,float,float>;
template class Graph<MyMesh, double,double,double>;
