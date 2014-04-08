#ifndef FACEPRIORITY_HH
#define FACEPRIORITY_HH

#include "MyMesh.hh"


template <typename M>
class FacePriorityT
{
public:
  typedef M MyMesh;
  float cost;
  int faceId;
  FacePriorityT(float c, int f){ cost = c; faceId = f; }
};

template <typename M>
class FacePriorityCompareT
{
public:
  typedef M MyMesh;
  bool operator() (FacePriorityT<M>* fp1, FacePriorityT<M>* fp2)
  {
    return (fp2->cost < fp1->cost);
  }
};

#endif
