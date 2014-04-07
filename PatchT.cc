#ifndef Patch_CC
#define Patch_CC

#include "PatchT.hh"

template <typename M>
PatchT<M>::PatchT(Vec c, Vec n, int handle)
{
  centroid = c;
  normal = n;
  faceHandles.push_back(handle);
  patchColour = QColor(rand() % 255, rand() % 255, rand() % 255);
}

template <typename M>
void 
PatchT<M>::addVertex(Vec v)
{
  verticies.push_back(v);
}


template <typename M>
void 
PatchT<M>::updateNormal()
{

}

template <typename M>
void 
PatchT<M>::updateCentroid()
{

}

#endif
