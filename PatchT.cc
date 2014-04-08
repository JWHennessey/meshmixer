#ifndef Patch_CC
#define Patch_CC

#include "PatchT.hh"

template <typename M>
PatchT<M>::PatchT(Vec c, Vec n, int handle)
{
  centroid = c;
  normal = n;
  patchColour = QColor(rand() % 255, rand() % 255, rand() % 255);
  //addFaceHandle(handle);
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
  centroid = Vec(0,0,0);
  for(int i=0; i < verticies.size(); i++)
    centroid += verticies[i];

  centroid = centroid / verticies.size();
}

template <typename M>
void
PatchT<M>::addFaceHandle(int handle)
{
  faceHandles.push_back(handle);
}

template <typename M>
void 
PatchT<M>::clear()
{
  verticies.clear();
  faceHandles.clear();
}

#endif
