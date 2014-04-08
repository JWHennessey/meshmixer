#ifndef PATCH_HH
#define PATCH_HH

#include "MyMesh.hh"
#include <eigen3/Eigen/Dense>

template <typename M>
class PatchT
{
public:
  typedef M MyMesh;
  PatchT(Vec c, Vec n, int handle);
  void updateCentroid();
  void updateNormal();
  void addVertex(Vec v);
  void addFaceHandle(int handle);
  void clear();
  QColor patchColour;
  std::vector<int> faceHandles;
  Vec centroid;
  Vec normal;
  std::vector<Vec> verticies;
};


#if defined(OM_INCLUDE_TEMPLATES) && !defined(PATCH_CC)
#  define SCENE_TEMPLATES
#  include "PatchT.cc"
#endif


#endif
