#ifndef GEOTREE_HH
#define GEOTREE_HH

#include "MyMesh.hh"
#include <QVector3D>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <QVector>
#include "PatchT.hh"
#include <queue>

template <typename M>
class GeoTreeT
{
public:
    typedef M MyMesh;
    M *mesh;
    GeoTreeT(M *m, int k);
    ~GeoTreeT();

private:
    int noPatches;
    std::vector<PatchT<M>*> patches;
    std::vector<std::priority_queue<float, std::vector<int> > > priorityQueues;
    void initPatches();
    Vec getFaceCentroid(typename M::FaceHandle fh);
    void updatePatchColours();
    void createPriorityQueues();

};


#if defined(OM_INCLUDE_TEMPLATES) && !defined(GEOTREE_CC)
#  define SCENE_TEMPLATES
#  include "GeoTreeT.cc"
#endif

#endif
