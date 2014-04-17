#ifndef GEOTREE_HH
#define GEOTREE_HH

#include "MyMesh.hh"
#include <QVector3D>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <QVector>
#include "PatchT.hh"
#include <queue>
#include "FacePriorityT.hh"

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
    //std::vector<std::priority_queue<FacePriorityT<M>*, std::vector<FacePriorityT<M>* >, FacePriorityCompareT<M> >* > priorityQueues;
    std::priority_queue<FacePriorityT<M>*, std::vector<FacePriorityT<M>* >, FacePriorityCompareT<M> > priorityQueues[100];
    void initPatches();
    void createPatches();
    void updatePatchColours();
    void updatePriorityQueues();
    void assignFaces();
    void updateCentroids();
    void clearPatches();
    PointMatrix getCentroids();
    double convergenceTest(PointMatrix old);
    Vec getFaceCentroid(typename M::FaceHandle fh);
};


#if defined(OM_INCLUDE_TEMPLATES) && !defined(GEOTREE_CC)
#  define SCENE_TEMPLATES
#  include "GeoTreeT.cc"
#endif

#endif
