#ifndef GEOTREE_CC
#define GEOTREE_CC

#include "GeoTreeT.hh"
#include <stdlib.h>

template <typename M>
GeoTreeT<M>::GeoTreeT(M *m, int k)
{
  mesh = m;
  noPatches = k;
  std::cout << "GeoTree Init" << "\n";
  QColor modelColor(255, 0, 0);
  mesh->request_face_colors();
  mesh->request_vertex_colors();
  initPatches();
  updatePatchColours();
}

template <typename M>
void
GeoTreeT<M>::initPatches()
{
  std::vector<int> initSeeds;
  for(int i = 0; i<noPatches; i++)
  {
    bool notUnique = true;
    while(notUnique)
    {
      int seed = rand() % mesh->n_faces();
      if(std::find(initSeeds.begin(), initSeeds.end(), seed) == initSeeds.end()) {
        typename M::FaceHandle face = mesh->face_handle(seed);
        Vec centroid = getFaceCentroid(face);
        Vec normal = Vec(mesh->normal(face)[0], mesh->normal(face)[1], mesh->normal(face)[2]);
        patches.push_back(new PatchT<M>(centroid, normal, seed));
        for (typename M::FaceVertexIter vf_it=mesh->fv_iter(face); vf_it; ++vf_it)
        {
          patches.back()->addVertex(Vec(mesh->point(*vf_it)[0], mesh->point(*vf_it)[1], mesh->point(*vf_it)[2]));
        }
        notUnique = false;
      }
    }
  }
}

template <typename M>
Vec
GeoTreeT<M>::getFaceCentroid(typename M::FaceHandle fh)
{
  Vec centroid = Vec(0,0,0);
  for (typename M::FaceVertexIter vf_it=mesh->fv_iter(fh); vf_it; ++vf_it)
  {
    centroid = centroid + Vec(mesh->point(*vf_it)[0], mesh->point(*vf_it)[1], mesh->point(*vf_it)[2]);
  }
  return centroid / 3;
}

template <typename M>
void
GeoTreeT<M>::updatePatchColours()
{
  for(int i = 0; i<noPatches; i++)
  {
    for(int j = 0; j < patches[i]->faceHandles.size(); j++)
    {
      typename M::FaceHandle face = mesh->face_handle(patches[i]->faceHandles[j]);
      QColor modelColor = patches[i]->patchColour;
      mesh->set_color(face,  OpenMesh::Vec3f(modelColor.redF(), modelColor.blueF(), modelColor.greenF()));
    }
  }
}


template <typename M>
void
GeoTreeT<M>::createPriorityQueues()
{
  for(int i = 0; i<noPatches; i++)
  {
    priorityQueues.push_back(std::priority_queue<float, std::vector<int> > ());
    for (int f_it = 0; f_it < mesh->n_faces(); f_it++)
    {
      typename M::FaceHandle face = mesh->face_handle(f_it);
      Vec centroid = getFaceCentroid(face);
      Vec normal = Vec(mesh->normal(face)[0], mesh->normal(face)[1], mesh->normal(face)[2]);
      
    }
  }
}

#endif
