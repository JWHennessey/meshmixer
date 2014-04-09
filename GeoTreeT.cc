#ifndef GEOTREE_CC
#define GEOTREE_CC

#include "GeoTreeT.hh"
#include <stdlib.h>
#include "FacePriorityT.hh"





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
  for(int i = 0; i<noPatches; i++)
  {
    std::priority_queue<FacePriorityT<M>*, std::vector<FacePriorityT<M>* >, FacePriorityCompareT<M> > queue;
    priorityQueues[i] = queue;

  }
  std::cout << "Seeds Chosen" << "\n";
  createPatches();
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
        notUnique = false;
      }
    }
    //std::cout << "No. Faces on patch " << i << " - " << patches[i]->faceHandles.size() << "\n";
    //std::cout << "No. Vertices on patch " << i << " - " << patches[i]->verticies.size() << "\n";
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
    //std::cout << i << "\n";
    for(int j = 0; j < patches[i]->faceHandles.size(); j++)
    {
      //std::cout << j << "\n";
      typename M::FaceHandle face = mesh->face_handle(patches[i]->faceHandles[j]);
      QColor modelColor = patches[i]->patchColour;
      mesh->set_color(face,  OpenMesh::Vec3f(modelColor.redF(), modelColor.blueF(), modelColor.greenF()));
    }
  }
}

template <typename M>
void
GeoTreeT<M>::updatePriorityQueues()
{

  for(int i = 0; i<noPatches; i++)
  {
    for (int f_it = 0; f_it < mesh->n_faces(); f_it++)
    {
      typename M::FaceHandle face = mesh->face_handle(f_it);
      Vec centroid = getFaceCentroid(face);
      Vec normal = Vec(mesh->normal(face)[0], mesh->normal(face)[1], mesh->normal(face)[2]);
      float dist = (centroid - patches[i]->centroid).norm();
      float dot = normal.dot(patches[i]->normal);
      normal.normalize();
      patches[i]->normal.normalize();
      float normDist = (normal - patches[i]->normal).norm() / 10.0f;
      std::cout << "Normal Dist " << normDist << "\n";
      std::cout << "Eucld Dist " << dist << "\n";

      //float normAngle = dot / prodNorms;
      //normAngle = sqrt(normAngle * normAngle);
      //This might not be right but I think parallel vectors will have
      //a value of 1 and the others will be 1.
      //float paraDist;
      //std::cout << normAngle << "\n";
      float cost = (dist); //+ normDist; //+ (1 - normAngle);// + (normAngle);
      FacePriorityT<M>* fp = new FacePriorityT<M>(cost, f_it);
      priorityQueues[i].push(fp);
    }
    std::cout << "PriorityQueue " << i << "  size after created " << priorityQueues[i].size() << "\n";
  }
}

template <typename M>
void
GeoTreeT<M>::createPatches()
{
  std::cout << "No Verticies on Mesh " << mesh->n_vertices() << "\n";
  std::cout << "No faces on Mesh " << mesh->n_faces() << "\n";
  bool notConverged = true;
  double diff = 100.0;
  while(notConverged)
  {
    PointMatrix oldCentroids = getCentroids();
    clearPatches();
    updatePriorityQueues();
    std::cout << "Priority Queues Done" << "\n";
    assignFaces();
    std::cout << "Assign Faces Done" << "\n";
    updateCentroids();
    std::cout << "Update Centroids Done" << "\n";
    double newDiff = convergenceTest(oldCentroids);
    std::cout << "Difference " << newDiff << "\n";
    if(newDiff > diff || newDiff < 0.001)
      notConverged = false;
    else
      diff = newDiff;
  }
}

template <typename M>
PointMatrix
GeoTreeT<M>::getCentroids()
{
  PointMatrix centroids(noPatches, 3);
  for(int i = 0; i<noPatches; i++)
  {
    centroids.row(i) = patches[i]->centroid;
  }
  return centroids;
}

template <typename M>
double
GeoTreeT<M>::convergenceTest(PointMatrix old)
{
  PointMatrix n = getCentroids();
  return (n - old).norm();
}

template <typename M>
void
GeoTreeT<M>::clearPatches()
{
  for(int i = 0; i<noPatches; i++)
  {
    patches[i]->clear();
  }
}

template <typename M>
void
GeoTreeT<M>::assignFaces()
{
  std::vector<int> assignedFaces;
  bool allFacesAssigned = false;
  while(!allFacesAssigned)
  {
    for(int i = 0; i<noPatches; i++)
    {
      bool itemTaken = false;
      while(!itemTaken)
      {
        if(priorityQueues[i].empty())
        {
          itemTaken = true;
        }
        else
        {
          FacePriorityT<M>* fp = priorityQueues[i].top();
          // if face isn't already assigned
          if(std::find(assignedFaces.begin(), assignedFaces.end(), fp->faceId) == assignedFaces.end()) 
          {
            // std::cout << "Face assigned" << "\n";
            patches[i]->addFaceHandle(fp->faceId);
            assignedFaces.push_back(fp->faceId);
            typename M::FaceHandle face = mesh->face_handle(fp->faceId);
            for (typename M::FaceVertexIter vf_it=mesh->fv_iter(face); vf_it; ++vf_it)
            {
              patches[i]->addVertex(Vec(mesh->point(*vf_it)[0], mesh->point(*vf_it)[1], mesh->point(*vf_it)[2]));
            }
            itemTaken = true;
          }
          priorityQueues[i].pop();
        }
      }
    }
    allFacesAssigned = assignedFaces.size() >= mesh->n_faces();
  }
  for(int i = 0; i<noPatches; i++)
  {
    std::cout << "PriorityQueue " << i << " size after faces assigned " << priorityQueues[i].size() << "\n";
    std::cout << "Patch " << i << " no. face handles " << patches[i]->faceHandles.size() << "\n";
    std::cout << "Patch " << i << " no. verticies " << patches[i]->verticies.size() << "\n";
  }
}

template <typename M>
void
GeoTreeT<M>::updateCentroids()
{
  for(int i = 0; i<noPatches; i++)
  {
    patches[i]->updateCentroid();
    patches[i]->updateNormal();
  }
}

#endif
