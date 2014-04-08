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
    priorityQueues.push_back(&queue);

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
    //std::priority_queue<FacePriorityT<M>*, std::vector<FacePriorityT<M>* >, FacePriorityCompareT<M> > queue;
    //priorityQueues.push_back(&queue);
    for (int f_it = 0; f_it < mesh->n_faces(); f_it++)
    {
      typename M::FaceHandle face = mesh->face_handle(f_it);
      Vec centroid = getFaceCentroid(face);
      Vec normal = Vec(mesh->normal(face)[0], mesh->normal(face)[1], mesh->normal(face)[2]);
      float dist = (centroid - patches[i]->centroid).norm();
      float normAngle = normal.dot(patches[i]->normal);
      //This might not be right but I think parallel vectors will have
      //a value of 1 and the others will be 1.
      float cost = dist * normAngle;
      FacePriorityT<M>* fp = new FacePriorityT<M>(cost, f_it);
      priorityQueues.back()->push(fp);
      //std::cout << "PQ Size" << priorityQueues[i]->size() << "\n";
      //std::cout << "PQ Top Cost" << priorityQueues[i]->top()->cost << "\n";
      //std::cout << priorityQueues.back()->size() << "\n";
    }
  }
}

template <typename M>
void
GeoTreeT<M>::createPatches()
{
  for(int i=0; i<2; i++)
  {
    clearPatches();
    updatePriorityQueues();
    std::cout << "Priority Queues Done" << "\n";
    assignFaces();
    std::cout << "Assign Faces Done" << "\n";
    updateCentroids();
    std::cout << "Update Centroids Done" << "\n";
  }
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
  //for(int i = 0; i<noPatches; i++)
  //{
    //std::cout << i << "\n";
    //std::cout << priorityQueues[i]->size() << "\n";
    //while(!priorityQueues[i]->empty())
    //{
      //FacePriorityT<M>* fp = priorityQueues[i]->top();
      //std::cout << fp->faceId << "  " << fp->cost << " \n";
      //priorityQueues[i]->pop();
    //}
  //}
  int counter = 0;
  std::vector<int> assignedFaces;
  bool allFacesAssigned = false;
  while(!allFacesAssigned)
  {
    for(int i = 0; i<noPatches; i++)
    {
      bool itemTaken = false;
      while(!itemTaken)
      {
        if(priorityQueues[i]->empty())
        {
          //std::cout << "Queue Empty" << "\n";
          //std::cout << "No. Faces " << mesh->n_faces() << "\n";
          //std::cout << "Assigned Faces " << assignedFaces.size() << "\n";
          ////for(int j = 0; j<assignedFaces.size(); j++)
          ////{
            ////std::cout << assignedFaces[j] << "\n";
          ////}
          itemTaken = true;
        }
        else
        {
          FacePriorityT<M>* fp = priorityQueues[i]->top();
          //std::cout << "Face Selected" << "\n";
          //if face isn't already assigned
          if(std::find(assignedFaces.begin(), assignedFaces.end(), fp->faceId) != assignedFaces.end()) 
          {
            //std::cout << "Face already assigned" << "\n";
          }
          else
          {
            //std::cout << "Face assigned" << "\n";
            patches[i]->addFaceHandle(fp->faceId);
            assignedFaces.push_back(fp->faceId);
            typename M::FaceHandle face = mesh->face_handle(fp->faceId);
            for (typename M::FaceVertexIter vf_it=mesh->fv_iter(face); vf_it; ++vf_it)
            {
              patches.back()->addVertex(Vec(mesh->point(*vf_it)[0], mesh->point(*vf_it)[1], mesh->point(*vf_it)[2]));
            }
            itemTaken = true;

          }
          priorityQueues[i]->pop();
        }
      }
    }
    //std::cout << mesh->n_faces() << "\n";
    //std::cout << assignedFaces.size() << "\n";
    //std::cout << mesh->n_faces() - assignedFaces.size() << "\n";
    //allFacesAssigned = true;
    allFacesAssigned = assignedFaces.size() >= mesh->n_faces();
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
