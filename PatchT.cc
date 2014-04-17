#ifndef Patch_CC
#define Patch_CC

#include "PatchT.hh"

template <typename M>
PatchT<M>::PatchT(Vec c, Vec n, int handle)
{
  centroid = c;
  normal = n;
  patchColour = QColor(rand() % 255, rand() % 255, rand() % 255);
  addFaceHandle(handle);
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
  PointMatrix mat(verticies.size(), 3);
  std::cout << "No. Faces on patch " << faceHandles.size() << "\n";
  std::cout << "No. Vertices on patch " << verticies.size() << "\n";
  for(int i = 0; i < verticies.size(); i++)
  {
    mat(i, 0) = verticies[i][0];
    mat(i, 1) = verticies[i][1];
    mat(i, 2) = verticies[i][2];
  }
  Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
  Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(mat.rows());
  Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
  Eigen::VectorXcd v = es.eigenvectors().row(1);
  normal = Vec(std::real(v[0]), std::real(v[1]), std::real(v[2]));
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
