#ifndef MyMesh_HH
#define MyMesh_HH
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <eigen3/Eigen/Dense>

struct MyTraits : public OpenMesh::DefaultTraits
{
  typedef OpenMesh::Vec3f Color;
  VertexAttributes( OpenMesh::Attributes::Normal |
                    OpenMesh::Attributes::Color  |
                    OpenMesh::Attributes::Status   );

  FaceAttributes( OpenMesh::Attributes::Normal  |
                  OpenMesh::Attributes::Color |  
                  OpenMesh::Attributes::Status );

  EdgeAttributes(OpenMesh::Attributes::Status);

  VertexTraits
  {
  private:
    float gauss_;
  public:
    
    VertexT() : gauss_(0.0f) { }
    
    const float& gauss() const { return gauss_; }
    void set_gauss(const float& _gauss) { gauss_ = _gauss; }
  };


  HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyMesh;
//typedef MyMesh::ConstFaceVertexIter FVI;
typedef MyMesh::Point Point;
typedef Eigen::MatrixX3d PointMatrix;
typedef Eigen::Vector3d Vec;

#endif
