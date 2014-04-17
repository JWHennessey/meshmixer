#ifndef MODEL_CC
#define MODEL_CC


#include "QtModelT.hh"
#include <math.h> 
#include <QFile>
#include <QTextStream>
#include <QVarLengthArray>
#include <QtOpenGL>
#include <random>
#include <cmath>
#include <math.h> 
#include <stdlib.h>


template <typename M>
QtModelT<M>::QtModelT(M& m)
  : modelColor(0, 0, 0)
  , vertical(0.0f)
  , horizontal(0.0f)
  , depth(0.0f)
  , deg2Rad(0.0174532925)
  , zAxis(0.0f)
{
  mesh = m;
  double min_x, max_x, min_y, max_y, min_z, max_z;
  bool first = true;
  /*
  OpenMesh::VPropHandleT< double > gauss;
  if(!mesh.get_property_handle(gauss, "Gaussian Curvature"))
    mesh.add_property(gauss, "Gaussian Curvature" );
  */
  // bounding box
  
  Vec3f bbMin, bbMax;
  bbMin = bbMax = OpenMesh::vector_cast<Vec3f>(mesh.point(*mesh.vertices_begin()));
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
  {
    if(first){
      min_x = mesh.point(*v_it)[0];
      max_x = mesh.point(*v_it)[0];
      min_y = mesh.point(*v_it)[1];
      max_y = mesh.point(*v_it)[1];
      min_z = mesh.point(*v_it)[2];
      max_z = mesh.point(*v_it)[2];
      first = false;
    }
    
    if(mesh.point(*v_it)[0] < min_x )
      min_x = mesh.point(*v_it)[0];
    else if(mesh.point(*v_it)[0] > max_x )
      max_x = mesh.point(*v_it)[0];

    if(mesh.point(*v_it)[1] < min_y )
      min_y = mesh.point(*v_it)[1];
    else if(mesh.point(*v_it)[0] > max_y )
      max_y = mesh.point(*v_it)[1];

    if(mesh.point(*v_it)[2] < min_z )
      min_z = mesh.point(*v_it)[2];
    else if(mesh.point(*v_it)[2] > max_z )
      max_z = mesh.point(*v_it)[2];
    
    bbMin.minimize( OpenMesh::vector_cast<Vec3f>(mesh.point(*v_it)));
    bbMax.maximize( OpenMesh::vector_cast<Vec3f>(mesh.point(*v_it)));
    mesh.data(*v_it).set_gauss(gauss_curvature(*v_it));
    //mesh.property(gauss,*v_it) = gauss_curvature(*v_it);
  }

  double diff, min;
  double diffX = max_x - min_x;
  double diffY = max_y - min_x;
  double diffZ = max_z - min_x;

  if(diffX > diffY && diffX > diffZ)
  {
    diff = diffX;
    min = min_x;
  }
  else if(diffY > diffX && diffY > diffZ)
  {
    diff = diffY;
    min = min_y;
  }
  else
  {
    diff = diffZ;
    min = min_z;
  }


  
  // set center and radius
  center = (bbMin+bbMax)*0.5;
  //horizontal = -center[0];
  //vertical = -center[1];
  //zAxis = -center[2];
  
  applyTransformations();
  //findBoundaryVertices();
  /*
  typedef typename M::Point Point;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    mesh.set_point( *v_it, Point(
          2.0*(mesh.point(*v_it)[0]-min)/(diff) - 1.0,
          2.0*(mesh.point(*v_it)[1]-min)/(diff) - 1.0,
          2.0*(mesh.point(*v_it)[2]-min)/(diff) - 1.0)
    );
  }
   */
  updateColour();
  calcNormals();
}


template <typename M>
std::vector<VertexHandle>
QtModelT<M>::findBoundaryRing(VertexHandle point){
  std::vector<VertexHandle> ring;
  ring.push_back(point);
  for (typename M::VertexVertexIter vv_it=mesh.vv_iter(point); vv_it; ++vv_it)
  {
    if (mesh.is_boundary(*vv_it)) {
      ring.push_back(*vv_it);
    }
  }
  size_t i = 0;
  while (i != ring.size()){
  i = ring.size();
  for (typename M::VertexVertexIter vv_it=mesh.vv_iter(ring[i-1]); vv_it; ++vv_it)
  {
    if (!std::find(ring.begin(), ring.end(), *vv_it)){
      if (mesh.is_boundary(*vv_it)) {
        ring.push_back(*vv_it);
      }
    }
  }
  }
  return ring;
}

template <typename M>
void
QtModelT<M>::findBoundaryVertices(){
  boundaryPoints.clear();
  boundaryPoints.reserve(mesh.n_vertices()/2);
  int i = 0;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
  {
    if (mesh.is_boundary(*v_it)) {
      boundaryPoints.push_back(*v_it);
      colourFaceFromVertexIndex(i);
    }
    i++;
  }
  boundaryMatrix.resize(boundaryPoints.size(), 3);
  for (int count = 0; count < boundaryPoints.size(); ++count)
  {
    VertexHandle v_it = boundaryPoints[count];
    boundaryMatrix(count, 0) = mesh.point(v_it)[0];
    boundaryMatrix(count, 1) = mesh.point(v_it)[1];
    boundaryMatrix(count, 2) = mesh.point(v_it)[2];
    count ++;
  }
}

template <typename M>
QtModelT<M>::~QtModelT()
{

}

template <typename M>
void
QtModelT<M>::select(int faceNumber){
  if (faceNumber > 2){
  typename M::FaceHandle face = mesh.face_handle(faceNumber-3);
  mesh.set_color(face, typename M::Color(0, 255, 255));
  }
}

template <typename M>
void
QtModelT<M>::colourFaceFromVertexIndex(int vertexNumber){
  typename M::VertexHandle point = mesh.vertex_handle(vertexNumber);
  for (typename M::VertexFaceIter vf_it=mesh.vf_begin(point); vf_it!=mesh.vf_end(point); ++vf_it)
  {
    mesh.set_color(*vf_it, typename M::Color(255, 0, 255));
  }
}

template <typename M>
void
QtModelT<M>::render()
{
    typename M::ConstFaceIter    fIt(mesh.faces_begin()),
                                 fEnd(mesh.faces_end());

    typename M::ConstFaceVertexIter fvIt;
    unsigned int index = 0;
    //std::cout << "Render" << "\n";
    //std::cout << horizontal << "\n";
    //std::cout << vertical << "\n";
    float matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    //glPopMatrix();
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(horizontal, vertical, zAxis);
    glMultMatrixf(matrix);
  
  //glEnable(GL_LIGHTING);
  //glShadeModel(GL_FLAT);
  
  glRotatef(meshRotation.x(), 1, 0, 0);
  glRotatef(meshRotation.y(), 0, 1, 0);
  glRotatef(meshRotation.z(), 0, 0, 1);
  
  glEnable(GL_DEPTH_TEST);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, mesh.points());
  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_FLOAT, 0, mesh.vertex_normals());
  glEnable(GL_NORMALIZE);
  glPushMatrix();
  glScalef(1, 1, 1);
  glLineWidth(10);
  glLoadName(index);
  glBegin(GL_LINES);
  glColor3f(255,0,0);
  glVertex3f(0, 0, 0);
  glVertex3f(1.5, 0, 0);
  glEnd();
  index++;
  glLoadName(index);
  glBegin(GL_LINES);
  glColor3f(0,255,0);
  glVertex3f(0, 0, 0);
  glVertex3f(0, 1.5, 0);
  glEnd();
  index++;
  glLoadName(index);
  glBegin(GL_LINES);
  glColor3f(0,0,255);
  glVertex3f(0, 0, 0);
  glVertex3f(0, 0, 1.5);
  index++;
  glEnd();
  for (; fIt!=fEnd; ++fIt)
  {
    glLoadName(index);
    //glPushName(index);
    glBegin(GL_TRIANGLES);
    glColor3fv(&mesh.color(*fIt)[0]);
    fvIt = mesh.cfv_iter(*fIt);
    glArrayElement(fvIt->idx());
    ++fvIt;
    glArrayElement(fvIt->idx());
    ++fvIt;
    glArrayElement(fvIt->idx());
    glEnd();
    index++;
  }
  glPopMatrix();
  glDisable(GL_NORMALIZE);

  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glPopMatrix();
  //glLoadMatrixf(matrix);
  //glPushMatrix();

}

template <typename M>
void
QtModelT<M>::applyTransformations()
{
  typedef typename M::Point Point;
  meshRotation = meshRotation * deg2Rad;
  Eigen::AngleAxis<float> aax(meshRotation.x(), Eigen::Vector3f(1, 0, 0));
  Eigen::AngleAxis<float> aay(meshRotation.y(), Eigen::Vector3f(0, 1, 0));
  Eigen::AngleAxis<float> aaz(meshRotation.z(), Eigen::Vector3f(0, 0, 1));
  Eigen::Quaternion<float> rotation = aax * aay * aaz;

  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    Eigen::Vector3f p = Eigen::Vector3f(mesh.point(*v_it)[0], mesh.point(*v_it)[1], mesh.point(*v_it)[2]);
    p = rotation * p;
    mesh.set_point( *v_it, Point(p[0], p[1], p[2]) );
    mesh.set_point( *v_it, mesh.point(*v_it) + Point(horizontal, vertical, depth) );
  }
  findBoundaryVertices();
  horizontal = 0.0f;
  vertical = 0.0f;
  depth = 0.0f;
  meshRotation.setX(0.0f);
  meshRotation.setY(0.0f);
  meshRotation.setZ(0.0f);
}

template <typename M>
PointMatrix//flann::Matrix<float>
QtModelT<M>::buildMatrix()
{
  PointMatrix m(mesh.n_vertices(), 3);
  int count = 0;

  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    m(count, 0) = mesh.point(*v_it)[0];
    m(count, 1) = mesh.point(*v_it)[1];
    m(count, 2) = mesh.point(*v_it)[2];
    count += 1;
  }
  return m;
}


template<typename M>
double
QtModelT<M>::gauss_curvature(VertexHandle _vh) {
  //if (mesh.status(_vh).deleted()) return 0.0;
  
  double gauss_curv = 2.0 * M_PI;
  
  /*
   
   TODO: Check the boundary case.
   
   If the vertex is a boundary vertex
   if ( _mesh.is_boundary(_vh) )
   gauss_curv = M_PI;
   
   */
  
  const Point p0 = mesh.point(_vh);
  
  typename M::CVOHIter voh_it(mesh.cvoh_iter(_vh));
  typename M::CVOHIter n_voh_it = voh_it;
  
  if ( ! voh_it->is_valid() )
    return 0.0;
  
  // move to next
  ++n_voh_it;
  
  for(; voh_it.is_valid(); ++voh_it, ++n_voh_it)
  {
    Point p1 = mesh.point(mesh.to_vertex_handle(   *voh_it));
    Point p2 = mesh.point(mesh.to_vertex_handle( *n_voh_it));
    
    gauss_curv -= acos(OpenMesh::sane_aarg( ((p1-p0).normalize() | (p2-p0).normalize()) ));
  }
  
  return gauss_curv;
}

template <typename M>
PointMatrix
QtModelT<M>::buildSampledMatrix()
{
  int noSamples = 5000;
  PointMatrix allMat = buildMatrix();
  PointMatrix randMat(noSamples, 3);
  for (int i = 0; i < noSamples; ++i )
  { 
    float ind = float(rand()) / RAND_MAX;

    randMat.row(i) = allMat.row(floor(ind * mesh.n_vertices() ) );
  }

  return randMat;

}

template <typename M>
PointMatrix//flann::Matrix<float>
QtModelT<M>::buildNormalMatrix()
{
  PointMatrix m(mesh.n_vertices(), 3);
  int count = 0;

  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    m(count, 0) = mesh.normal(*v_it)[0];
    m(count, 1) = mesh.normal(*v_it)[1];
    m(count, 2) = mesh.normal(*v_it)[2];
    count += 1;
  }
  return m;
}


template <typename M>
void
QtModelT<M>::updateTransformations(Matrix<double, 3, 3>& R, double x, double y, double z)
{
  typedef typename M::Point Point;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    Eigen::Vector3d p = Eigen::Vector3d(mesh.point(*v_it)[0], mesh.point(*v_it)[1], mesh.point(*v_it)[2]);
    p = R * p;
    mesh.set_point( *v_it, Point(p[0], p[1], p[2]) );
    mesh.set_point( *v_it, mesh.point(*v_it) - Point(x, y, z) );
  }
  render();
}


template <typename M>
int
QtModelT<M>::getNoVerticies()
{
  return mesh.n_vertices();
}

template <typename M>
void
QtModelT<M>::updateRotation(QVector3D& rotationVec)
{
  /*
  double x = 0.0, y = 0.0, z = 0.0;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    x += mesh.point(*v_it)[0];
    y += mesh.point(*v_it)[1];
    z += mesh.point(*v_it)[2];
  }
  x = x / mesh.n_vertices();
  y = y / mesh.n_vertices();
  z = z / mesh.n_vertices();

  Eigen::Vector3f t = Eigen::Vector3f(x, y, z);

  typedef typename M::Point Point;
  QVector3D modelRotation = rotationVec;
  modelRotation = modelRotation * 0.0174532925;
  Eigen::AngleAxis<float> aax(modelRotation.x(), Eigen::Vector3f(1, 0, 0));
  Eigen::AngleAxis<float> aay(modelRotation.y(), Eigen::Vector3f(0, 1, 0));
  Eigen::AngleAxis<float> aaz(modelRotation.z(), Eigen::Vector3f(0, 0, 1));
  Eigen::Quaternion<float> rotation = aax * aay * aaz;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    Eigen::Vector3f p = Eigen::Vector3f(mesh.point(*v_it)[0], mesh.point(*v_it)[1], mesh.point(*v_it)[2]);
    p = p - t;
    p = rotation * p;
    p = p + t;
    mesh.set_point( *v_it, Point(p[0], p[1], p[2]) );
  }
   */
  meshRotation += rotationVec;
}

template <typename M>
void
QtModelT<M>::updateHorizontal(float x)
{
  horizontal += x;
}

template <typename M>
void
QtModelT<M>::updateVertical(float x)
{
  vertical -= x;
}

template <typename M>
void
QtModelT<M>::updateZAxis(float x)
{
  zAxis += x;
}

template <typename M>
void
QtModelT<M>::updateColour()
{
  mesh.request_face_colors();
  mesh.request_vertex_colors();
  for (typename M::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it)
  {
    mesh.set_color(*f_it, OpenMesh::Vec3f(modelColor.redF(), modelColor.blueF(), modelColor.greenF()));
  }
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
  {
    mesh.set_color(*v_it, OpenMesh::Vec3f(modelColor.redF(), modelColor.blueF(), modelColor.greenF()));
  }
}

template <typename M>
void
QtModelT<M>::clearColour()
{
  mesh.request_vertex_colors();
  modelColor.setRgb(10, 10, 10);
  updateColour();
}

template <typename M>
void
QtModelT<M>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
  /*
  std::cout << "paint" << "\n";
  painter->drawRect(boundingRect());
  painter->beginNativePainting();
  glPushMatrix();
  glTranslatef(horizontal, vertical, 0);
  glRotatef(modelRotation.x(), 1, 0, 0);
  glRotatef(modelRotation.y(), 0, 1, 0);
  glRotatef(modelRotation.z(), 0, 0, 1);

  glEnable(GL_DEPTH_TEST);
  glEnableClientState(GL_VERTEX_ARRAY);

  if ( mesh.has_vertex_colors() )
  {
    glEnableClientState( GL_COLOR_ARRAY );
    glColorPointer(3, GL_FLOAT, 0, mesh.vertex_colors());
  }
  //glColor4f(modelColor.redF(), modelColor.greenF(), modelColor.blueF(), 1.0f);
  glVertexPointer(3, GL_FLOAT, 0, mesh.points());
  glDrawArrays( GL_POINTS, 0, static_cast<GLsizei>(mesh.n_vertices()) );

  glDisableClientState(GL_VERTEX_ARRAY);
  glDisable(GL_DEPTH_TEST);

  glDisableClientState(GL_COLOR_ARRAY);  

  glPopMatrix();

  painter->endNativePainting();
*/
}


template <typename M>
QRectF 
QtModelT<M>::boundingRect() const
{
  return QRectF(0,0, 1024, 768);
}

template <typename M>
void
QtModelT<M>::calcNormals()
{

  //std::cout << "calcNormals()" << "\n";
  OpenMesh::IO::Options opt;
  mesh.request_vertex_normals();
  // Add face normals as default property
  mesh.request_face_normals();

  // If the file did not provide vertex normals, then calculate them
  if ( !opt.check( OpenMesh::IO::Options::VertexNormal ) &&
       mesh.has_face_normals() && mesh.has_vertex_normals() )
  {
    // let the mesh update the normals
    mesh.update_face_normals();
    mesh.update_vertex_normals();
  }

}

template <typename M>
void
QtModelT<M>::nearestNeighbours(double radius, MapTable* resultTable)
{
  //build kd tree
  PointMatrix pAll = buildMatrix();

  std::cout << pAll << "\n";

  typedef nanoflann::KDTreeEigenMatrixAdaptor<PointMatrix>  kd_tree_t;
  kd_tree_t mat_index(3, pAll, 10);
  mat_index.index->buildIndex();

  //find neighbourhood for each point
  int i = 0;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
  {
    typename M::Point p = mesh.point(*v_it);
    std::vector<double> query_pt(3);
    query_pt[0] = p[0];
    query_pt[1] = p[1];
    query_pt[2] = p[2];

    std::cout << p << "\n";

    std::vector< std::pair< size_t, double > > resultPairs;
    resultPairs.reserve(mesh.n_vertices());

    size_t count = mat_index.index->radiusSearch(&query_pt[0], radius, resultPairs, nanoflann::SearchParams(true));
    std::cout << resultPairs.size() << "\n";
    resultTable->push_back(resultPairs);
    i++;
  }
  std::cout << resultTable->size() << "\n";
}

template<typename M>
void
QtModelT<M>::scale(float alpha)
{
  typedef typename M::Point Point;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    mesh.set_point( *v_it, Point(alpha*mesh.point(*v_it)));
  }

}

#endif
