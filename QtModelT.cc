#ifndef MODEL_CC
#define MODEL_CC

#include <stdio.h>
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
#include <stdio.h>
#include "maxflow/graph.h"
//#include <omp.h>

template <typename M>
QtModelT<M>::QtModelT(M& m)
  : modelColor(0, 0, 0)
  , vertical(0.0f)
  , horizontal(0.0f)
  , depth(0.0f)
  , deg2Rad(0.0174532925)
  , zAxis(0.0f)
  , dest(-1)
  , showFuzzy(false)
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

    if (mesh.is_boundary(*v_it)) boundaryPoints.push_back(*v_it);


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
      colourFaceFromVertexIndex(i,Point(255,0,255));
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
QtModelT<M>::addToStroke(int f){
  stroke.push_back(f);
  fuzzyRegion.insert(f);
  typename M::FaceHandle face = mesh.face_handle(f);

  Vec n = Vec(mesh.normal(face)[0], mesh.normal(face)[1], mesh.normal(face)[2]);
  if(stroke.size() == 1){
    firstStrokeNorm = n;
  }else{
    lastStrokeNorm = n;
  }
  mesh.set_color(face, typename M::Color(0, 255, 255));
}

//template <typename M>
//void
//QtModelT<M>::select(int faceNumber){
  //if (faceNumber > 2){
    //typename M::FaceHandle face = mesh.face_handle(faceNumber-3);
    //mesh.set_color(face, typename M::Color(0, 255, 255));
  //}
//}

template <typename M>
void
QtModelT<M>::colourFaceFromVertexIndex(int vertexNumber, Point col){
  typename M::VertexHandle point = mesh.vertex_handle(vertexNumber);
  for (typename M::VertexFaceIter vf_it=mesh.vf_begin(point); vf_it!=mesh.vf_end(point); ++vf_it)
  {
    mesh.set_color(*vf_it, typename M::Color(col[0],col[1],col[2]));
  }
}

template <typename M>
void
QtModelT<M>::addToFuzzyRegion(int f){

  if(std::find(stroke.begin(), stroke.end(), f) == stroke.end()) {
    fuzzyRegion.insert(f);
    typename M::FaceHandle face = mesh.face_handle(f);
    //mesh.set_color(face, typename M::Color(0, 255, 0));
  }
}

template <typename M>
void
QtModelT<M>::select(int faceNumber){
    if(stroke.empty())
      addToStroke(faceNumber);
    else if(facesConnected(stroke.back(), faceNumber))
    { 
      //std::cout << "Connected" << "\n";
      addToStroke(faceNumber);
      typename M::FaceHandle fh1 = mesh.face_handle(faceNumber);
      for (typename M::FaceFaceIter ff_it1=mesh.ff_iter(fh1); ff_it1; ++ff_it1)
      {
        for (typename M::FaceFaceIter ff_it2=mesh.ff_iter(ff_it1.handle()); ff_it2; ++ff_it2)
        {
          //for (typename M::FaceFaceIter ff_it3=mesh.ff_iter(ff_it2.handle()); ff_it3; ++ff_it3)
          //{
            //addToFuzzyRegion(ff_it3.handle().idx());
          //}
          addToFuzzyRegion(ff_it2.handle().idx());
        }
        addToFuzzyRegion(ff_it1.handle().idx());
      }
    }
    else
    {
      //std::cout << "Not Connected" << "\n";
      stroke.clear();
      strokeVertices.clear();
      clearColour();
      addToStroke(faceNumber);
    }
}


template <typename M>
void
QtModelT<M>::addSink(int faceNumber){
  if(inRegion(faceNumber)) sourceRegion.erase(faceNumber);

  sinkRegion.insert(faceNumber);
  typename M::FaceHandle fh = mesh.face_handle(faceNumber);
  mesh.set_color(fh, typename M::Color(255, 0, 0));

}

template <typename M>
void
QtModelT<M>::addSource(int faceNumber){
  if(inRegion(faceNumber)) sinkRegion.erase(faceNumber);

  sourceRegion.insert(faceNumber);
  typename M::FaceHandle fh = mesh.face_handle(faceNumber);
  mesh.set_color(fh, typename M::Color(0, 0, 255));

}

//template <typename M>
//void
//QtModelT<M>::removeFromRegion(int faceNumber){

//}

template <typename M>
bool
QtModelT<M>::facesConnected(int f1, int f2){
  //typename M::FaceHandle fh1 = mesh.face_handle(f1);
  //typename M::FaceHandle fh2 = mesh.face_handle(f2);
  //bool ret = false;
  //for (typename M::FaceVertexIter vf_it1=mesh.fv_iter(fh1); vf_it1; ++vf_it1)
  //{
    //for (typename M::FaceVertexIter vf_it2=mesh.fv_iter(fh2); vf_it2; ++vf_it2)
    //{
      //strokeVertices.push_back(mesh.point(*vf_it1));
      ////std::cout << mesh.point(*vf_it1) << "\n";
      ////std::cout << mesh.point(*vf_it2) << "\n";
      //if(mesh.point(*vf_it1) == mesh.point(*vf_it2))
      //{
        //ret = true;
      //}
    //}
  //}
  //return ret;
  return true;
}

template <typename M>
std::vector<int>
QtModelT<M>::getStroke()
{
  return stroke;
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
    //glLoadIdentity();
    glTranslatef(horizontal, vertical, zAxis);
    //glMultMatrixf(matrix);
  
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
  index = 0;
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
  QVector3D localRotation = rotationVec;
  localRotation = localRotation * 0.0174532925;
  Eigen::AngleAxis<float> aax(localRotation.y(), Eigen::Vector3f(1, 0, 0));
  Eigen::AngleAxis<float> aay(localRotation.x(), Eigen::Vector3f(0, 1, 0));
  Eigen::AngleAxis<float> aaz(localRotation.z(), Eigen::Vector3f(0, 0, 1));
  Eigen::Quaternion<float> rotation = aax * aay * aaz;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    Eigen::Vector3f p = Eigen::Vector3f(mesh.point(*v_it)[0], mesh.point(*v_it)[1], mesh.point(*v_it)[2]);
    p = p - t;
    p = rotation * p;
    p = p + t;
    mesh.set_point( *v_it, Point(p[0], p[1], p[2]) );
  }

  //modelRotation += rotationVec;

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
  dest = -1;
  fuzzyRegion.clear();
  prev.clear();
  stroke.clear();
  strokeVertices.clear();
  sourceRegion.clear();
  sinkRegion.clear();
  showFuzzy = false; 
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

template<typename M>
void
QtModelT<M>::createGeoTree(int k)
{
  std::cout << "Create GeoTree with " << k << " patches" << "\n";
  geoTree = new GeoTreeT<M>(&mesh, k);
}

template<typename M>
M
QtModelT<M>::cut()
{
  std::cout << "Cut" << "\n";
  M m;
  for ( auto it = sinkRegion.begin(); it != sinkRegion.end(); ++it )
  {
    typename M::FaceHandle fh = mesh.face_handle(*it);
    std::vector<typename M::VertexHandle>  face_vhandles;
    for (typename M::FaceVertexIter vf_it=mesh.fv_iter(fh); vf_it; ++vf_it)
    {
      face_vhandles.push_back(m.add_vertex(mesh.point(*vf_it)));
    }
    m.add_face(face_vhandles);
  }
  deleteSink();
  typedef typename M::Point Point;
  QVector3D tempR = meshRotation * deg2Rad;
  Eigen::AngleAxis<float> aax(tempR.x(), Eigen::Vector3f(1, 0, 0));
  Eigen::AngleAxis<float> aay(tempR.y(), Eigen::Vector3f(0, 1, 0));
  Eigen::AngleAxis<float> aaz(tempR.z(), Eigen::Vector3f(0, 0, 1));
  Eigen::Quaternion<float> rotation = aax * aay * aaz;

  for (typename M::VertexIter v_it=m.vertices_begin(); v_it!=m.vertices_end(); ++v_it) 
  {
    Eigen::Vector3f p = Eigen::Vector3f(m.point(*v_it)[0], m.point(*v_it)[1], m.point(*v_it)[2]);
    p = rotation * p;
    m.set_point( *v_it, Point(p[0], p[1], p[2]) );
    m.set_point( *v_it, m.point(*v_it) + Point(horizontal, vertical, depth) );
  }
  
  return m;
}

template<typename M>
double
QtModelT<M>::cost(int u, int v)
{
  //std::cout << u << "\n";
  typename M::VertexHandle u_vh = mesh.vertex_handle(u);
  typename M::VertexHandle v_vh = mesh.vertex_handle(v);
  double dist = (mesh.point(v_vh) - mesh.point(u_vh)).norm();
  dist += inverseGeodesic(v);
  //std::cout << "Before Normal Dist" << "\n";
  dist += normalDistance(v);
  //std::cout << "After Normal Dist" << "\n ... \n";
  return dist;
}

template<typename M>
double
QtModelT<M>::normalDistance(int vertex)
{
 //return 0.0;
 typename M::VertexHandle vh = mesh.vertex_handle(vertex);
 Vec n = Vec(mesh.normal(vh)[0], mesh.normal(vh)[1], mesh.normal(vh)[2] );
 double cosOpeningAngle = cos(firstStrokeNorm.dot(lastStrokeNorm));
 double nDotv = n.dot(strokeNormal);
 if(cosOpeningAngle <= nDotv)
   return 1.0;
 else
   return (nDotv + 1.0 / cosOpeningAngle + 1.0);
}

//Perhaps this could be the Geodesic distance to the vertex closest to the centroid.
template<typename M>
double
QtModelT<M>::inverseGeodesic(int vertex)
{
  typename M::VertexHandle vh = mesh.vertex_handle(vertex);
  Vec p = Vec(mesh.point(vh)[0], mesh.point(vh)[1], mesh.point(vh)[2] );
  return 1.0 / (strokeCentroid - p).norm();
}


template<typename M>
void
QtModelT<M>::calcStrokeProxies()
{
  PointMatrix mat(stroke.size(), 3);
  int i = 0;
  for ( auto it = stroke.begin(); it != stroke.end(); ++it )
   {
    typename M::FaceHandle fh = mesh.face_handle(*it);
    Vec v = getFaceCentroid(fh);
    mat(i, 0) = v[0];
    mat(i, 1) = v[1];
    mat(i, 2) = v[2];
    i++;
  }
  Eigen::MatrixXd meanMat = mat.colwise().mean();
  strokeCentroid = Vec(meanMat(0, 0), meanMat(0, 1), meanMat(0, 2));
  Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
  Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(mat.rows());
  Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
  Eigen::VectorXcd v = es.eigenvectors().row(1);
  strokeNormal = Vec(std::real(v[0]), std::real(v[1]), std::real(v[2]));
}


template <typename M>
Vec
QtModelT<M>::getFaceCentroid(typename M::FaceHandle fh)
{
  Vec centroid = Vec(0,0,0);
  for (typename M::FaceVertexIter vf_it=mesh.fv_iter(fh); vf_it; ++vf_it)
  {
    centroid = centroid + Vec(mesh.point(*vf_it)[0], mesh.point(*vf_it)[1], mesh.point(*vf_it)[2]);
  }
  return centroid / 3;
}

template<typename M>
bool
QtModelT<M>::createSourceAndSink()
{
  sourceRegion.clear();
  sinkRegion.clear();
  std::cout << "Source and Sink" << "\n";
  bool notUnique = true;
  int assignedFaces = fuzzyRegion.size();
  bool completed = false;
  while(notUnique)
  {
    std::cout << "Not unique source" << "\n";
    std::cout << "No Faces" << mesh.n_faces() << "\n";
    //for(int f = 1; f < mesh.n_faces(); f++)
    //{
      //std::cout << f << "\n";
      //if(!inRegion(f))
      //{
        //regionGrow(f, &sourceRegion, 0);
        //notUnique = false;
      //}
    //}
    int f = rand() % mesh.n_faces();
    if(!inRegion(f))
    {
      regionGrow(f, &sourceRegion, 0);
      notUnique = false;
    }
  }
  notUnique = true;
  //while(notUnique)
  //{
    std::cout << "Not unique sink" << "\n";

    for(int f = 1; f < mesh.n_faces(); f++)
    {
      //std::cout << f << "\n";
      if(!inRegion(f))
      {
        regionGrow(f, &sinkRegion, 1);
        completed = true;
      }
    }
  //}
  //assignedFaces = fuzzyRegion.size() + sinkRegion.size() + sourceRegion.size();
  //std::cout << "In Fuzzy Region " << fuzzyRegion.size() << " / " << mesh.n_faces() << "\n";
  //std::cout << "In Source Region " << sourceRegion.size() << " / " << mesh.n_faces() << "\n";
  //std::cout << "In Sink Region " << sinkRegion.size() << " / " << mesh.n_faces() << "\n";
  return completed;
  //if(completed){
    //sinkRegion.clear();
    //sourceRegion.clear();
    //fuzzyRegion.clear();
    //stroke.clear();
    //std::cout << "Bad Path Selected" << "\n";
  //}
}

template<typename M>
void
QtModelT<M>::regionGrow(int f, std::unordered_set<int>* region, int type)
{
  //std::cout << f << " Region Grow " << type << "\n";
  region->insert(f);
  
  typename M::FaceHandle fh = mesh.face_handle(f);
  if(type == 1)
    mesh.set_color(fh, typename M::Color(255, 0, 0));
  else
    mesh.set_color(fh, typename M::Color(0, 0, 255));


  for (typename M::FaceFaceIter ff_it=mesh.ff_iter(fh); ff_it; ++ff_it)
  {
    f = ff_it.handle().idx();
    //std::cout << "FFI " << f << "\n";
    if(!inRegion(f))
    {
      //std::cout << "Not in region" << "\n";
      regionGrow(f, region, type);
    }
    else
    {
      //std::cout << "In region" << "\n";
    }
  }
}

template<typename M>
bool
QtModelT<M>::inRegion(int f)
{
 std::unordered_set<int>::const_iterator gotFR = fuzzyRegion.find (f);
 if(gotFR != fuzzyRegion.end())
 {
  //std::cout << "In Fuzzy Region " << fuzzyRegion.size() << " / " << mesh.n_faces() << "\n";
   return true;
 }
 std::unordered_set<int>::const_iterator gotSoR = sourceRegion.find (f);
 if(gotSoR != sourceRegion.end())
 {
  //std::cout << "In Source Region" << sourceRegion.size() << " / " << mesh.n_faces() << "\n";;
  return true;
 }
 std::unordered_set<int>::const_iterator gotSiR = sinkRegion.find (f);
 if(gotSiR != sinkRegion.end())
 {
  //std::cout << "In Sink Region" << sinkRegion.size() << " / " << mesh.n_faces() << "\n";;
  return true;
 }
 return false;
}

//template<typename M>
//void
//QtModelT<M>::sourceSinkWeightThread(Graph<M, double,double,double> *g, int *mapping[std::size_t], const int start, const int end)
//{
  //for(int i = start; i<end; i++)
  //{
    //std::cout << "t-weight" << "\n";
    //g->add_tweights(i, distToSource(mapping[i]), distToSink(mapping[i]) );
  //}

//}

template<typename M>
void
QtModelT<M>::graphCut()
{
  std::cout << "Graph Cut" << "\n";
  
  typedef Graph<M, double,double,double> GraphType;
  //estimated nunber of nodes and estimated number of edges
  GraphType *g = new GraphType(fuzzyRegion.size(), fuzzyRegion.size()*4);

  int mapping[fuzzyRegion.size()];
  int index = 0;
  for ( auto it = fuzzyRegion.begin(); it != fuzzyRegion.end(); ++it )
  {
    g->add_node(); 
    mapping[index++] = *it;
  }
  std::cout << "Mapping Created" << "\n";
  
  //std::thread first(&QtModelT<M>::sourceSinkWeightThread(g, &mapping, 0, fuzzyRegion.size()), this);


  //first.join();
  //#pragma omp parallel for
  for(int i = 0; i<fuzzyRegion.size(); i++)
  {
    std::cout << "t-weight" << "\n";
    g->add_tweights(i, distToSource(mapping[i]), distToSink(mapping[i]) );
  }
  std::cout << "t-weights done" << "\n";

  for(int i = 0; i<fuzzyRegion.size(); i++)
  {
    std::cout << fuzzyRegion.size() - i << "\n";
    typename M::FaceHandle fh = mesh.face_handle(mapping[i]);
    for (typename M::FaceFaceIter ff_it=mesh.ff_iter(fh); ff_it; ++ff_it)
    {
      int otherFaceId = -1;
      for(int j = 0; j<fuzzyRegion.size(); j++){
        if(mapping[j] == ff_it.handle().idx())
        {
          otherFaceId = j;
        }
      }
      if(otherFaceId != -1)
      {
        double dist = faceDist(mapping[i], ff_it.handle().idx());
        g->add_edge(i, otherFaceId,  dist, dist );
      }
    }
  }


  int flow = g->maxflow();

  printf("Flow = %d\n", flow);
  printf("Minimum cut:\n");

  for(int i = 0; i<fuzzyRegion.size(); i++)
  {
    if (g->what_segment(i) == GraphType::SOURCE){
      sinkRegion.insert(mapping[i]);
      std::cout << i << " is in the SOURCE set\n";
      typename M::FaceHandle fh = mesh.face_handle(mapping[i]);
      mesh.set_color(fh, typename M::Color(255, 0, 0));
    }
    else
    {
      sourceRegion.insert(mapping[i]);
      std::cout << i << " is in the SINK set\n";
      typename M::FaceHandle fh = mesh.face_handle(mapping[i]);
      mesh.set_color(fh, typename M::Color(0, 0, 255));
    }
  }

  //typename GraphType::arc_id a;
  //a = g->get_first_arc();
  //typename GraphType::node_id u, v;
  //printf("Arcos do corte mínimo:\n");
  //for (int i = 0; i < g->get_arc_num(); i++) {
    //g->get_arc_ends(a, u, v);
    //if (g->what_segment(u) != g->what_segment(v)) {
      //printf("%d -> %d\n", u, v);
    //}
    //a = g->get_next_arc(a);
  //}

  delete g;

}


template<typename M>
double
QtModelT<M>::distToSource(int fId)
{
  double dist = 1000.0;
  typename M::FaceHandle fh1 = mesh.face_handle(fId);
  Vec f1 = getFaceCentroid(fh1);
  for ( auto it = sourceRegion.begin(); it != sourceRegion.end(); ++it )
  {
    typename M::FaceHandle fh2 = mesh.face_handle(*it);
    Vec f2 = getFaceCentroid(fh2);
    double d = (f1 - f2).norm();
    if(d < dist)
      dist = d;
  }
  std::cout << "Dist To Source " << dist << "\n";
  return dist;
}

template<typename M>
double
QtModelT<M>::distToSink(int fId)
{
  double dist = 1000.0;
  typename M::FaceHandle fh1 = mesh.face_handle(fId);
  Vec f1 = getFaceCentroid(fh1);
  for ( auto it = sinkRegion.begin(); it != sinkRegion.end(); ++it )
  {
    typename M::FaceHandle fh2 = mesh.face_handle(*it);
    Vec f2 = getFaceCentroid(fh2);
    double d = (f1 - f2).norm();
    if(d < dist)
      dist = d;
  }
  std::cout << "Dist To Sink " << dist << "\n";
  return dist;
}

template<typename M>
double
QtModelT<M>::faceDist(int fId1, int fId2)
{

  double dist = 0.0;
  //Barycenter Euclidian Distance
  typename M::FaceHandle fh1 = mesh.face_handle(fId1);
  Vec f1 = getFaceCentroid(fh1);
  typename M::FaceHandle fh2 = mesh.face_handle(fId2);
  Vec f2 = getFaceCentroid(fh2);
  dist += (f1 - f2).norm();

  //Normal Difference
  Vec n1 = Vec(mesh.normal(fh1)[0], mesh.normal(fh1)[1], mesh.normal(fh1)[2]);
  Vec n2 = Vec(mesh.normal(fh2)[0], mesh.normal(fh2)[1], mesh.normal(fh2)[2]);
  dist += (n1.dot(n2) * 0.1);
  std::cout << "Face Dist " << dist << "\n";
  if(dist < 0.0) dist = 0.001;
  return dist;
}


template<typename M>
void
QtModelT<M>::deleteSink()
{
  for ( auto it = sinkRegion.begin(); it != sinkRegion.end(); ++it )
  {
    typename M::FaceHandle fh = mesh.face_handle(*it);
    mesh.delete_face(fh, false);
  }
  mesh.garbage_collection();
}

template<typename M>
void
QtModelT<M>::toggleFuzzy()
{
  for ( auto it = fuzzyRegion.begin(); it != fuzzyRegion.end(); ++it )
  {
    typename M::FaceHandle fh = mesh.face_handle(*it);
    if(!showFuzzy)
      mesh.set_color(fh, typename M::Color(0, 255, 0));
    else{
     std::unordered_set<int>::const_iterator gotSoR = sourceRegion.find(*it);
     std::unordered_set<int>::const_iterator gotSiR = sinkRegion.find(*it);

     if(gotSoR != sourceRegion.end())
     {
      mesh.set_color(fh, typename M::Color(0, 0, 255));
     }
     else if(gotSiR != sinkRegion.end())
     {
      mesh.set_color(fh, typename M::Color(255, 0, 0));

     }
     else
     {
      mesh.set_color(fh, OpenMesh::Vec3f(modelColor.redF(), modelColor.blueF(), modelColor.greenF()));
     }
    }
 }
  for (size_t i = 0;  i < stroke.size(); i++ )
  {
    typename M::FaceHandle fh = mesh.face_handle(stroke[i]);
    mesh.set_color(fh, typename M::Color(0, 255, 255));

  }
  showFuzzy = !showFuzzy;
}

template<typename M>
void
QtModelT<M>::autoSelect()
{
  
  if(sourceRegion.size() != 0 && sinkRegion.size()!=0)
  {
  std::cout << "Auto Select" << "\n";
  dest = -1;
  prev.clear();
  calcStrokeProxies();
  typename M::FaceHandle face = mesh.face_handle(stroke.front());
  typename M::FaceVertexIter fv_it = mesh.fv_iter(face);
  source = fv_it.handle().idx();
  double distances[mesh.n_vertices()];
  int previous[mesh.n_vertices()];
  std::unordered_set<int> q;
  distances[source] = 0;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    int v = v_it.handle().idx();
    if(v != source)
    {
      distances[v] = 1000.0;
      previous[v] = NULL;
    }
    q.insert(v);
  }
  while(!q.empty())
  {
    double min = 1001.0;
    int u;
    for ( auto i = q.begin(); i != q.end(); ++i )
    {
      if(distances[*i] < min)
      {
        min = distances[*i];
        u = *i;
      }
    }
    min = 1001.0;
    //std::cout << "Q Size Before " << q.size() << "\n";
    //std::cout << "Removed " << u << "\n";
    q.erase(u);
    //std::cout << "Q Size After" << q.size() << "\n";
    typename M::VertexHandle vh = mesh.vertex_handle(u);
    for (typename M::VertexVertexIter vv_it=mesh.vv_iter(vh); vv_it; ++vv_it)
    {
      int v = vv_it.handle().idx();
      float alt = distances[u] + cost(u, v);
      if(alt < distances[v])
      {
        distances[v] = alt;
        previous[v] = u;
      }
    }
  }
  //std::cout << "Q Empty" << "\n";
  typename M::FaceHandle destFace = mesh.face_handle(stroke.back());
  typename M::FaceVertexIter dest_fv_it = mesh.fv_iter(destFace);
  dest = dest_fv_it.handle().idx();
  for(int i = 0; i<mesh.n_vertices(); i++)
  {
    prev.push_back(previous[i]);
    //std::cout << previous[i] << "\n";
    //std::cout << distances[i] << "\n";
  }
  //std::cout << "Prev built" << "\n";
  int to;
  int from = dest;
  int counter = 0;
  while(from != source)
  {
    //std::cout << "Add to Fuzzy" << counter++ << "\n";
    //addToFuzzyRegion(from);
    typename M::VertexHandle vh1 = mesh.vertex_handle(from);
    for (typename M::VertexFaceIter ff_it1=mesh.vf_iter(vh1); ff_it1; ++ff_it1)
    {
        for (typename M::FaceFaceIter ff_it2=mesh.ff_iter(ff_it1.handle()); ff_it2; ++ff_it2)
        {
          //for (typename M::FaceFaceIter ff_it3=mesh.ff_iter(ff_it2.handle()); ff_it3; ++ff_it3)
          //{
            //addToFuzzyRegion(ff_it3.handle().idx());
          //}
          addToFuzzyRegion(ff_it2.handle().idx());
        }
        addToFuzzyRegion(ff_it1.handle().idx());
    }
    //std::cout << "Source " << source << "\n";
    //std::cout << "From " << from << "\n";
    to = prev[from];
    from = to;
  }
  if(createSourceAndSink())
  {
    showFuzzy = true;
    toggleFuzzy();
    graphCut();
  }
  else{
    sinkRegion.clear();
    sourceRegion.clear();
    fuzzyRegion.clear();
    stroke.clear();
    clearColour();
    std::cout << "Bad Path Selected" << "\n"; 
  }
  }
}

template <typename M>
void
QtModelT<M>::flipRegions()
{
  std::unordered_set<int> temp = sourceRegion;
  sourceRegion = sinkRegion;
  sinkRegion = temp;

  for ( auto it = sourceRegion.begin(); it != sourceRegion.end(); ++it )
  {
    typename M::FaceHandle fh = mesh.face_handle(*it); 
    mesh.set_color(fh, typename M::Color(0, 0, 255));
  }
  for ( auto it = sinkRegion.begin(); it != sinkRegion.end(); ++it )
 {
    typename M::FaceHandle fh = mesh.face_handle(*it); 
    mesh.set_color(fh, typename M::Color(255, 0, 0));
  }

}

template <typename M>
M
QtModelT<M>::copy()
{
  std::cout << "Copy" << "\n";
  M m;
  for ( auto it = sinkRegion.begin(); it != sinkRegion.end(); ++it )
  {
    typename M::FaceHandle fh = mesh.face_handle(*it);
    std::vector<typename M::VertexHandle>  face_vhandles;
    for (typename M::FaceVertexIter vf_it=mesh.fv_iter(fh); vf_it; ++vf_it)
    {
      face_vhandles.push_back(m.add_vertex(mesh.point(*vf_it)));
    }
    m.add_face(face_vhandles);
  }
  typedef typename M::Point Point;
  QVector3D tempR = meshRotation * deg2Rad;
  Eigen::AngleAxis<float> aax(tempR.x(), Eigen::Vector3f(1, 0, 0));
  Eigen::AngleAxis<float> aay(tempR.y(), Eigen::Vector3f(0, 1, 0));
  Eigen::AngleAxis<float> aaz(tempR.z(), Eigen::Vector3f(0, 0, 1));
  Eigen::Quaternion<float> rotation = aax * aay * aaz;

  for (typename M::VertexIter v_it=m.vertices_begin(); v_it!=m.vertices_end(); ++v_it) 
  {
    Eigen::Vector3f p = Eigen::Vector3f(m.point(*v_it)[0], m.point(*v_it)[1], m.point(*v_it)[2]);
    p = rotation * p;
    m.set_point( *v_it, Point(p[0], p[1], p[2]) );
    m.set_point( *v_it, m.point(*v_it) + Point(horizontal, vertical, depth) );
  }
  return m;
}

#endif
