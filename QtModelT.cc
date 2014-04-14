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
  , dest(-1)
{
  mesh = m;

  double min_x, max_x, min_y, max_y, min_z, max_z;
  bool first = true;
  boundaryPoints.reserve(mesh.n_vertices()/2);
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


  typedef typename M::Point Point;
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    mesh.set_point( *v_it, Point(
          2.0*(mesh.point(*v_it)[0]-min)/(diff) - 1.0,
          2.0*(mesh.point(*v_it)[1]-min)/(diff) - 1.0,
          2.0*(mesh.point(*v_it)[2]-min)/(diff) - 1.0)
    );
  }

  updateColour();
  calcNormals();
}

template <typename M>
void
QtModelT<M>::findBoundaryVertices(){
  boundaryPoints.clear();
  boundaryPoints.reserve(mesh.n_vertices()/2);
  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
  {
    if (mesh.is_boundary(*v_it)) boundaryPoints.push_back(*v_it);
  }
  PointMatrix m(boundaryPoints.size(), 3);
  for (int count = 0; count < boundaryPoints.size(); ++count)
  {
    VertexHandle v_it = boundaryPoints[count];
    m(count, 0) = mesh.point(v_it)[0];
    m(count, 1) = mesh.point(v_it)[1];
    m(count, 2) = mesh.point(v_it)[2];
    count += 1;
  }
  boundaryMatrix = m;
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

template <typename M>
void
QtModelT<M>::addToFuzzyRegion(int f){

  if(std::find(stroke.begin(), stroke.end(), f) == stroke.end()) {
    fuzzyRegion.insert(f);
    typename M::FaceHandle face = mesh.face_handle(f);
    mesh.set_color(face, typename M::Color(0, 255, 0));
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
    glPushMatrix();
    glTranslatef(horizontal, vertical, zAxis);
    //glRotatef(modelRotation.x(), 1, 0, 0);
    //glRotatef(modelRotation.y(), 0, 1, 0);
    //glRotatef(modelRotation.z(), 0, 0, 1);

  //glEnable(GL_LIGHTING);
  //glShadeModel(GL_FLAT);
  glEnable(GL_DEPTH_TEST);
  
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, mesh.points());
  
  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_FLOAT, 0, mesh.vertex_normals());
  

  for (; fIt!=fEnd; ++fIt)
  {
    glLoadName(index);
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

  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);


  //if(dest >= 0)
  //{
    //glBegin(GL_LINES);
    //glLineWidth(5.0f);
    //glColor3b (0, 255, 0);
    //int to;
    //int from = dest;
    //while(from != source)
    //{
      //to = prev[from];
      //typename M::VertexHandle to_vh = mesh.vertex_handle(to);
      //typename M::VertexHandle from_vh = mesh.vertex_handle(from);
      //glVertex3f(mesh.point(from_vh)[0], mesh.point(from_vh)[1], mesh.point(from_vh)[2]);
      //glVertex3f(mesh.point(to_vh)[0], mesh.point(to_vh)[1], mesh.point(to_vh)[2]);
      //from = to;
    //}
    ////for(int i = 0; i<strokeVertices.size()-1; i++)
    ////{
        ////glVertex3f(strokeVertices[i][0], strokeVertices[i][1], strokeVertices[i][2]);
        ////glVertex3f(strokeVertices[i+1][0], strokeVertices[i+1][1], strokeVertices[i+1][2]);
    ////}
    //glEnd();
  //}
  /*
    glEnable(GL_LIGHTING);
    glShadeModel(GL_FLAT);

    glEnable(GL_DEPTH_TEST);
    glEnableClientState(GL_VERTEX_ARRAY);
  

    for (; fIt!=fEnd; ++fIt)
    {
        glLoadName(index);
        glBegin(GL_TRIANGLES);
        glNormal3fv( &mesh.normal(*fIt)[0] );
        fvIt = mesh.cfv_iter(*fIt);
        glVertex3fv( &mesh.point(*fvIt)[0] );
        ++fvIt;
        glVertex3fv( &mesh.point(*fvIt)[0] );
        ++fvIt;
        glVertex3fv( &mesh.point(*fvIt)[0] );
        glEnd();
        index++;
     }
*/
    //glBegin(GL_LINES);
    //glLineWidth(2.0f);
    //for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
    //{
      //glColor3b (255, 255, 255);
      //glVertex3f(mesh.point(*v_it)[0], mesh.point(*v_it)[1], mesh.point(*v_it)[2]);
      //glVertex3f(mesh.point(*v_it)[0]+mesh.normal(*v_it)[0], mesh.point(*v_it)[1]+mesh.normal(*v_it)[1], mesh.point(*v_it)[2]+mesh.normal(*v_it)[2]);

    //}
    //glEnd();

    glPopMatrix();

}

template <typename M>
void
QtModelT<M>::applyTransformations()
{
  typedef typename M::Point Point;
  modelRotation = modelRotation * deg2Rad;
  Eigen::AngleAxis<float> aax(modelRotation.x(), Eigen::Vector3f(1, 0, 0));
  Eigen::AngleAxis<float> aay(modelRotation.y(), Eigen::Vector3f(0, 1, 0));
  Eigen::AngleAxis<float> aaz(modelRotation.z(), Eigen::Vector3f(0, 0, 1));
  Eigen::Quaternion<float> rotation = aax * aay * aaz;

  for (typename M::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) 
  {
    Eigen::Vector3f p = Eigen::Vector3f(mesh.point(*v_it)[0], mesh.point(*v_it)[1], mesh.point(*v_it)[2]);
    p = rotation * p;
    mesh.set_point( *v_it, Point(p[0], p[1], p[2]) );
    mesh.set_point( *v_it, mesh.point(*v_it) + Point(horizontal, vertical, depth) );
  }
  horizontal = 0.0f;
  vertical = 0.0f;
  depth = 0.0f;
  modelRotation.setX(0.0f);
  modelRotation.setY(0.0f);
  modelRotation.setZ(0.0f);
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


template <typename M>
PointMatrix
QtModelT<M>::buildSampledMatrix()
{
  int noSamples = 5000;
  PointMatrix allMat = buildMatrix();
  PointMatrix randMat(noSamples, 3);
  for ( unsigned i = 0U; i < noSamples; ++i )
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
  double x, y, z = 0.0;
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
  //modelRotation += rotationVec;
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
}

template <typename M>
void
QtModelT<M>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
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
void
QtModelT<M>::cut()
{
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
          for (typename M::FaceFaceIter ff_it3=mesh.ff_iter(ff_it2.handle()); ff_it3; ++ff_it3)
          {
            addToFuzzyRegion(ff_it3.handle().idx());
          }
          addToFuzzyRegion(ff_it2.handle().idx());
        }
        addToFuzzyRegion(ff_it1.handle().idx());
    }
    //std::cout << "Source " << source << "\n";
    //std::cout << "From " << from << "\n";
    to = prev[from];
    from = to;
  }
  createSourceAndSink();
  graphCut();
  fuzzyRegion.clear();
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
void
QtModelT<M>::createSourceAndSink()
{
  sourceRegion.clear();
  sinkRegion.clear();
  std::cout << "Source and Sink" << "\n";
  bool notUnique = true;
  int assignedFaces = fuzzyRegion.size();
  while(notUnique)
  {
    int f = rand() % mesh.n_faces();
    if(!inRegion(f))
    {
      regionGrow(f, &sourceRegion, 0);
      notUnique = false;
    }
  }
  notUnique = true;
  while(notUnique)
  {
    int f = rand() % mesh.n_faces();
    if(!inRegion(f))
    {
      regionGrow(f, &sinkRegion, 1);
      notUnique = false;
    }
  }
  assignedFaces = fuzzyRegion.size() + sinkRegion.size() + sourceRegion.size();
  //std::cout << "In Fuzzy Region " << fuzzyRegion.size() << " / " << mesh.n_faces() << "\n";
  //std::cout << "In Source Region " << sourceRegion.size() << " / " << mesh.n_faces() << "\n";
  //std::cout << "In Sink Region " << sinkRegion.size() << " / " << mesh.n_faces() << "\n";
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

template<typename M>
void
QtModelT<M>::graphCut()
{
  std::cout << "Graph Cut" << "\n";
}


#endif
