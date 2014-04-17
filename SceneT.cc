#ifndef SCENE_CC
#define SCENE_CC

#include <iostream>
#include "Scene.hh"
#include <QVector3D>
#include "GLUT/glut.h"
#include <QtGui>
//#include <QtOpenGL>
#include <QDoubleSpinBox>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Tools/Utils/Timer.hh>
#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif


template <typename M>
QDialog
*SceneT<M>::createDialog(const QString &windowTitle) const
{
  QDialog *dialog = new QDialog(0, Qt::CustomizeWindowHint | Qt::WindowTitleHint);

  dialog->setWindowOpacity(0.8);
  dialog->setWindowTitle(windowTitle);
  dialog->setLayout(new QVBoxLayout);

  return dialog;
}

template <typename M>
SceneT<M>::SceneT()
: m_backgroundColor(0.0f, 0.0f, 0.0f)
, m_distance(4.5f)
, m_vertical(0.0f)
, m_horizontal(0.0f)
, TANSLATE_SPEED(0.01f)
, deg2Rad(0.0174532925)
, inPaintingMode(false)
, mouseClicked(false)
{
  modelCount = 0;
  QWidget *controls = createDialog(tr("Controls"));

  m_modelButton = new QPushButton(tr("Load model"));
  controls->layout()->addWidget(m_modelButton);

  groupBox = new QGroupBox(tr("Select Mesh"));
  radio1 = new QRadioButton(tr("All"));
  radio2 = new QRadioButton(tr("M1"));
  radio3 = new QRadioButton(tr("M2"));
  radio4 = new QRadioButton(tr("M3"));
  radio5 = new QRadioButton(tr("M4"));
  radio6 = new QRadioButton(tr("M5"));
  radio7 = new QRadioButton(tr("M6"));
  radio8 = new QRadioButton(tr("M7"));
  radio9 = new QRadioButton(tr("M8"));
  radio10 = new QRadioButton(tr("M9"));
  radio11 = new QRadioButton(tr("M10"));

  groupBox->setHidden(true);
  groupBox->setFocusPolicy(Qt::StrongFocus);
  radio3->setHidden(true);
  radio4->setHidden(true);
  radio5->setHidden(true);
  radio6->setHidden(true);
  radio7->setHidden(true);
  radio8->setHidden(true);
  radio9->setHidden(true);
  radio10->setHidden(true);
  radio11->setHidden(true);

  radio1->setChecked(true);
  QVBoxLayout *vbox = new QVBoxLayout;
  vbox->addWidget(radio1);
  vbox->addWidget(radio2);
  vbox->addWidget(radio3);
  vbox->addWidget(radio4);
  vbox->addWidget(radio5);
  vbox->addWidget(radio6);
  vbox->addWidget(radio7);
  vbox->addWidget(radio8);
  vbox->addWidget(radio9);
  vbox->addWidget(radio10);
  vbox->addWidget(radio11);
  vbox->addStretch(1);
  groupBox->setLayout(vbox);
  controls->layout()->addWidget(groupBox);

  removeModelButton = new QPushButton(tr("Remove model"));
  controls->layout()->addWidget(removeModelButton);
  removeModelButton->setHidden(true);

  mouseControlBox = new QGroupBox(tr("Mouse Control"));
  mouseControlBox->setHidden(true);
  translateRadio = new QRadioButton(tr("Translate"));
  rotateRadio = new QRadioButton(tr("Rotate"));
  paintFacesRadio = new QRadioButton(tr("Paint Faces (with Alt)"));
  translateRadio->setChecked(true);
  QVBoxLayout *vb = new QVBoxLayout;
  vb->addWidget(translateRadio);
  vb->addWidget(rotateRadio);
  vb->addWidget(paintFacesRadio);
  vb->addStretch(1);
  mouseControlBox->setLayout(vb);
  controls->layout()->addWidget(mouseControlBox);

  removeModelButton = new QPushButton(tr("Remove model"));
  controls->layout()->addWidget(removeModelButton);
  removeModelButton->setHidden(true);

  clearFacesButton = new QPushButton(tr("Clear Painted Faces"));
  controls->layout()->addWidget(clearFacesButton);
  clearFacesButton->setHidden(true);

  cutButton = new QPushButton(tr("Cut"));
  controls->layout()->addWidget(cutButton);
  cutButton->setHidden(true);

  copyButton = new QPushButton(tr("Copy"));
  controls->layout()->addWidget(copyButton);
  copyButton->setHidden(true);


  pasteButton = new QPushButton(tr("Paste"));
  controls->layout()->addWidget(pasteButton);
  pasteButton->setHidden(true);

  deleteButton = new QPushButton(tr("Delete"));
  controls->layout()->addWidget(deleteButton);
  deleteButton->setHidden(true);


  geoTreeSpinBox = new QSpinBox();
  geoTreeSpinBox->setMinimum(1);
  geoTreeSpinBox->setMaximum(20);
  geoTreeSpinBox->setPrefix("k: ");
  controls->layout()->addWidget(geoTreeSpinBox);
  geoTreeSpinBox->setHidden(true);

  geoTreeButton = new QPushButton(tr("GeoTree"));
  controls->layout()->addWidget(geoTreeButton);
  geoTreeButton->setHidden(true);


  //QWidget *widgets[] = { meshes, controls, examples  };
  QWidget *widgets[] = { controls };

  for (uint i = 0; i < sizeof(widgets) / sizeof(*widgets); ++i) {
    QGraphicsProxyWidget *proxy = new QGraphicsProxyWidget(0, Qt::Dialog);
    proxy->setWidget(widgets[i]);
    addItem(proxy);
  }

  QPointF pos(10, 10);
  foreach (QGraphicsItem *item, items()) {
    item->setFlag(QGraphicsItem::ItemIsMovable);
    item->setCacheMode(QGraphicsItem::DeviceCoordinateCache);
    const QRectF rect = item->boundingRect();
    item->setPos(pos.x() - rect.x(), pos.y() - rect.y());
    pos += QPointF(0, 10 + rect.height());
  }
}

template <typename M>
void
SceneT<M>::setDefaultMaterial(void)
{
  GLfloat mat_a[] = {0.1, 0.1, 0.1, 1.0};
  GLfloat mat_d[] = {0.7, 0.7, 0.5, 1.0};
  GLfloat mat_s[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat shine[] = {120.0};
  
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_a);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_d);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_s);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
}


//----------------------------------------------------------------------------

template <typename M>
void
SceneT<M>::setDefaultLight(void)
{
  GLfloat pos1[] = { 0.1,  0.1, -0.02, 0.0};
  GLfloat pos2[] = {-0.1,  0.1, -0.02, 0.0};
  GLfloat pos3[] = { 0.0,  0.0,  0.1,  0.0};
  GLfloat col1[] = { 1.0,  1.0,  1.0,  0.3};
  GLfloat col2[] = { 1.0,  1.0,  1.0,  0.3};
  GLfloat col3[] = { 1.0,  1.0,  1.0,  0.3};
 
  glEnable(GL_LIGHT0);    
  glLightfv(GL_LIGHT0,GL_POSITION, pos1);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,  col1);
  glLightfv(GL_LIGHT0,GL_SPECULAR, col1);
  
  glEnable(GL_LIGHT1);  
  glLightfv(GL_LIGHT1,GL_POSITION, pos2);
  glLightfv(GL_LIGHT1,GL_DIFFUSE,  col2);
  glLightfv(GL_LIGHT1,GL_SPECULAR, col2);
  
  glEnable(GL_LIGHT2);  
  glLightfv(GL_LIGHT2,GL_POSITION, pos3);
  glLightfv(GL_LIGHT2,GL_DIFFUSE,  col3);
  glLightfv(GL_LIGHT2,GL_SPECULAR, col3);
}



template <typename M>
void
SceneT<M>::drawForeground(QPainter *painter, const QRectF &rect)
{
  //std::cout << "Draw Foreground" << "\n";
  painter->beginNativePainting();

  if(modelCount > 0){
    glSelectBuffer(65535, PickBuffer);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluPerspective(70, width() / height(), 0.01, 1000);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glEnable(GL_LIGHTING);
    glShadeModel(GL_FLAT);
    glTranslatef(m_horizontal, m_vertical, -m_distance);
    glRotatef(m_rotation.x(), 1, 0, 0);
    glRotatef(m_rotation.y(), 0, 1, 0);
    glRotatef(m_rotation.z(), 0, 0, 1);
    glBegin(GL_LINES);
    glLineWidth(2);
    glColor3f(255,255,255);
    glVertex3f(0, 0, 0);
    glVertex3f(2, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 2, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 2);
    glEnd();
    
    glEnable(GL_MULTISAMPLE);
    glColorMaterial ( GL_FRONT_AND_BACK, GL_EMISSION ) ;
    glEnable( GL_COLOR_MATERIAL ) ;
    for (int i = 0; i != modelCount; i++) {
      if(models[i] != NULL) models[i]->render();
    }
    setDefaultMaterial();
    setDefaultLight();
    glDisable(GL_MULTISAMPLE);
    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
  }
  //painter->endNativePainting();
}

template <typename M>
void
SceneT<M>::softICP(QtModelT<M>* m1, QtModelT<M>* m2)
{
  typedef nanoflann::KDTreeEigenMatrixAdaptor< PointMatrix > my_kd_tree_t;
  m1->applyTransformations();
  m2->applyTransformations();
  double sizeM1 = computeSnapRegionSize(m1,m2);
  double sizeM2 = computeSnapRegionSize(m2,m1);
  double snapSize = sizeM1 < sizeM2 ? sizeM1 : sizeM2;
  snapSize = sizeM1;
  std::vector<size_t> m1SnapRegion = computeSnapRegion(m1, snapSize);
  std::vector<size_t> m2SnapRegion = computeSnapRegion(m2, snapSize);
  std::cout << snapSize << " " << sizeM1 << sizeM2 << " snap\n";

  //OpenMesh::VPropHandleT< double > gauss;
  int iterations = 1;
  int iterCount = 0;
  int vcount = m1->mesh.n_vertices();
  while (iterCount < iterations){
    float scaling[vcount];
    std::vector< Matrix<double, 3, 3> > rotations;
    rotations.reserve(m1->mesh.n_vertices());
    Matrix<double, Dynamic, 3> translations(vcount,3);
    
    //Find correspondence φ of SA to SB
    std::vector<size_t> correspondVector = findCorrespondence(m1,m2,m1SnapRegion,m2SnapRegion);
    int elasticity = 1;
    my_kd_tree_t b_mat_index(3, m1->boundaryMatrix, 10);
    b_mat_index.index->buildIndex();
    PointMatrix m1Matrix = m1->buildMatrix();

    //For each point pi in MA Find the local neighborhood N(pi) Calculate the transformation Ti based on φ|N(pi)
    for (size_t vi = 0; vi < vcount; vi++){

      std::vector<double> query_pt(3);
      query_pt[0] =  m1Matrix(vi,0);
      query_pt[1] =  m1Matrix(vi,1);
      query_pt[2] =  m1Matrix(vi,2);
      
      //calculate distance to snapping region
      std::vector<size_t> ret_index(1);
      std::vector<double> out_dist_sqr(1);
      nanoflann::KNNResultSet<double> resultSet(1);
      resultSet.init(&ret_index[0], &out_dist_sqr[0]);
      
      b_mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(0));
      double dist = pow(out_dist_sqr[0],0.5) - snapSize;
      
      if (dist == 0) dist = DBL_MIN;
      long indices = pow((iterCount*elasticity)/dist,2);
      std::cout << exp(-indices) << " exp\n";
      double radius = snapSize * exp(-indices);
      
      std::cout << radius << " radius\n";
      //for points outside the snapping region offset with their distance to the snapping region
      if (dist > 0) {
        radius += dist;
        std::cout << radius << " radius2\n";
      }
      std::vector<size_t> neighbours = findLocalNeighbourhood(query_pt, m1Matrix, radius);
      std::vector<size_t> intersection;
      intersection.reserve(m1SnapRegion.size());
      
      std::sort(neighbours.begin(), neighbours.end());
      std::sort(m1SnapRegion.begin(), m1SnapRegion.end());
      std::set_intersection(neighbours.begin(),neighbours.end(),m1SnapRegion.begin(),m1SnapRegion.end(),std::back_inserter(intersection));
      std::cout << intersection.size() << " " << neighbours.size() << " " << m1SnapRegion.size() << " inter\n";

      //bounding boxes
      Vec3f lbbMin, lbbMax, cbbMin, cbbMax;
      //lbbMin = lbbMax = OpenMesh::vector_cast<Vec3f>(m1->mesh.point(m1->mesh.vertex_handle(intersection[0])));
      
      PointMatrix localMatrix(intersection.size(),3);
      PointMatrix correspondMatrix(intersection.size(),3);
      int x = 0;
      for (size_t i = 0; i < m1SnapRegion.size(); i++){
        for (size_t j = 0; j < intersection.size(); j++){
          if (m1SnapRegion[i] == intersection[j]){
            VertexHandle vh = m1->mesh.vertex_handle(intersection[j]);
            Point p = m1->mesh.point(vh);
            lbbMin.minimize( OpenMesh::vector_cast<Vec3f>(p));
            lbbMax.maximize( OpenMesh::vector_cast<Vec3f>(p));
            localMatrix(x,0) = p[0];
            localMatrix(x,1) = p[1];
            localMatrix(x,2) = p[2];
            
            VertexHandle c_vh = m2->mesh.vertex_handle(correspondVector[i]);
            Point q = m2->mesh.point(c_vh);
            cbbMin.minimize( OpenMesh::vector_cast<Vec3f>(q));
            cbbMax.maximize( OpenMesh::vector_cast<Vec3f>(q));
            correspondMatrix(x,0) = q[0];
            correspondMatrix(x,1) = q[1];
            correspondMatrix(x,2) = q[2];
            x++;
          }
        }
      }
      //std::cout << correspondMatrix(0,0) << "corMat \n";
      float scale = (lbbMax - lbbMin).length() / (cbbMax - cbbMax).length();
      Matrix<double, 1, 3> localMean = localMatrix.colwise().mean();
      Matrix<double, 1, 3> coreMean = correspondMatrix.colwise().mean();

      PointMatrix qHat(correspondMatrix.rows(), 3);
      PointMatrix pHat(localMatrix.rows(), 3);
      
      calcBaryCenteredPoints(qHat, correspondMatrix);
      calcBaryCenteredPoints(pHat, localMatrix);
      
      Matrix<double, 3, 3> A(3,3);
      generateAMatrix(A, pHat, qHat);
      
      JacobiSVD< MatrixXd > svd(A, ComputeThinU | ComputeThinV);
      
      Matrix<double, 3, 3> R = svd.matrixU() * svd.matrixV().transpose();
      Matrix<double, 1, 3> temp1 = (R * localMean.transpose());
      Matrix<double, 1, 3> temp2 = coreMean;
      Matrix<double, 1, 3> t = temp1 - temp2;
      
      float li = iterCount/iterations;
      int it = 0;
      for (typename M::VertexIter v_it=m1->mesh.vertices_begin(); v_it!=m1->mesh.vertices_end(); ++v_it)
      {
        Matrix<double, 3, 3> R = rotations[it];
        float S = scaling[it];
        Matrix<double, 1, 3> T = translations.row(it);
        R = R * li;
        S = S * li;
        T = T * li;
        //scaling
        Point vertex = m1->mesh.point(*v_it);
        Eigen::Vector3d p = Eigen::Vector3d(vertex[0], vertex[1], vertex[2]);
        p = R * p;
        m1->mesh.set_point( *v_it, Point(p[0], p[1], p[2]) +  Point(T[0], T[1], T[2]));
        it++;
      }
    }
    std::cout << "Iteration " << iterCount << "\n";
    iterCount++;
  }
}

template <typename M>
void
runSoftICP(){
  
}

template <typename M>
std::vector<size_t>
findCorrespondence(QtModelT<M>* m1, QtModelT<M>* m2, std::vector<size_t> m1SnapRegion, std::vector<size_t> m2SnapRegion){
  std::vector<size_t> correspondVector;
  correspondVector.reserve(m1SnapRegion.size());
  for (size_t i = 0; i < m1SnapRegion.size(); i++){
    VertexHandle m1Handle = m1->mesh.vertex_handle(m1SnapRegion[i]);
    typename M::Normal m1Norm = m1->mesh.normal(m1Handle);
    Point p = m1->mesh.point(m1Handle);
    float closest = 999;
    size_t correspond = 0;
    for (size_t j = 0; j < m2SnapRegion.size(); j++){
      VertexHandle m2Handle = m2->mesh.vertex_handle(m2SnapRegion[j]);
      typename M::Normal m2Norm = m2->mesh.normal(m2Handle);
      Point q = m2->mesh.point(m2Handle);
      float euclidist = (q-p).length();
      float normdist = acosf(m1Norm|m2Norm);
      float gaussdist = m1->mesh.data(m1Handle).gauss() - m2->mesh.data(m2Handle).gauss();//set_gaussian//m1->mesh.property(gauss,m1Handle);
      float wdist = 0.6 * euclidist + 0.2 * normdist + 0.2 * gaussdist;
      if (wdist < closest){
        closest = wdist;
        correspond = j;
      }
    }
    correspondVector.push_back(correspond);
  }
  return correspondVector;
}

template <typename M>
std::vector<size_t>
SceneT<M>::findLocalNeighbourhood(std::vector<double> query_pt, PointMatrix m1Matrix, double radius){
  typedef nanoflann::KDTreeEigenMatrixAdaptor< PointMatrix >  my_kd_tree_t;
  my_kd_tree_t mat_index(3, m1Matrix, 10);
  mat_index.index->buildIndex();
  std::vector<size_t> neighbours;
  std::vector< std::pair< size_t, double > > resultPairs;
  resultPairs.reserve(m1Matrix.rows());
  size_t count = mat_index.index->radiusSearch(&query_pt[0], radius, resultPairs, nanoflann::SearchParams(true));
  neighbours.reserve(count);
  for (size_t j = 0; j < count; j++){
    neighbours.push_back(resultPairs[j].first);
  }
  return neighbours;
}

template <typename M>
double
computeSnapRegionSize(QtModelT<M>* m1, QtModelT<M>* m2)
{
  typedef nanoflann::KDTreeEigenMatrixAdaptor< PointMatrix >  my_kd_tree_t;
  
  const size_t num_results = 1;
  PointMatrix m1Matrix = m1->buildMatrix();
  PointMatrix m2Matrix = m2->buildMatrix();
  //We find the set of closest points to the other shape boundary loop
  
  //m1 find closest points on m2 boundary
  my_kd_tree_t mat_index(3, m1Matrix, 10);
  mat_index.index->buildIndex();
  
  int m2bSize = m2->boundaryPoints.size();
  Matrix<double, Dynamic, 1> mappings(m2bSize, 1);
  Matrix<double, Dynamic, 1> distances(m2bSize, 1);
  
  //for each point in m2 boundary loop find closest point in m1
  for(int i=0; i < m2bSize; i++)
  {
    std::vector<double> query_pt(3);
    query_pt[0] =  m2->boundaryMatrix(i,0);
    query_pt[1] =  m2->boundaryMatrix(i,1);
    query_pt[2] =  m2->boundaryMatrix(i,2);
    
    std::vector<size_t> ret_index(num_results);
    std::vector<double> out_dist_sqr(num_results);
    
    nanoflann::KNNResultSet<double> resultSet(num_results);
    
    resultSet.init(&ret_index[0], &out_dist_sqr[0] );
    mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(0));
    
    
    glVertex3f(query_pt[0], query_pt[1], query_pt[2]);
    glVertex3f(m1Matrix(ret_index[0],0),m1Matrix(ret_index[0],1),m1Matrix(ret_index[0],2));
    
    mappings(i, 0) = ret_index[0];
    distances(i, 0) = out_dist_sqr[0];
  }
  
  //Next, we compute the geodesic distance for each point in the set to their own boundary loop.
  my_kd_tree_t b_mat_index(3, m1->boundaryMatrix, 10);
  b_mat_index.index->buildIndex();
  double snapMax = 0;
  
  for(int i=0; i < m2bSize; i++)
  {
    int setIndex = mappings(i,0);
    std::vector<double> query_pt(3);
    query_pt[0] =  m1Matrix(setIndex,0);
    query_pt[1] =  m1Matrix(setIndex,1);
    query_pt[2] =  m1Matrix(setIndex,2);
    //m1->colourFaceFromVertexIndex(m1->boundaryPoints[i]);
    std::vector<size_t> ret_index(num_results);
    std::vector<double> out_dist_sqr(num_results);
    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_index[0], &out_dist_sqr[0]);
    
    b_mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(0));
    //We define the snapping region size R as the maximum of the geodesic distances.
    double dist = pow(out_dist_sqr[0],0.5);
    if (dist > snapMax) {
      std::cout << dist << " > " << snapMax << "\n";
      snapMax = dist;
    }
  }
  std::cout << snapMax << " maxdist \n";
  return snapMax;
}

template <typename M>
std::vector<size_t>
SceneT<M>::computeSnapRegion(QtModelT<M>* m1, float snapMax)
{
  typedef nanoflann::KDTreeEigenMatrixAdaptor< PointMatrix >  my_kd_tree_t;
  PointMatrix m1Matrix = m1->buildMatrix();
  my_kd_tree_t mat_index(3, m1Matrix, 10);
  mat_index.index->buildIndex();
  std::vector<size_t> snapRegion;
  //Thus, the sub-mesh within R geodesic distance from the boundary loop is defined as the snapping region
  //for each boudary point find all points within radius
  for(size_t i=0; i < m1->boundaryPoints.size(); i++)
  {
    std::vector<double> query_pt(3);
    query_pt[0] = m1->boundaryMatrix(i,0);
    query_pt[1] = m1->boundaryMatrix(i,1);
    query_pt[2] = m1->boundaryMatrix(i,2);
    
    std::vector< std::pair< size_t, double > > resultPairs;
    resultPairs.reserve(m1->getNoVerticies());
    size_t count = mat_index.index->radiusSearch(&query_pt[0], snapMax, resultPairs, nanoflann::SearchParams(true));
    snapRegion.reserve(count);
    for (size_t j = 0; j < count; j++){
      m1->colourFaceFromVertexIndex(resultPairs[j].first);
      snapRegion.push_back(resultPairs[j].first);
    }
  }
  //std::cout << snapRegion.size() << " ";
  std::sort(snapRegion.begin(), snapRegion.end());
  snapRegion.erase(std::unique(snapRegion.begin(), snapRegion.end()), snapRegion.end());
  //std::cout << snapRegion.size() << " unique\n";
  return snapRegion;
}

template <typename M>
void
SceneT<M>::generateAMatrix(Matrix<double, 3, 3>  &A, const PointMatrix &qHat, const PointMatrix &pHat)
{
  A = (pHat.transpose()*qHat)*(qHat.transpose()*qHat).inverse();
}

template <typename M>
void
SceneT<M>::calcBaryCenteredPoints(PointMatrix &matHat, const PointMatrix &mat)
{
  Matrix<double, 1, 3> mean = mat.colwise().mean();
  matHat = mat.rowwise() - mean;
}

/*
template <typename M>
bool
runICPTwoMeshes(QtModelT<M>* m1, QtModelT<M>* m2, int iterations)
{
  //q = still = m1, p = mover = m2;

  bool notConverged = true;
  int iterCount = 0;
  while(notConverged)
  {
    const float max_range = 0.1;
    
    const size_t num_results = 1;
    PointMatrix qAll = m1->buildMatrix();
    PointMatrix pAll = m2->buildMatrix();
    
    typedef nanoflann::KDTreeEigenMatrixAdaptor< PointMatrix >  my_kd_tree_t;
    my_kd_tree_t mat_index(3, qAll, 10);
    mat_index.index->buildIndex();
    
    int rowCount = m2->getNoVerticies();
    std::cout << rowCount << "rows \n";
    Matrix<double, Dynamic, 1> mappings(rowCount, 1);
    Matrix<double, Dynamic, 1> distances(rowCount, 1);
    
    for(int i=0; i < rowCount; i++)
    {
      std::vector<double> query_pt(3);
      query_pt[0] = pAll(i, 0);
      query_pt[1] = pAll(i, 1);
      query_pt[2] = pAll(i, 2);
      
      std::vector<size_t>   ret_index(num_results);
      std::vector<double>    out_dist_sqr(num_results);
      
      nanoflann::KNNResultSet<double> resultSet(num_results);
      
      resultSet.init(&ret_index[0], &out_dist_sqr[0] );
      mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(0));
      
      mappings(i, 0) = ret_index[0];
      distances(i, 0) = out_dist_sqr[0];
    }
    
    double distThreshold = (max_range < 4 * distances.mean()) ? max_range : 4 * distances.mean();
    std::vector<int> goodPairs;
    goodPairs.reserve(pAll.rows());
    for(int i = 0; i < rowCount; i++)
    {
      if(distances(i, 0) <= distThreshold)
        goodPairs.push_back(i);
    }
    
    PointMatrix q(goodPairs.size(), 3);
    PointMatrix p(goodPairs.size(), 3);
    for(int c = 0; c < (int)goodPairs.size(); c++)
    {
      p.row(c) = pAll.row(goodPairs[c]);
      q.row(c) = qAll.row(mappings(goodPairs[c], 0));
      //std::cout << "Mapping" << p.row(c) << "\t" << q.row(c) << "\t" << (p.row(c)-q.row(c)).norm()  << "\n";
    }
    
    PointMatrix qHat(q.rows(), 3);
    PointMatrix pHat(p.rows(), 3);
    
    getBaryCenteredPoints(qHat, q);
    getBaryCenteredPoints(pHat, p);
    
    Matrix<double, 3, 3> A(3,3);
    generateA(A, pHat, qHat);
    
    JacobiSVD< MatrixXd > svd(A, ComputeThinU | ComputeThinV);
    
    Matrix<double, 3, 3> R = svd.matrixU() * svd.matrixV().transpose();
    Matrix<double, 1, 3> temp1 = (R * p.colwise().mean().transpose());
    Matrix<double, 1, 3> temp2 = q.colwise().mean();
    Matrix<double, 1, 3> t = temp1 - temp2;
    std::cout << "Iteration " << iterCount << "\n";
    
    //std::cout << R << "\n";
    //std::cout << t << "\n";
    
    if (isnan(R(0,0) + R(1,1) + R(2,2)) || isinf(t.norm()) ) {
      std::cout << "threw error " << "\n\n";
      std::cout << "Sum of diagonal of R " << (R(0,0) + R(1,1) + R(2,2)) << "\n";
      std::cout << "Norm of t " << t.norm() << "\n\n";
      return false;
    }
    m2->updateTransformations(R, t(0, 0), t(0,1), t(0, 2));
    
    std::cout << "Sum of diagonal of R " << (R(0,0) + R(1,1) + R(2,2)) << "\n";
    std::cout << "Norm of t " << t.norm() << "\n\n";
    
    if(t.norm() < 0.00001 && (R(0,0) + R(1,1) + R(2,2)) > 2.999)
    {
      std::cout << "Converged in " << iterCount << " iterations"  << "\n";
      notConverged = false;
      return true;
    } else if(iterCount >= iterations-1)
    {
      notConverged = false;
      std::cout << "No convergence after " << iterations << " iterations"   << "\n";
      return false;
    }
    iterCount++;
  }
  return false;
}
*/


template <typename M>
void
SceneT<M>::drawBackground(QPainter *painter, const QRectF &)
{
  if (painter->paintEngine()->type() != QPaintEngine::OpenGL
      && painter->paintEngine()->type() != QPaintEngine::OpenGL2)
  {
    qWarning("OpenGLScene: drawBackground needs a QGLWidget to be set as viewport on the graphics view");
    return;
  }
  painter->beginNativePainting();
  glClearColor(m_backgroundColor.redF(), m_backgroundColor.greenF(), m_backgroundColor.blueF(), 1.0f);
  //glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  painter->endNativePainting();
}

template <typename M>
void
SceneT<M>::loadMesh(const QString filePath)
{
  if (filePath.isEmpty())
    return;

  std::cout << filePath.toStdString() << "\n";
  m_modelButton->setEnabled(false);
  QApplication::setOverrideCursor(Qt::BusyCursor);
  if(OpenMesh::IO::read_mesh(m_mymesh, filePath.toStdString(), _options))
  {
    if(modelCount > 9)
    {
      std::cout << "Too Many Models. Please Remove All Models" << "\n";
    }
    else
    {
      models[modelCount] = new QtModelT<M>(m_mymesh);
      models[modelCount]->modelnumber = modelCount;
      modelCount++;
      switch(modelCount)
      {
        case 1:
          removeModelButton->setHidden(false);
          groupBox->setHidden(false);
          mouseControlBox->setHidden(false);
          clearFacesButton->setHidden(false);
          cutButton->setHidden(false);
          deleteButton->setHidden(false);
          pasteButton->setHidden(false);
          copyButton->setHidden(false);
          geoTreeButton->setHidden(false);
          geoTreeSpinBox->setHidden(false);
          break;
        case 2:
          radio3->setHidden(false);
          break;
        case 3:
          radio4->setHidden(false);
          break;
        case 4:
          radio5->setHidden(false);
          break;
        case 5:
          radio6->setHidden(false);
          break;
        case 6:
          radio7->setHidden(false);
          break;
        case 7:
          radio8->setHidden(false);
          break;
        case 8:
          radio9->setHidden(false);
          break;
        case 9:
          radio10->setHidden(false);
          break;
        case 10:
          radio11->setHidden(false);
          break;
      }
    }
    clickRadioButton(modelCount+1);

    //std::clog << m_mymesh.n_vertices() << " vertices, "
    //<< m_mymesh.n_edges()    << " edge, "
    //<< m_mymesh.n_faces()    << " faces\n";
  }
  else
  {
    std::cout << "Error Loading Mesh" << "\n";
  }
  m_modelButton->setEnabled(true);
  QApplication::restoreOverrideCursor();

}


template <typename M>
void
SceneT<M>::wheelEvent(QGraphicsSceneWheelEvent *event)
{
  QGraphicsScene::wheelEvent(event);
  if (event->isAccepted())
    return;

  const int radioId = whichRadioButton();
  if(radioId == 1)
  {
    m_distance *= qPow(1.2, -event->delta() / 120.);
    //std::cout << m_distance << "\n";
  }
  else
  {
    models[radioId-2]->scale(qPow(1.2, -event->delta() / 120.));
  }
  event->accept();
  update();
}

template <typename M>
void
SceneT<M>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  QGraphicsScene::mouseMoveEvent(event);
 if(inPaintingMode && mouseClicked) paintFaces(event); 
  
  if (event->isAccepted())
    return;
  if (event->buttons() & Qt::LeftButton) {
    QPointF delta = event->scenePos() - event->lastScenePos();
    const int radioId = whichRadioButton();
    QVector3D angularImpulse = QVector3D(delta.y(), delta.x(), 0) * 0.1;
    if (radioId == 0){
      
    }
    if(mouseRadioSelected() == 1)
    {
      if(radioId  == 1){
        //std::cout << m_distance << "\n";
        m_vertical -= delta.y() * (TANSLATE_SPEED);
        m_horizontal += delta.x() * (TANSLATE_SPEED);
      }
      else//if(radioId  != 1)
      {

        moveMeshInOneAxis(event);
        models[radioId-2]->updateHorizontal(delta.x() * TANSLATE_SPEED);
        models[radioId-2]->updateVertical(delta.y() * TANSLATE_SPEED);
        //models[radioId-2]->updateHorizontal(p[0] * TANSLATE_SPEED);
        //models[radioId-2]->updateVertical(p[1] * TANSLATE_SPEED);
        //models[radioId-2]->updateZAxis(p[2] * TANSLATE_SPEED);

        //models[radioId-2]->updateHorizontal(delta.x() * TANSLATE_SPEED);
        //models[radioId-2]->updateVertical(delta.y() * TANSLATE_SPEED);
        //models[radioId-2]->updateZAxis(0 * TANSLATE_SPEED);

      }
    }
    if(mouseRadioSelected() == 2)
    {
      
      if(radioId  == 1){
        m_rotation += angularImpulse;
        ////for (int i = 0; i != modelCount; i++) {
          ////if(models[i] != NULL) models[i]->updateRotation(angularImpulse);
        ////}

      }
      if(radioId  != 1)
      {
        models[radioId-2]->updateRotation(angularImpulse);
      }
    }
    for (int i = 0; i != modelCount; i++) {
      if(models[i] != NULL) models[i]->render();
    }
    event->accept();
    update();
  }
}

template <typename M>
void
SceneT<M>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  //int selected = whichRadioButton() - 2;
  //if(modelCount > 0 && event->button() == LeftButton && selected >= 0 && mouseRadioSelected() == 3){
    //paintFaces(event);
    //inPaintingMode = true;
  //}else {

  mouseClicked = true;

  //int index = getClickedMeshIndex(event);
  //if (index >= 0)
    //clickRadioButton(index+2);
  //else
    //clickRadioButton(1);

  QGraphicsScene::mousePressEvent(event);
  if (event->isAccepted())
    return;
  m_mouseEventTime = m_time.elapsed();
  event->accept();
  update();
 // }
}

template <typename M>
void
SceneT<M>::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{

  mouseClicked = false;
  QGraphicsScene::mouseReleaseEvent(event);
  if (event->isAccepted())
    return;
  const int delta = m_time.elapsed() - m_mouseEventTime;
  event->accept();
  update();
}

template <typename M>
void
SceneT<M>::keyReleaseEvent( QKeyEvent* event)
{
  inPaintingMode = false;
}

template <typename M>
void
SceneT<M>::keyPressEvent( QKeyEvent* event)
{
  const int radioId = whichRadioButton();
  if(radioId  == 1)
  {
    switch(event->key())
    {
      case Key_Up:
        m_vertical += TANSLATE_SPEED;
        break;
      case Key_Down:
        m_vertical -= TANSLATE_SPEED;
        break;
      case Key_Right:
        m_horizontal += TANSLATE_SPEED;
        break;
      case Key_Left:
        m_horizontal -= TANSLATE_SPEED;
        break;
    }
  }
  else
  {
    switch(event->key())
    {
      case Key_Alt:
        if(mouseRadioSelected() == 3)
          inPaintingMode = true;
        break;
      case Key_Up:
        models[radioId-2]->updateVertical(-TANSLATE_SPEED);
        break;
      case Key_Down:
        models[radioId-2]->updateVertical(TANSLATE_SPEED);
        break;
      case Key_Right:
        models[radioId-2]->updateHorizontal(TANSLATE_SPEED);
        break;
      case Key_Left:
        models[radioId-2]->updateHorizontal(-TANSLATE_SPEED);
        break;
      case Key_S:
        if (models[0] != NULL && models[1] != NULL) softICP(models[0], models[1]);
    }
  }
  event->accept();
  update();
}

template <typename M>
int
SceneT<M>::getClickedMeshIndex(QGraphicsSceneMouseEvent *event){
  int selected = 0;
  GLuint Nhits = 0;
  while (models[selected] != NULL && Nhits == 0){
  glRenderMode(GL_SELECT);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  gluPickMatrix(event->scenePos().x(), (GLdouble)(viewport[3]-event->scenePos().y()), 0.01, 0.1, viewport);
  gluPerspective(70, width() / height(), 0.0001, 1000);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslatef(m_horizontal, m_vertical, -m_distance);
  glRotatef(m_rotation.x(), 1, 0, 0);
  glRotatef(m_rotation.y(), 0, 1, 0);
  glRotatef(m_rotation.z(), 0, 0, 1);
  
  glEnable(GL_MULTISAMPLE);
  glInitNames();
  glPushName( 0xffffffff );
  if (models[selected] != NULL) models[selected]->render();
  glDisable(GL_MULTISAMPLE);
  
  glPopMatrix();
  
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  
  Nhits = glRenderMode(GL_RENDER);
  selected++;
  }
  if (models[selected-1] != NULL){
    return selected-1;
  }else {
    return -1;
  }

}

template <typename M>
void
SceneT<M>::paintFaces(QGraphicsSceneMouseEvent *event)
{
    int selected = whichRadioButton() - 2;
    glRenderMode(GL_SELECT);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    gluPickMatrix(event->scenePos().x(), (GLdouble)(viewport[3]-event->scenePos().y()), 0.01, 0.1, viewport);
    gluPerspective(70, width() / height(), 0.01, 1000);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(m_horizontal, m_vertical, -m_distance);
    glRotatef(m_rotation.x(), 1, 0, 0);
    glRotatef(m_rotation.y(), 0, 1, 0);
    glRotatef(m_rotation.z(), 0, 0, 1);
    
    glEnable(GL_MULTISAMPLE);
    glInitNames();
    glPushName( 0xffffffff );
    if (models[selected] != NULL) models[selected]->render();

    //std::cout << selected << "q\n";

    glDisable(GL_MULTISAMPLE);
    
    glPopMatrix();
    
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
  
    GLuint Nhits = glRenderMode(GL_RENDER);
    clicked == false;

    //std::cout << "Nhits " << Nhits << "\n";

    //std::cout << Nhits << " hits\n";

    if (Nhits > 0){
      GLuint item;
      GLuint front;
      GLuint min = 0;
      for(size_t i = 0, index = 0; i < Nhits; i++ )
      {
        GLuint nitems = PickBuffer[index++];
        //std::cout << "nitems " << nitems << "\n";
        GLuint zmin = PickBuffer[index++];
        GLuint zmax = PickBuffer[index++];
        //std::cout << "zmin " << zmin << "\n";
        //std::cout << "zmax " << zmax << "\n";
        if(nitems == 1)
        {
          if(min == 0){
            min = zmin;
            item = PickBuffer[index++];
          }
          else if(min > zmin)
          {
            min = zmin;
            item = PickBuffer[index++];
          }
        }
      }
      models[selected]->select(item);
    }
}

template <typename M>
void
SceneT<M>::moveMeshInOneAxis(QGraphicsSceneMouseEvent *event)
{
  int selected = whichRadioButton() - 2;
  glRenderMode(GL_SELECT);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  gluPickMatrix(event->scenePos().x(), (GLdouble)(viewport[3]-event->scenePos().y()), 0.01, 0.1, viewport);
  gluPerspective(70, width() / height(), 0.01, 1000);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslatef(m_horizontal, m_vertical, -m_distance);
  glRotatef(m_rotation.x(), 1, 0, 0);
  glRotatef(m_rotation.y(), 0, 1, 0);
  glRotatef(m_rotation.z(), 0, 0, 1);
  
  glEnable(GL_MULTISAMPLE);
  glInitNames();
  glPushName( 0xffffffff );
  if (models[selected] != NULL) models[selected]->render();
  glDisable(GL_MULTISAMPLE);
  
  glPopMatrix();
  
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  
  GLuint Nhits = glRenderMode(GL_RENDER);
  clicked == false;
  std::cout << Nhits << " hits\n";
  if (Nhits > 0){
    GLuint item;
    GLuint front;
    for(size_t i = 0, index = 0; i < Nhits; i++ )
    {
      GLuint nitems = PickBuffer[index++];
      std::cout << index << " x" << index+1 << " \n";
      index+= 2;
      for(size_t j = 0; j < nitems; j++ )
      {
        item = PickBuffer[index++];
      }
      QPointF delta = event->scenePos() - event->lastScenePos();
      //typedef typename M::Point Point;
      QVector3D modelRotation = models[selected]->meshRotation;
      if (item == 0){
        models[selected]->updateHorizontal(modelRotation.x() * TANSLATE_SPEED);
      }
    }
  }
  
  
}
template <typename M>
void
SceneT<M>::clickRadioButton(int index){
  switch (index) {
    case 1:
      radio1->click();
      break;
    case 2:
      radio2->click();
      break;
    case 3:
      radio3->click();
      break;
    case 4:
      radio4->click();
      break;
    case 5:
      radio5->click();
      break;
    case 6:
      radio6->click();
      break;
    case 7:
      radio7->click();
      break;
    default:
      break;
  }
}

template <typename M>
int
SceneT<M>::whichRadioButton()
{
  if(radio1->isChecked())
    return 1;
  else if(radio2->isChecked())
    return 2;
  else if(radio3->isChecked())
    return 3;
  else if(radio4->isChecked())
    return 4;
  else if(radio5->isChecked())
    return 5;
  else if(radio6->isChecked())
    return 6;
  else if(radio7->isChecked())
    return 7;
  else if(radio8->isChecked())
    return 8;
  else if(radio9->isChecked())
    return 9;
  else if(radio10->isChecked())
    return 10;
  else if(radio11->isChecked())
    return 11;
  else
    return 0;
}

template <typename M>
void
SceneT<M>::removeMesh()
{
  const int radioId = whichRadioButton();
  if(radioId  == 1)
  {

    for(int i = 0; i != modelCount; i++) {
        models[i] = NULL;
    }
    modelCount = 0;
    removeModelButton->setHidden(true);
    groupBox->setHidden(true);
    mouseControlBox->setHidden(true);
    radio2->setHidden(false);
    radio3->setHidden(true);
    radio4->setHidden(true);
    radio5->setHidden(true);
    radio6->setHidden(true);
    radio7->setHidden(true);
    radio8->setHidden(true);
    radio9->setHidden(true);
    radio10->setHidden(true);
    radio11->setHidden(true);
    clearFacesButton->setHidden(true);
    cutButton->setHidden(true);
    deleteButton->setHidden(true);
    pasteButton->setHidden(true);
    copyButton->setHidden(true);
    geoTreeButton->setHidden(true);
    geoTreeSpinBox->setHidden(true);
  }
  else
  {
    models[radioId-2] = NULL;
    removeRadio(radioId);
    bool clear = true;
    for(int i = 0; i != modelCount; i++) {
      if(models[i] != NULL)
        clear = false;
    }

    if(clear){
      removeModelButton->setHidden(true);
      groupBox->setHidden(true);
      mouseControlBox->setHidden(true);
      radio2->setHidden(false);
      clearFacesButton->setHidden(true);
      cutButton->setHidden(true);
      deleteButton->setHidden(true);
      pasteButton->setHidden(true);
      copyButton->setHidden(true);
      geoTreeButton->setHidden(true);
      geoTreeSpinBox->setHidden(true);
      modelCount = 0;
    }
  }

}

template <typename M>
void
SceneT<M>::removeRadio(int radioId)
{
  switch(radioId)
  {
    case 2:
      radio2->setHidden(true);
      break;
    case 3:
      radio3->setHidden(true);
      break;
    case 4:
      radio4->setHidden(true);
      break;
    case 5:
      radio5->setHidden(true);
      break;
    case 6:
      radio6->setHidden(true);
      break;
    case 7:
      radio7->setHidden(true);
      break;
    case 8:
      radio8->setHidden(true);
      break;
    case 9:
      radio9->setHidden(true);
      break;
    case 10:
      radio10->setHidden(true);
      break;
    case 11:
      radio11->setHidden(true);
      break;

  }
  radio1->setChecked(true);

}


template <typename M>
int
SceneT<M>::mouseRadioSelected()
{
  if(translateRadio->isChecked())
    return 1;
  else if(rotateRadio->isChecked())
    return 2;
  else if(paintFacesRadio->isChecked())
    return 3;
  else
    return 0;
}


template <typename M>
void
SceneT<M>::cut()
{
  std::cout << "Cut Pressesed" << "\n";
  const int radioId = whichRadioButton();
  if(radioId != 1)
  {
    models[radioId-2]->cut();
  }
}

template <typename M>
void
SceneT<M>::paste()
{
  std::cout << "Paste Pressesed" << "\n";
}

template <typename M>
void
SceneT<M>::copy()
{
  std::cout << "Copy Pressesed" << "\n";
}

template <typename M>
void
SceneT<M>::clearFaces()
{
  std::cout << "Clear Faces Pressesed" << "\n";
  const int radioId = whichRadioButton();
  if(radioId != 1 && models[radioId-2] != NULL){
    models[radioId-2]->clearColour();
  }
}


template <typename M>
void
SceneT<M>::deleteSection()
{
  std::cout << "Delete Pressesed" << "\n";
  const int radioId = whichRadioButton();
  if(radioId != 1 && models[radioId-2] != NULL){
    models[radioId-2]->deleteSink();
  }
}

template <typename M>
void
SceneT<M>::geoTree()
{
  std::cout << "GeoTree Pressesed" << "\n";
  const int radioId = whichRadioButton();
  if(radioId != 1 && models[radioId-2] != NULL){
    models[radioId-2]->createGeoTree(geoTreeSpinBox->value());
  }

}
#endif
