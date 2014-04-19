#ifndef MODEL_HH
#define MODEL_HH

#include <iostream>
#include <fstream>
#include <QString>
#include <QVector>
#include <QPainter>
#include <QGraphicsItem>
#include "Scene.hh"
#include "MyMesh.hh"
#include <QVector3D>
#include <eigen3/Eigen/Dense>
#include <nanoflann.hpp>
#include "GeoTreeT.hh"
#include <queue>
#include <unordered_set>
#include "maxflow/graph.h"
//#include <flann/io/hdf5.h>

using namespace Qt;
using namespace OpenMesh;
using namespace Eigen;

template <typename M>
class QtModelT : public QGraphicsItem
{
public:
    typedef M MyMesh;
    typedef std::vector<std::vector< std::pair< size_t, double > > > MapTable;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    QRectF boundingRect() const;
    int modelnumber;
    void calcNormals();
    std::vector<VertexHandle> boundaryPoints;
    PointMatrix boundaryMatrix;
  
public:
    M mesh;
    QtModelT(M& m);
    ~QtModelT();
    void render();
    void select(int faceNumber);
    void addSink(int faceNumber);
    void addSource(int faceNumber);
    void updateColour();
    void updateRotation(QVector3D& rotationVec);
    void updateHorizontal(float x);
    void updateVertical(float x);
    void updateZAxis(float x);
    void applyTransformations();
    M* getMesh(){ return &mesh; }
    PointMatrix buildMatrix();
    PointMatrix buildSampledMatrix();
    PointMatrix buildNormalMatrix();
    int getNoVerticies();
    void updateTransformations(Matrix<double, 3, 3>& R, double x, double y, double z);
    void nearestNeighbours(double radius, MapTable* resultTable);
    void scale(float alpha);
    void clearColour();
    void colourFaceFromVertexIndex(int vertexNumber, Point col);
    void createGeoTree(int k);
    std::vector<int> getStroke();
    void deleteSink();
    QVector3D meshRotation;
    void toggleFuzzy();
    void autoSelect();
    void flipRegions();
    M cut();
    M copy();
    void exportMesh();
private:
    bool facesConnected(int f1, int f2);
    void addToStroke(int f);
    void addToFuzzyRegion(int f);
    QVector3D modelRotation;
    //void sourceSinkWeightThread(Graph<M, double,double,double>  *g, int  *mapping[150000000], const int start, const int end); 
    
    QColor modelColor;
    GLfloat vertical;
    GLfloat horizontal;
    GLfloat depth;
    GLfloat zAxis;
    GeoTreeT<M> *geoTree;
    std::unordered_set<int> fuzzyRegion;
    std::unordered_set<int> sourceRegion;
    std::unordered_set<int> sinkRegion;
    std::vector<int> stroke;
    std::vector<Point> strokeVertices;
    Vec getFaceCentroid(typename M::FaceHandle fh);
    Vec strokeNormal;
    Vec strokeCentroid;
    Vec firstStrokeNorm;
    Vec lastStrokeNorm;
    double strokeOpeningAngle;
    void calcStrokeProxies();
    double cost(int u, int v);
    int dest;
    int source;
    std::vector<int> prev;
    double normalDistance(int vertex);
    double inverseGeodesic(int vertex);
    bool createSourceAndSink();
    void graphCut();
    bool inRegion(int f);
    void regionGrow(int f, std::unordered_set<int>* region, int type);
    double distToSource(int fId);
    double distToSink(int fId);
    double faceDist(int fId1, int fId2);
    double gauss_curvature(VertexHandle _vh);
    const float deg2Rad;
    void findBoundaryVertices();
    std::vector<VertexHandle> findBoundaryRing(VertexHandle point);
    Vec3f center;
    bool showFuzzy;
    //static int SIZE;
};

struct Dist
{
    float dist;
    int id;

    Dist(int i, float d) : id(i), dist(d)
    {
      //
    }

    bool operator<(const struct Dist& other) const
    {
        //Your priority logic goes here
        return dist > other.dist;
    }
};


#if defined(OM_INCLUDE_TEMPLATES) && !defined(MODEL_CC)
#  define SCENE_TEMPLATES
#  include "QtModelT.cc"
#endif

#endif
