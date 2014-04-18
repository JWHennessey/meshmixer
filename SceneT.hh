#ifndef SCENET_HH
#define SCENET_HH

#include <QGraphicsScene>
#include "GLUT/glut.h"
#include <QLabel>
#include <QTime>
//#include <QtOpenGL>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "QtModelT.hh"
#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Utils/GenProg.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
#include <nanoflann.hpp>
#include <eigen3/Eigen/Dense>

using namespace Qt;
using namespace OpenMesh;
using namespace Eigen;

template <typename M>
class SceneT : public QGraphicsScene
{
public:
  typedef M MyMesh;
  QtModelT<M>* models [10];
public:
    SceneT();
    void drawBackground(QPainter *painter, const QRectF &rect);
    void drawForeground(QPainter *painter, const QRectF &rect);

protected:
    void loadMesh(const QString filePath);

protected:
    QWidget *m_modelButton;
    QWidget *removeModelButton;
    QWidget *m_ex1Button;
    QWidget *m_ex2Button;
    QWidget *m_ex3Button;
    QWidget *clearFacesButton;
    QWidget *autoSelectButton;
    QWidget *toggleFuzzyButton;
    QWidget *cutButton;
    QWidget *deleteButton;
    QWidget *pasteButton;
    QWidget *copyButton;
    QWidget *flipRegionsButton;
    QWidget *geoTreeButton;
    QSpinBox *geoTreeSpinBox;
    void wheelEvent(QGraphicsSceneWheelEvent * wheelEvent);
    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
    void keyPressEvent( QKeyEvent* event);
    void keyReleaseEvent( QKeyEvent* event);
    GLuint PickBuffer[65535];
    bool clicked;
    QPointF clickLocation;
    void removeMesh();
    void clearFaces();
    void paste();
    void cut();
    void copy();
    void deleteSection();
    void geoTree();
    void softICP(QtModelT<M>* m1, QtModelT<M>* m2);
    std::vector<size_t> computeSnapRegion(QtModelT<M>* m1, QtModelT<M>* m2);
    void clickRadioButton(int index);
    int getClickedMeshIndex(QGraphicsSceneMouseEvent *event);
    void toggleFuzzy();
    void autoSelect();
    void flipRegions();
    void addMesh(MyMesh m_mymesh);

private:
    int modelCount;
    QDialog *createDialog(const QString &windowTitle) const;
    int whichRadioButton();
    void setDefaultMaterial();
    void setDefaultLight();
    void updateGTDistances();
    void removeRadio(int radioId);
    int mouseRadioSelected();
    MyMesh m_mymesh;
    QColor m_backgroundColor;
    OpenMesh::IO::Options _options;
    bool mouseClicked;
    QTime m_time;
    int m_mouseEventTime;

    float m_distance;
    float m_vertical;
    float m_horizontal;
    float radius;
    QVector3D m_rotation;

    QGraphicsRectItem *m_lightItem;
    const float TANSLATE_SPEED;
    const float deg2Rad;

    QWidget *meshes;
    QGroupBox* groupBox;
    QRadioButton* radio1;
    QRadioButton* radio2;
    QRadioButton* radio3;
    QRadioButton* radio4;
    QRadioButton* radio5;
    QRadioButton* radio6;
    QRadioButton* radio7;
    QRadioButton* radio8;
    QRadioButton* radio9;
    QRadioButton* radio10;
    QRadioButton* radio11;
    QGroupBox* mouseControlBox;
    QRadioButton* translateRadio;
    QRadioButton* rotateRadio;
    QRadioButton* paintStrokeRadio;
    QRadioButton* paintSinkRadio;
    QRadioButton* paintSourceRadio;
    void paintFaces(QGraphicsSceneMouseEvent *event);
    bool inPaintingMode;
    bool inRotatingMode;
    void moveMeshInOneAxis(QGraphicsSceneMouseEvent *event);
  //void softICP(QtModelT<M>* m1, QtModelT<M>* m2);
    std::vector<size_t> computeSnapRegion(QtModelT<M>* m1, float snapMax);
    void calcBaryCenteredPoints(PointMatrix &matHat, const PointMatrix &mat);
    void generateAMatrix(Matrix<double, 3, 3>  &A, const PointMatrix &m1Hat, const PointMatrix &m2Hat);
    std::vector<size_t> findLocalNeighbourhood(std::vector<double> query_pt, PointMatrix m1Matrix, double radius);
};

#if defined(OM_INCLUDE_TEMPLATES) && !defined(SCENE_CC)
#  define SCENE_TEMPLATES
#  include "SceneT.cc"
#endif

#endif
