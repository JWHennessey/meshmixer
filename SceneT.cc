#ifndef SCENE_CC
#define SCENE_CC

#include <iostream>
#include "Scene.hh"
#include <QVector3D>
#include "GLUT/glut.h"
#include <QtGui>
#include <QtOpenGL>
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
, m_vertical(-0.1f)
, m_horizontal(0.0f)
, TANSLATE_SPEED(0.01f)
, deg2Rad(0.0174532925)
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
  translateRadio->setChecked(true);
  QVBoxLayout *vb = new QVBoxLayout;
  vb->addWidget(translateRadio);
  vb->addWidget(rotateRadio);
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
  GLfloat col1[] = { 0.7,  0.7,  0.8,  1.0};
  GLfloat col2[] = { 0.8,  0.7,  0.7,  1.0};
  GLfloat col3[] = { 1.0,  1.0,  1.0,  1.0};
 
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
    setDefaultMaterial();

    setDefaultLight();  
    glLoadIdentity();

    glEnable(GL_LIGHTING);
    glShadeModel(GL_FLAT);
    //glutSolidTeapot(0.5);
    
    glTranslatef(m_horizontal, m_vertical, -m_distance);
    glRotatef(m_rotation.x(), 1, 0, 0);
    glRotatef(m_rotation.y(), 0, 1, 0);
    glRotatef(m_rotation.z(), 0, 0, 1);

    glEnable(GL_MULTISAMPLE);

    for (int i = 0; i != modelCount; i++) {
      if(models[i] != NULL) models[i]->render();
    }
    glDisable(GL_MULTISAMPLE);

    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
  }
  //painter->endNativePainting();
}

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
  if (event->isAccepted())
    return;
  if (event->buttons() & Qt::LeftButton) {
    QPointF delta = event->scenePos() - event->lastScenePos();
    const int radioId = whichRadioButton();
    QVector3D angularImpulse = QVector3D(delta.y(), delta.x(), 0) * 0.1;
    if(mouseTranslate())
    {
      if(radioId  == 1){
        //std::cout << m_distance << "\n";
        m_vertical -= delta.y() * (TANSLATE_SPEED);
        m_horizontal += delta.x() * (TANSLATE_SPEED);
      }
      else
      {
        typedef typename M::Point Point;
        QVector3D modelRotation = m_rotation;
        modelRotation = modelRotation * deg2Rad;
        Eigen::AngleAxis<float> aax(modelRotation.x(), Eigen::Vector3f(1, 0, 0));
        Eigen::AngleAxis<float> aay(modelRotation.y(), Eigen::Vector3f(0, 1, 0));
        Eigen::AngleAxis<float> aaz(modelRotation.z(), Eigen::Vector3f(0, 0, 1));
        Eigen::Quaternion<float> rotation = aax * aay * aaz;

        Eigen::Vector3f p = Eigen::Vector3f(delta.x(), delta.y(), 0);
        p = rotation * p;
        //delta = R * delta;
        models[radioId-2]->updateHorizontal(p[0] * TANSLATE_SPEED);
        models[radioId-2]->updateVertical(p[1] * TANSLATE_SPEED);
        models[radioId-2]->updateZAxis(p[2] * TANSLATE_SPEED);
        //models[radioId-2]->updateHorizontal(delta.x() * TANSLATE_SPEED);
        //models[radioId-2]->updateVertical(delta.y() * TANSLATE_SPEED);
        //models[radioId-2]->updateZAxis(0 * TANSLATE_SPEED);

      }
    }
    else
    {
      
      if(radioId  == 1){
        m_rotation += angularImpulse;
        //for (int i = 0; i != modelCount; i++) {
          //if(models[i] != NULL) models[i]->updateRotation(angularImpulse);
        //}

      }
      else
      {
        models[radioId-2]->updateRotation(angularImpulse);
      }
    }
    event->accept();
    update();
  }
}

template <typename M>
void
SceneT<M>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  clickLocation = event->scenePos();
  clicked = true;
  QGraphicsScene::mousePressEvent(event);
  if (event->isAccepted())
    return;
  m_mouseEventTime = m_time.elapsed();
  event->accept();
  update();
}

template <typename M>
void
SceneT<M>::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{

  //glInitNames();
  //glPushName(0);
  if(modelCount > 0){
    glRenderMode(GL_SELECT);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    gluPickMatrix(event->scenePos().x(), (GLdouble)(viewport[3]-event->scenePos().y()), 1, 1, viewport);
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
    for (int i = 0; i != modelCount; i++) {
      if (models[i] != NULL) models[i]->render();
    }
    
    glDisable(GL_MULTISAMPLE);
    
    glPopMatrix();
    
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
  }
  GLuint Nhits = glRenderMode(GL_RENDER);
  clicked == false;
  std::cout << Nhits << " y\n";
  if (Nhits > 0){
    GLuint item;
    GLuint front;
    for(size_t i = 0, index = 0; i < Nhits; i++ )
    {
      GLuint nitems = PickBuffer[index++];
      index+= 2;
      for(size_t j = 0; j < nitems; j++ )
      {
        item = PickBuffer[index++];
        std::cout << Nhits << " z" << item << " \n";
      }
      models[0]->select(item);
    }
  } else {
  QGraphicsScene::mouseReleaseEvent(event);
  if (event->isAccepted())
    return;
  const int delta = m_time.elapsed() - m_mouseEventTime;
  event->accept();
  update();
  }
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
    }
  }
  event->accept();
  update();
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

  }
  else
  {
    models[radioId-2] = NULL;
    removeRadio(radioId);
    bool clear = true;
    for(typename std::vector<QtModelT<M>*>::size_type i = 0; i != modelCount; i++) {
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
bool
SceneT<M>::mouseTranslate()
{
  if(translateRadio->isChecked())
    return true;
  else
    return false;
}


template <typename M>
void
SceneT<M>::cut()
{
  std::cout << "Cut Pressesed" << "\n";
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
}


template <typename M>
void
SceneT<M>::deleteSection()
{
  std::cout << "Delete Pressesed" << "\n";
}

#endif
