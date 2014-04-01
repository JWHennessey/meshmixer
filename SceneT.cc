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
{

  QWidget *controls = createDialog(tr("Controls"));
  
  groupBox = new QGroupBox(tr("Select Mesh"));
  radio1 = new QRadioButton(tr("All"));
  radio2 = new QRadioButton(tr("M1"));
  radio3 = new QRadioButton(tr("M2"));
  radio4 = new QRadioButton(tr("M3"));
  radio5 = new QRadioButton(tr("M4"));
  radio6 = new QRadioButton(tr("M5"));
  groupBox->setHidden(true);
  groupBox->setFocusPolicy(Qt::StrongFocus);
  radio3->setHidden(true);
  radio4->setHidden(true);
  radio5->setHidden(true);
  radio6->setHidden(true);

  radio1->setChecked(true);
  QVBoxLayout *vbox = new QVBoxLayout;
  vbox->addWidget(radio1);
  vbox->addWidget(radio2);
  vbox->addWidget(radio3);
  vbox->addWidget(radio4);
  vbox->addWidget(radio5);
  vbox->addWidget(radio6);

  vbox->addStretch(1);
  groupBox->setLayout(vbox);
  controls->layout()->addWidget(groupBox);

  m_modelButton = new QPushButton(tr("Load model"));
  controls->layout()->addWidget(m_modelButton);

  removeModelButton = new QPushButton(tr("Remove model"));
  controls->layout()->addWidget(removeModelButton);
  removeModelButton->setHidden(true);
  
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
  /*
  if (clicked){
    glSelectBuffer(65535, PickBuffer);
    glRenderMode(GL_SELECT);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    gluPickMatrix(clickLocation.x(), (GLdouble)(viewport[3]-clickLocation.y()), 10, 10, viewport);
    std::cout << clickLocation.x() << "\n";
  } else {
    glRenderMode(GL_RENDER);
  }
   */
  if(models.size() > 0){
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
      for(typename std::vector<QtModelT<M>*>::size_type i = 0; i != models.size(); i++) {
        models[i]->render();
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

    models.push_back(new QtModelT<M>(m_mymesh));

    switch(models.size())
    {
      case 1:
        removeModelButton->setHidden(false);
        groupBox->setHidden(false);
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

  m_distance *= qPow(1.2, -event->delta() / 120.);
  //std::cout << m_distance << "\n";
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
    const QPointF delta = event->scenePos() - event->lastScenePos();
    QVector3D angularImpulse = QVector3D(delta.y(), delta.x(), 0) * 0.1;
    
    const int radioId = whichRadioButton();
    if(radioId  == 1){
      m_rotation += angularImpulse;
    }
    else
    {
      models[radioId-2]->updateRotation(angularImpulse);
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
  if(models.size() > 0){
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
    for(typename std::vector<QtModelT<M>*>::size_type i = 0; i != models.size(); i++) {
      models[i]->render();
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
        models[radioId-2]->updateVertical(TANSLATE_SPEED);
        break;
      case Key_Down:
        models[radioId-2]->updateVertical(-TANSLATE_SPEED);
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
  else
    return 0;
}




#endif
