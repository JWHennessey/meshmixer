#ifndef SCENE_HH
#define SCENE_HH

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <SceneT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>
#include "MyMesh.hh"

using namespace Qt;
using namespace OpenMesh;


class Scene : public SceneT<MyMesh>
{
  Q_OBJECT
public:

    /// default constructor
    Scene() : SceneT<MyMesh>()
    {
      connect(m_modelButton, SIGNAL(clicked()), this, SLOT(selectMesh()));
      connect(removeModelButton, SIGNAL(clicked()), this, SLOT(removeMeshSlot()));
      connect(clearFacesButton, SIGNAL(clicked()), this, SLOT(clearFacesSlot()));
      connect(cutButton, SIGNAL(clicked()), this, SLOT(cutSlot()));
      connect(pasteButton, SIGNAL(clicked()), this, SLOT(pasteSlot()));
      connect(copyButton, SIGNAL(clicked()), this, SLOT(copySlot()));
      connect(deleteButton, SIGNAL(clicked()), this, SLOT(deleteSlot()));

    }
public slots:
    void selectMesh()
    {
      QString selfilter = tr("Meshes (*.stl *.obj)");
      loadMesh(QFileDialog::getOpenFileName(0, tr("Choose mesh"), QString(), tr("All files (*.*);;Meshes (*.stl *.obj)" ), &selfilter));
    }
    void removeMeshSlot(){
      removeMesh();
    }
    void clearFacesSlot(){
      clearFaces();
    }
    void cutSlot(){
      cut();
    }
    void pasteSlot(){
      paste();
    }
    void copySlot(){
      copy();
    }
    void deleteSlot(){
      deleteSection();
    }

};
#endif
