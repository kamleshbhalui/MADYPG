#ifndef __TESTMESHSIMULATION__H__
#define __TESTMESHSIMULATION__H__

#include "AbstractMeshProvider.h"


#include "meshio.h"

class TestMeshSimulation : public AbstractMeshProvider
{
public:
  TestMeshSimulation() {
    m_indicesDirty = true;

    // std::string fname = "presim/cube10cm.obj";
    std::string fname = "presim/suzanne_t.obj";
    // std::string fname = "presim/0100_00.obj";
    m_mesh = load_obj_mesh(fname);

    // m_mesh.X.resize(4, 3);
    // m_mesh.X << 0.0f, 0.0f, 0.0f,
    //     1.0f, 0.0f, 0.0f,
    //     1.0f, 1.0f, 0.0f,
    //     0.0f, 1.0f, 0.0f;
    // m_mesh.X *= 0.1;
    // m_mesh.U = m_mesh.X;
    // m_mesh.F.resize(2, 3);
    // m_mesh.F << 0, 1, 2,
    //     0, 2, 3;
  }
  ~TestMeshSimulation() {}

  // update/step the simulation/replay/animation
  virtual void update() {
    // m_mesh.X *= 1.01f;

    m_indicesDirty = false;
  }
};

#endif // __TESTMESHSIMULATION__H__