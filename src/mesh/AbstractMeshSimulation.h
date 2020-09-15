#ifndef __ABSTRACTMESHSIMULATION__H__
#define __ABSTRACTMESHSIMULATION__H__

#include "Mesh.h"
class AbstractMeshSimulation
{
public:
  AbstractMeshSimulation() {}
  virtual ~AbstractMeshSimulation() {}

  // update/step the simulation/replay/animation
  virtual void update() = 0;

  // mesh topology has changed in the last update
  bool meshIndicesDirty() { return m_indicesDirty; } // for rendering

  // material space has changed (coords or topology) in the last update
  // currently assuming no coord-changes
  bool materialSpaceChanged() { return m_indicesDirty; }

  Mesh &getMesh() { return m_mesh; }
  const Mesh &getMesh() const { return m_mesh; }

protected:
  bool m_indicesDirty = false;
  Mesh m_mesh;
};

#endif // __ABSTRACTMESHSIMULATION__H__