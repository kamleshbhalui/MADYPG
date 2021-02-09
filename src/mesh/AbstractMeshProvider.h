#ifndef __ABSTRACTMESHPROVIDER__H__
#define __ABSTRACTMESHPROVIDER__H__

#include "Mesh.h"
#include "Trafo.h"

// abstraction of mesh animation, to provide updating meshes to the YarnMapper
// specialized for obj-file sequences and PBD-simulation
class AbstractMeshProvider {
 public:
  struct Obstacle {
    Trafo transformation;
    Mesh mesh;
  };

  AbstractMeshProvider() {}
  virtual ~AbstractMeshProvider() {}

  // update/step the simulation/replay/animation
  virtual void update() = 0;

  // mesh topology has changed in the last update
  bool meshIndicesDirty() { return m_indicesDirty; }  // for rendering

  // material space has changed (coords or topology) in the last update
  // currently assuming no coord-changes
  bool materialSpaceChanged() { return m_indicesDirty; }

  Mesh &getMesh() { return m_mesh; }
  const Mesh &getMesh() const { return m_mesh; }
  const std::vector<Obstacle> &getObstacles() const { return m_obstacles; }

  // only implemented for pbd
  virtual void applyForce(float fx, float fy, float fz) {}

 protected:
  std::vector<Obstacle> m_obstacles;
  bool m_indicesDirty = false;
  Mesh m_mesh;
};

#endif  // __ABSTRACTMESHPROVIDER__H__