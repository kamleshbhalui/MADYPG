#ifndef __PBDSIMULATION__H__
#define __PBDSIMULATION__H__

#include "AbstractMeshProvider.h"
#include <Simulation/SimulationModel.h>

class PBDSimulation : public AbstractMeshProvider {
 public:
  struct Settings {
    int simulationMethod = 3;
    int bendingMethod = 2;
    int substeps = 4;
    float timestep = 0.001;
    float density = 0.3; // doesn't do anything?
  } m_settings;

  PBDSimulation(const Settings& settings);

  ~PBDSimulation() {}

  void loadMesh();
  void update();

 private:

  std::shared_ptr<PBD::SimulationModel> m_model;
  // ...
};

#endif  // __PBDSIMULATION__H__