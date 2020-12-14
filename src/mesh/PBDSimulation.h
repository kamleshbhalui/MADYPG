#ifndef __PBDSIMULATION__H__
#define __PBDSIMULATION__H__

#include "AbstractMeshProvider.h"
#include <Simulation/SimulationModel.h>
#include <Simulation/DistanceFieldCollisionDetection.h>
#include <Simulation/CubicSDFCollisionDetection.h>

class PBDSimulation : public AbstractMeshProvider {
 public:
  struct Settings {
    int simulationMethod = 4;
    int bendingMethod = 2;
    int substeps = 1;
    float timestep = 0.001;
    float density = 0.3; // doesn't do anything?
  } m_settings;

  PBDSimulation(const Settings& settings);

  ~PBDSimulation();

  void update();

  void applyForce(float fx , float fy, float fz);

 private:
  std::shared_ptr<PBD::SimulationModel> m_model;
  PBD::DistanceFieldCollisionDetection m_cd; 
  PBD::CubicSDFCollisionDetection m_cd2; 
  std::vector<std::shared_ptr<PBD::CubicSDFCollisionDetection::Grid>> m_rbsdfs;
  // ...

  bool loadClothMesh(const std::string& filepath);
  bool loadBoxObstacle(const std::string& filepath);
  bool loadSDFObstacle(const std::string& filepath);
};

#endif  // __PBDSIMULATION__H__