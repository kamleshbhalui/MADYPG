#ifndef __YARNMAPPER__H__
#define __YARNMAPPER__H__

#include "../utils/debug_timer.h"
// #include "../mesh/AbstractMeshProvider.h"
#include <memory>

#include "../mesh/ObjSeqAnimation.h"
#include "Grid.h"
#include "Model.h"
#include "YarnSoup.h"

class YarnMapper {
 public:
  struct Settings {
    std::string modelfolder = "";
    bool flat_normals       = false;
    bool flat_strains       = false;
    bool shepard_weights    = true;
    bool deform_reference   = true;
    bool shell_map          = true;
    bool default_same_tri   = false;
    bool repeat_frame       = false;
    ObjSeqAnimation::Settings objseq_settings;
    // TODO XPBDMeshProvider::Settings xpbd_settings;
  } m_settings;

  YarnMapper() : m_initialized(false) {}

  void initialize() { step(); }

  ~YarnMapper() {}

  void step();
  void deform_reference(const Mesh& mesh, bool flat_strains = false);
  void shell_map(const Mesh& mesh, bool flat_normals = false);

  const std::vector<uint32_t>& getIndices() const {
    return m_soup.getIndices();
  }
  const MatrixGLf& getVertexData() const { return m_soup.get_Xws(); }
  const std::shared_ptr<AbstractMeshProvider> getMeshSimulation() {
    return m_meshProvider;
  }
  float getRadius() const { return m_model->getPYP().r; }

  bool isInitialized() const { return m_initialized; }

  // timing
  Debug::MovingAverageTimer<10, std::chrono::microseconds> m_timer;

 private:
  bool m_initialized;
  std::unique_ptr<Model> m_model;
  Grid m_grid;
  YarnSoup m_soup;

  // simulation mesh stuff
  std::shared_ptr<AbstractMeshProvider> m_meshProvider;

  // yarn stuff
  std::vector<uint32_t> I;
};

#endif  // __YARNMAPPER__H__