#ifndef __YARNMAPPER__H__
#define __YARNMAPPER__H__

#include "../utils/debug_timer.h"
// #include "../mesh/AbstractMeshProvider.h"
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/TimeQuery.h>
#include "shaders/ShellMapShader.h"
#include "shaders/DeformShader.h"

#include <memory>

#include "../mesh/BinSeqAnimation.h"
#include "../mesh/ObjSeqAnimation.h"
#include "../mesh/PBDSimulation.h"
#include "Grid.h"
// #include "ModelV1.h"
#include "Model.h"
#include "Model4D.h"
#include "YarnSoup.h"

class YarnMapper {
 public:
  struct Settings {
    std::string modelfolder = "";
    float min_yarn_length_per_r = 16.0f;
    bool flat_normals       = false;
    bool flat_strains       = false;
    bool shepard_weights    = true;
    float deform_reference  = 1.0f;
    float linearized_bending = 1.0f;
    // bool use_linearized_bending = false;
    bool shell_map          = true;
    bool default_same_tri   = true;
    bool repeat_frame       = false;
    bool gpu_compute        = true;
    float phong_deformation = 0.5f;
    float svdclamp = 0.8;
    enum Provider {
      ObjSeq = 0,
      BinSeq = 1,
      PBD   = 2,
      COUNT
    } provider_type = Provider::BinSeq;
    ObjSeqAnimation::Settings objseq_settings;
    BinSeqAnimation::Settings binseq_settings;
    PBDSimulation::Settings pbd_settings;
  } m_settings;

  YarnMapper() : m_initialized(false) {}

  void initialize() { step(); }

  ~YarnMapper() {}

  void step();
  void deform_reference(const Mesh& mesh, bool flat_strains = false);
  void shell_map(const Mesh& mesh, bool flat_normals = false);

  // const std::vector<uint32_t>& getIndices() const {
  //   return m_soup.getIndices();
  // }

  const std::shared_ptr<AbstractMeshProvider> getMeshSimulation() {
    return m_meshProvider;
  }
  float getRadius() const { return m_model->getPYP().r; }

  bool isInitialized() const { return m_initialized; }

  VectorBuffer<VertexWSData>& getVertexBuffer() { return m_soup.get_Xws(); }
  VectorBuffer<uint32_t>& getIndexBuffer() { return m_soup.getIndexBuffer(); }

  void reloadShaders() {
    m_deformShader = DeformShader();
    m_shellMapShader = ShellMapShader();
  }

  void applyForce(float fx , float fy, float fz) {
    m_meshProvider->applyForce(fx, fy, fz);
  }

  // timing
  Debug::MovingAverageTimer<1000, std::chrono::microseconds> m_timer;
  // Debug::MovingAverageTimer<10, std::chrono::microseconds> m_timer;

  // export the current yarn geometry to an FBX file
  bool export2fbx(const std::string& filename);
  bool export2fbx_cloth(const std::string& filename);

  void dbg_compare_bending();
  void dbg_compare_nsamples();
  bool dbg_compare_keep_uv=false;

  // #define DO_DEBUG_STATS
  struct DebugSettings {
    std::vector<float> strain_toggle = std::vector<float>{1,1,1,1,1,1};
    int hist_nbins = 40;
    float hist_min = -0.8f;
    float hist_max = 2.0f;
    float hist_bend = 100.0f;
    std::vector<std::vector<float>> hist_counts;
    int hist_stepcount = 0;

    bool toggle = false;
  } m_dbg;

 private:
  bool m_initialized;
  // std::unique_ptr<ModelV0> m_model;
  std::unique_ptr<Model> m_model;
  std::unique_ptr<Model4D> m_model4D;
  Grid m_grid;
  YarnSoup m_soup;

  // simulation mesh stuff
  std::shared_ptr<AbstractMeshProvider> m_meshProvider;

  // yarn stuff
  Magnum::GL::Buffer m_buf_tri;
  Magnum::GL::Buffer m_buf_bary;

  // GPU compute
  DeformShader m_deformShader;
  ShellMapShader m_shellMapShader;
  Magnum::GL::TimeQuery m_glq1{Magnum::GL::TimeQuery::Target::TimeElapsed},
              m_glq2{Magnum::GL::TimeQuery::Target::TimeElapsed},
              m_glq3{Magnum::GL::TimeQuery::Target::TimeElapsed};
};

#endif  // __YARNMAPPER__H__