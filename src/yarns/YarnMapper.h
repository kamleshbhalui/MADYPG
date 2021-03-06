#ifndef __YARNMAPPER__H__
#define __YARNMAPPER__H__

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/TimeQuery.h>

#include <memory>

#include "../mesh/BinSeqAnimation.h"
#include "../mesh/ObjSeqAnimation.h"
#include "../mesh/PBDSimulation.h"
#include "../utils/debug_timer.h"
#include "Grid.h"
#include "Model.h"
#include "Model4D.h"
#include "YarnSoup.h"
#include "shaders/DeformShader.h"
#include "shaders/ShellMapShader.h"

class YarnMapper {
 public:
  struct Settings {
    std::string modelfolder     = "";
    float min_yarn_length_per_r = 16.0f;
    bool shepard_weights        = true;
    float deform_reference      = 1.0f;
    float linearized_bending    = 1.0f;
    bool shell_map              = true;
    bool repeat_frame           = false;
    bool gpu_compute            = true;
    float phong_deformation     = 0.5f;
    float svdclamp              = 0.8;
    enum Provider {
      ObjSeq = 0,
      BinSeq = 1,
      PBD    = 2,
      COUNT
    } provider_type = Provider::BinSeq;
    ObjSeqAnimation::Settings objseq_settings;
    BinSeqAnimation::Settings binseq_settings;
    PBDSimulation::Settings pbd_settings;
  } m_settings;

  YarnMapper();
  ~YarnMapper() {}

  // initialize everything with a first step (get first mesh frame etc.)
  void initialize() { step(); }

  // update the mesh animation, and animate yarns accordingly
  // local mechanics-aware deformation and embedded mapping
  void step();

  // cpu implementation of local deformation
  void deform_reference(const Mesh& mesh);
  // cpu implementation of embedded deformation
  void shell_map(const Mesh& mesh);

  const std::shared_ptr<AbstractMeshProvider> getMeshSimulation() {
    return m_meshProvider;
  }
  float getRadius() const { return m_model->getPYP().r; }
  bool isInitialized() const { return m_initialized; }
  VectorBuffer<VertexWSData>& getVertexBuffer() { return m_soup.get_Xws(); }
  VectorBuffer<uint32_t>& getIndexBuffer() { return m_soup.getIndexBuffer(); }
  int getNumVertices() const { return m_soup.numVertices(); }

  void reloadShaders() {
    m_deformShader   = DeformShader();
    m_shellMapShader = ShellMapShader();
  }

  // interface to apply force to PBD sim, ignored for other mesh animations
  void applyForce(float fx, float fy, float fz) {
    m_meshProvider->applyForce(fx, fy, fz);
  }

  // timing
  // Debug::MovingAverageTimer<1000, std::chrono::microseconds> m_timer;
  Debug::MovingAverageTimer<30, std::chrono::microseconds> m_timer;

  // export the current yarn centerlines to an NPY file
  bool export2npy(const std::string& filename_X, const std::string& filename_I, bool xyz_only=false);
  // export the current yarn geometry to an FBX file
  bool export2fbx(const std::string& filename);
  // export the current cloth mesh to an FBX file
  bool export2fbx_cloth(const std::string& filename);

  // debugging methods used in generation error/comparison figures
  void dbg_compare_bending();
  void dbg_compare_nsamples();
  bool dbg_compare_keep_uv = false;  // for overwriting uv for error coloring

 private:
  bool m_initialized;
  std::unique_ptr<Model> m_model;      // yarn displacement data
  std::unique_ptr<Model4D> m_model4D;  // cpu-only 4D-bending model
  Grid m_grid;                         // background UV grid
  YarnSoup m_soup;                     // tiled yarn vertex soup
  std::shared_ptr<AbstractMeshProvider> m_meshProvider;

  // GPU compute
  DeformShader m_deformShader;
  ShellMapShader m_shellMapShader;
  Magnum::GL::TimeQuery m_glq0{Magnum::GL::TimeQuery::Target::TimeElapsed},
      m_glq1{Magnum::GL::TimeQuery::Target::TimeElapsed},
      m_glq2{Magnum::GL::TimeQuery::Target::TimeElapsed},
      m_glq3{Magnum::GL::TimeQuery::Target::TimeElapsed};
};

#endif  // __YARNMAPPER__H__