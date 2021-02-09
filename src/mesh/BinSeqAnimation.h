#ifndef __BINSEQANIMATION__H__
#define __BINSEQANIMATION__H__

#include "AbstractMeshProvider.h"

// ...
#include "../io/framesio.h"
#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

// Binary version of an obj-sequence animation, using pre-baked binary 'Frame's
// for fast loading
class BinSeqAnimation : public AbstractMeshProvider {
 public:
  struct Settings {
    std::string filepath;
    bool repeat = true;
    float scale = 1.0;
  } m_settings;

  BinSeqAnimation(const Settings& settings) : m_settings(settings) {
    m_indicesDirty = true;
    m_iter         = 0;

    if (!deserialize_frames(m_settings.filepath, m_frames)) {
      return;
    }

    bool rescale = std::abs(m_settings.scale - 1.0f) > 0.001f;
    if (rescale) {
      for (auto& frame : m_frames) {
        frame.cloth_U *= m_settings.scale;
        frame.cloth_V *= m_settings.scale;
        for (size_t i = 0; i < frame.obs_V.size(); ++i)
          frame.obs_V[i] *= m_settings.scale;
      }
    }

    update();  // load first frame

    // initial flags
    m_indicesDirty = true;
  }

  ~BinSeqAnimation() {}

  // load next frame: update cloth mesh / obstacle transformation
  void update() {
    if (m_frames.empty()) {
      m_indicesDirty = false;
      return;
    }

    // repeat
    if (m_iter >= m_frames.size() && m_settings.repeat) {
      m_iter = 0;
    }

    // load next
    if (m_iter < m_frames.size()) {
      const auto& frame = m_frames[m_iter];

      // update cloth
      // m_mesh.X.cpu() = frame.cloth_V;
      copy_to(frame.cloth_V, m_mesh.X.cpu());
      if (frame.cloth_U.size() > 0 || frame.cloth_F.size() > 0 ||
          frame.cloth_Fms.size() > 0) {
        copy_to(frame.cloth_U, m_mesh.U.cpu());
        copy_to(frame.cloth_F, m_mesh.F.cpu());
        copy_to(frame.cloth_Fms, m_mesh.Fms.cpu());
        m_indicesDirty = true;
      } else {
        m_indicesDirty = false;
      }

      // update obstacles
      int n_obstacles =
          std::max(std::max(frame.obs_V.size(), frame.obs_F.size()),
                   frame.obs_trafo.size());
      if (int(m_obstacles.size()) < n_obstacles)
        m_obstacles.resize(n_obstacles);
      for (size_t i = 0; i < frame.obs_V.size(); ++i)
        copy_to(frame.obs_V[i], m_obstacles[i].mesh.X.cpu());
      for (size_t i = 0; i < frame.obs_F.size(); ++i)
        copy_to(frame.obs_F[i], m_obstacles[i].mesh.F.cpu());
      for (size_t i = 0; i < frame.obs_trafo.size(); ++i)
        m_obstacles[i].transformation = frame.obs_trafo[i];

      ++m_iter;
    }
  }

 private:
  std::vector<Frame> m_frames;
  size_t m_iter;

  template <typename MatrixT, typename T>
  void copy_to(const MatrixT& A, std::vector<T>& B) {
    B.clear();
    B.reserve(A.rows());
    for (size_t i = 0; i < size_t(A.rows()); i++) {
      B.emplace_back();
      B.back().map() = A.row(i);
    }

    // parallel: it appears that this is (negligibly but still)
    // slower than serial copy! let this be a lesson to default to non-parallel
    // implementation for trivial operations and when parallelizing always
    // measure!
    // B.resize(A.rows());
    // tbb::parallel_for(size_t(0),size_t(A.rows()),[&](size_t i){
    //   B[i].map() = A.row(i);
    // });
  }
};

#endif  // __BINSEQANIMATION__H__