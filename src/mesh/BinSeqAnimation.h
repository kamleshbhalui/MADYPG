#ifndef __BINSEQANIMATION__H__
#define __BINSEQANIMATION__H__

#include "AbstractMeshProvider.h"

// TODO CPP includes

#include "../io/framesio.h"
#include "../utils/debug_logging.h"

class BinSeqAnimation : public AbstractMeshProvider {
 public:
  struct Settings {
    std::string filepath;
    bool repeat = true;
  } m_settings;

  BinSeqAnimation(const Settings& settings) : m_settings(settings) {
    m_indicesDirty = true;
    m_iter         = 0;

    if (!deserialize_frames(m_settings.filepath, m_frames)) {
      return;
    }

    update();  // load first frame

    // initial flags
    m_indicesDirty = true;
  }

  ~BinSeqAnimation() {}

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
      // NOTE optionally could do a parallel copy, but likely not a bottleneck
      m_mesh.X = frame.cloth_V;
      if (frame.cloth_U.size() > 0 || frame.cloth_F.size() > 0 ||
          frame.cloth_Fms.size() > 0) {
        m_mesh.U       = frame.cloth_U;
        m_mesh.F       = frame.cloth_F;
        m_mesh.Fms     = frame.cloth_Fms;
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
        m_obstacles[i].mesh.X = frame.obs_V[i];
      for (size_t i = 0; i < frame.obs_F.size(); ++i)
        m_obstacles[i].mesh.F = frame.obs_F[i];
      for (size_t i = 0; i < frame.obs_trafo.size(); ++i)
        m_obstacles[i].transformation = frame.obs_trafo[i];

      ++m_iter;
    }
  }

 private:
  std::vector<Frame> m_frames;
  size_t m_iter;
};

#endif  // __BINSEQANIMATION__H__