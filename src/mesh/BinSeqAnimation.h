#ifndef __BINSEQANIMATION__H__
#define __BINSEQANIMATION__H__

#include "AbstractMeshProvider.h"

// TODO CPP includes

#include "../io/framesio.h"
#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

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

  // template <typename T, typename _scalar, int _nrows, int _ncols, int _options>
  // void copy_matrix_to_buffer(VectorBuffer<T> &buf, _scalar TODOMAT) {
    
  //   buf.cpu().resize(M.rows());

  //   // FUCK ITS ALL GOING DOWNHILL BAHHHHHHHH
  // }

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


  template<typename MatrixT, typename T>
  void copy_to(const MatrixT& A, std::vector<T> &B) {
    B.clear();
    B.reserve(A.rows());
    for (size_t i = 0; i < size_t(A.rows()); i++) {
      B.emplace_back();
      B.back().map() = A.row(i);
    }

    // parallel: it appears that this is (kind of negligibly but observably) slower than serial copy! let this be a lesson to default to non-parallel stuff for trivial operations and when parallelizing always measure! 
    // B.resize(A.rows());
    // tbb::parallel_for(size_t(0),size_t(A.rows()),[&](size_t i){
    //   B[i].map() = A.row(i);
    // });
  }
};

#endif  // __BINSEQANIMATION__H__