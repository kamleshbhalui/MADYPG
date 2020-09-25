#ifndef __MODEL__H__
#define __MODEL__H__

#include "PeriodicYarnPattern.h"

// CPP INCLUDES

#include <filesystem>  // requires C++17, g++ >= 8 & linking to stdc++fs)
#include <set>
#include <sstream>
namespace fs = std::filesystem;
#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

void deserialize_matrix(const std::string& filename, MatrixGLf& M);

class Model {
 public:
  Model(const std::string& folder) {
    m_initialized = false;
    if (!fs::is_directory(folder)) {
      Debug::error("Could not load yarn model directory:", folder);
      return;
    }
    const auto pypfile = fs::path(folder) / "pyp";
    if (!fs::exists(pypfile) || fs::is_directory(pypfile)) {
      Debug::error("No 'pyp' file found in folder:", folder);
      return;
    }

    // Load model / pyp
    m_pyp.deserialize(pypfile);  // DEBUG hardcoded file
    m_pyp.rectangulize();  // TODO potentially assume that this is true for new
    // pyp, but it should take only a millisecond anyway
    m_pyp.compute_parametric();  // TODO potentially cache // actually ideally
                                 // already store with python model creation

    // TODO NOTE: most of this stuff is bad code/design and needs to be revamped
    // once the model format is more clear.
    // NOTE currently (unnecessarily) using
    //      vix -> y,t for deform(vix)
    //      pyarnixs[y],pyarn[T] for model creation

    // make Tex2D per yarn per strain combo
    int i0 = 0, i1 = 2;  // DEBUG fixed sx sy
    {
      std::set<scalar> S0_set, S1_set;
      auto subfolder = fs::path(folder) / (strain_names[i0] + strain_names[i1]);
      for (auto& p : fs::directory_iterator(subfolder)) {
        // split "GT_s0_s1"
        std::stringstream ss(p.path().filename());
        std::string dummy;
        scalar s0, s1;
        std::getline(ss, dummy, '_');
        std::getline(ss, dummy, '_');
        s0 = std::stod(dummy);
        std::getline(ss, dummy, '_');
        s1 = std::stod(dummy);
        S0_set.insert(s0);
        S1_set.insert(s1);
      }
      std::vector<scalar> S0, S1;
      std::copy(S0_set.begin(), S0_set.end(), std::back_inserter(S0));
      std::sort(S0.begin(), S0.end());
      std::copy(S1_set.begin(), S1_set.end(), std::back_inserter(S1));
      std::sort(S1.begin(), S1.end());

      // preload all the GT data into memory
      int n0 = int(S0.size()), n1 = int(S1.size());
      std::vector<MatrixGLf> GT_preload(n0 * n1);
      static constexpr int C = 8;
      char buf0[C], buf1[C];
      for (size_t i = 0; i < S0.size(); i++) {
        std::snprintf(buf0, C, "%.3f", double(S0[i]));
        for (size_t j = 0; j < S1.size(); j++) {
          std::snprintf(buf1, C, "%.3f", double(S1[j]));
          GT_preload[i * n1 + j].resize(m_pyp.Q.rows(), 4);

          deserialize_matrix(
              subfolder / (std::string("GT_") + buf0 + std::string("_") + buf1),
              GT_preload[i * n1 + j]);
        }
      }

      m_tex2Ds_sxsy.resize(m_pyp.param_y2v.size());
      threadutils::parallel_for(
          size_t(0), m_pyp.param_y2v.size(), [&](size_t y) {
            auto& tex       = m_tex2Ds_sxsy[y];
            const auto& y2v = m_pyp.param_y2v[y];
            // store per-yarn T reference
            tex.T = &m_pyp.param_y2t[y];  // TODO: ptr
            // reserve matrix 2d-grid
            tex.X_arr.resize(n0 * n1);

            for (size_t i = 0; i < S0.size(); i++) {
              for (size_t j = 0; j < S1.size(); j++) {
                size_t k = i * n1 + j;
                tex.X_arr[k].resize(y2v.size(), 4);  // alloc matrix
                // fill matrix
                // NOTE: TODO assuming that we are loading
                // "GT" = [Sinv Rt utilde, dtheta]
                for (size_t h = 0; h < y2v.size(); h++)
                  tex.X_arr[k].row(h)
                      << m_pyp.Q.row(y2v[h]) + GT_preload[k].row(y2v[h]);
              }
            }

            tex.S0 = S0;  // copy..
            tex.S1 = S1;
            tex.init();
          });
    }
    m_initialized = true;
  }

  const std::tuple<Vector4s, scalar, scalar> deformation(const Vector6s& strain,
                                                         int vix) {
    Vector4s g = Vector4s::Zero();

    // TODO not pass vix?
    int y         = m_pyp.param_v2y[vix];
    scalar t      = m_pyp.param_v2t[vix];
    Vector4s xref = m_pyp.Q.row(vix);

    int i0 = 0, i1 = 2;
    {
      Vector4s x;
      scalar dbg0, dbg1;
      std::tie(x, dbg0, dbg1) =
          m_tex2Ds_sxsy[y].blerp(strain[i0], strain[i1], t);
      g += x - xref;
      return std::make_tuple(g, dbg0, dbg1);
    }
  }

  const PeriodicYarnPattern& getPYP() const { return m_pyp; }
  bool isInitialized() const { return m_initialized; }

 private:
  bool m_initialized;
  PeriodicYarnPattern m_pyp;
  static const std::vector<std::string> strain_names;

  // Tex2D[y]: S0, S1, T[y], X: ij-> GT[ixs[y]]
  //  precomp: invlen0, invlen1
  //  blerp: get cell, in cell: Gxx = pwise(T,Xij,t), use
  struct Tex2D {
    void init() {
      inv0.resize(int(S0.size()) - 1);
      inv1.resize(int(S1.size()) - 1);
      for (size_t i = 0; i < inv0.size(); i++)
        inv0[i] = 1 / (S0[i + 1] - S0[i]);
      for (size_t i = 0; i < inv1.size(); i++)
        inv1[i] = 1 / (S1[i + 1] - S1[i]);
    }

    Vector4s pwiselerp(const std::vector<scalar>& T, const MatrixGLf& X,
                       scalar t) {
      int i = std::distance(T.begin(),
                            std::upper_bound(T.begin() + 1, T.end() - 1, t)) -
              1;
      scalar a = (t - T[i]) / (T[i + 1] - T[i]);  // TODO precompute inv
      a        = std::min(std::max(scalar(0), a), scalar(1));
      // TODO periodic ! if i == last, interpolate between last and first, and
      // also use t % RL or something
      return (1 - a) * X.row(i) + a * X.row(i + 1);
    }

    std::tuple<Vector4s, scalar, scalar> blerp(scalar s0, scalar s1, scalar t) {
      // s0=0;s1=0;

      // find closest/containing cell
      int i0 = std::distance(S0.begin(),
                             std::upper_bound(S0.begin() + 1, S0.end() - 1,
                                              s0)) -
               1;  // first cell (i,i+1) s.t. s0 < S0[i+1] else last = size-2
      int i1 = std::distance(S1.begin(), std::upper_bound(S1.begin() + 1,
                                                          S1.end() - 1, s1)) -
               1;

      // Debug::log(i0,S0.size(),i1,S1.size());

      // coeff
      scalar a = (s0 - S0[i0]) * inv0[i0];
      scalar b = (s1 - S1[i1]) * inv1[i1];

      // Debug::log("   ",S0[i0],s0,S0[i0+1],"_____",(1-a)*S0[i0]+a*S0[i0+1]);

      // clamp extrapolation
      a = std::min(std::max(scalar(0), a), scalar(1));
      b = std::min(std::max(scalar(0), b), scalar(1));

      // bilinear
      Vector4s G00 = pwiselerp(*T, X_arr[i0 * S1.size() + i1], t);
      Vector4s G01 = pwiselerp(*T, X_arr[i0 * S1.size() + i1 + 1], t);
      Vector4s G10 = pwiselerp(*T, X_arr[(i0 + 1) * S1.size() + i1], t);
      Vector4s G11 = pwiselerp(*T, X_arr[(i0 + 1) * S1.size() + i1 + 1], t);

      Vector4s g = (1 - a) * (1 - b) * G00 + a * (1 - b) * G10 +
                   (1 - a) * b * G01 + a * b * G11;

      return std::make_tuple(g, a, b);
    }

    std::vector<scalar> inv0, inv1;
    std::vector<scalar> S0, S1;
    std::vector<MatrixGLf> X_arr;
    std::vector<scalar>* T;  // DEBUG as pointer
  };
  std::vector<Tex2D> m_tex2Ds_sxsy;

  /*
  def robust_eigenstuff(II):
    # input: II = IIxx, IIxy, IIyy
    A = 0.5 * (II[0] + II[2])
    B = 0.5 * (II[0] - II[2])
    eps = 1e-8
    S = np.sqrt(B**2 + II[1]**2 + eps)
    k = -1 if (II[0] - II[2]) < 0 else 1
    lam1 = A + S
    lam2 = A - S
    c2 = 0.5 + k * (0.5 - II[1]**2 / ((B + k * S)**2 + II[1]**2))
    return lam1, lam2, c2
  */
};

#endif  // __MODEL__H__