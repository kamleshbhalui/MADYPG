#ifndef __MODEL4D__H__
#define __MODEL4D__H__

#include "PeriodicYarnPattern.h"
#include "YarnSoup.h" // vertexmsdata

// CPP INCLUDES

#include <cnpy/cnpy.h>
#include <filesystem>  // requires C++17, g++ >= 8 & linking to stdc++fs)
// #include <set>
#include <fstream>
#include <sstream>
#include <regex>
namespace fs = std::filesystem;
#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

// compile time integer pow https://stackoverflow.com/a/1506856
template<int X, int P>
struct Pow { enum { result = X*Pow<X,P-1>::result }; };
template<int X>
struct Pow<X,0> { enum { result = 1 }; };
template<int X>
struct Pow<X,1> { enum { result = X }; }; 

class Model4D {

 public:
  Model4D(const std::string& folder) {
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
    m_pyp.recompute_VE_table();
    m_pyp.rectangulize();

    // for bend4Dx, bend4Dy: load axes info and data 
    {
      m_tex_SIIx.getAxes().resize(4);
      if (!load_axs(fs::path(folder) / "bend4Dx" / "axes.txt", m_tex_SIIx.getAxes()))
        return;
      m_tex_SIIx.initAxes();
      
      cnpy::NpyArray npyarr =
          cnpy::npy_load(fs::path(folder) / "bend4Dx" / "data.npy");
      auto& data = m_tex_SIIx.getData();
      data.resize(npyarr.shape[0]);
      float* npydata = npyarr.data<float>();
      threadutils::parallel_for(size_t(0), npyarr.shape[0], [&](size_t i) {
        auto g_i = data[i].map();
        for (size_t j = 0; j < 4; j++) {
          g_i[j] = npydata[4 * i + j];
        }
      });
    }

    {
      m_tex_SIIy.getAxes().resize(4);
      if (!load_axs(fs::path(folder) / "bend4Dy" / "axes.txt", m_tex_SIIy.getAxes()))
        return;
      m_tex_SIIy.initAxes();
      
      cnpy::NpyArray npyarr =
          cnpy::npy_load(fs::path(folder) / "bend4Dy" / "data.npy");
      auto& data = m_tex_SIIy.getData();
      data.resize(npyarr.shape[0]);
      float* npydata = npyarr.data<float>();
      threadutils::parallel_for(size_t(0), npyarr.shape[0], [&](size_t i) {
        auto g_i = data[i].map();
        for (size_t j = 0; j < 4; j++) {
          g_i[j] = npydata[4 * i + j];
        }
      });
    }

 
    m_initialized = true;
  }


  const std::tuple<Vector4s, float, float> deformation(const Vector6s& strain,
                                                         int pix) {
    Vector4s g = Vector4s::Zero();

    float l1, l2, c2;
    std::tie(l1,l2,c2) = robust_eigenstuff(strain[3],strain[4],strain[5]);

    // g += c2 g_x(S, l1) + (1-c2) g_x(S, l2)
    // g += c2 g_y(S, l2) + (1-c2) g_y(S, l1)
    // g += g_x(S,0)


    g += c2 * (m_tex_SIIx.sample({strain[0],strain[1],strain[2],l1}, pix) + m_tex_SIIy.sample({strain[0],strain[1],strain[2],l2}, pix));
    g += (1 - c2) * (m_tex_SIIx.sample({strain[0],strain[1],strain[2],l2}, pix) + m_tex_SIIy.sample({strain[0],strain[1],strain[2],l1}, pix));

    // TODO INCORRECT GROUND TRUTH MODEL: counts any sxsasy defo twice..
    g -= m_tex_SIIx.sample({strain[0],strain[1],strain[2],0}, pix); // hijacking one of the texs for subtracting the doubly counted pure inplane defo...

    return std::make_tuple(g, 0, 0);
  }

  const PeriodicYarnPattern& getPYP() const { return m_pyp; }
  bool isInitialized() const { return m_initialized; }

 private:
  bool m_initialized;
  PeriodicYarnPattern m_pyp;

  bool load_axs(const std::string& filepath, std::vector<std::vector<float>>& axs) {
    std::fstream ifs(filepath.c_str(), std::ios::in);
    if (!ifs) {
      // std::cerr << "Error: failed to open file " << filepath << "\n";
      return false;
    }

    size_t axi = 0;
    std::string line;
    float val;
    while (!ifs.eof() && axi < axs.size()) {
      assert(axi < axs.size());
      auto& ax = axs[axi];
      ax.clear();
      std::getline(ifs, line);
      std::stringstream ss(line);
      while ((ss >> val)) {
        ax.push_back(val);
      }
      ++axi;
    }

    return true;
  }
  
  struct DeformationEntry {
    float x, y, z, th;
    Eigen::Map<Eigen::Matrix<float, 4, 1, Eigen::ColMajor>, Eigen::Unaligned>
    map() {
      return Eigen::Map<Eigen::Matrix<float, 4, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<float*>(this), 4);
    }
    Eigen::Map<const Eigen::Matrix<float, 4, 1, Eigen::ColMajor>,
               Eigen::Unaligned>
    map() const {
      return Eigen::Map<const Eigen::Matrix<float, 4, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<const float*>(this),
                                          4);
    }
  };

  std::tuple<float, float, float> robust_eigenstuff(float IIxx, float IIxy, float IIyy) {
    float A = 0.5f * (IIxx + IIyy);
    float B = 0.5f * (IIxx - IIyy);
    float eps = 1e-8;
    float IIxy2 = IIxy*IIxy;
    float S = std::sqrt(B*B + IIxy2 + eps);
    int k = (IIxx - IIyy) < 0 ? -1 : 1;
    float BkS = B + k * S;
    float lam1 = A + S;
    float lam2 = A - S;
    float c2 = 0.5f + k * (0.5f - IIxy2 / (BkS*BkS + IIxy2));
    return std::make_tuple(lam1, lam2, c2);
  }

  template <int N>
  class TexND {
    public:

    void initAxes() {
      // precompute axes inverse
      m_axesinv.resize(m_axes.size());
      m_accsize.reserve(m_axes.size()+1);
      for (size_t i = 0; i < N; ++i) {
        auto& ax = m_axes[i];
        auto& axinv = m_axesinv[i];
        axinv.resize(int(ax.size()) - 1);
        for (size_t j = 0; j < axinv.size(); ++j)
          axinv[j] = 1 / (ax[j + 1] - ax[j]);
        
        if (i == 0)
          m_accsize.push_back(1);
        else
          m_accsize.push_back(m_accsize.back() * m_axes[i-1].size());
      }
      m_accsize.push_back(m_accsize.back() * m_axes.back().size());
    }

    // std::tuple<Vector4s, float, float> sample (const std::vector<float>& in) {
    Vector4s sample (const std::vector<float>& in, int pix) {
      std::vector<float> A;
      std::vector<int> C;
      A.reserve(N);
      C.reserve(N);
      // assert(int(in.size())==N && "Incorrect texture sample input length.");
      // Debug::msgassert("Incorrect texture sample input length.", int(in.size())==N);

      for (int i = 0; i < N; i++) {
        float val = in[i];
        const auto& ax = m_axes[i];
        int c = std::distance(ax.begin(),
                              std::upper_bound(ax.begin() + 1, ax.end() - 1,
                                                val)) -
                1;  // first cell (c,c+1) s.t. val < ax[i+1] else last = size-2
        float a = (val - ax[c]) * m_axesinv[i][c];
        a = std::min(std::max(float(0), a), float(1)); // clamp extrapolation
        A.push_back(a);
        C.push_back(c);
      }

      AlignedVector<Vector4s> samples;
      samples.reserve(Pow<2,N>::result);

      for (int i = 0; i < Pow<2,N>::result; i++) {
        // for halfing procedure order if N=4:
        // c0   c1 c2 c3
        // c0+1 c1 c2 c3
        // ...
        // slowest change in last axis, so halfing happens over that one
        // ie halfing inverse axes order

        int dix = data_index(C, i, pix);
        samples.push_back(m_data[dix].map());
      }

      ND_lerp_half(samples, A);
      return samples[0];
    }

    std::vector<std::vector<scalar>>& getAxes(){return m_axes;}
    std::vector<DeformationEntry>& getData(){return m_data;}
    
    private:

    int data_index (const std::vector<int>& axixs, int binary_offset, int pix) {
      // python:
      // dix = isx + (len(SX)) * isa + (len(SX) * len(SA)) * isy + (len(SX) * len(SA) * len(SY)) * ibend + (len(SX) * len(SA) * len(SY) * len(bendarr)) * vix
      // axes order: SX SA SY bend
      
      // binary offset offsets individual ixs for hypercube cell lookup
      int dix = 0;
      for (int i = 0; i < N; i++) {
        // e.g. binary=0010 means using axix[1]+1, and axix[i]+0 otherwise
        dix += m_accsize[i] * (axixs[i] + ((binary_offset >> i) & 0b1));
      }
      return dix + m_accsize.back() * pix;
    }

    void ND_lerp_half(AlignedVector<Vector4s> & samples, const std::vector<float>& alpha) {
      size_t sz = samples.size();
      for (int i = 0; i < N; i++) {
        // samples[ < len(samples)/2] = (1-a[i]) samples[ < len(samples)/2] + (a[i]) samples[ >= len(samples)/2]
        // samples.resize(len(samples)/2)
        size_t sz2 = sz/2;
        float a = alpha[N-1-i]; // inverse order bc sample layout such that last axis changes slowest
        for (size_t j = 0; j < sz2; j++)
          samples[j] = (1-a) * samples[j] + a * samples[j + sz2]; // lerp into first half
        sz = sz2;
      }
      assert(sz == 1);
    }

    std::vector<DeformationEntry> m_data;
    std::vector<std::vector<float>> m_axes;
    std::vector<std::vector<float>> m_axesinv;
    std::vector<int> m_accsize;
  };

  TexND<4> m_tex_SIIx, m_tex_SIIy;

};

#endif  // __MODEL4D__H__