#ifndef __MODELV1__H__
#define __MODELV1__H__

#include "PeriodicYarnPattern.h"
#include "YarnSoup.h" // vertexmsdata

// CPP INCLUDES

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

class ModelV1 {

 public:
  ModelV1(const std::string& folder) {
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
    // m_pyp.rectangulize();  // TODO potentially assume that this is true for new
    // pyp, but it should take only a millisecond anyway
    // TODO ACTUALlY RECTANGULIZING HERE MIGHT BE BAD
    // BECAUSE IF IM USING THE DATA TO LOOK UP G = Xdef - Xref, and Xref is rectangulized while Xdef is deformed version of unrectangulized
    // THEN THIS IS INCONSISTENT AND WILL PRODUCE ERRORS, SPURIOUS DEFORMATION 
    // m_pyp.compute_parametric();  // TODO potentially cache // actually ideally
                                 // already store with python model creation

    //  load v2t and v2y
    v2t.reserve(m_pyp.Q.rows());
    v2y.reserve(m_pyp.Q.rows());
    if (!load_vector(fs::path(folder) / "v2t", v2t)) {
      Debug::error("No 'v2t' file found in folder:", folder);
      return;
    }
    if (!load_vector(fs::path(folder) / "v2y", v2y)) {
      Debug::error("No 'v2y' file found in folder:", folder);
      return;
    }


    // for file sxsasy/y%02d: load tex 

    static const auto regex_y02d = std::regex("y(\\d{2})");
    {
      const auto subfolder = fs::path(folder) / "sxsasy";
      if (!fs::is_directory(subfolder)) {
        Debug::error("Could not load yarn model directory:", subfolder);
        return;
      }
      std::vector<std::pair<std::string, int>> files;
      for (auto& p : fs::directory_iterator(subfolder)) {
        if (!fs::is_regular_file(p))
          continue;
        const auto& str = p.path().filename().string();
        std::smatch match;
        if (std::regex_search(str, match, regex_y02d)) {
          int id = std::stoi(match[1].str());
          if (id >= int(m_tex_sxsasy.size()))
            m_tex_sxsasy.resize(id + 1);
          files.push_back(std::make_pair(p.path().string(), id));
        }
      }
      threadutils::parallel_for_each(files, [&](const std::pair<std::string, int>& f) {
        m_tex_sxsasy[f.second].load(f.first);
      });
    }

    // bending
    {
      const auto subfolder = fs::path(folder) / "IIx";
      if (!fs::is_directory(subfolder)) {
        Debug::error("Could not load yarn model directory:", subfolder);
        return;
      }
      std::vector<std::pair<std::string, int>> files;
      for (auto& p : fs::directory_iterator(subfolder)) {
        if (!fs::is_regular_file(p))
          continue;
        const auto& str = p.path().filename().string();
        std::smatch match;
        if (std::regex_search(str, match, regex_y02d)) {
          int id = std::stoi(match[1].str());
          if (id >= int(m_tex_IIx.size()))
            m_tex_IIx.resize(id + 1);
          files.push_back(std::make_pair(p.path().string(), id));
        }
      }
      threadutils::parallel_for_each(files, [&](const std::pair<std::string, int>& f) {
        m_tex_IIx[f.second].load(f.first);
      });
    }
    {
      const auto subfolder = fs::path(folder) / "IIy";
      if (!fs::is_directory(subfolder)) {
        Debug::error("Could not load yarn model directory:", subfolder);
        return;
      }
      std::vector<std::pair<std::string, int>> files;
      for (auto& p : fs::directory_iterator(subfolder)) {
        if (!fs::is_regular_file(p))
          continue;
        const auto& str = p.path().filename().string();
        std::smatch match;
        if (std::regex_search(str, match, regex_y02d)) {
          int id = std::stoi(match[1].str());
          if (id >= int(m_tex_IIy.size()))
            m_tex_IIy.resize(id + 1);
          files.push_back(std::make_pair(p.path().string(), id));
        }
      }
      threadutils::parallel_for_each(files, [&](const std::pair<std::string, int>& f) {
        m_tex_IIy[f.second].load(f.first);
      });
    }
 
    m_initialized = true;
  }

  // const void set_parametric(std::vector<VertexMSData> & Xms, const std::vector<int>& pypix) {
  //   threadutils::parallel_for(size_t(0), Xms.size(), [&](size_t i) {
  //     auto& v = Xms[i];
  //     int pix = pypix[i];
  //     v.py = v2y[pix];
  //     v.pt = v2t[pix];
  //   });
  // }

  const std::tuple<Vector4s, scalar, scalar> deformation(const Vector6s& strain,
                                                         int pix, bool dbg) {
                                                        //  int vix, bool dbg) {
    Vector4s g = Vector4s::Zero();

    int y         = v2y[pix]; 
    scalar t      = v2t[pix];
    // Vector4s xref = m_pyp.Q.row(vix);
    Vector4s xref = m_pyp.Q.row(pix);
    {
      assert(y < int(m_tex_sxsasy.size()));
      scalar dbg0=0, dbg1=0;
      Vector4s x;
      // std::tie(x,dbg0,dbg1) = m_tex_sxsasy[y].sample({t,strain[0],strain[1],strain[2]});
      x = m_tex_sxsasy[y].sample({t,strain[0],strain[1],strain[2]});
      g += x - xref;

      if (!dbg) { // bending
        float l1, l2, c2;
        std::tie(l1,l2,c2) = robust_eigenstuff(strain[3],strain[4],strain[5]);

        // g += c2 g_x(l1) + (1-c2) g_x(l2)
        // g += c2 g_y(l2) + (1-c2) g_y(l1)

        g += c2 * (m_tex_IIx[y].sample({t,l1}) + m_tex_IIy[y].sample({t,l2}));
        g += (1 - c2) * (m_tex_IIx[y].sample({t,l2}) + m_tex_IIy[y].sample({t,l1}));
        g -= 2 * xref;

      }
      
      return std::make_tuple(g, dbg0, dbg1);
    }
  }

  const PeriodicYarnPattern& getPYP() const { return m_pyp; }
  bool isInitialized() const { return m_initialized; }

 private:
  bool m_initialized;
  PeriodicYarnPattern m_pyp;
  // static const std::vector<std::string> strain_names;

  // vix -> (t,y) until this is part of the main program
  std::vector<scalar> v2t;
  std::vector<int> v2y;

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

  template <typename T>
  T str2num(const std::string& str);

  template<typename T>
  bool load_vector(const std::string& filename, std::vector<T>& vec) {

    std::fstream ifs(filename.c_str(), std::ios::in);
    if (!ifs) {
      std::cerr << "Error: failed to open file " << filename << "\n";
      return false;
    }

    std::string line;
    while (std::getline(ifs, line))
      vec.push_back(str2num<T>(line));

    return true;
  }

  template <int N>
  class TexND {
    public:
    bool load(const std::string& fname) {
      Debug::log("loading",fname);
      
      std::fstream ifs(fname.c_str(), std::ios::in);
      if (!ifs) {
        std::cerr << "Error: failed to open file " << fname << "\n";
        return false;
      }

      std::string line;
      
      
      m_axes.resize(N);
      for (size_t i = 0; i < N; ++i) {
        auto& ax = m_axes[i];
        if (!std::getline(ifs, line))
          return false;
        {
          std::istringstream iss(line);
          int n = std::distance(std::istream_iterator<std::string>(iss),
                      std::istream_iterator<std::string>());
          if (i > 0)
            n -= 1; // ignore strain index
          ax.resize(n);
        }
        {
          std::stringstream linestream(line);
          scalar value;
          if (i > 0)
            linestream >> value; // get rid of dummy first integer
          for (size_t j = 0; j < ax.size(); j++)
            linestream >> ax[j];
        }
      }

      // precompute axes inverse
      m_axesinv.resize(m_axes.size());
      m_accsize.reserve(m_axes.size());
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
      
      // allocate data and load it
      int totalsize = 1;
      for (const auto& ax : m_axes) {
        totalsize *= ax.size();
      }
      m_data.reserve(totalsize);
      for (int i = 0; i < totalsize; ++i) {
        if (!std::getline(ifs, line))
          return false;
        std::stringstream linestream(line);
        for (size_t j = 0; j < 4; j++)
          linestream >> m_data[i][j];
      }

      return true;
    }
    // std::tuple<Vector4s, scalar, scalar> sample (const std::vector<scalar>& in) {
    Vector4s sample (const std::vector<scalar>& in) {
      std::vector<scalar> A;
      std::vector<int> C;
      A.reserve(N);
      C.reserve(N);
      // assert(int(in.size())==N && "Incorrect texture sample input length.");
      // Debug::msgassert("Incorrect texture sample input length.", int(in.size())==N);

      for (int i = 0; i < N; i++) {
        scalar val = in[i];
        const auto& ax = m_axes[i];
        int c = std::distance(ax.begin(),
                              std::upper_bound(ax.begin() + 1, ax.end() - 1,
                                                val)) -
                1;  // first cell (c,c+1) s.t. val < ax[i+1] else last = size-2
        scalar a = (val - ax[c]) * m_axesinv[i][c];
        a = std::min(std::max(scalar(0), a), scalar(1)); // clamp extrapolation
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

        int dix = data_index(C, i); // using i as binary offset means we are creating data
        samples.push_back(m_data[dix]);
      }

      ND_lerp_half(samples, A);
      // return samples[0];
      // int K = 0;
      // return std::make_tuple(samples[0], C[K] * 1.0f/m_axes[K].size(), A[K]);
      return samples[0];
    }
    
    private:

    int data_index (const std::vector<int>& axixs, int binary_offset = 0) {
      // binary offset offsets individual ixs for hypercube cell lookup
      int dix = 0;
      for (int i = 0; i < N; i++)
        dix += m_accsize[i] * (axixs[i] + ((binary_offset >> i) & 0b1));
      return dix;
    }

    void ND_lerp_half(AlignedVector<Vector4s> & samples, const std::vector<scalar>& alpha) {
      size_t sz = samples.size();
      for (int i = 0; i < N; i++) {
        // samples[ < len(samples)/2] = (1-a[i]) samples[ < len(samples)/2] + (a[i]) samples[ >= len(samples)/2]
        // samples.resize(len(samples)/2)
        size_t sz2 = sz/2;
        scalar a = alpha[N-1-i]; // inverse order bc sample layout such that last axis changes slowest
        for (size_t j = 0; j < sz2; j++)
          samples[j] = (1-a) * samples[j] + a * samples[j + sz2]; // lerp into first half
        sz = sz2;
      }
      assert(sz == 1);
    }

    AlignedVector<Vector4s> m_data;
    std::vector<std::vector<scalar>> m_axes; // e.g. ax0 = [t0, t1, t2, t3, ...]
    std::vector<std::vector<scalar>> m_axesinv;
    std::vector<int> m_accsize;
  };
  std::vector<TexND<4>> m_tex_sxsasy;
  std::vector<TexND<2>> m_tex_IIx, m_tex_IIy;

};

#endif  // __MODELV1__H__