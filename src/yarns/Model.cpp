#include "Model.h"

#include <cnpy/cnpy.h>

#include <filesystem>  // requires C++17, g++ >= 8 & linking to stdc++fs)
// #include <set>
#include <fstream>
#include <sstream>
namespace fs = std::filesystem;
#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

Model::Model(const std::string& folder) {
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
  // m_pyp.rectangulize();  // TODO potentially assume that this is true for
  // new pyp, but it should take only a millisecond anyway
  // TODO ACTUALlY RECTANGULIZING HERE MIGHT BE BAD
  // BECAUSE IF IM USING THE DATA TO LOOK UP G = Xdef - Xref, and Xref is
  // rectangulized while Xdef is deformed version of unrectangulized THEN THIS
  // IS INCONSISTENT AND WILL PRODUCE ERRORS, SPURIOUS DEFORMATION
  // m_pyp.compute_parametric();  // TODO potentially cache // actually
  // ideally already store with python model creation

  {  // sxsasy texture
     // load sxsy axes.txt
    m_tex_sxsasy_axes.cpu().resize(1);
    AxesInfo& axinf = m_tex_sxsasy_axes.cpu()[0];
    std::vector<axwrapper> axs{axwrapper(axinf.lenSX, axinf.SX, axinf.invSX),
                               axwrapper(axinf.lenSA, axinf.SA, axinf.invSA),
                               axwrapper(axinf.lenSY, axinf.SY, axinf.invSY)};
    if (!load_axes(fs::path(folder) / "sxsasy" / "axes.txt", axs)) {
      Debug::error("No 'axes.txt' file found in sxsasy subfolder of:", folder);
      return;
    }

    // load data.npy
    cnpy::NpyArray npyarr =
        cnpy::npy_load(fs::path(folder) / "sxsasy" / "data.npy");
    auto& data = m_tex_sxsasy_data.cpu();
    Debug::msgassert("data.npy is not 2D", npyarr.shape.size() == 2);
    Debug::msgassert(
        "data.npy is not (N*len(Q)) x 4",
        npyarr.shape[1] == 4 && (npyarr.shape[0] % m_pyp.Q.rows() == 0));
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

const std::tuple<Vector4s, scalar, scalar> Model::deformation(
    const Vector6s& strain, uint32_t pix) {
  scalar dbg0 = 0, dbg1 = 0;
  Vector4s g = sample3D(strain.head<3>(), pix);
  return std::make_tuple(g, dbg0, dbg1);
}

bool Model::load_axes(const std::string& filepath,
                      std::vector<axwrapper>& axes) {
  std::fstream ifs(filepath.c_str(), std::ios::in);
  if (!ifs) {
    std::cerr << "Error: failed to open file " << filepath << "\n";
    return false;
  }

  size_t axi = 0;
  std::string line;
  while (!ifs.eof() && axi < axes.size()) {
    assert(axi < axes.size());
    auto& ax = axes[axi];
    std::getline(ifs, line);
    std::stringstream ss(line);
    size_t i = 0;
    while ((ss >> ax.data[i]) && i < AXES_MAX_LENGTH) {
      ++i;
    }
    ax.len = uint32_t(i);

    for (size_t j = 0; j < ax.len - 1; ++j) {
      ax.invdata[j] = 1.0f / (ax.data[j + 1] - ax.data[j]);
    }
    ++axi;
  }

  return true;
}

// // DEBUG
//   uint upper_bound(uint istart,uint iend, uint which_array, float value) {
//     auto axes = m_tex_sxsasy_axes.cpu()[0];
//   uint i,step,count;
//   count = iend-istart;
//   float arrvalue;
//   while (count > 0) {
//     step = count/2;
//     i = istart + step;

//     if (which_array == 0)
//       arrvalue = axes.SX[i];
//     else if (which_array == 1)
//       arrvalue = axes.SA[i];
//     else
//       arrvalue = axes.SY[i];

//     if (value < arrvalue) {
//       count = step;
//     } else {
//       istart = i+1;
//       count -= step + 1;
//     }
//   }
//   return istart;
// }

Vector4s Model::sample3D(Vector3s strain, uint32_t pix) {
  const auto& axes = m_tex_sxsasy_axes.cpu()[0];
  const auto& data = m_tex_sxsasy_data.cpu();

  auto sample_at = [&](int i_sx, int i_sa, int i_sy, int pix) {
    // from python:
    // isx + (len(SX)) * isa + (len(SX) * len(SA)) * isy + (len(SX) * len(SA)
    // * len(SY)) * vix
    int loc =
        i_sx + axes.lenSX * (i_sa + axes.lenSA * (i_sy + axes.lenSY * pix));
    return Vector4s(data[loc].map());
  };
  // if(pix == 10)
  // Debug::log(2 + axes.lenSX * (3 + axes.lenSA * (4 + axes.lenSY * 5)));

  // bad code block here is finding the cell/sample-index the lookup strain
  // falls into, as well as its lerp-weight/fraction
  float a_sx, a_sa, a_sy;
  int i_sx, i_sa, i_sy;
  for (int i = 0; i < 3; i++) {
    scalar val = strain[i];
    float& a   = i == 0 ? a_sx : (i == 1 ? a_sa : a_sy);
    int& c     = i == 0 ? i_sx : (i == 1 ? i_sa : i_sy);
    const uint32_t& len =
        i == 0 ? axes.lenSX : (i == 1 ? axes.lenSA : axes.lenSY);
    const auto& ax = i == 0 ? axes.SX : (i == 1 ? axes.SA : axes.SY);
    const auto& invax =
        i == 0 ? axes.invSX : (i == 1 ? axes.invSA : axes.invSY);
    // first cell (c,c+1) s.t. val < ax[i+1] else last = size-2
    c = std::distance(ax, std::upper_bound(ax + 1, ax + len - 1, val)) - 1;
    a = (val - ax[c]) * invax[c];
    a = std::min(std::max(scalar(0), a), scalar(1));  // clamp extrapolation

    // int check =  upper_bound(1, len - 1, i, val) - 1;
    // if (check != c) {
    //   Debug::log("DIFF",c, check );
    // }
  }

  Vector4s g = Vector4s::Zero();
  g += (1 - a_sx) * (1 - a_sa) * (1 - a_sy) * sample_at(i_sx, i_sa, i_sy, pix);
  g += (1 - a_sx) * (1 - a_sa) * a_sy * sample_at(i_sx, i_sa, i_sy + 1, pix);
  g += (1 - a_sx) * a_sa * (1 - a_sy) * sample_at(i_sx, i_sa + 1, i_sy, pix);
  g += (1 - a_sx) * a_sa * a_sy * sample_at(i_sx, i_sa + 1, i_sy + 1, pix);
  g += a_sx * (1 - a_sa) * (1 - a_sy) * sample_at(i_sx + 1, i_sa, i_sy, pix);
  g += a_sx * (1 - a_sa) * a_sy * sample_at(i_sx + 1, i_sa, i_sy + 1, pix);
  g += a_sx * a_sa * (1 - a_sy) * sample_at(i_sx + 1, i_sa + 1, i_sy, pix);
  g += a_sx * a_sa * a_sy * sample_at(i_sx + 1, i_sa + 1, i_sy + 1, pix);

  return g;
}