#include "meshio.h"

#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

int count_words(const std::string &str);
int count_words(const std::string &str) {
  std::istringstream iss(str);
  return std::distance(std::istream_iterator<std::string>(iss),
                       std::istream_iterator<std::string>());
}

// int count_remaining_words(std::stringstream &ss)
// {
//   auto spos = ss.tellg();
//   auto state = ss.rdstate();
//   int count = std::distance(std::istream_iterator<std::string>(ss),
//   std::istream_iterator<std::string>()); ss.seekg(spos); ss.setstate(state);
//   std::cout<<state<<spos<<"\n";
//   return count;
// }

Mesh load_obj_mesh(const std::string &filename, float scale) {
  Mesh mesh;
  std::deque<VectorGLi> f;      // polygonal faces (world-space)
  std::deque<VectorGLi> fms;    // polygonal faces (material-space)
  AlignedDeque<VectorGL3f> v;   // world coordinates
  AlignedDeque<VectorGL2f> vt;  // texture coordinates / material space

  std::fstream ifs(filename.c_str(), std::ios::in);
  if (!ifs) {
    std::cerr << "Error: failed to open file " << filename << "\n";
    return mesh;
  }

  std::string line;
  std::string kw;
  while (std::getline(ifs, line)) {
    // trim(line);
    if (line.empty())
      continue;

    std::stringstream linestream(line);
    linestream >> kw;

    if (kw == "vt") {
      VectorGL2f vec;
      linestream >> vec[0] >> vec[1];
      vt.push_back(vec * scale);
    } else if (kw == "v") {
      VectorGL3f vec;
      linestream >> vec[0] >> vec[1] >> vec[2];
      v.push_back(vec * scale);
    } else if (kw == "f") {
      int nprimverts = count_words(line) - 1;
      // Debug::msgassert("OBJ: face has unexpected number of indices",
      //                  nprimverts == 3);
      VectorGLi ixs_ms(nprimverts);
      VectorGLi ixs_ws(nprimverts);
      std::string w;
      int j = 0;
      while (linestream >> w) {
        std::stringstream wstream(w);
        int msix, wsix;
        char c;
        wstream >> wsix >> c >> msix;
        if (wstream.fail()) {  // couldn't read format 'ws/ms'. trying instead
                               // just 'ws', ie only one set of face indices /
                               // no seam vertices
          wstream.clear();
          wstream.str(w);
          wstream >> wsix;
          msix = wsix;
          Debug::msgassert(
              "OBJ: Couldn't read arcsim face indices: '" + line + "'",
              !wstream.fail());
        }
        ixs_ms[j] = msix - 1;  // NOTE: obj files are 1-indexed
        ixs_ws[j] = wsix - 1;
        ++j;
      }
      assert(j == nprimverts);
      f.push_back(ixs_ws);
      fms.push_back(ixs_ms);
    }
  }

  // allocate
  mesh.X.resize(v.size(), 3);
  mesh.U.resize(vt.size(), 2);
  if (f.size() > 0) {
    mesh.F.resize(f.size(), f.back().size());
    mesh.Fms.resize(fms.size(), f.back().size());
  }

  // fill
  threadutils::parallel_for(size_t(0), v.size(),
                            [&](size_t i) { mesh.X.row(i) = v[i]; });
  threadutils::parallel_for(size_t(0), vt.size(),
                            [&](size_t i) { mesh.U.row(i) = vt[i]; });
  threadutils::parallel_for(size_t(0), f.size(),
                            [&](size_t i) { mesh.F.row(i) = f[i]; });
  threadutils::parallel_for(size_t(0), fms.size(),
                            [&](size_t i) { mesh.Fms.row(i) = fms[i]; });


  // TODO triangulate quad faces using a mesh.method


  // copy ms into ws // DEBUG
  mesh.F = mesh.Fms;
  mesh.X.resize(mesh.U.rows(), 3);
  mesh.X.col(0) = mesh.U.col(0);
  mesh.X.col(1) = mesh.U.col(1);
  mesh.X.col(2).setZero();

  // mesh.X.resize(4, 3);
  // mesh.X << 0.0f, 0.0f, 0.0f,
  //     1.0f, 0.0f, 0.0f,
  //     1.0f, 1.0f, 0.0f,
  //     0.0f, 1.0f, 0.0f;
  // mesh.X *= 0.1;
  // mesh.U = mesh.X;
  // mesh.F.resize(2, 3);
  // mesh.F << 0, 1, 2,
  //     0, 2, 3;

  Debug::msgassert(
      "OBJ: v and vt data size inconsistent!",
      v.size() <=
          vt.size());  // NOTE: one worldspace vertex can have multiple uvs

  // TODO optionally: assert each index in faces is in [0,nvertices)

  return mesh;
}