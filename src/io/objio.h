#ifndef __OBJIO_H__
#define __OBJIO_H__

#include <assert.h> /* assert */

#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "../utils/threadutils.h"

int count_words(const std::string &str);

// NOTE: due to circumstances this file contains two version of obj loading
// one with templates/Eigen matrices, and one with vector<structs>
// we use the first one in obj2bin.cpp to create binary sequences, and the 
// latter otherwise..

template <typename Mf, typename Mi>
bool load_obj(const std::string &path, Mf &V, Mi &F, Mf &U, Mi &Fms,
              bool with_faces, bool with_uvs) {
  // NOTE: assuming that the obj file contains: v, vt, f
  //  with f either just v0 v1 v2, or v0/uv0 v1/uv1 v2/uv2

  typedef typename Mf::Scalar scalar_type;
  typedef typename Mi::Scalar index_type;

  std::fstream ifs(path.c_str(), std::ios::in);
  if (!ifs) {
    std::cerr << "Error: failed to open file " << path << "\n";
    return false;
  }

  std::deque<std::string> str_f;
  std::deque<std::string> str_v;
  std::deque<std::string> str_vt;

  // TODO maybe go through file once to get number of verts/faces and reserve
  // the str things. if going through file twice and resetting the ifs is faster
  // than constantly reallocating the strings this should help.

  std::string line, kw;
  while (std::getline(ifs, line)) {
    if (line.empty())
      continue;

    std::stringstream linestream(line);
    linestream >> kw;

    if (kw == "v") {
      str_v.push_back(line);
    }
    if (kw == "vt" && with_uvs) {
      str_vt.push_back(line);
    }
    if (kw == "f" && with_faces) {
      str_f.push_back(line);
    }
  }

  {
    V.resize(str_v.size(), 3);
    threadutils::parallel_for(size_t(0), str_v.size(), [&](size_t i) {
      std::string kw;
      std::stringstream linestream(str_v[i]);
      linestream >> kw >> V(i, 0) >> V(i, 1) >> V(i, 2);
    });
  }
  if (with_uvs) {
    U.resize(str_vt.size(), 2);
    threadutils::parallel_for(size_t(0), str_vt.size(), [&](size_t i) {
      std::string kw;
      std::stringstream linestream(str_vt[i]);
      linestream >> kw >> U(i, 0) >> U(i, 1);
    });
  }

  if (with_faces) {
    // serial because of potential quad triangulation
    std::vector<std::vector<index_type>> tmp_f;
    tmp_f.reserve(str_f.size());

    for (size_t i = 0; i < str_f.size(); ++i) {
      std::string kw, w;
      std::stringstream linestream(str_f[i]);
      int N = count_words(str_f[i]) - 1;  // number of vertices in primitive
      assert((N == 3 || N == 4) && "OBJ: expecting triangles or quads");

      linestream >> kw;
      std::vector<index_type> ixs_ms(N), ixs_ws(N);
      size_t j = 0;
      while (linestream >> w) {
        std::stringstream wstream(w);
        index_type msix, wsix;
        char c;
        wstream >> wsix >> c >> msix;
        if (wstream.fail()) {  // couldn't read format 'ws/ms'. trying instead
                               // just 'ws', ie only one set of face indices /
                               // no seam vertices
          wstream.clear();
          wstream.str(w);
          wstream >> wsix;
          msix = wsix;
          assert(!wstream.fail() &&
                 std::string("OBJ: Couldn't read face indices: '" + line + "'")
                     .c_str());
        }
        ixs_ms[j] = msix - 1;  // NOTE: obj files are 1-indexed
        ixs_ws[j] = wsix - 1;
        ++j;
      }

      // push back triangulated
      if (N == 3) {
        tmp_f.push_back(
            {ixs_ms[0], ixs_ms[1], ixs_ms[2], ixs_ws[0], ixs_ws[1], ixs_ws[2]});
      } else if (N == 4) {
        Mf &X =
            (U.size() > 0) ? U : V;  // use uvs for triangulation if it exists
        auto &ixs          = (U.size() > 0) ? ixs_ms : ixs_ws;
        scalar_type diag02 = 0, diag13 = 0;
        diag02 = (X.row(ixs[0]) - X.row(ixs[2])).squaredNorm();
        diag13 = (X.row(ixs[1]) - X.row(ixs[3])).squaredNorm();
        if (diag02 <= diag13) {
          tmp_f.push_back({ixs_ms[0], ixs_ms[1], ixs_ms[2], ixs_ws[0],
                           ixs_ws[1], ixs_ws[2]});
          tmp_f.push_back({ixs_ms[0], ixs_ms[2], ixs_ms[3], ixs_ws[0],
                           ixs_ws[2], ixs_ws[3]});
        } else {
          tmp_f.push_back({ixs_ms[0], ixs_ms[1], ixs_ms[3], ixs_ws[0],
                           ixs_ws[1], ixs_ws[3]});
          tmp_f.push_back({ixs_ms[1], ixs_ms[2], ixs_ms[3], ixs_ws[1],
                           ixs_ws[2], ixs_ws[3]});
        }
      }
    }

    F.resize(tmp_f.size(), 3);
    if (with_uvs)
      Fms.resize(tmp_f.size(), 3);

    threadutils::parallel_for(size_t(0), tmp_f.size(), [&](size_t i) {
      const auto &frow = tmp_f[i];
      F.row(i) << frow[3], frow[4], frow[5];
      if (with_uvs)
        Fms.row(i) << frow[0], frow[1], frow[2];
    });

  }  // with_faces

  return true;
}

template <typename Mf, typename Mi>
bool load_obj(const std::string &path,
              Mf &V) {  // only load updated vertex data
  Mf U;
  Mi F, Fms;
  return load_obj(path, V, F, U, Fms, false, false);
}

template <typename Mf, typename Mi>
bool load_obj(const std::string &path, Mf &V,
              Mi &F) {  // only load world space data, ignore uvs
  Mf U;
  Mi Fms;
  return load_obj(path, V, F, U, Fms, true, false);
}

template <typename Mf, typename Mi>
bool load_obj(const std::string &path, Mf &V, Mi &F, Mf &U,
              Mi &Fms) {  // load full mesh data
  return load_obj(path, V, F, U, Fms, true, true);
}

#include "../mesh/Mesh.h"
// NOTE: same as above but now for the vector<struct> data types ...

bool load_obj(const std::string &path, std::vector<Mesh::WSVertex> &X,
              std::vector<Mesh::Face> &F, std::vector<Mesh::MSVertex> &U,
              std::vector<Mesh::Face> &Fms, bool with_faces, bool with_uvs);
bool load_obj(const std::string &path, std::vector<Mesh::WSVertex> &X,
              std::vector<Mesh::Face> &F, std::vector<Mesh::MSVertex> &U,
              std::vector<Mesh::Face> &Fms);
bool load_obj(const std::string &path, std::vector<Mesh::WSVertex> &X,
              std::vector<Mesh::Face> &F);
bool load_obj(const std::string &path, std::vector<Mesh::WSVertex> &X);

template <typename Mf, typename Mi>
bool load_obj(const std::string &path, Mf &V, Mi &F, Mf &U, Mi &Fms,
              bool with_faces, bool with_uvs);

#endif  //__OBJIO_H__
