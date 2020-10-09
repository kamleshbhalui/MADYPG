#include "objio.h"

int count_words(const std::string &str) {
  std::istringstream iss(str);
  return std::distance(std::istream_iterator<std::string>(iss),
                       std::istream_iterator<std::string>());
}

bool load_obj(const std::string &path, std::vector<Mesh::WSVertex> &X,
              std::vector<Mesh::Face> &F, std::vector<Mesh::MSVertex> &U,
              std::vector<Mesh::Face> &Fms, bool with_faces, bool with_uvs) {
  // NOTE: assuming that the obj file contains: v, vt, f
  //  with f either just v0 v1 v2, or v0/uv0 v1/uv1 v2/uv2

  typedef uint32_t index_type;
  typedef float scalar_type;

  std::fstream ifs(path.c_str(), std::ios::in);
  if (!ifs) {
    std::cerr << "Error: failed to open file " << path << "\n";
    return false;
  }

  std::deque<std::string> str_f;
  std::deque<std::string> str_v;
  std::deque<std::string> str_vt;

  // TODO maybe go through file once to get number of verts/faces and reserve the str things. if going through file twice and resetting the ifs is faster than constantly reallocating the strings this should help.

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
    X.resize(str_v.size());
    threadutils::parallel_for(size_t(0), str_v.size(), [&](size_t i) {
      std::string kw;
      std::stringstream linestream(str_v[i]);
      linestream >> kw >> X[i].x >> X[i].y >> X[i].z;
    });
  }
  if (with_uvs) {
    U.resize(str_vt.size());
    threadutils::parallel_for(size_t(0), str_vt.size(), [&](size_t i) {
      std::string kw;
      std::stringstream linestream(str_vt[i]);
      linestream >> kw >> U[i].u >> U[i].v;
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
        scalar_type diag02 = 0, diag13 = 0;
        if (U.size() > 0) { // use uvs for triangulation if it exists
          diag02 = (U[ixs_ms[0]].map() - U[ixs_ms[2]].map()).squaredNorm();
          diag13 = (U[ixs_ms[1]].map() - U[ixs_ms[3]].map()).squaredNorm();
        } else {
          diag02 = (X[ixs_ws[0]].map() - X[ixs_ws[2]].map()).squaredNorm();
          diag13 = (X[ixs_ws[1]].map() - X[ixs_ws[3]].map()).squaredNorm();
        }
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

    F.resize(tmp_f.size());
    if (with_uvs)
      Fms.resize(tmp_f.size());

    threadutils::parallel_for(size_t(0), tmp_f.size(), [&](size_t i) {
      const auto &frow = tmp_f[i];
      F[i].map() << frow[3], frow[4], frow[5];
      if (with_uvs)
        Fms[i].map() << frow[0], frow[1], frow[2];
    });

  }  // with_faces

  return true;
}

bool load_obj(const std::string &path, std::vector<Mesh::WSVertex> &X,
              std::vector<Mesh::Face> &F, std::vector<Mesh::MSVertex> &U,
              std::vector<Mesh::Face> &Fms) {
  return load_obj(path, X, F, U, Fms, true, true);
}

bool load_obj(const std::string &path, std::vector<Mesh::WSVertex> &X,
              std::vector<Mesh::Face> &F) {
  std::vector<Mesh::MSVertex> U;
  std::vector<Mesh::Face> Fms;
  return load_obj(path, X, F, U, Fms, true, false);
}

bool load_obj(const std::string &path, std::vector<Mesh::WSVertex> &X) {
  std::vector<Mesh::MSVertex> U;
  std::vector<Mesh::Face> F, Fms;
  return load_obj(path, X, F, U, Fms, true, false);
}