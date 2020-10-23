
#include "PeriodicYarnPattern.h"

#include <deque>
#include <fstream>
#include <iostream>
#include <limits>

#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

// trim string from left and right (inc. comments)

static inline void trim(std::string &, const std::string &, const std::string &,
                        const std::string &);
static inline void trim(std::string &str,
                        const std::string &lchars  = "\t\n\v\f\r ",
                        const std::string &rchars  = "\t\n\v\f\r ",
                        const std::string &comment = "#") {
  if (!comment.empty()) {
    auto cpos = str.find(comment);
    if (cpos != str.npos)
      str.erase(cpos);
  }
  str.erase(str.find_last_not_of(rchars) + 1);
  str.erase(0, str.find_first_not_of(lchars));
}

// load a file consisting of blocks separated by empty lines into deque(blocks)
// of deques(block) of strings(line)
std::deque<std::deque<std::string>> preload_chunks(const std::string &);
std::deque<std::deque<std::string>> preload_chunks(
    const std::string &filename) {
  std::deque<std::deque<std::string>> chunks;

  std::ifstream ifs(filename);
  if (!ifs) {
    Debug::msgassert("Couldn't load pyp file!", bool(ifs));
    return chunks;
  }

  std::string line;
  while (std::getline(ifs, line)) {
    trim(line);
    if (!line.empty()) {
      chunks.push_back({});
      chunks.back().push_back(line);
      while (std::getline(ifs, line)) {
        trim(line);
        if (line.empty())
          break;
        chunks.back().push_back(line);
      }
    }
  }
  return chunks;
}

template <typename value_type, int ncols, int options>
void deserialize_matrix_chunk(
    const std::deque<std::string> &chunk,
    Eigen::Matrix<value_type, Eigen::Dynamic, ncols, options> &M) {
  assert(M.cols() > 0);  // need to set size of M before calling this
  for (int i = 0; i < int(chunk.size()) - 1; i++) {
    std::stringstream ss(chunk[i + 1]);
    for (int j = 0; j < M.cols(); j++) {
      ss >> M(i, j);
    }
  }
}

template void deserialize_matrix_chunk<float, Eigen::Dynamic, Eigen::RowMajor>(
    const std::deque<std::string> &, MatrixGLf &);
template void deserialize_matrix_chunk<int, Eigen::Dynamic, Eigen::RowMajor>(
    const std::deque<std::string> &, MatrixXXRMi &);
void deserialize_vector_chunk(const std::deque<std::string> &,
                              std::vector<scalar> &);
void deserialize_vector_chunk(const std::deque<std::string> &chunk,
                              std::vector<scalar> &V) {
  for (int i = 0; i < int(chunk.size()) - 1; i++) {
    std::stringstream ss(chunk[i + 1]);
    ss >> V[i];
  }
}

void PYP::deserialize(const std::string &filename) {
  auto chunks = preload_chunks(filename);

  for (auto &chunk : chunks) {
    std::string &header = chunk[0];  // assume trimmed
    trim(header);
    if (header == "*pattern") {
      std::stringstream ss(chunk[1]);
      int n;
      ss >> this->px >> this->py >> this->r;
      ss >> n;
      this->Q.resize(n, 4);
      this->E.resize(n, 4);
      this->RL.resize(n);
    } else if (header == "*Q") {
      deserialize_matrix_chunk(chunk, this->Q);
    } else if (header == "*E") {
      deserialize_matrix_chunk(chunk, this->E);
    } else if (header == "*RL") {
      deserialize_vector_chunk(chunk, this->RL);
    }
  }

  // set min Q because used a couple of times
  Qmin = Q.leftCols<2>().colwise().minCoeff();
}

void PYP::rectangulize()  // NOTE: untested, currently all the pyp files should
                          // satisfy this by default
{
  if (VE.rows() != Q.rows())
    recompute_VE_table();

  // for each vert find di dj such that its position + di*py dj px is within
  // [minx,minx+px) etc apply its new position and apply di dj depending on
  // order to its incident edges
  Vector2s minxy = Qmin;
  // minxy << 0,0;

  for (int i = 0; i < int(Q.rows()); i++) {
    Vector2s xy = Q.row(i).head<2>();
    int dx      = int(std::floor((xy[0] - minxy[0]) / px));
    int dy      = int(std::floor((xy[1] - minxy[1]) / py));

    if (dx != 0 || dy != 0) {
      Q.row(i).head<2>() -= MakeVec(px * dx, py * dy);
      E.row(VE(i, 0)).tail<2>() += MakeVec(dx, dy);  // modify prev edge
      E.row(VE(i, 1)).tail<2>() -= MakeVec(dx, dy);  // modify next edge
    }
  }

  assert("edge crossing multiple boundaries" && [&]() {
    for (int i = 0; i < int(E.rows()); i++)
      if (abs(E(i, 2)) + abs(E(i, 3)) > 1)
        return false;
    return true;
  }());
}

bool PYP::isPeriodicEdge(int eix) {
  return !(E(eix, 2) == 0 && E(eix, 3) == 0);
}

void PYP::recompute_VE_table() {
  VE = MatrixXXRMi::Constant(
      Q.rows(), 2, -1);  // vertex edge table: vix -> [eix_prev, eix_next]
  threadutils::parallel_for(0, int(E.rows()), [&](int i) {
    VE(E(i, 0), 1) = i;  // set as edge as next of its first vertex
    VE(E(i, 1), 0) = i;  // and as prev of its second vertex
  });
}

std::vector<uint32_t> PYP::compute_simple_yarns() {
  if (VE.rows() != Q.rows())
    recompute_VE_table();

  std::vector<uint32_t> ixs;
  ixs.reserve(E.rows());

  for (int i = 0; i < int(E.rows()); ++i) {
    if (isPeriodicEdge(i)) {
      int eix = i;
      int vix = E(eix, 1);
      ixs.push_back(vix);
      eix = VE(vix, 1);

      while (!isPeriodicEdge(eix)) {
        vix = E(eix, 1);
        ixs.push_back(vix);
        eix = VE(vix, 1);
      }
      ixs.push_back(std::numeric_limits<uint32_t>::max());
    }
  }

  return ixs;
}

// yarn index and restlength along-yarn-parameter t ϵ [0, L)
void PYP::compute_parametric() {
  // TODO DEPRECATE, INSTEAD USE FROM MODEL
  if (VE.rows() != Q.rows())
    recompute_VE_table();

  std::deque<int> starts;
  std::vector<bool> visited;
  visited.resize(Q.rows(), false);

  param_v2t.resize(Q.rows());
  param_v2y.resize(Q.rows());
  int y = 0;
  for (int e = 0; e < int(E.rows()); ++e) {
    if (isPeriodicEdge(e)) {
      int vix = E(e, 1);
      if (!visited[vix]) {
        param_y2v.push_back({});
        auto &y2v = param_y2v.back();
        int eix  = VE(vix, 1);
        scalar L = 0;
        do {
          visited[vix]    = true;
          y2v.push_back(vix);
          // y2t.push_back(L);
          param_v2y[vix] = y;
          param_v2t[vix]    = L;
          L += RL[vix]; // NOTE using vix because RL in HYLC computed indexed per vertex for its outgoing edge
          vix = E(eix, 1);
          eix = VE(vix, 1);
        } while (!visited[vix]);
        ++y;
      }
    }
  }

  param_y2t.resize(param_y2v.size());
  threadutils::parallel_for(size_t(0),param_y2v.size(), [&](size_t i){
      auto &y2t = param_y2t[i];
      auto &y2v = param_y2v[i];
      y2t.reserve(y2v.size());
      for (int vix : y2v) {
        y2t.push_back(param_v2t[vix]);
      }
  });

  Debug::log("Y2V 0:\n",param_y2v[0],"\n");
  Debug::log("Y2V 1:\n",param_y2v[1],"\n");

  Debug::log("Y2T 0:\n",param_y2t[0],"\n");
  Debug::log("Y2T 1:\n",param_y2t[1],"\n");
}