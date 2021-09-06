#ifndef __YARNSOUP__H__
#define __YARNSOUP__H__

#include "Grid.h"
#include "PeriodicYarnPattern.h"
#include "../VectorBuffer.h"

// barycentric coordinates and assigned triangle index
struct VertexBaryData {
  float a, b, c;
  int32_t tri;  // can be -1 if unassigned
};

// material space yarn vertex data
struct VertexMSData {
  // NOTE: instead could just define float uvht[4];
  float u, v, h, t;  // material-space vertexedge
  float nx, ny, nz;  // edge normal
  float a;           // arc length parameter along garment yarn
  uint32_t pix;      // vertex index as in pyp
  template <int From, int To> // like python slice arr[from:to]
  auto mapSlice() {
    constexpr int Len = To - From;
    return Eigen::Map<Eigen::Matrix<float, Len, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this) + From, Len);
  }
  auto mapXT() { return mapSlice<0,4>(); }
  auto mapN() { return mapSlice<4,7>(); }
};

// world space yarn vertex data
struct VertexWSData {
  float x, y, z, t;  // world-space vertexedge
  float a;           // arc length parameter
  float nx, ny, nz;  // edge normal
  float u, v;        // mesh-space ref coords
  float r = 1;       // radius scale // NOTE: deprecated

  template <int From, int To> // like python slice arr[from:to]
  auto mapSlice() {
    constexpr int Len = To - From;
    return Eigen::Map<Eigen::Matrix<float, Len, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this) + From, Len);
  }

  auto mapX() { return mapSlice<0,3>(); }
  auto mapN() { return mapSlice<5,8>(); }
  auto map() { return mapSlice<0,11>(); }
};

// Class that owns and creates yarns by tiling over a mesh
class YarnSoup {
 public:
  // tile the periodic yarn pattern (pyp) over the filled cells of a grid
  void fill_from_grid(const PYP& pyp, const Grid& grid);
  // find the underlying triangle for each yarn vertex, compute bary. coords
  void assign_triangles(const Grid& grid, const Mesh& mesh);
  // mark edges of vertices without assigned triangles as deleted
  void cut_outside();
  // recompute barycentric coords (and optionally re-find triangle)
  void reassign_triangles(const Grid& grid, const Mesh& mesh,
                          bool default_same = false);
  // compute arc lengths of yarns, and generate (concatenated) index lists for yarns that are longer than the min. length
  void generate_index_list(const std::vector<scalar> & pypRL, scalar min_length=0);

  auto& getIndexBuffer() { return m_indices; }
  auto& get_Xms() { return X_ms; }
  const auto& get_Xms() const { return X_ms; }
  auto& get_Xws() { return X_ws; }
  const auto& get_Xws() const { return X_ws; }
  int vertexArraySize() const { return X_ms.getCPUSize(); } // vertex array size, including 'deleted' vertices that are not part of any yarn (ie not in any index list, will be ignored in computation and will not be rendered)

  const VectorBuffer<VertexBaryData>& get_B0() const { return B0; }
  VectorBuffer<VertexBaryData>& get_B0() { return B0; }
  const VectorBuffer<VertexBaryData>& get_B() const { return B; }
  VectorBuffer<VertexBaryData>& get_B() { return B; }

  // YarnSoup() {}
  // YarnSoup(YarnSoup && rhs) = default;
  // YarnSoup(const YarnSoup&) = default;
  // void reset() {
  //   X_ms.cpu().clear();
  //   X_ws.cpu().clear();
  //   m_indices.cpu().clear();
  //   B.cpu().clear();
  //   B0.cpu().clear();
  //   E.resize(0,0);
  // }

  int numVertices() const { return m_nvertices; }
  int numYarns() const { return m_nyarns; }
  
 private:
  VectorBuffer<VertexMSData> X_ms;                   // [u v h t a ... ] undeformed material space
  VectorBuffer<VertexWSData> X_ws;  // [x y z t ...] deformed world space
  MatrixXXRMi E;                    // [v0 v1] edges

  VectorBuffer<VertexBaryData> B, B0;  // [abc,tri] for WS and MS

  VectorBuffer<uint32_t>
      m_indices;  // [ y0v0 y0v1 y0v2 ... delim y1v0 y1v1 ... ]

  int m_nvertices;
  int m_nyarns;

};

#endif  // __YARNSOUP__H__