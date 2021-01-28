#ifndef __YARNSOUP__H__
#define __YARNSOUP__H__

#include "Grid.h"
#include "PeriodicYarnPattern.h"
// #include "../EigenDefinitions.h"
#include "../VectorBuffer.h"

struct VertexBaryData {
  float a, b, c;
  int32_t tri;  // can be -1 if unassigned
};
struct VertexMSData {
  // NOTE: instead could just define float uvht[4];
  float u, v, h, t;  // material-space vertexedge
  float bx,by,bz;    // edge binormal
  float a;           // arc length parameter along garment yarn
  uint32_t pix;       // vertex index as in pyp
  template <int From, int To> // like python slice arr[from:to]
  auto mapSlice() {
    constexpr int Len = To - From;
    return Eigen::Map<Eigen::Matrix<float, Len, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this) + From, Len);
  }
  auto mapXT() { return mapSlice<0,4>(); }
  auto mapB() { return mapSlice<4,7>(); }
};
struct VertexWSData {
  float x, y, z, t;  // world-space vertexedge
  float a;           // arc length parameter
  float dx, dy, dz;  // edge director
  float u, v;        // mesh-space ref coords
  float r = 1;       // radius scale
  // float pad0;
  // float yu,yv; // arclength & radial texture coord

  template <int From, int To> // like python slice arr[from:to]
  auto mapSlice() {
    constexpr int Len = To - From;
    return Eigen::Map<Eigen::Matrix<float, Len, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this) + From, Len);
  }

  auto mapX() { return mapSlice<0,3>(); }
  auto mapD() { return mapSlice<5,8>(); }
};
// struct EdgeTBData {
//   // float tx,ty,tz; // tangent
//   float bx,by,bz;  // binormal
//   // auto mapT() {
//   //   return Eigen::Map<Eigen::Matrix<float, 3, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this), 3);
//   // }
//   auto mapB() {
//     // return Eigen::Map<Eigen::Matrix<float, 3, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this)+3, 3);
//     return Eigen::Map<Eigen::Matrix<float, 3, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this), 3);
//   }
// };

class YarnSoup {
 public:
  void fill_from_grid(const PYP& pyp, const Grid& grid);
  void cut_outside();
  void assign_triangles(const Grid& grid, const Mesh& mesh);
  void reassign_triangles(const Grid& grid, const Mesh& mesh,
                          bool default_same = false);
  void generate_index_list(const std::vector<scalar> & pypRL, scalar min_length=0);
  void compute_arclengths(const std::vector<scalar> & pypRL);
  auto& getIndexBuffer() { return m_indices; }
  auto& get_Xms() { return X_ms; }
  const auto& get_Xms() const { return X_ms; }
  auto& get_Xws() { return X_ws; }
  const auto& get_Xws() const { return X_ws; }
  int num_vertices() const { return X_ms.getCPUSize(); }

  const VectorBuffer<VertexBaryData>& get_B0() const { return B0; }
  VectorBuffer<VertexBaryData>& get_B0() { return B0; }
  const VectorBuffer<VertexBaryData>& get_B() const { return B; }
  VectorBuffer<VertexBaryData>& get_B() { return B; }
  // VectorBuffer<EdgeTBData>& get_TB() { return TB; }

  // const auto& getPYPixmap(int vix) { return m_pypix; }
  // int getParametric(int vix) { return m_pypix[vix]; }
  // int getParametric(int vix) { return X_ms[vix].pix; } // maybe reenable this for V1 model if using that to show bending

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

  int numVertices() { return m_nvertices; }
  int numYarns() { return m_nyarns; }
  
 private:
  VectorBuffer<VertexMSData> X_ms;                   // [u v h t a ... ] undeformed material space
  VectorBuffer<VertexWSData> X_ws;  // [x y z t ...] deformed world space
  MatrixXXRMi E;                    // [v0 v1] edges

  VectorBuffer<VertexBaryData> B, B0;  // [abc,tri] for WS and MS

  // VectorBuffer<EdgeTBData> TB;  // edge tangent and binormal, actually just refd1 now, TODO might merge into Xms

  VectorBuffer<uint32_t>
      m_indices;  // [ y0v0 y0v1 y0v2 ... delim y1v0 y1v1 ... ]
  // std::vector<int>
  //     m_pypix;  // index in pyp // TODO replace with something parametric: y,t

  int m_nvertices;
  int m_nyarns;

};

#endif  // __YARNSOUP__H__