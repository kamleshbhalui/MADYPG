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
  float x, y, z, t;  // material-space vertexedge
  uint32_t paramy;   // yarn id
  float paramt;      // arc length param
  // float rl;
  // int32_t next;
  // uint32_t paramy; // yarn id
  // float paramt;
};
struct VertexWSData {
  float x, y, z, t;  // world-space vertexedge
  float u, v;        // mesh-space ref coords
  float r = 1;       // radius scale
  // float yu,yv; // arclength & radial texture coord
};

class YarnSoup {
 public:
  void fill_from_grid(const PYP& pyp, const Grid& grid);
  void cut_outside();
  void assign_triangles(const Grid& grid, const Mesh& mesh);
  void reassign_triangles(const Grid& grid, const Mesh& mesh,
                          bool default_same = false);
  void generate_index_list();
  // const std::vector<uint32_t>& getIndices() const { return m_indices; } //
  // TODO DELETE
  auto& getIndexBuffer() { return m_indices; }
  MatrixGLf& get_Xms() { return X_ms; }
  const MatrixGLf& get_Xms() const { return X_ms; }
  auto& get_Xws() { return X_ws; }
  const auto& get_Xws() const { return X_ws; }
  int num_vertices() const { return X_ms.rows(); }

  // int get_tri(int vix) const { return m_v2tri[vix]; }
  // const Vector3s& get_bary(int vix) const { return m_vbary[vix]; }
  // Vector3s get_bary(int vix) const { return m_vbary.row(vix); }

  // const std::vector<int>& get_tri_array() const { return m_v2tri; }
  // auto& get_bary_array() const { return B; }
  const VectorBuffer<VertexBaryData>& get_B0() const { return B0; }
  VectorBuffer<VertexBaryData>& get_B0() { return B0; }
  const VectorBuffer<VertexBaryData>& get_B() const { return B; }
  VectorBuffer<VertexBaryData>& get_B() { return B; }

  int getParametric(int vix) { return m_pypix[vix]; }
  int getNext(int vix) const {  // DEBUG for fancy radius
    if (VE(vix, 1) < 0)
      return -1;
    return E(VE(vix, 1), 1);
  }

 private:
  MatrixGLf X_ms;                   // [u v h t] undeformed material space
  VectorBuffer<VertexWSData> X_ws;  // [x y z t] deformed world space
  MatrixXXRMi E;                    // [v0 v1] edges

  VectorBuffer<VertexBaryData> B, B0;  // [abc,tri] for WS and MS

  MatrixXXRMi
      VE;  // [eix_prev, eix_next] incident edges, ONLY NEEDED FOR FANCY RADIUS

  VectorBuffer<uint32_t>
      m_indices;  // [ y0v0 y0v1 y0v2 ... delim y1v0 y1v1 ... ]
  std::vector<int>
      m_pypix;  // index in pyp // TODO replace with something parametric: y,t
};

#endif  // __YARNSOUP__H__