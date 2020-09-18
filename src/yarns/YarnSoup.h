#ifndef __YARNSOUP__H__
#define __YARNSOUP__H__

#include "../EigenDefinitions.h"
#include "PeriodicYarnPattern.h"
#include "Grid.h"

class YarnSoup {
 public:
  void fill_from_grid(const PYP& pyp, const Grid& grid);
  void cut_outside();
  void assign_triangles(const Grid& grid, const Mesh& mesh);
  void generate_index_list();
  const std::vector<uint32_t>& getIndices() const { return m_indices; }
  MatrixGLf& get_Xms() { return X_ms; }
  const MatrixGLf& get_Xms() const { return X_ms; }
  MatrixGLf& get_Xws() { return X_ws; }
  const MatrixGLf& get_Xws() const { return X_ws; }
  int num_vertices() const { return X_ms.rows(); }
  int get_tri(int vix) const { return m_v2tri[vix]; }
  const Vector3s& get_bary(int vix) const { return m_vbary[vix]; }


 private:
  MatrixGLf X_ms;  // [u v h t] undeformed material space
  MatrixGLf X_ws;  // [x y z t] deformed world space
  MatrixXXRMi E;   // [v0 v1] edges

  std::vector<uint32_t> m_indices; // [ y0v0 y0v1 y0v2 ... delim y1v0 y1v1 ... ]
  std::vector<int> m_v2tri; // triangle association per vertex
  AlignedVector<Vector3s> m_vbary; // tri barycentric coordinates per vertex
};

#endif // __YARNSOUP__H__