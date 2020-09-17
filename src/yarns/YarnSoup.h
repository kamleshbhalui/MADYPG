#ifndef __YARNSOUP__H__
#define __YARNSOUP__H__

#include "../EigenDefinitions.h"
#include "PeriodicYarnPattern.h"
#include "Grid.h"

class YarnSoup {
 public:
  void fill_from_grid(const PYP& pyp, const Grid& grid);
  void assign_triangles(const Grid& grid, const Mesh& mesh);
  void generate_index_list();
  const std::vector<uint32_t>& getIndices() const { return m_indices; }
  const MatrixGLf& getVertexData() const { return X_ms; }
  // const MatrixGLf& getVertexData() const { return X_ws; }


 private:
  MatrixGLf X_ms;  // [u v h t] undeformed material space
  MatrixGLf X_ws;  // [x y z t] deformed world space
  MatrixXXRMi E;   // [v0 v1] edges

  std::vector<uint32_t> m_indices; // [ y0v0 y0v1 y0v2 ... delim y1v0 y1v1 ... ]
  std::vector<int> m_v2tri; // triangle association per vertex
};

#endif // __YARNSOUP__H__