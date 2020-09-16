#ifndef __GRID__H__
#define __GRID__H__

#include <utility>

#include "../EigenDefinitions.h"
#include "../mesh/Mesh.h"
#include "PeriodicYarnPattern.h"

class Grid {
 public:
  // create grid to cover mesh uv bounds with cell size (pyp.px, pyp.py)
  void fromTiling(const Mesh& mesh, const PYP& pyp);
  void overlap_triangles(const Mesh& mesh, float eps = 1e-3);

 private:
  int nx, ny;       // number of cells
  float cx, cy;     // cell size
  Vector2s offset;  // position lower left corner of cell(0,0)

  std::vector<std::vector<std::pair<int, int>>>
      tri2cells;  // tri -> cells [(i0,j0),(i1,j1)]

  template <typename Vec>
  std::pair<int, int> getIndex(const Vec& point) {
    return std::make_pair(int(std::floor((point[1] - offset[1])/cy)),
                          int(std::floor((point[0] - offset[0])/cx)));
  }
  Vector2s lowerLeft(int i, int j);
  void overlap_triangle(std::vector<std::pair<int, int>>& tri2cell,
                        const Eigen::Matrix<scalar, 3, 2, Eigen::RowMajor>& coords, float eps);
};

#endif  // __GRID__H__