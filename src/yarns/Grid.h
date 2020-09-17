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

  const std::vector<std::vector<std::pair<int, int>>>& getTri2cell() const {
    return tri2cells;
  }
  int getNx() const { return nx; }
  int getNy() const { return ny; }
  bool inside(int i, int j) const {
    return !(i < 0 || j < 0) && i < ny && j < nx;
  }
  Vector2s lowerLeft(int i, int j) const;
  const Vector2s& getPivot() const { return pivot; }

 private:
  int nx, ny;       // number of cells
  float cx, cy;     // cell size
  Vector2s offset;  // position lower left corner of cell(0,0)
  Vector2s pivot;   // offset to cell corners (e.g. to align with pyp)

  std::vector<std::vector<std::pair<int, int>>>
      tri2cells;  // tri -> cells [(i0,j0),(i1,j1)]

  template <typename Vec>
  std::pair<int, int> getIndex(const Vec& point) {
    return getIndex(point[0], point[1]);
  }
  std::pair<int, int> getIndex(scalar x, scalar y) {
    return std::make_pair(int(std::floor((y - offset[1]) / cy)),
                          int(std::floor((x - offset[0]) / cx)));
  }
  void overlap_triangle(
      std::vector<std::pair<int, int>>& tri2cell,
      const Eigen::Matrix<scalar, 3, 2, Eigen::RowMajor>& coords, float eps);
};

#endif  // __GRID__H__