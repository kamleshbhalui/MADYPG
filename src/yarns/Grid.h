#ifndef __GRID__H__
#define __GRID__H__

#include <utility>

#include "../EigenDefinitions.h"
#include "../mesh/Mesh.h"
#include "PeriodicYarnPattern.h"

// (sparse) background grid over mesh material space,
// storing triangles overlapping cells.
// cell size equals the periodic pattern size.
class Grid {
 public:
  // create grid to cover mesh uv bounds with cell size (pyp.px, pyp.py)
  void fromTiling(const Mesh& mesh, const PYP& pyp);
  // compute which cells overlap which triangles and store the results
  void overlap_triangles(const Mesh& mesh, float eps = 1e-3);

  int getNx() const { return nx; }
  int getNy() const { return ny; }
  bool inside(int i, int j) const {
    return !(i < 0 || j < 0) && i < ny && j < nx;
  }
  Vector2s lowerLeft(int i, int j) const;
  const Vector2s& getPivot() const { return pivot; }

  template <typename Vec>
  std::pair<int, int> getIndex(const Vec& point) const {
    return getIndex(point[0], point[1]);
  }
  std::pair<int, int> getIndex(scalar x, scalar y) const {
    return std::make_pair(int(std::floor((y - offset[1]) / cy)),
                          int(std::floor((x - offset[0]) / cx)));
  }

  // sparse grid index mapping functions (full i,j -> sparse k)
  const MatrixXXi& get_ij2k() const { return ij2k; }
  const std::vector<std::pair<int, int>>& get_k2ij() const { return k2ij; }
  bool filled(int i, int j) const { return ij2k(i, j) >= 0; }
  int numFilled() const { return n_filled; }
  const std::deque<int>& cell2tris(int i, int j) const {
    return cellk2tris[ij2k(i, j)];
  }

  void print();

 private:
  int nx, ny;       // number of cells
  float cx, cy;     // cell size
  Vector2s offset;  // position lower left corner of cell(0,0)
  Vector2s pivot;   // offset to cell corners (e.g. to align with pyp)

  MatrixXXi ij2k;                         // map cell i,j -> filled cell index k
  std::vector<std::pair<int, int>> k2ij;  // map cell index k -> filled cell i,j
  int n_filled;
  std::vector<std::deque<int>> cellk2tris;

  // for one triangle with coordinates 'coords' check which cells it overlaps
  void overlap_triangle(
      std::vector<std::pair<int, int>>& tri2cell,
      const Eigen::Matrix<scalar, 3, 2, Eigen::RowMajor>& coords, float eps);
};

#endif  // __GRID__H__