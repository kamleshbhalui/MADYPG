#include "Grid.h"

#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

// def rotate(edge):
//     return np.array([-edge[1], edge[0]])

// 2D version of https://stackoverflow.com/a/17503268
bool intersect_box_tri(const Matrix2s& box,
                       const Eigen::Matrix<scalar, 3, 2, Eigen::RowMajor>& tri,
                       bool check_aabb);
bool intersect_box_tri(const Matrix2s& box,
                       const Eigen::Matrix<scalar, 3, 2, Eigen::RowMajor>& tri,
                       bool check_aabb) {
  // box: [[xmin, xmax],[ymin,ymax]]
  // tri: [xyz xyz xyz]

  if (check_aabb) {
    for (int i = 0; i < 2; ++i) {
      scalar t0 = tri.col(i).minCoeff();
      scalar t1 = tri.col(i).maxCoeff();
      if (t0 > box(i, 1) || t1 < box(i, 0))
        return false;
    }
  }

  // for each triangle edge
  // project onto its normal a: the two samey verts of the edge, the
  // opposing vert -> [t0, t1] # project box half-lengths to compute a
  // "radius" fast and agnostic of position # hx = xlength/2, hy = ylength/2
  // r = hx * |ax| + hy * |ay|
  // also project box center c = mid . a
  // if min(t0,t1) >  c + r or max(t0,t1) < c - r: return false; not
  // intersecting
  Vector2s m, e, a;
  scalar t0, t1, r, c;

  scalar hx = scalar(0.5) * (box(0, 1) - box(0, 0));
  scalar hy = scalar(0.5) * (box(1, 1) - box(1, 0));
  m << box(0, 0) + hx, box(1, 0) + hy;
  for (int i = 0; i < 3; ++i) {
    e = tri.row((i + 1) % 3) - tri.row(i);
    a << -e(1), e(0);                  // edge normal
    t0 = tri.row((i + 2) % 3).dot(a);  // projected opposing vertex
    t1 = tri.row(i).dot(a);  // projected incident vertex (same for both)
    if (t0 > t1)
      std::swap(t0, t1);
    r = hx * std::abs(a[0]) + hy * std::abs(a[1]);  // projected box radius
    c = m.dot(a);                                   // projected box center
    if (t0 > c + r || t1 < c - r)
      return false;
  }

  return true;
}

Vector2s Grid::lowerLeft(int i, int j) const {
  Vector2s ll;
  ll << cx * j, cy * i;
  ll += offset;
  return ll;
}

void Grid::fromTiling(const Mesh& mesh, const PYP& pyp) {
  // compute uv bounds PARALLELIZABLE (reduce both min and max simulatenously
  // or separately)
  Vector2s uv_min = mesh.U.colwise().minCoeff();
  Vector2s uv_max = mesh.U.colwise().maxCoeff();
  // grow bounds for added robustness
  uv_min -= 0.1f * MakeVec(pyp.px, pyp.py);
  uv_max += 0.1f * MakeVec(pyp.px, pyp.py);

  pivot = pyp.Qmin; // offset to cell corners (e.g. to align with pyp)

  cx = pyp.px;
  cy = pyp.py;

  // compute offset and number
  int i_startx = int(std::floor((uv_min[0] - pivot[0]) / cx));
  int i_endx   = int(std::floor((uv_max[0] - pivot[0]) / cx));
  nx           = i_endx - i_startx + 1;
  int i_starty = int(std::floor((uv_min[1] - pivot[1]) / cy));
  int i_endy   = int(std::floor((uv_max[1] - pivot[1]) / cy));
  ny           = i_endy - i_starty + 1;

  offset[0] = i_startx * cx + pivot[0];
  offset[1] = i_starty * cy + pivot[1];

  Debug::logf("Grid: [%d x %d]\n", ny, nx);
  // Debug::log("grid from to",offset.transpose(), (offset +
  // MakeVec(nx*cx,ny*cy)).transpose());
  // Debug::log("uv from to", uv_min.transpose(),
  // uv_max.transpose());
}

void Grid::overlap_triangles(const Mesh& mesh, float eps) {
  eps *= (cx + cy) * 0.5f;  // tolerance relative to cell size

  int n_tris = mesh.Fms.rows();

  // NOTE for now tri2cells is temporary, and used to construct its inverse
  std::vector<std::vector<std::pair<int, int>>>
      tri2cells;  // tri -> cells [(i0,j0),(i1,j1)]
  tri2cells.resize(n_tris);

  threadutils::parallel_for(0, n_tris, [&](int tri) {
    auto ixs = mesh.Fms.row(tri);
    Eigen::Matrix<scalar, 3, 2, Eigen::RowMajor> coords;
    coords << mesh.U.row(ixs[0]), mesh.U.row(ixs[1]), mesh.U.row(ixs[2]);

    overlap_triangle(tri2cells[tri], coords, eps);
  });

  // mark cells as filled
  ij2k = MatrixXXi::Constant(ny, nx, -1);
  k2ij.reserve(nx * ny);  // conservative / naive upper bound
  int k = 0;
  for (const auto& cells : tri2cells) { // cell list of triangle
    for (const auto& cell : cells) { // cell in list
      if (ij2k(cell.first, cell.second) < 0) {
        ij2k(cell.first, cell.second) = k++;
        k2ij.push_back(cell);
      }
    }
  }
  n_filled = k;

  // finally invert to cell2tris
  cellk2tris.resize(n_filled);
  for (int tri = 0; tri < n_tris; ++tri) {
    for (const auto& cell : tri2cells[tri]) {
      int k = ij2k(cell.first, cell.second);
      // Debug::log(cell.first, cell.second,"-->",k,"/",n_filled);
      cellk2tris[k].push_back(tri);
    }
  }
  
  // Debug::log("----------------------");
  // for (int i = ny-1; i >= 0; i--) {
  //   for (int j = 0; j < nx; j++) {
  //     std::cout << (filled(i, j) ? "F" : ".") << " ";
  //   }
  //   std::cout << "\n";
  // }
  // Debug::log("----------------------");

  Debug::logf("Grid filled %d/%d\n", n_filled, ny * nx);
}

void Grid::overlap_triangle(
    std::vector<std::pair<int, int>>& tri2cell,
    const Eigen::Matrix<scalar, 3, 2, Eigen::RowMajor>& coords, float eps) {
  tri2cell.clear();
  tri2cell.reserve(
      9);  // arbitrary allocation: assumption triangle size ~ cell size
  // alternative: #(ijmax-ijmin)/2 for axis aligned equilateral tri

  // potentially overlapping cell indices from triangle bounds
  auto tmin = coords.colwise().minCoeff();
  auto tmax = coords.colwise().maxCoeff();
  // Vector2s tmin = coords.colwise().minCoeff();
  // Vector2s tmax = coords.colwise().maxCoeff();

  auto ijmin = getIndex(tmin);
  auto ijmax = getIndex(tmax);

  for (int i = ijmin.first; i <= ijmax.first; ++i) {
    if (i < 0 || i >= ny) {
      assert(false && "triangle outside of grid!");
      continue;
    }
    for (int j = ijmin.second; j <= ijmax.second; ++j) {
      if (j < 0 || j >= nx) {
        assert(false && "triangle outside of grid!");
        continue;
      }

      auto ll = lowerLeft(i, j);
      Matrix2s cell_box;  // [[minx,maxx],
                          //  [miny,maxy]]
      cell_box.row(0) << ll[0] - eps, ll[0] + cx + eps;
      cell_box.row(1) << ll[1] - eps, ll[1] + cy + eps;

      if (intersect_box_tri(cell_box, coords, false))
        tri2cell.push_back(std::make_pair(i, j));
    }
  }
}