#include "Grid.h"

#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

// def rotate(edge):
//     return np.array([-edge[1], edge[0]])

// 2D version of https://stackoverflow.com/a/17503268
bool intersect_box_tri(const Matrix2s& box,
                       const Eigen::Matrix<scalar, 3, 2, Eigen::RowMajor>& tri,
                       bool check_aabb, float eps) {
  // box: [[xmin, xmax],[ymin,ymax]]
  // tri: [xyz xyz xyz]

  if (check_aabb) {
    for (int i = 0; i < 2; ++i) {
      scalar t0 = tri.col(i).minCoeff();
      scalar t1 = tri.col(i).maxCoeff();
      if (t0 > box(i, 1) + eps || t1 < box(i, 0) - eps)
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
    r = hx * std::abs(a[0]) + hy * std::abs(a[1]) +
        eps;       // projected box radius
    c = m.dot(a);  // projected box center
    if (t0 > c + r || t1 < c - r)
      return false;
  }

  return true;
}

Vector2s Grid::lowerLeft(int i, int j) {
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

  Vector2s pivot;  // offset to cell centers
  pivot << 0, 0;

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
}

void Grid::overlap_triangles(const Mesh& mesh, float eps) {
  eps *= (cx + cy) * 0.5f;  // tolerance relative to cell size

  int n_tris = mesh.Fms.rows();

  // tri2cells.clear();
  tri2cells.resize(n_tris);

  // TODO parallel
  threadutils::parallel_for(0, n_tris, [&](int tri) {
    // for (int tri = 0; tri < n_tris; ++tri) {
    auto ixs = mesh.Fms.row(tri);
    // Vector2s uv0 = mesh.U.row(ixs[0]);
    // Vector2s uv1 = mesh.U.row(ixs[1]);
    // Vector2s uv2 = mesh.U.row(ixs[2]);
    Eigen::Matrix<scalar, 3, 2, Eigen::RowMajor> coords;
    coords << mesh.U.row(ixs[0]), mesh.U.row(ixs[1]), mesh.U.row(ixs[2]);

    overlap_triangle(tri2cells[tri], coords, eps);

    // }
  });

  // for (size_t i = 0; i < tri2cells[tri].size(); i++) {
  //   Debug::log(tri, ": ", tri2cells[tri]);
  // }

  // TODO CONTINUE HERE: create filled[i,j] array: dense!
  // grid.allocate_array("filled", dtype=np.bool, fill_value=False)
  // filled = grid["filled"]
  // for tri, cells in enumerate(tri2cells):
  //     for i, j in cells:
  //         filled[i, j] = True

  // TODO PRINT nfilled cells / total cellsnxny

  Debug::log("IDK");
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

  auto ijmin = getIndex(tmin);
  auto ijmax = getIndex(tmax);

  for (int i = ijmin.first; i <= ijmax.first; ++i) {
    if (i < 0 || i >= ny) {
      assert(false && "triangle out side of grid!");
      continue;
    }
    for (int j = ijmin.second; j <= ijmax.second; ++j) {
      if (j < 0 || j >= nx) {
        assert(false && "triangle out side of grid!");
        continue;
      }

      auto ll = lowerLeft(i, j);
      Matrix2s cell_box;  // [[minx,maxx],
                          //  [miny,maxy]]
      cell_box.row(0) << ll[0], ll[0] + cx;
      cell_box.row(1) << ll[1], ll[1] + cy;

      if (intersect_box_tri(cell_box, coords, false, eps))
        tri2cell.push_back(std::make_pair(i, j));
    }
  }
}