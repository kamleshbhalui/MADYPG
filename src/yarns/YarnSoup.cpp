#include "YarnSoup.h"

#include <limits>

#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

void YarnSoup::fill_from_grid(const PYP& pyp, const Grid& grid) {
  // NOTE potentially parallelizable in some way using atomics

  // generate mapping from cell i,j to incrementing k for filled cells
  // (overlapping any triangle) which will be used to tile the yarn pattern
  // and define its global memory layout
  const auto& t2c  = grid.getTri2cell();
  int nx           = grid.getNx();
  int ny           = grid.getNy();
  MatrixXXi cellix = MatrixXXi::Constant(ny, nx, -1);
  std::vector<std::pair<int, int>>
      ix2ij;               // inverse for efficient parallel iteration
  ix2ij.reserve(nx * ny);  // conservative / naive upper bound
  int k = 0;
  for (const auto& cells : t2c) {
    for (const auto& cell : cells) {
      if (cellix(cell.first, cell.second) < 0) {
        cellix(cell.first, cell.second) = k++;
        ix2ij.push_back(cell);
      }
    }
  }
  int n_tiles = k;
  Debug::logf("Grid filled %d/%d\n", n_tiles, ny * nx);

  int n_verts_tile = pyp.Q.rows();
  int n_edges_tile = pyp.E.rows();  // NOTE: should be the same as Q
  int n_verts      = n_tiles * n_verts_tile;
  // bc periodic edges are counted once between two neighboring tiles, and we
  // ignore(/special case) periodic edges into nonexisting tiles
  int n_edges = n_tiles * n_edges_tile;

  auto filled = [&](int i, int j) { return cellix(i, j) >= 0; };

  // Debug::log("----------------------");
  // for (int i = ny-1; i >= 0; i--) {
  //   for (int j = 0; j < nx; j++) {
  //     std::cout << (filled(i, j) ? "F" : ".") << " ";
  //   }
  //   std::cout << "\n";
  // }
  // Debug::log("----------------------");

  // maps from local pyp index to global index at cell i j
  auto global_vix = [&](int local_vix, int i, int j) {
    int cix = cellix(i, j);
    assert(cix >= 0);
    return local_vix + n_verts_tile * cix;
  };
  // auto global_eix = [&](int local_eix, int i, int j) {
  //   int cix = cellix(i, j);
  //   assert(cix >= 0);
  //   return local_eix + n_edges_tile * cix;
  // };

  // allocate memory for vertex data, vertex/edge topology
  X_ms.resize(n_verts, 4);
  X_ws.resize(n_verts, 4);
  E = MatrixXXRMi::Constant(n_edges, 2,
                            -1);  // NOTE parallelizable? (or just set -1
                                  // explicitly in edge copy for borders)

  // TODO IGNORE UNTIL USED
  // self.v_pp_id = np.full(
  //     (n_verts,), -1, dtype=np.int)  # periodic patch id
  // self.v_tri_id = np.full((n_verts,), -1, dtype=np.int)
  // self.v_tri_bary = np.empty((n_verts, 3), dtype=np.float)
  // # incident edges: eix_prev eix_next
  // self.v_edges = np.full((n_verts, 2), -1, dtype=np.int)

  // parallel over filled grid cells
  threadutils::parallel_for(0, n_tiles, [&](int cix) {
    int i, j;
    std::tie(i, j) = ix2ij[cix];

    // copy and shift vertices
    int vix_shift       = n_verts_tile * cix;
    Vector2s vpos_shift = grid.lowerLeft(i, j) - grid.getPivot();
    // grid has been shifted to be aligned with pyp using pivot, but to copy its
    // geometry we have to undo the pivot
    for (int lvix = 0; lvix < n_verts_tile; ++lvix) {
      X_ms.row(vix_shift + lvix) = pyp.Q.row(lvix);
      X_ms.row(vix_shift + lvix).head<2>() += vpos_shift;
      // self.v_pp_id[vix_shift + lvix] = lvix TODO IGNORE UNTIL USED
    }

    int eix_shift = n_edges_tile * cix;
    int lv0, lv1, dj, di, gvix0, gvix1;
    bool neighbor_exists;  // NOTE: includes the same-tile case di=dj=0
    for (int leix = 0; leix < n_edges_tile; ++leix) {
      auto edge = pyp.E.row(leix);
      lv0       = edge(0);
      lv1       = edge(1);
      dj        = edge(2);
      di        = edge(3);

      // edge in bounds ?
      neighbor_exists = grid.inside(i + di, j + dj);

      // edge to filled tile ?
      if (neighbor_exists)
        neighbor_exists = filled(i + di, j + dj);

      if (neighbor_exists) {
        gvix0 = global_vix(lv0, i, j);
        gvix1 = global_vix(lv1, i + di, j + dj);
        // NOTE arbitrarily choosing the v0 -> v1+dij direction (compared to
        // v0-dij -> v1) so that only one of the two cells creates the edge in
        // case of tile crossing
        E.row(eix_shift + leix) << gvix0, gvix1;
      } else {
        assert(!(di == 0 && dj == 0));  // else own tile does not exist
        // periodic edge across cell boundary with no adjacent filled cell,
        // leaving edge as -1 to be removed, i.e. uv boundary
      }
    }
  });
}

void YarnSoup::generate_index_list() {
  auto delim = std::numeric_limits<uint32_t>::max();
  // [ y0v0 y0v1 y0v2 ... delim y1v0 y1v1 ... ]

  // construct temporary vertex-edge table
  MatrixXXRMi VE = MatrixXXRMi::Constant(
      X_ms.rows(), 2, -1);  // vertex edge table: vix -> [eix_prev, eix_next]
  threadutils::parallel_for(0, int(E.rows()), [&](int i) {
    if (E(i, 0) >= 0 && E(i, 1) >= 0) {
      VE(E(i, 0), 1) = i;  // set as edge as next of its first vertex
      VE(E(i, 1), 0) = i;  // and as prev of its second vertex
    }
  });

  auto hasPrev = [&](int vix) { return VE(vix, 0) >= 0; };
  auto hasNext = [&](int vix) { return VE(vix, 1) >= 0; };
  auto getNext = [&](int vix) { return E(VE(vix, 1), 1); };

  // mark as starts vertices with no prev edge (< 0) but a next edge (>= 0)
  std::deque<int> starts;
  for (int i = 0; i < int(VE.rows()); ++i) {
    if (!hasPrev(i) && hasNext(i))
      starts.push_back(i);
  }
  std::vector<float> lengths(starts.size());  // #verts per yarn

  // for each start (parallel ?)
  //   walk until end (nextedge<0) to count number of verts
  //   store count in lengths

  threadutils::parallel_for(0, int(starts.size()), [&](int i) {
    int vix = starts[i];
    int c   = 1;
    while (hasNext(vix)) {
      vix = getNext(vix);
      // TODO  assert(self.v_tri_id[vw.vix] >= 0)
      ++c;
      if (c >= VE.rows() + 1) {
        Debug::msgassert("infinite loop in yarnsoup", false);
        break;
      }
    }
    lengths[i] = c;
  });

  // sum up lengths and at the same time overwrite array to encode accumulated
  // lengths (inc delims)
  int total = 0;
  for (size_t i = 0; i < starts.size(); ++i) {
    int len    = lengths[i];
    lengths[i] = total;
    total += len + 1;  // inc delim
  }

  // use accum. lengths to set indices per yarn in parallel
  m_indices.resize(total);
  threadutils::parallel_for(0, int(starts.size()), [&](int i) {
    int offset            = lengths[i];
    int vix               = starts[i];
    m_indices[offset + 0] = vix;
    int c                 = 1;
    while (hasNext(vix)) {
      vix                   = getNext(vix);
      m_indices[offset + c] = vix;
      ++c;
      if (c >= VE.rows() + 1) {
        Debug::msgassert("infinite loop in yarnsoup", false);
        break;
      }
    }
    m_indices[offset + c] = delim;
  });
}