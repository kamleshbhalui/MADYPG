#include "YarnSoup.h"

#include <limits>

#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"

void YarnSoup::fill_from_grid(const PYP& pyp, const Grid& grid) {
  // NOTE potentially parallelizable in some way using atomics

  // generate mapping from cell i,j to incrementing k for filled cells
  // (overlapping any triangle) which will be used to tile the yarn pattern
  // and define its global memory layout
  int n_tiles      = grid.numFilled();
  const auto& ij2k = grid.get_ij2k();
  const auto& k2ij = grid.get_k2ij();

  int n_verts_tile = pyp.Q.rows();
  int n_edges_tile = pyp.E.rows();  // NOTE: should be the same as Q
  int n_verts      = n_tiles * n_verts_tile;
  // bc periodic edges are counted once between two neighboring tiles, and we
  // ignore(/special case) periodic edges into nonexisting tiles
  int n_edges = n_tiles * n_edges_tile;

  // maps from local pyp index to global index at cell i j
  auto global_vix = [&](int local_vix, int i, int j) {
    int cix = ij2k(i, j);
    assert(cix >= 0);
    return local_vix + n_verts_tile * cix;
  };
  // auto global_eix = [&](int local_eix, int i, int j) {
  //   int cix = ij2k(i, j);
  //   assert(cix >= 0);
  //   return local_eix + n_edges_tile * cix;
  // };



  // precompute tangents
  AlignedVector<Vector3s> tangents;
  tangents.reserve(pyp.Q.rows());
  for (int i = 0; i < pyp.Q.rows(); ++i) {
    Vector3s x0 = pyp.Q.block<1,3>(i,0);
    auto edge = pyp.E.row(pyp.VE(i,1));
    Vector3s x1 = pyp.Q.block<1,3>(i,0);
    Vector3s shift;
    shift << edge(2) * pyp.px, edge(3) * pyp.py, 0;  // account for periodic jump
    tangents[i] = (x1 + shift - x0).normalized();
  }

  // allocate memory for vertex data, vertex/edge topology
  auto& Xms = X_ms.cpu();
  Xms.resize(n_verts);
  // m_pypix.resize(n_verts);  // TODO parametric fake
  E = MatrixXXRMi::Constant(n_edges, 2,
                            -1);  // NOTE parallelizable? (or just set -1
                                  // explicitly in edge copy for borders)

  // parallel over filled grid cells
  threadutils::parallel_for(0, n_tiles, [&](int cix) {
    int i, j;
    std::tie(i, j) = k2ij[cix];

    // copy and shift vertices
    int vix_shift       = n_verts_tile * cix;
    Vector2s vpos_shift = grid.lowerLeft(i, j) - grid.getPivot();
    // grid has been shifted to be aligned with pyp using pivot, but to copy its
    // geometry we have to undo the pivot
    for (int lvix = 0; lvix < n_verts_tile; ++lvix) {
      auto& X = Xms[vix_shift + lvix];
      auto Xvec = X.mapXT();
      Xvec = pyp.Q.row(lvix);
      X.a = 0;
      Xvec.head<2>() += vpos_shift;

      X.pix = lvix;  // TODO parametric fake
      // m_pypix[vix_shift + lvix] = lvix;  // TODO parametric fake

      // optionally check if due to floating point precision ij match the
      // shifted position, and if not try a heuristic of shrinking the copy away
      // from its cell boundaries if even this would fail, then the vertex
      // should get no triangle assigned and thus be marked for cutting
      auto testij = grid.getIndex(Xvec.head<2>());
      if (testij.first != i || testij.second != j) {
        Xvec = pyp.Q.row(lvix) * scalar(0.99);
        Xvec.head<2>() += vpos_shift;
      }
    }

    for (int lvix = 0; lvix < n_verts_tile; ++lvix) {
      auto Bi = Xms[vix_shift + lvix].mapB();
      Vector3s b = pyp.RefD1.row(lvix);
      const Vector3s& t = tangents[lvix];
      Bi = (b - b.dot(t) * t).normalized(); // orthonormalize just in case
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
        neighbor_exists = grid.filled(i + di, j + dj);

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

void YarnSoup::assign_triangles(const Grid& grid, const Mesh& mesh) {
  auto& bary0 = B0.cpu();
  auto& Xms = X_ms.cpu();

  bary0.resize(X_ms.getCPUSize());

  threadutils::parallel_for(0, int(X_ms.getCPUSize()), [&](int vix) {
    bary0[vix].tri = -1;

    int i, j;
    Vector2s p     = Xms[vix].mapXT().head<2>();
    std::tie(i, j) = grid.getIndex(p);

    if (!grid.inside(
            i, j)) {  // NOTE this might happen after deforming reference
                      // (vertex might get pushed outside of previous bounds)
                      // actually not using this after deform anymore, but for
                      // rib it still ends up here, why?
                      // grid not allocated large enough? it might also be because of nonrectangulized pyp ie not all vertices of pyp are inside. then would have to make grid larger by tolerance s.t. those verts are respected 
      Debug::log("WARNING NODE OUTSIDE");
      Debug::log(i, grid.getNy(), j, grid.getNx());
      // TODO fall back to closest grid node and check with its triangles!
      i = std::min(std::max(0, i), grid.getNy() - 1);
      j = std::min(std::max(0, j), grid.getNx() - 1);
    }

    // - due to floating point precision, it might be that the translated tiled
    // vertex is numerically shifted to an adjacent grid cell
    // - this can be fixed by either finding those cases and adapting positions
    // to be inside the prescribed cell, or by cutting them (set tri=-1), since
    // this would happen only at boundaries where there is no neighbor or the
    // neighbor is not filled
    // - this might also happen after deforming ms, in which case there should
    // be an option to use the previously assigned cell (TODO other method
    // "reassign_triangles"?)
    if (!grid.filled(i, j))
      return;  // fallback: ignore/cut vertex

    const auto& tris = grid.cell2tris(i, j);

    for (int tri : tris) {
      Vector3s abc = mesh.barycentric_ms(tri, p);
      if (barycentric_inside(abc)) {
        bary0[vix].tri = tri;
        bary0[vix].a = abc[0];
        bary0[vix].b = abc[1];
        bary0[vix].c = abc[2];
        break;
      }
    }

    // NOTE: for now no fallback closest triangle, or other stuff like
    // 'use_previous' and 'default_previous' ...
  });
}

void YarnSoup::reassign_triangles(const Grid& grid, const Mesh& mesh,
                                  bool default_same) {
  // for reassigning triangles of yarn vertices (that where not cut, i.e. at the
  // beginning they where inside some triangle) the following issues can happen
  // a) the vertex has moved such that it moved outside of the mesh bounds, and
  // is not inside any triangle any more b) if a) then it could have been pushed
  // into a non-filled grid cell or even outside of the grid for those cases the
  // ideal solution would be to find the closest triangle by distance however
  // for now we simply implemented a fallback to the previously assigned
  // triangle in that case, which makes sense, since its deformation didn't push
  // it into some other existing one, and assuming that the deformation is not
  // too large to jump across multiple triangles. if there are artifacts near
  // mesh boundaries, then looking into closest-triangle checks makes sense

  // NOTE: using X_ws here as predeformed ref coords!

  auto& bary0 = B0.cpu();
  auto& bary = B.cpu();
  bary.resize(bary0.size());


  if (default_same) {  // just use previous triangle
    threadutils::parallel_for(size_t(0), X_ws.getCPUSize(), [&](size_t vix) {
      int tri = bary0[vix].tri;
      if (tri < 0)
        return;  // skip unassigned
      Vector2s p   = X_ws.row<float,2>(vix);
      Vector3s abc = mesh.barycentric_ms(tri, p);
      bary[vix].tri = tri;
      bary[vix].a = abc[0];
      bary[vix].b = abc[1];
      bary[vix].c = abc[2];
    });
    return;
  }

  threadutils::parallel_for(size_t(0), X_ws.getCPUSize(), [&](size_t vix) {
    int tri = bary0[vix].tri;
    if (tri < 0)
      return;  // skip unassigned

    Vector2s p   = X_ws.row<float,2>(vix);

    Vector3s abc;
    int i, j;
    std::tie(i, j) = grid.getIndex(p);

    bool empty_cell;
    if (grid.inside(i, j))
      empty_cell = !grid.filled(i, j);
    else
      empty_cell = true;

    if (empty_cell) {  // fall back to previous, see notes above
      abc = mesh.barycentric_ms(tri, p);
    } else {
      const auto& tris = grid.cell2tris(i, j);
      bool found       = false;
      for (int tri_ : tris) {
        abc = mesh.barycentric_ms(tri_, p);
        if (barycentric_inside(abc)) {
          found = true;
          tri   = tri_;
          break;
        }
      }
      if (!found) {  // didnt hit any triangle, see notes above
        abc = mesh.barycentric_ms(tri, p);
      }
    }
    bary[vix].tri = tri;
    bary[vix].a = abc[0];
    bary[vix].b = abc[1];
    bary[vix].c = abc[2];
  });
}

void YarnSoup::cut_outside() {
  auto& bary0 = B0.cpu();
  // remove edges if any of its vertices is not assigned to some triangle
  threadutils::parallel_for(0, int(E.rows()), [&](int eix) {
    int v0 = E(eix, 0);
    int v1 = E(eix, 1);
    bool cut;
    if (v0 < 0 || v1 < 0) {  // bad edge (already marked for deletion somehow)
      cut = true;
    } else {
      cut = bary0[v0].tri < 0 || bary0[v1].tri < 0;  // any vertex unassigned
    }
    if (cut)
      E.row(eix) << -1, -1;  // "delete"
  });
}

void YarnSoup::generate_index_list(const std::vector<scalar> & pypRL, scalar min_length) {
  auto delim = std::numeric_limits<uint32_t>::max();
  // [ y0v0 y0v1 y0v2 ... delim y1v0 y1v1 ... ]

  auto& Xms = X_ms.cpu();

  // construct temporary vertex-edge table
  MatrixXXRMi VE = MatrixXXRMi::Constant(
      X_ms.getCPUSize(), 2, -1);  // vertex edge table: vix -> [eix_prev, eix_next]
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

  std::vector<int> counts(starts.size());  // #verts per yarn
  std::vector<scalar> arclens(starts.size());  // arclengths of yarn

  // for each start (parallel ?)
  //   walk until end (nextedge<0) to count number of verts
  //   store count in counts
  //   .. also accumulate rest lengths as arclength

  threadutils::parallel_for(0, int(starts.size()), [&](int i) {
    int vix = starts[i];
    int c   = 1;
    float arc = 0;
    Xms[vix].a = arc;
    while (hasNext(vix)) {
      arc += pypRL[Xms[vix].pix];
      vix = getNext(vix);
      Xms[vix].a = arc;
      // TODO  assert(self.v_tri_id[vw.vix] >= 0)
      ++c;
      if (c >= VE.rows() + 1) {
        Debug::msgassert("infinite loop in yarnsoup", false);
        break;
      }
    }
    counts[i] = c;
    arclens[i] = arc;
  });

  std::vector<int> skip_yarn(starts.size()); // int because bool is not thread-safe
  threadutils::parallel_for(size_t(0), starts.size(), [&](size_t i) {
    skip_yarn[i] = (counts[i] < 4 || arclens[i] < min_length) ? 1 : 0;
    // NOTE: 4 due to rendering using 4 adjacent verts
  });

  // sum up counts and at the same time overwrite array to encode accumulated counts (inc delims)
  int total = 0;
  int nyarns = 0;
  for (size_t i = 0; i < starts.size(); ++i) {
    if (skip_yarn[i]) {
      counts[i] = 0;
      continue;
    }
    int len    = counts[i];
    counts[i] = total; // replace with offset
    total += len + 1;  // inc delim
    ++nyarns;
  }

  int n_vertices_total = total - nyarns; // - #delims

  this->m_nvertices = n_vertices_total;
  this->m_nyarns = nyarns;
  Debug::log("# yarn vertices:",
             Debug::format_locale(n_vertices_total, "en_US.UTF-8"));
  Debug::log("# yarns:",
             Debug::format_locale(nyarns, "en_US.UTF-8"));

  // use accum. counts to set indices per yarn in parallel
  auto& indices = m_indices.cpu();
  indices.resize(total);
  threadutils::parallel_for(0, int(starts.size()), [&](int i) {
    if (skip_yarn[i])
      return;
    int offset            = counts[i];
    int vix               = starts[i];
    indices[offset + 0] = vix;
    int c                 = 1;
    while (hasNext(vix)) {
      vix                   = getNext(vix);
      indices[offset + c] = vix;
      ++c;
      if (c >= VE.rows() + 1) {
        Debug::msgassert("infinite loop in yarnsoup", false);
        break;
      }
    }
    indices[offset + c] = delim;
  });

  
  // // compute arc length
  // threadutils::parallel_for(0, int(starts.size()), [&](int i) {
  //   float arc = 0;
  //   int vix               = starts[i];
  //   Xms[vix].a = arc;
  //   while (hasNext(vix)) {
  //     // arc += pypRL[m_pypix[vix]];
  //     arc += pypRL[Xms[vix].pix];
  //     vix = getNext(vix);
  //     Xms[vix].a = arc;
  //   }
  // });
  // float arc = 0;
  // for (size_t i = 0; i < indices.size(); i++)
  // {
  //   int vix               = indices[i];
  //   if (vix != delim) {
  //     X_ms(vix,4) = arc;
  //     arc += pypRL[m_pypix[vix]];
  //     Debug::log(" ",X_ms(vix,4));
  //   }
  //   else {
  //     arc = 0;
  //     Debug::log(" ------ ");
  //   }
  // }

  // for (size_t i = 0; i < indices.size(); i++)
  // {
  //   int vix               = indices[i];
  //   if (vix != delim)
  //     Debug::log(" ",X_ms(vix,4));
  //   else
  //     Debug::log(" ------ ");
  // }
}