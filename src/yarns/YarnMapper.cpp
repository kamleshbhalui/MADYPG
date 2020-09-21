#include "YarnMapper.h"

#include "../utils/debug_includes.h"
#include "../utils/threadutils.h"

void YarnMapper::step() {
  m_timer.tick();

  if (m_initialized) {
    // step the mesh provider
    m_meshProvider->update();
    m_timer.tock("mesh provider update");
  } else {
    // set up mesh provider
    m_meshProvider = std::static_pointer_cast<AbstractMeshProvider>(
        std::make_shared<ObjSeqAnimation>(m_settings.objseq_settings));

    Debug::log("# mesh vertices:",
               Debug::format_locale(m_meshProvider->getMesh().Fms.rows(),
                                    "en_US.UTF-8"));

    m_model = std::make_unique<Model>(m_settings.modelfolder);

    m_timer.tock("mesh provider & pyp/model init");
  }

  Mesh& mesh = m_meshProvider->getMesh();
  // if (false)  //
  // {
  //   // copy ms into ws // DEBUG
  //   mesh.F = mesh.Fms;
  //   mesh.X.resize(mesh.U.rows(), 3);
  //   mesh.X.col(0) = mesh.U.col(0);
  //   mesh.X.col(1) = mesh.U.col(1);
  //   mesh.X.col(2).setZero();
  // }

  if (m_meshProvider->materialSpaceChanged() || !m_initialized) {
    m_grid.fromTiling(mesh, m_model->getPYP());
    m_grid.overlap_triangles(mesh);

    m_timer.tock("grid setup");

    if (!m_initialized) {
      // ... soup
      // NOTE about soup matrices: keep as 'good' matrix only the things
      // needed for gpu, the rest of the stuff could be joint into a vector of
      // structs of vertexdata ? dep on how its used ? . shader might want per
      // vertex:
      //    x y z (t)
      //    u_mesh v_mesh ; for meshspace texturing (but from undeformed ref!)
      //    parametric_t(=acc.rest length) for yarnspace texturing/twist
      //    geom shader promotes and additionally produces parametric_c
      //    (around circle)

      m_soup.fill_from_grid(m_model->getPYP(), m_grid);

      // fancy: do some random uv displacement with multi-level 3D
      // displacement noise (sth like n levels with n strength values) ?
      // actually might be better to do in shader using a texture:
      // meshuv->noise, actually no. still want to do that in uv space (and
      // somehow account for yarns being pushed outside of uvmesh: tribary
      // choose closest tri)

      m_timer.tock("soup tiling");
    }

    mesh.compute_invDm();
    m_soup.assign_triangles(m_grid, mesh);

    m_timer.tock("mesh invdm & soup triangles");

    if (!m_initialized) {
      m_soup.cut_outside();  // 'delete' unassigned vertices

      // NOTE/TODO: currently skipping deleting by restlength and pruning
      // arrays

      m_soup.generate_index_list();  // TODO prune by length while assembling!
                                     // @ first parallel count: also sum up RL
                                     // and prune (by not adding their indices
                                     // and also deleting their edges), or
                                     // just prune before assembly but that
                                     // means I need to pass through yarns
                                     // before.. (i guess if this happens once
                                     // it doesnt matter, and its nicer
                                     // however the code is more readable!!)

      m_timer.tock("soup cut & index list");
    }

    mesh.compute_v2f_map(m_settings.shepard_weights);

    m_timer.tock("mesh v2f");
  }  // end of MS change

  // @ WORLD SPACE MESH CHANGES

  mesh.compute_face_normals();  // TODO make optional in case obj file has
                                // normals
  mesh.compute_vertex_normals();

  // DEBUG in lieu of masm for now just copy soup ms to ws
  // TODO load model (model class <--- pyp) / ws change / ... / masm / ...
  // TODO mesh.strains_face
  // TODO mesh.face2vertex(strains)

  if (false) {
    // TODO // face_strains, vertex_strains @ mesh
    deform_reference(mesh, m_settings.flat_strains);

    // TODO soup.assign_triangles again (at least redo bary)
    // ... maybe soup.reassign_triangles (default_prev, find_closest, ....) 
  } else {
    int n     = m_soup.num_vertices();
    auto& Xms = m_soup.get_Xms();
    auto& Xws = m_soup.get_Xws();
    Xws.resize(Xms.rows(), Xms.cols());
    threadutils::parallel_for(0, n, [&](int i) { Xws.row(i) = Xms.row(i); });
  }

  shell_map(mesh, m_settings.flat_normals);

  m_timer.tock("ws change");
  m_initialized = true;
}

void YarnMapper::shell_map(const Mesh& mesh, bool flat_normals) {
  // x = phi(xi1, xi2) + h n(xi1, xi2); where xi1 xi2 h are pre-deformed
  // reference coordinates

  int nverts = m_soup.num_vertices();
  threadutils::parallel_for(0, nverts, [&](int vix) {
    int tri = m_soup.get_tri(vix);
    if (tri < 0)  // skip unassigned vertex
      return;
    // NOTE: assuming Xws are (deformed or copied) ms coords
    //   and abc are bary coords _after_ ms-deformation
    auto x_ws       = m_soup.get_Xws().row(vix);
    const auto& abc = m_soup.get_bary(vix);
    scalar h        = x_ws(2);

    Vector3s n;
    if (flat_normals) {
      n = mesh.normals[tri];
    } else {
      auto ms_ixs = mesh.Fms.row(tri);
      n           = mesh.vertex_normals[ms_ixs[0]] * abc[0] +
          mesh.vertex_normals[ms_ixs[1]] * abc[1] +
          mesh.vertex_normals[ms_ixs[2]] * abc[2];
      n.normalize();
    }

    auto ws_ixs  = mesh.F.row(tri);
    Vector3s phi = mesh.X.row(ws_ixs[0]) * abc[0] +
                   mesh.X.row(ws_ixs[1]) * abc[1] +
                   mesh.X.row(ws_ixs[2]) * abc[2];
    ;

    x_ws.head<3>() = phi + h * n;
  });
}
