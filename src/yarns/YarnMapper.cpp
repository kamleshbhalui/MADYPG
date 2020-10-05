#include "YarnMapper.h"

#include "../utils/debug_includes.h"
#include "../utils/threadutils.h"

void YarnMapper::step() {
  m_timer.tick();

  if (!m_initialized) {
    // set up mesh provider
    switch(m_settings.provider_type) {
      default:
      case Settings::XPBD:
        Debug::log("XPBD not implemented, defaulting to ObjSeq.");
      case Settings::ObjSeq:
        m_meshProvider = std::static_pointer_cast<AbstractMeshProvider>(
            std::make_shared<ObjSeqAnimation>(m_settings.objseq_settings));
        break;
      case Settings::BinSeq:
        m_meshProvider = std::static_pointer_cast<AbstractMeshProvider>(
            std::make_shared<BinSeqAnimation>(m_settings.binseq_settings));
        break;
    }

    Debug::log("# mesh faces:",
               Debug::format_locale(m_meshProvider->getMesh().Fms.rows(),
                                    "en_US.UTF-8"));

    m_model = std::make_unique<Model>(m_settings.modelfolder);

    m_timer.tock("mesh provider & pyp/model init");
  }

  Mesh& mesh = m_meshProvider->getMesh();
  if (mesh.empty() ||
      !m_model->isInitialized()) {  // couldn't load model or no mesh
    m_initialized = true;
    return;
  }

  if (m_initialized) {
    if (!m_settings.repeat_frame) {
      // step the mesh provider
      m_meshProvider->update();
      m_timer.tock("mesh provider update");
    }
  }

  if ((m_meshProvider->materialSpaceChanged() && !m_settings.repeat_frame) ||
      !m_initialized) {
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
    mesh.compute_v2f_map(m_settings.shepard_weights);
    mesh.compute_face_adjacency();  // TODO cache in obj file / or binary cache
    m_timer.tock("mesh invdm v2f adjacency");

    m_soup.assign_triangles(m_grid, mesh);

    m_timer.tock("soup tri bary");

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
  }  // end of MS change

  // @ WORLD SPACE MESH CHANGES

  mesh.compute_face_data();  // TODO make normals optional in case obj file has
                             // normals ? actually might not matter since
                             // strains need area(=normal)
  if (!m_settings.flat_normals)
    mesh.compute_vertex_normals();
  if (!m_settings.flat_strains)
    mesh.compute_vertex_strains();

  m_timer.tock("mesh normals & strains");

  if (m_settings.deform_reference > scalar(1e-10)) {
    deform_reference(mesh, m_settings.flat_strains);
    m_timer.tock("deform");
    m_soup.reassign_triangles(m_grid, mesh, m_settings.default_same_tri);
    m_timer.tock("bary");
  } else {
    int n     = m_soup.num_vertices();
    auto& Xms = m_soup.get_Xms();
    auto& Xws = m_soup.get_Xws();
    Xws.resize(Xms.rows(), Xms.cols() + 2 + 1);
    threadutils::parallel_for(0, n, [&](int i) {
      Xws.row(i) << Xms.row(i), Xms.block<1, 2>(i, 0), 1;
    });
    m_timer.tock("deform");
    m_soup.reassign_triangles(m_grid, mesh, m_settings.default_same_tri);
    m_timer.tock("bary");
  }


  if (m_settings.shell_map)
    shell_map(mesh, m_settings.flat_normals);
  else {
    int nverts = m_soup.num_vertices();
    threadutils::parallel_for(0, nverts, [&](int vix) {
      auto x = m_soup.get_Xws().row(vix);
      x << x(0), x(2), -x(1), x(3), x(0), x(1), 1;
    });
  }

  m_timer.tock("shell map");
  m_initialized = true;
}

void YarnMapper::deform_reference(const Mesh& mesh, bool flat_strains) {
  int nverts = m_soup.num_vertices();
  auto& Xms  = m_soup.get_Xms();
  auto& Xws  = m_soup.get_Xws();
  Xws.resize(Xms.rows(), Xms.cols() + 2 + 1);
  threadutils::parallel_for(0, nverts, [&](int vix) {
    int tri = m_soup.get_tri(vix);
    if (tri < 0)  // skip unassigned vertex
      return;

    const auto& abc = m_soup.get_bary(vix);

    Vector6s s;
    if (flat_strains) {
      s = mesh.strains[tri];
    } else {
      auto ms_ixs = mesh.Fms.row(tri);
      s           = mesh.vertex_strains[ms_ixs[0]] * abc[0] +
          mesh.vertex_strains[ms_ixs[1]] * abc[1] +
          mesh.vertex_strains[ms_ixs[2]] * abc[2];
    }

    Vector4s g;
    float dbg0, dbg1;

    std::tie(g, dbg0, dbg1) =
        m_model->deformation(s, m_soup.getParametric(vix));
    g *= m_settings.deform_reference;
    // store deformed ms coordinates intm. in ws coords
    Xws.row(vix) << Xms.row(vix) + g.transpose(), Xms.block<1, 2>(vix, 0),1;
    // Xws.row(vix) << Xms.row(vix) + g.transpose(), s[0],s[2],1; // DEBUG:
    // visualize strain
    // Xws.row(vix) << Xms.row(vix) + g.transpose(), dbg0, dbg1,
        // 1;  // DEBUG: visualize other
  });
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

  const auto& pyp = m_model->getPYP();

  auto& Xws = m_soup.get_Xws();
  threadutils::parallel_for(0, nverts, [&](int vix) {
    int lix      = m_soup.getParametric(vix);
    int next_vix = m_soup.getNext(vix);
    scalar ratio;
    if (next_vix >= 0)
      ratio = pyp.RL[lix] /
              (Xws.block<1, 3>(next_vix, 0) - Xws.block<1, 3>(vix, 0)).norm();
    else
      ratio = 1;

    Xws(vix, 6) = std::min(std::max(scalar(0.8), ratio), scalar(1.2));
  });
}
