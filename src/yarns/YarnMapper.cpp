#include "YarnMapper.h"

#include "../utils/debug_includes.h"
#include "../utils/threadutils.h"

void YarnMapper::step() {
  m_timer.tick();

  if (!m_initialized) {
    // set up mesh provider
    switch (m_settings.provider_type) {
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
               Debug::format_locale(m_meshProvider->getMesh().Fms.cpu().size(),
                                    "en_US.UTF-8"));

    // m_model = std::make_unique<ModelV0>(m_settings.modelfolder);
    m_model = std::make_unique<Model>(m_settings.modelfolder);

    m_timer.tock("@init: mesh provider & pyp/model init");
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
    m_grid.fromTiling(mesh, m_model->getPYP()); // maybe only once unless fast
    m_grid.overlap_triangles(mesh);

    m_timer.tock("@uv: grid setup");

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

      m_timer.tock("@init: soup tiling");
    }


    mesh.compute_invDm();
    mesh.compute_v2f_map(m_settings.shepard_weights);
    mesh.compute_face_adjacency();  // TODO cache in obj file / or binary cache
    m_timer.tock("@uv: mesh invdm v2f adjacency");

    m_soup.assign_triangles(m_grid, mesh);

    m_timer.tock("@uv: soup tri bary");

    m_soup.get_B0().bufferData(Magnum::GL::BufferUsage::StaticDraw);
    mesh.invDmU.bufferData(Magnum::GL::BufferUsage::StaticDraw);
    mesh.U.bufferData(Magnum::GL::BufferUsage::StaticDraw);
    mesh.Fms.bufferData(Magnum::GL::BufferUsage::StaticDraw);
    mesh.F.bufferData(
        Magnum::GL::BufferUsage::StaticDraw);  // push to gpu // current hint
                                               // static bc expecting no change
                                               // in matspace
    m_timer.tock("@uv: gpu buffers");

    if (!m_initialized) {
      m_soup.cut_outside();  // 'delete' unassigned vertices

      // NOTE/TODO: currently skipping deleting by restlength and pruning
      // arrays


      m_soup.generate_index_list(m_model->getPYP().RL);  // TODO prune by length while assembling!
                                     // @ first parallel count: also sum up RL
                                     // and prune (by not adding their indices
                                     // and also deleting their edges), or
                                     // just prune before assembly but that
                                     // means I need to pass through yarns
                                     // before.. (i guess if this happens once
                                     // it doesnt matter, and its nicer
                                     // however the code is more readable!!)

      m_soup.getIndexBuffer().bufferData(Magnum::GL::BufferUsage::StaticDraw);

      m_soup.get_TB().bufferData(); // edge reference binormal to gpu

      m_timer.tock("@init: soup cut & index list & buffer");
    }
  }  // end of MS change

  // @ WORLD SPACE MESH CHANGES

  mesh.X.bufferData();  // push to gpu

  mesh.compute_face_data();
  if (!m_settings.flat_normals)
    mesh.compute_vertex_normals();

  { // face/vertex-deformation gradients to cpu&gpu, for phong deformation and edge binormal transformation
    mesh.compute_vertex_defF();
    mesh.defF.bufferData();
    mesh.vertex_defF.bufferData();
  }

  if (!m_settings.flat_strains)
    mesh.compute_vertex_strains();
  
  #ifdef DO_DEBUG_STATS
  m_dbg.hist_stepcount++;
  float invscale = 1.0f / mesh.strains.size();
  float invrge = 1.0f/(m_dbg.hist_max - m_dbg.hist_min);
  m_dbg.hist_counts.resize(3);
  for (auto& counts : m_dbg.hist_counts)
    counts.resize(m_dbg.hist_nbins, 0);
  for (size_t i = 0; i < mesh.strains.size(); i++)
  {
    const auto& s = mesh.strains[i];
    for (size_t j = 0; j < 3; j++)
    {
      auto& counts = m_dbg.hist_counts[j];
      int bin = std::max(0,std::min(int((s[j] - m_dbg.hist_min) * invrge * (m_dbg.hist_nbins - 1)),m_dbg.hist_nbins-1));

      // ++counts[bin];
      // c = (c*prevn+ 1 )/newn;
      // counts[bin] = (counts[bin]*(m_dbg.hist_stepcount-1) + invscale) / m_dbg.hist_stepcount;
      counts[bin] += invscale;
    }
  }
  #endif

  m_timer.tock("mesh normals & strains");

  if (m_settings.deform_reference > scalar(1e-10)) {
    deform_reference(mesh, m_settings.flat_strains);
    m_timer.tock("deform");
  } else {
    int n     = m_soup.num_vertices();
    auto& Xms = m_soup.get_Xms();
    auto& Xws = m_soup.get_Xws();
    Xws.cpu().resize(Xms.rows());
    // Xws.resize(Xms.rows(), Xms.cols() + 2 + 1);
    threadutils::parallel_for(0, n, [&](int i) {
      Xws.row<float, 11>(i) << Xms.row(i).transpose(), 
      m_soup.get_TB().cpu()[i].mapB(),
      // 0,0,0,
          Xms.block<1, 2>(i, 0).transpose(), 1;
    });
    m_timer.tock("deform");
  }

// {  auto delim = std::numeric_limits<uint32_t>::max();
//   auto& indices = m_soup.getIndexBuffer().cpu();
//   auto& X_ws = m_soup.get_Xws().cpu();
//   for (size_t i = 0; i < indices.size(); i++)
//   {
//     int vix               = indices[i];
//     if (vix != delim) {
//       Debug::log(" ",X_ws[vix].a);
//     }
//     else {
//       Debug::log(" ------ ");
//     }
//   }}


  if (m_settings.gpu_compute) {

    m_soup.get_Xws().bufferData(Magnum::GL::BufferUsage::StreamDraw);


    if (m_settings.flat_normals)
      mesh.normals.bufferData();
    else {
      mesh.vertex_normals.bufferData();
    }
   
    m_ssshader.compute(m_soup.get_Xws().getGPUSize(), m_soup.get_Xws().gpu(),
                   m_soup.get_B0().gpu(), mesh.invDmU.gpu(),
                   m_settings.flat_normals ? mesh.normals.gpu()
                                           : mesh.vertex_normals.gpu(),
                   mesh.X.gpu(), mesh.F.gpu(), mesh.Fms.gpu(),
                   m_settings.shell_map, m_settings.flat_normals);

    m_timer.tock("bary & shell map & buf");
  } else {
    m_soup.reassign_triangles(m_grid, mesh, m_settings.default_same_tri);

    if (m_settings.shell_map)
      shell_map(mesh, m_settings.flat_normals);
    else {
      int nverts = m_soup.num_vertices();
      threadutils::parallel_for(0, nverts, [&](int vix) {
        auto x = m_soup.get_Xws().row<float, 4>(vix);
        x << x(0), x(2), -x(1), x(3), 0, 1,0,0, x(0), x(1), 1; // NOTE: not correct d1. TODO make consistent cpugpu debug functionality
      });
    }

    m_soup.get_Xws().bufferData(
        Magnum::GL::BufferUsage::StreamDraw);  // TODO buffer when needed, for
                                               // cpu-compute mode that is after
                                               // shell mapping, for gpu-compute
                                               // it is before deformref
    m_timer.tock("bary & shell map & buf");
  }


  // TODO CONTINUE HERE buffer TB after shellmap, shader compute TB update from Xws BARRIERS, bind and use in yarngeom



  m_timer.tock("rest");

  m_initialized = true;
}

void YarnMapper::deform_reference(const Mesh& mesh, bool flat_strains) {
  int nverts = m_soup.num_vertices();
  auto& Xms  = m_soup.get_Xms();
  auto& Xws  = m_soup.get_Xws();
  // TODO i sometimes need  both cpudata for access, and buf object for row need
  // to write that nicer TODO CONTINUE GETTING RID OF OLD FUNCS AND USE B0 AND B
  // ... maybe always use .cpu()[vix].. or make rowmap available to cpudata?  I
  // guess .cpu is a more verbose and nicer method
  auto& B0 = m_soup.get_B0().cpu();
  auto& TB = m_soup.get_TB().cpu();
  Xws.cpu().resize(Xms.rows());
  // Xws.resize(Xms.rows(), Xms.cols() + 2 + 1);
  threadutils::parallel_for(0, nverts, [&](int vix) {
    const auto& bary = B0[vix];
    int tri          = bary.tri;
    if (tri < 0)  // skip unassigned vertex
      return;

    Vector3s abc;
    abc << bary.a, bary.b, bary.c; // TODO CLAMP BARY HERE AND FOR NORMALS IN SHMAP (GPU AND CPU)

    Vector6s s;
    if (flat_strains) {
      s = mesh.strains[tri];
    } else {
      auto ms_ixs = mesh.Fms.cpu()[tri].map();
      s           = mesh.vertex_strains[ms_ixs[0]] * abc[0] +
          mesh.vertex_strains[ms_ixs[1]] * abc[1] +
          mesh.vertex_strains[ms_ixs[2]] * abc[2];
    }

    // DEBUG
    #ifdef DO_DEBUG_STATS
    {
      for (size_t i = 0; i < 6; i++)
      {
        s[i]*= m_dbg.strain_toggle[i];
      }
    }
    #endif

    Vector4s g;
    float dbg0, dbg1;

    std::tie(g, dbg0, dbg1) =
        m_model->deformation(s, m_soup.getParametric(vix));
    g *= m_settings.deform_reference;
    // store deformed ms coordinates intm. in ws coords
    // NOTE: buffer.row is a colvector so it expects colvector comma init
    // Xws.row<float, 7>(vix) << Xms.row(vix).transpose() + g,
    //     Xms.block<1, 2>(vix, 0).transpose(), 1;
    Xws.row<float, 11>(vix) << Xms.block<1,4>(vix,0).transpose() + g, Xms(vix,4),
        // 0,0,0,
        TB[vix].mapB(),
        // dbg0,dbg1, 1;
        // (s[1] + 0.5f),(s[2] + 0.5f), 1;
        Xms.block<1, 2>(vix, 0).transpose(), 1;

    // if (vix == 1200)
    //   Debug::log("  ",dbg0);
  });
}

void YarnMapper::shell_map(const Mesh& mesh, bool flat_normals) {
  // x = phi(xi1, xi2) + h n(xi1, xi2); where xi1 xi2 h are pre-deformed
  // reference coordinates

  auto& U    = mesh.U.cpu();
  auto& defF    = mesh.defF.cpu();
  auto& defFv    = mesh.vertex_defF.cpu();
  auto& N    = mesh.normals.cpu();
  auto& Nv   = mesh.vertex_normals.cpu();
  auto& B    = m_soup.get_B().cpu();
  int nverts = m_soup.num_vertices();
  threadutils::parallel_for(0, nverts, [&](int vix) {
    const auto& bary = B[vix];
    int tri          = bary.tri;
    if (tri < 0)  // skip unassigned vertex
      return;
    // NOTE: assuming Xws are (deformed or copied) ms coords
    //   and abc are bary coords _after_ ms-deformation
    auto x_ws = m_soup.get_Xws().row<float, 3>(vix);
    scalar h  = x_ws(2);
    Vector3s abc;
    abc << bary.a, bary.b, bary.c;

    Vector3s n;
    if (flat_normals) {
      n = N[tri].map();
    } else {
      auto ms_ixs = mesh.Fms.cpu()[tri].map();
      // auto ms_ixs = mesh.Fms.row(tri);
      n = Nv[ms_ixs[0]].map() * abc[0] + Nv[ms_ixs[1]].map() * abc[1] +
          Nv[ms_ixs[2]].map() * abc[2];
      // n           = mesh.vertex_normals[ms_ixs[0]] * abc[0] +
      //     mesh.vertex_normals[ms_ixs[1]] * abc[1] +
      //     mesh.vertex_normals[ms_ixs[2]] * abc[2];
      n.normalize();
    }

    auto ws_ixs = mesh.F.cpu()[tri].map();
    // auto ws_ixs  = mesh.F.row(tri);
    Vector3s phi = mesh.X.row<float, 3>(ws_ixs[0]) * abc[0] +
                   mesh.X.row<float, 3>(ws_ixs[1]) * abc[1] +
                   mesh.X.row<float, 3>(ws_ixs[2]) * abc[2];
    ;
    
    
    // phong deformation:
    // phi = sum_i b_i x_i + alpha * sum_i b_i F_i (X-Xi)
    // gradphi = F_tri + alpha * sum_i (b_i F_i + F_i(X-Xi) gradb_i^T)
    if (m_settings.phong_deformation > 0) {
      auto ms_ixs = mesh.Fms.cpu()[tri].map(); // TODO get from top
      for (size_t i = 0; i < 3; ++i) {
        // Debug::log(ms_ixs[i],defFv.size());
        phi += m_settings.phong_deformation * (
          abc[i] * (defFv[ms_ixs[i]].map() * (x_ws.head<2>() - U[ms_ixs[i]].map()))
        );
      }
    }

    {
      MatrixNMs<3,3> Q; // transformation from uvh to xyz
      Q.setZero();
      Q.block<3,2>(0,0) = defF[tri].map();
      // TODO phong gradient

      Q.col(2) = n; // usually not orthogonal, would need to take transpose inverse of Q! TODO discuss

      // if (vix == 555) {
      //   Debug::log(Q);
      // }
      
      // Q.col(0).normalize();
      // Q.col(1).normalize();
      // Q.col(2) = Q.col(0).cross(Q.col(1)); // IF USING THIS AND SIMPLE Q THAN NEED TO ORTHOGONALIZE COL1!

      auto _xws = m_soup.get_Xws().row<float, 11>(vix); // TODO def above
      // _xws.block<3,1>(3,0) = Q * _xws.block<3,1>(3,0); // b = Q * b
      _xws.block<3,1>(5,0) = Q.inverse().transpose() * _xws.block<3,1>(5,0); // b = Q** * b
    }

    x_ws.head<3>() = phi + h * n;




    // TODO HERE AND GPU:
    // transform B with normalized gradient of shell map and interpolated normal. NOTE transforming using surface at vertex position although location in middle of edge would be more correct. assumption that surface varies little compared to segments, and that later reorthonormalization hides issues anyway. 
  });

  // const auto& pyp = m_model->getPYP();

  // auto& Xws = m_soup.get_Xws();
  // threadutils::parallel_for(0, nverts, [&](int vix) {
  //   int lix      = m_soup.getParametric(vix);
  //   int next_vix = m_soup.getNext(vix);
  //   scalar ratio;
  //   if (next_vix >= 0)
  //     ratio = pyp.RL[lix] /
  //             (Xws.row<float,3>(next_vix) - Xws.row<float,3>(vix)).norm();
  //   else
  //     ratio = 1;

  //   Xws.cpu()[vix].r = std::min(std::max(scalar(0.8), ratio), scalar(1.2));

  //   // TODO consider volumepreservation-based formula (better for paper)
  //   // TODO consider making a parameter
  //   // TODO this might also be a bit buggy? maybe some divisions by 0?
  // });
}
