#include "YarnMapper.h"

#include "../utils/debug_includes.h"
#include "../utils/threadutils.h"
#include "../io/export_fbx.h"

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
    // m_model = std::make_unique<ModelV1>(m_settings.modelfolder);
    m_model = std::make_unique<Model>(m_settings.modelfolder);

    // buffer model texture to gpu
    m_model->getTexAxes().bufferData();
    m_model->getTexData().bufferData();

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
      // ....

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

      m_soup.get_Xms().bufferData(Magnum::GL::BufferUsage::StaticDraw); // NOTE: important to do this here after assigning arc lengths. otherwise geometry shader will produce NaN due to arclength based weigths.

      m_soup.getIndexBuffer().bufferData(Magnum::GL::BufferUsage::StaticDraw);

      // m_soup.get_TB().bufferData(); // edge reference binormal to gpu

      m_timer.tock("@init: soup cut & index list & buffer");
    }
  }  // end of MS change


  // @ WORLD SPACE MESH CHANGES

  mesh.X.bufferData();  // push to gpu

  bool compute_bending_strains = false;
  mesh.compute_face_data(m_settings.svdclamp, compute_bending_strains);
  if (!m_settings.flat_normals)
    mesh.compute_vertex_normals();

  { // face/vertex-deformation gradients to cpu&gpu, for phong deformation and edge binormal transformation
    mesh.compute_vertex_defF();
    mesh.defF.bufferData();
    if (m_settings.phong_deformation > 0)
      mesh.vertex_defF.bufferData();
  }

  // DEBUG
  #ifdef DO_DEBUG_STATS
  {
    auto& strns = mesh.strains.cpu();
    threadutils::parallel_for(size_t(0),strns.size(),[&](size_t i){
      auto s = strns[i].map();
      for (size_t k = 0; k < 6; k++)
      {
        s[k]*= m_dbg.strain_toggle[std::min(k,size_t(3))]; // NOTE min 3 to mult all bend with same factor
      }
    });
  }
  #endif
  
  if (!m_settings.flat_strains)
    mesh.compute_vertex_strains();

  
  #ifdef DO_DEBUG_STATS
  m_dbg.hist_stepcount++;
  auto& strns = mesh.strains.cpu();
  float invscale = 1.0f / strns.size();
  float invrge = 1.0f/(m_dbg.hist_max - m_dbg.hist_min);
  float invrgeB = 1.0f/(m_dbg.hist_bend + m_dbg.hist_bend);
  m_dbg.hist_counts.resize(6);
  for (auto& counts : m_dbg.hist_counts)
    counts.resize(m_dbg.hist_nbins, 0);
  for (size_t i = 0; i < strns.size(); i++)
  {
    const auto s = strns[i].map();
    for (size_t j = 0; j < 6; j++)
    {
      auto& counts = m_dbg.hist_counts[j];
      int bin;
      if(j < 3)
        bin = std::max(0,std::min(int((s[j] - m_dbg.hist_min) * invrge * (m_dbg.hist_nbins - 1)),m_dbg.hist_nbins-1));
      else
        bin = std::max(0,std::min(int((s[j] + m_dbg.hist_bend) * invrgeB * (m_dbg.hist_nbins - 1)),m_dbg.hist_nbins-1));

      
      // ++counts[bin];
      // c = (c*prevn+ 1 )/newn;
      // counts[bin] = (counts[bin]*(m_dbg.hist_stepcount-1) + invscale) / m_dbg.hist_stepcount;
      counts[bin] += invscale;
    }
  }
  #endif

  // float s=0,b=0,o=0;
  // for (size_t i = 0; i < mesh.strains.size(); i++)
  // {
  //   const auto& st = mesh.strains[i];
  //   s += st(0);
  //   b += st(3);
  //   o = std::max(abs(st[1]),o);
  //   o = std::max(abs(st[2]),o);
  //   o = std::max(abs(st[4]),o);
  //   o = std::max(abs(st[5]),o);
  // }
  // s/=mesh.strains.size();
  // b/=mesh.strains.size();
  // Debug::log("STRAINS",s,b,o);

  


  m_timer.tock("mesh normals & strains");

  if (m_settings.gpu_compute) {
    // TODO option to switch between vertex or flat strains, here fore buffering and in deformshader for usage.
    mesh.vertex_strains.bufferData();
    
    if (m_settings.flat_normals)
      mesh.normals.bufferData();
    else {
      mesh.vertex_normals.bufferData();
    }

    // allocate Xws on gpu if necessary
    if (m_soup.get_Xws().getGPUSize() != m_soup.get_Xms().getGPUSize())
      m_soup.get_Xws().allocateGPU(m_soup.get_Xms().getCPUSize());

    // DEFORM REFERENCE
    m_deformShader.compute(m_soup.get_Xws().getGPUSize(), m_soup.get_Xws().gpu(), m_soup.get_Xms().gpu(), m_soup.get_B0().gpu(), mesh.vertex_strains.gpu(), mesh.Fms.gpu(), m_model->getTexAxes().gpu(), m_model->getTexData().gpu(), m_settings.deform_reference);

    m_timer.tock("deform");

    // SHELL MAP
    if(m_settings.shell_map)
      m_shellMapShader.compute(m_soup.get_Xws().getGPUSize(), m_soup.get_Xws().gpu(),
                   m_soup.get_B0().gpu(), mesh.invDmU.gpu(),
                   m_settings.flat_normals ? mesh.normals.gpu()
                                           : mesh.vertex_normals.gpu(),
                   mesh.X.gpu(), mesh.F.gpu(), mesh.Fms.gpu(),
                   mesh.defF.gpu(), mesh.vertex_defF.gpu(), mesh.U.gpu(),
                   m_settings.flat_normals, m_settings.phong_deformation);
    else
      Magnum::GL::Renderer::setMemoryBarrier(Magnum::GL::Renderer::MemoryBarrier::VertexAttributeArray);

    m_timer.tock("bary & shell map & buf");

  } else { // CPU compute
    // DEFORM REFERENCE
    if (m_settings.deform_reference > scalar(1e-10)) {
      deform_reference(mesh, m_settings.flat_strains);
    } else {
      int n     = m_soup.num_vertices();
      auto& Xms = m_soup.get_Xms().cpu();
      auto& Xws = m_soup.get_Xws();
      Xws.cpu().resize(Xms.size());
      // Xws.resize(Xms.rows(), Xms.cols() + 2 + 1);
      threadutils::parallel_for(0, n, [&](int i) {
        Xws.row<float, 11>(i) << Xms[i].mapXT(), Xms[i].a,
        Xms[i].mapB(),
        // m_soup.get_TB().cpu()[i].mapB(),
        // 0,0,0,
          Xms[i].mapXT().head<2>(), 1;
      });
    }
    m_timer.tock("deform");

    // SHELL MAP
    m_soup.reassign_triangles(m_grid, mesh, m_settings.default_same_tri);

    if (m_settings.shell_map)
      shell_map(mesh, m_settings.flat_normals);
    else {
      int nverts = m_soup.num_vertices();
      threadutils::parallel_for(0, nverts, [&](int vix) {
        auto& x = m_soup.get_Xws().cpu()[vix];
        float tmp = x.y;
        x.y = x.z;
        x.z = -tmp;
      });
    }

    // finally buffer Xws for yarnshader to draw
    m_soup.get_Xws().bufferData(
        Magnum::GL::BufferUsage::StreamDraw);
    glMemoryBarrier(GLbitfield(GL_ALL_BARRIER_BITS));
    m_timer.tock("bary & shell map & buf");
  }

  // { // DEBUG check for opengl errors.
  //   auto err = Magnum::GL::Renderer::error();
  //   if(err != Magnum::GL::Renderer::Error::NoError) {
  //     Debug::error("ERROREND", int(err));
  //   }
  // }
  // Debug::log("---------------");

  m_initialized = true;
}

void YarnMapper::deform_reference(const Mesh& mesh, bool flat_strains) {
  int nverts = m_soup.num_vertices();
  auto& Xms  = m_soup.get_Xms().cpu();
  auto& Xws  = m_soup.get_Xws();
  auto& S  = mesh.strains.cpu();
  auto& Sv  = mesh.vertex_strains.cpu();
  // TODO i sometimes need  both cpudata for access, and buf object for row need
  // to write that nicer TODO CONTINUE GETTING RID OF OLD FUNCS AND USE B0 AND B
  // ... maybe always use .cpu()[vix].. or make rowmap available to cpudata?  I
  // guess .cpu is a more verbose and nicer method
  auto& B0 = m_soup.get_B0().cpu();
  // auto& TB = m_soup.get_TB().cpu();
  Xws.cpu().resize(Xms.size());
  // Xws.resize(Xms.rows(), Xms.cols() + 2 + 1);
  threadutils::parallel_for(0, nverts, [&](int vix) {
    const auto& bary = B0[vix];
    int tri          = bary.tri;
    if (tri < 0)  // skip unassigned vertex
      return;

    Vector3s abc;
    abc << bary.a, bary.b, bary.c;
    // NOTE might want to clamp barycentric coords in case vertex was layed out
    // outside of a triangle,and so using naive bary coords would actually
    // extrapolate strain. or now assuming negligible issue

    Vector6s s;
    if (flat_strains) {
      s = S[tri].map();
    } else {
      auto ms_ixs = mesh.Fms.cpu()[tri].map();
      s           = Sv[ms_ixs[0]].map() * abc[0] +
          Sv[ms_ixs[1]].map() * abc[1] +
          Sv[ms_ixs[2]].map() * abc[2];
    }

    /*const*/ auto& X = Xms[vix];

    Vector4s g;
    float dbg0, dbg1;

    std::tie(g, dbg0, dbg1) =
        m_model->deformation(s, X.pix/*, m_dbg.toggle*/);
        // m_model->deformation(s, m_soup.getParametric(vix), m_dbg.toggle);
    g *= m_settings.deform_reference;
    // store deformed ms coordinates intm. in ws coords
    // NOTE: buffer.row is a colvector so it expects colvector comma init
    // Xws.row<float, 7>(vix) << Xms.row(vix).transpose() + g,
    //     Xms.block<1, 2>(vix, 0).transpose(), 1;
    Xws.row<float, 11>(vix) << X.mapXT() + g, X.a,
        // 0,0,0,
        X.mapB(),
        // TB[vix].mapB(),
        // dbg0,dbg1, 1;
        // (s[0] + 0.5f),(s[2] + 0.5f), 1;
        X.mapXT().head<2>(), 1;
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
      // TODO wrong? not including transformation from uvh0 to uvhdeformed? approximate. or using new barycoords works fine maybe..
      MatrixNMs<3,3> Q; // transformation from uvh to xyz
      Q.setZero();
      Q.block<3,2>(0,0) = defF[tri].map();
      // TODO phong gradient (without it it's just using the flat triangle & interpolated normal, which is probably a decent approximation (since we orthogonalize in the shader) until cheap-transformed refd1 becomes colinear with phongtransformed tangent).

      Q.col(2) = n; // usually not orthogonal, would need to take transpose inverse of Q! TODO discuss

      // if (vix == 555) {
      //   Debug::log(Q);
      // }
      
      // Q.col(0).normalize();
      // Q.col(1).normalize();
      // Q.col(2) = Q.col(0).cross(Q.col(1)); // IF USING THIS AND SIMPLE Q THEN NEED TO ORTHOGONALIZE COL1!

      auto _xws = m_soup.get_Xws().row<float, 11>(vix); // TODO def above
      // _xws.block<3,1>(3,0) = Q * _xws.block<3,1>(3,0); // b = Q * b
      _xws.block<3,1>(5,0) = Q.inverse().transpose() * _xws.block<3,1>(5,0); // b = Q** * b
    }

    x_ws.head<3>() = phi + h * n;

    // NOTE transforming using surface (~phong grad and normal) at vertex position although location in middle of edge would be more correct. assumption that surface varies little compared to segments, and that later reorthonormalization hides issues anyway. 
  });
}

// TODO this takes forever for big garments. probably better to instead prepare appropriate buffers and copy this stuff into a compute shader, though maybe the bottleneck is the disk write..
// also if doing this on cpu then have to get Xws back from gpu!
bool YarnMapper::export2fbx(const std::string& filename) {
  if (!m_initialized)
    return false;

  static std::vector<uint32_t> vcids; // id of circle of vertices, for all non-tip midline vertices
  static std::vector<uint32_t> fcids; // id of circle of faces, for all non-tips exluding also the second to last midline vertex

  auto LIM = std::numeric_limits<uint32_t>::max();
  const auto& I = m_soup.getIndexBuffer().cpu(); // midline vertex indices of LIM separated yarns
  
  // on time setup of vcids and fcids
  static uint32_t nvcid; // number of vertex ids
  static uint32_t nfcid; // number of face ids
  if (vcids.size() == 0) {
    nvcid = 0;
    nfcid = 0;
    vcids.reserve(I.size());
    fcids.reserve(I.size());
    for (size_t i = 0; i < I.size(); ++i) {
      if (i == 0 || i == I.size() - 1 || I[i] == LIM) { // ignore outside verts or tips
        vcids.push_back(LIM);
      } else {
        if (I[i-1] == LIM || I[i+1] == LIM) // LIM-delimited tip
          vcids.push_back(LIM);
        else
          vcids.push_back(nvcid++); // accepted location with existing neighbor left and right
      }
    
      // for faces same as above but one additional skip for second to last vert per yarn
      if (i == 0 || i >= I.size() - 2 || I[i] == LIM) { 
        fcids.push_back(LIM);
      } else {
        if (I[i-1] == LIM || I[i+1] == LIM || I[i+2] == LIM)
          fcids.push_back(LIM);
        else
          fcids.push_back(nfcid++); // accepted location with existing neighbor left and right
      }
    }
  }

  int NSEGS = 8; // number of cylinder sides
  float invSEG = 1.0f/(NSEGS);
  int NVERTICES = NSEGS + 1; // + 1 bc of radial seam
  float R = m_model->getPYP().r; // yarn radius
  float normalTwist      = 1.0f; // TODO normal map params from app. as func params here?
  float normalNum        = 4.0f;

  // allocate export data
  uint32_t num_total_verts = NVERTICES * nvcid;
  uint32_t num_total_faces = 2 * NSEGS * nfcid;
  FBXExportData data;
  data.vertices.resize(num_total_verts * 3);
  data.uv0.resize(num_total_verts * 2);
  data.uv1.resize(num_total_verts * 2);
  data.faces.resize(num_total_faces * 3);
  data.uvfaces.resize(num_total_faces * 3);
  
  auto& Xws = m_soup.get_Xws().cpu(); // xyzt | arc | d1x d1y d1z | u v | rlocal
  threadutils::parallel_for(size_t(0),I.size(),[&](size_t i) {
    uint32_t vcid = vcids[i];
    if (vcid == LIM)
      return;

    // [vert_0] ---edge_A--- [vert_1] ---edge_B--- [vert_2]

    //indices of left, self and right yarn segment
    uint32_t left = I[i-1];
    uint32_t self = I[i];
    uint32_t right = I[i+1];
    auto& x0 = Xws[left];
    auto& x1 = Xws[self];
    auto& x2 = Xws[right];

     // edge rest length
    float lA = x1.a - x0.a;
    float lB = x2.a - x1.a;

    Vector3s tA = x1.mapX() - x0.mapX();
    Vector3s tB = x2.mapX() - x1.mapX();
    
    Vector3s tv = (lA*tA + lB*tB); // vertex tangent
    Vector3s nv = (lA*x0.mapD() + lB*x1.mapD()); // vertex normal
    nv = (nv - tv * tv.dot(nv)/tv.squaredNorm()).normalized(); // orthonormalize
    Vector3s bv = tv.cross(nv).normalized(); // vertex binormal

    // TODO edge current length, and scale radius
    float rA = 1 * R; // TODO vol.preserve
    float rB = 1 * R; // TODO vol.preserve
    float r1 = (lA*rA + lB*rB)/(lA+lB);

    if (i == 66) {
      // Debug::log(x0.mapX());
      // Debug::log(x1.mapX());
      // Debug::log(x2.mapX());
      // Debug::log(tv);
      // Debug::log(nv);
      // Debug::log(bv);
    }


    for(int cs=0; cs<NSEGS+1; cs++) { // iterate circular vertices, including duplicate end
      float alpha = cs * invSEG; // 0 to 1, inclusive

      float a = alpha * 2.0f * float(M_PI);
      float ca = std::cos(a);
      float sa = std::sin(a);
      // TODO edge twist theta ~> get average theta and shift a ?

      Vector3s n = ca*nv + sa * bv;
      Vector3s p = x1.mapX() + r1 * n;  
      
      // if (i == 66) {
      //   std::cout<<p[0]<<", "<<p[1]<<", "<<p[2]<<"\n";
      //   // Debug::log(x0.mapX());
      //   // Debug::log(x1.mapX());
      //   // Debug::log(x2.mapX());
      //   // Debug::log(tv);
      //   // Debug::log(nv);
      //   // Debug::log(bv);
      // }

      uint32_t vertix = vcid * NVERTICES + cs;
      data.vertices[3*vertix + 0] = p[0];
      data.vertices[3*vertix + 1] = p[1];
      data.vertices[3*vertix + 2] = p[2];
      data.uv0[2*vertix + 0] = x1.u;
      data.uv0[2*vertix + 1] = x1.v;
      data.uv1[2*vertix + 0] = x1.a;
      data.uv1[2*vertix + 1] =  normalNum*(-alpha - 0.1591549f * normalTwist * x1.a / r1);
    }


    uint32_t fcid = fcids[i];
    if (fcid == LIM)
      return;

    uint32_t vcid_next = vcids[i+1]; // assert == vcid+1?


    for(int cs=0; cs<NSEGS; cs++) { // iterate circular segments
      // in 2D pseudo indices, with i radial and f the circle/crosssection index:
      // { [f,i], [f+1,i+1], [f+1,i] }
      // { [f,i], [f,i+1], [f+1,i+1] }
      // note that we model a radial seam explicitly, so cs+1 makes sense even for the cs==NSEGS-1

      uint32_t facedofA = 3*(fcid * 2*NSEGS + 2*cs + 0);
      uint32_t facedofB = 3*(fcid * 2*NSEGS + 2*cs + 1);
      uint32_t vertix0 = vcid * NVERTICES + cs; // "[f,cs]"
      uint32_t vertix1 = vcid_next * NVERTICES + cs; // "[f+1,cs]"

      auto F0 = std::vector<uint32_t>{vertix0,vertix1+1,vertix1};
      // auto F0 = std::vector<uint32_t>{0,1,2};
      auto F1 = std::vector<uint32_t>{vertix0,vertix0+1,vertix1+1};

      // last poly index negated with ~
      data.faces[facedofA+0] = F0[0];
      data.faces[facedofA+1] = F0[1];
      data.faces[facedofA+2] = ~int32_t(F0[2]);
      data.faces[facedofB+0] = F1[0];
      data.faces[facedofB+1] = F1[1];
      data.faces[facedofB+2] = ~int32_t(F1[2]);

      data.uvfaces[facedofA+0] = F0[0];
      data.uvfaces[facedofA+1] = F0[1];
      data.uvfaces[facedofA+2] = F0[2];
      data.uvfaces[facedofB+0] = F1[0];
      data.uvfaces[facedofB+1] = F1[1];
      data.uvfaces[facedofB+2] = F1[2];
    }
  });

  return export_fbx(filename, data);
}