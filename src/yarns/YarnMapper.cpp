#include "YarnMapper.h"

#include <Magnum/GL/Renderer.h>

#include "../io/export_fbx.h"
#include "../utils/debug_includes.h"
#include "../utils/threadutils.h"

#define MODEL4D "data/yarnmodels/model_stock_bend4D"
// #define MODEL4D "data/yarnmodels/model_rib_bend4D"

YarnMapper::YarnMapper() : m_initialized(false) {}


void YarnMapper::step() {
  m_glq3.end();
  // m_timer.tick();
  m_timer.tock("outside CPU");
  m_timer.tockDuration("outside GPU",
                      m_glq3.result<Magnum::UnsignedInt>() / 1000);

  if (!m_initialized) {
    // set up mesh provider
    switch (m_settings.provider_type) {
      default:
      case Settings::PBD:
        m_meshProvider = std::static_pointer_cast<AbstractMeshProvider>(
            std::make_shared<PBDSimulation>(m_settings.pbd_settings));
        break;
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

    #ifdef MODEL4D
    Debug::log("LOADING 4D MODEL.");
    m_model4D = std::make_unique<Model4D>(MODEL4D);
    #endif

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
      m_timer.tock("mesh: update");
    }
  }

  // m_initialized = false; // DEBUG
  // m_soup.reset();
  // m_grid = Grid();

  if ((m_meshProvider->materialSpaceChanged() && !m_settings.repeat_frame) ||
      !m_initialized) {
    // m_grid = Grid();
    if (!m_initialized)  // DEBUG doing this only once makes remeshed stuff
                         // flicker, and then still lose verts?
      m_grid.fromTiling(mesh,
                        m_model->getPYP());  // maybe only once unless fast or
                                             // unnecessary (fixed boundary~)
    m_grid.overlap_triangles(
        mesh);  // NOTE: for tiny patterns with small periodic patch size, the
                // grid is very large and this computation becomes an overhead.
                // if that is important, consider implementing this on the gpu.

    // m_grid.print();

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
    mesh.compute_face_adjacency();  // TODO cache in obj file / or binary cache?
    m_timer.tock("@uv: mesh invdm v2f adjacency");

    m_soup.assign_triangles(m_grid, mesh);

    // { // TODO DEBUGGING, remeshing losing yarn verts. number of nonneg
    // assigned tris remains same.. but still every _reset_/revisit of 0th frame
    // loses stuff.
    //   int c = 0;
    //   auto& B0 = m_soup.get_B0().cpu();
    //   for (size_t i = 0; i < B0.size(); i++)
    //   {
    //     c += B0[i].tri < 0;
    //   }
    //   Debug::logf("skipping %d/%d = %.2f%%\n",c,B0.size(),c*1.0/B0.size());
    // }

    m_timer.tock("@uv: soup tri bary");

    m_soup.get_B0().bufferData(Magnum::GL::BufferUsage::StaticDraw);
    mesh.invDmU.bufferData(Magnum::GL::BufferUsage::StaticDraw);
    mesh.U.bufferData(Magnum::GL::BufferUsage::StaticDraw);
    mesh.Fms.bufferData(Magnum::GL::BufferUsage::StaticDraw);
    mesh.F.bufferData(
        Magnum::GL::BufferUsage::StaticDraw);  // push to gpu // current hint
                                               // static bc expecting no change
                                               // in matspace
    m_timer.tock("@uv: gpu buffers");  // NOTE: this might not actually count
                                       // the gpu time spent?

    if (!m_initialized) {
      m_soup.cut_outside();  // 'delete' unassigned vertices

      // generate index list of yarns for rendering,
      // and also set arclengths along each yarn,
      // remove yarns that are too short
      m_soup.generate_index_list(
          m_model->getPYP().RL,
          m_settings.min_yarn_length_per_r * m_model->getPYP().r);
      // however the code is more readable!!)

      // NOTE/TODO: currently skipping pruning arrays (actually shrinking arrays
      // by removing "deleted" tri=-1 vertices)

      m_soup.get_Xms().bufferData(
          Magnum::GL::BufferUsage::
              StaticDraw);  // NOTE: important to do this here after assigning
                            // arc lengths. otherwise geometry shader will
                            // produce NaN due to arclength based weigths.

      m_soup.getIndexBuffer().bufferData(Magnum::GL::BufferUsage::StaticDraw);

      m_timer.tock("@init: soup cut & index list & buffer");
    }
  }  // end of MS change

  // @ WORLD SPACE MESH CHANGES

  mesh.X.bufferData();  // push to gpu

  // bool compute_bending_strains = m_settings.linearized_bending > 0.0001f;
  mesh.compute_face_data(m_settings.svdclamp /*, compute_bending_strains*/);
  mesh.compute_vertex_normals();

  // face/vertex-deformation gradients to cpu&gpu, for phong deformation and
  // edge binormal transformation
  mesh.compute_vertex_defF();
  mesh.defF.bufferData();
  if (m_settings.phong_deformation > 0)
    mesh.vertex_defF.bufferData();

  mesh.compute_vertex_strains();

  m_timer.tock("mesh: strains");

  if (m_settings.gpu_compute) {
    // and in deformshader for usage.
    mesh.vertex_strains.bufferData();
    mesh.vertex_normals.bufferData();

    m_glq1.begin();

    // allocate Xws on gpu if necessary
    if (m_soup.get_Xws().getGPUSize() != m_soup.get_Xms().getGPUSize())
      m_soup.get_Xws().allocateGPU(m_soup.get_Xms().getCPUSize());

    // DEFORM REFERENCE
    m_deformShader.compute(
        m_soup.get_Xws().getGPUSize(), m_soup.get_Xws().gpu(),
        m_soup.get_Xms().gpu(), m_soup.get_B0().gpu(),
        mesh.vertex_strains.gpu(), mesh.Fms.gpu(), m_model->getTexAxes().gpu(),
        m_model->getTexData().gpu(), m_settings.deform_reference,
        m_settings.linearized_bending, m_settings.svdclamp);

    // memory barrier so that the next compute shader gets updated values
    // Magnum::GL::Renderer::setMemoryBarrier(Magnum::GL::Renderer::MemoryBarrier::ShaderStorage);
    // // the version of Magnum used in this project defines the wrong constant
    // here, defaulting to direct opengl instead
    glMemoryBarrier(GLbitfield(GL_SHADER_STORAGE_BARRIER_BIT));

    m_glq1.end();
    // m_timer.tock("deform");
    m_timer.tockDuration("yarns: deform", m_glq1.result<Magnum::UnsignedInt>() / 1000);

    // SHELL MAP
    m_glq2.begin();
    if (m_settings.shell_map)
      m_shellMapShader.compute(
          m_soup.get_Xws().getGPUSize(), m_soup.get_Xws().gpu(),
          m_soup.get_B0().gpu(), mesh.invDmU.gpu(),
          mesh.vertex_normals.gpu(),
          mesh.X.gpu(), mesh.F.gpu(), mesh.Fms.gpu(), mesh.defF.gpu(),
          mesh.vertex_defF.gpu(), mesh.U.gpu(),
          m_settings.phong_deformation);

    // memory barrier so that vertex shader gets updated values
    Magnum::GL::Renderer::setMemoryBarrier(
        Magnum::GL::Renderer::MemoryBarrier::VertexAttributeArray);

    m_glq2.end();
    // m_timer.tock("bary & shell map & buf");
    m_timer.tockDuration("yarns: map",
                        m_glq2.result<Magnum::UnsignedInt>() / 1000);

  } else {  // CPU compute
    // DEFORM REFERENCE
    if (m_settings.deform_reference > scalar(1e-10)) {
      deform_reference(mesh);
    } else {
      int n     = m_soup.num_vertices();
      auto& Xms = m_soup.get_Xms().cpu();
      auto& Xws = m_soup.get_Xws();
      Xws.cpu().resize(Xms.size());
      // Xws.resize(Xms.rows(), Xms.cols() + 2 + 1);
      threadutils::parallel_for(0, n, [&](int i) {
        Xws.row<float, 11>(i) << Xms[i].mapXT(), Xms[i].a, Xms[i].mapB(),
            // m_soup.get_TB().cpu()[i].mapB(),
            // 0,0,0,
            Xms[i].mapXT().head<2>(), 1;
      });
    }
    m_timer.tock("yarns: deform");

    // SHELL MAP
    constexpr bool same_triangle = true;
    // note: actually just recomputing barycentric coords in assumed same triangle
    m_soup.reassign_triangles(m_grid, mesh, same_triangle);

    if (m_settings.shell_map)
      shell_map(mesh);
    else {
      int nverts = m_soup.num_vertices();
      threadutils::parallel_for(0, nverts, [&](int vix) {
        auto& x   = m_soup.get_Xws().cpu()[vix];
        float tmp = x.y;
        x.y       = x.z;
        x.z       = -tmp;
      });
    }

    // finally buffer Xws for yarnshader to draw
    m_soup.get_Xws().bufferData(Magnum::GL::BufferUsage::StreamDraw);
    glMemoryBarrier(GLbitfield(GL_ALL_BARRIER_BITS));
    m_timer.tock("yarns: map");
  }

  // { // DEBUG check for opengl errors.
  //   auto err = Magnum::GL::Renderer::error();
  //   if(err != Magnum::GL::Renderer::Error::NoError) {
  //     Debug::error("ERROREND", int(err));
  //   }
  // }
  // Debug::log("---------------");

  m_initialized = true;
  m_glq3.begin();
  m_timer.tick();
}

void hermitian_eig_clamp(Vector6s& m, float mineigval);
void hermitian_eig_clamp(Vector6s& m, float mineigval) {
  // following:
  // Closed-form expressions of the eigen decomposition of 2x2 and 3x3 Hermitian matrices
  // https://hal.archives-ouvertes.fr/hal-01501221/document

  // mat C = a,c;c,b  (vec m ~ a,c,b,_,_,_)
  float ab = (m[0] - m[2]); // a - b
  float d = 4*m[1]*m[1] + ab*ab; // d^2
  if (d < 0.0001f) {
    // d == 0 ---> a == b && c == 0 ---> l1 = l2 = a = b
    float l12 = m[0]; 
    l12 = std::max(l12, mineigval);
    m[0] = l12;
    m[2] = l12;
    // m[1] = 0; 
  } else {
    d = sqrt(d);
    float l1 = (m[0]+m[2]-d)*0.5f;
    float l2 = (m[0]+m[2]+d)*0.5f;
    l1 = std::max(l1, mineigval);
    l2 = std::max(l2, mineigval);
    float A = (l2-l1);
    float B = d*(l1+l2);
    float C = A*ab;
    float inv2d = 0.5f/d;
    m[0] = (B+C)*inv2d;
    m[2] = (B-C)*inv2d;
    m[1] *= 2*A*inv2d; 
  }
}

void YarnMapper::deform_reference(const Mesh& mesh) {
  int nverts = m_soup.num_vertices();
  auto& Xms  = m_soup.get_Xms().cpu();
  auto& Xws  = m_soup.get_Xws();
  auto& S    = mesh.strains.cpu();
  auto& Sv   = mesh.vertex_strains.cpu();
  
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
    auto ms_ixs = mesh.Fms.cpu()[tri].map();
    s = Sv[ms_ixs[0]].map() * abc[0] + Sv[ms_ixs[1]].map() * abc[1] +
        Sv[ms_ixs[2]].map() * abc[2];

    /*const*/ auto& X = Xms[vix];


  #define LINEARIZED_BENDING
  #ifdef LINEARIZED_BENDING
    if (!m_dbg.toggle) // for 4dbend comparisong
      if (m_settings.linearized_bending > 0.001f) {
        s.head<3>() = s.head<3>() - m_settings.linearized_bending * 2 * X.h * s.tail<3>();
        s.tail<3>().setZero();
      }
    if(m_settings.svdclamp > 0)
      hermitian_eig_clamp(s, m_settings.svdclamp);
    float sx = std::sqrt(s[0]) - 1;
    float sy = std::sqrt(s[2]) - 1;
    float sa = s[1]/((sx+1)*(sy+1));
    s.head<3>() << sx,sa,sy;
  #else
    //...
  #endif

    Vector4s g;

    #ifdef MODEL4D
    if (m_dbg.toggle) {
      g = m_model4D->deformation(s, X.pix);
    } else {
      g = m_model->deformation(s, X.pix);
    }
    #else
    g = m_model->deformation(s, X.pix);
    #endif

    g *= m_settings.deform_reference;
    // store deformed ms coordinates intm. in ws coords
    // NOTE: buffer.row is a colvector so it expects colvector comma init
    // Xws.row<float, 7>(vix) << Xms.row(vix).transpose() + g,
    //     Xms.block<1, 2>(vix, 0).transpose(), 1;

    if(!dbg_compare_keep_uv)
      Xws.row<float, 11>(vix) << X.mapXT() + g, X.a,
          // 0,0,0,
          X.mapB(),
          // TB[vix].mapB(),
          // dbg0,dbg1, 1;
          // (s[0] + 0.5f),(s[2] + 0.5f), 1;
          X.mapXT().head<2>(), 1;
    else
      Xws.row<float, 11>(vix) << X.mapXT() + g, X.a, X.mapB(), Xws.cpu()[vix].u, Xws.cpu()[vix].v, 1;

        
  });
}

void YarnMapper::shell_map(const Mesh& mesh) {
  // x = phi(xi1, xi2) + h n(xi1, xi2); where xi1 xi2 h are pre-deformed
  // reference coordinates

  auto& U     = mesh.U.cpu();
  auto& F     = mesh.F.cpu();
  auto& Fms     = mesh.Fms.cpu();
  auto& defF  = mesh.defF.cpu();
  auto& defFv = mesh.vertex_defF.cpu();
  auto& N     = mesh.normals.cpu();
  auto& Nv    = mesh.vertex_normals.cpu();
  auto& B     = m_soup.get_B().cpu();
  int nverts  = m_soup.num_vertices();
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
    auto ms_ixs = mesh.Fms.cpu()[tri].map();
    n = Nv[ms_ixs[0]].map() * abc[0] + Nv[ms_ixs[1]].map() * abc[1] +
        Nv[ms_ixs[2]].map() * abc[2];
    n.normalize();

    auto ws_ixs = F[tri].map();
    Vector3s phi = mesh.X.row<float, 3>(ws_ixs[0]) * abc[0] +
                   mesh.X.row<float, 3>(ws_ixs[1]) * abc[1] +
                   mesh.X.row<float, 3>(ws_ixs[2]) * abc[2];
    ;

    // phong deformation:
    // phi = sum_i b_i x_i + alpha * sum_i b_i F_i (X-Xi)
    // gradphi = F_tri + alpha * sum_i (b_i F_i + F_i(X-Xi) gradb_i^T)
    if (m_settings.phong_deformation > 0) {
      for (size_t i = 0; i < 3; ++i) {
        // Debug::log(ms_ixs[i],defFv.size());
        phi += m_settings.phong_deformation *
               (abc[i] * (defFv[ms_ixs[i]].map() *
                          (x_ws.head<2>() - U[ms_ixs[i]].map())));
      }
    }

    {
      // TODO wrong? not including transformation from uvh0 to uvhdeformed?
      // approximate. or using new barycoords works fine maybe..
      MatrixNMs<3, 3> Q;  // transformation from uvh to xyz
      Q.setZero();
      Q.block<3, 2>(0, 0) = defF[tri].map();
      // TODO phong gradient (without it it's just using the flat triangle &
      // interpolated normal, which is probably a decent approximation (since we
      // orthogonalize in the shader) until cheap-transformed refd1 becomes
      // colinear with phongtransformed tangent).

      Q.col(2) = n;  // usually not orthogonal, would need to take transpose
                     // inverse of Q! TODO discuss

      auto _xws = m_soup.get_Xws().row<float, 11>(vix);  // TODO def above
      _xws.block<3, 1>(5, 0) =
          Q.inverse().transpose() * _xws.block<3, 1>(5, 0);  // b = Q** * b
    }

    x_ws.head<3>() = phi + h * n;

    // NOTE transforming using surface (~phong grad and normal) at vertex
    // position although location in middle of edge would be more correct.
    // assumption that surface varies little compared to segments, and that
    // later reorthonormalization hides issues anyway.
  });
}

// TODO this takes forever for big garments. probably better to instead prepare
// appropriate buffers and copy this stuff into a compute shader, though maybe
// the bottleneck is the disk write.. also if doing this on cpu then have to get
// Xws back from gpu!
bool YarnMapper::export2fbx(const std::string& filename) {
  if (!m_initialized)
    return false;

  if (m_settings.gpu_compute) {
    Debug::warning("yarn fbx export only works in cpu-mode.");
    return false;
  }


  static std::vector<uint32_t>
      vcids;  // id of circle of vertices, for all non-tip midline vertices
  static std::vector<uint32_t>
      fcids;  // id of circle of faces, for all non-tips exluding also the
              // second to last midline vertex

  auto LIM      = std::numeric_limits<uint32_t>::max();
  const auto& I = m_soup.getIndexBuffer()
                      .cpu();  // midline vertex indices of LIM separated yarns

  // on time setup of vcids and fcids
  static uint32_t nvcid;  // number of vertex ids
  static uint32_t nfcid;  // number of face ids
  if (vcids.size() == 0) {
    nvcid = 0;
    nfcid = 0;
    vcids.reserve(I.size());
    fcids.reserve(I.size());
    for (size_t i = 0; i < I.size(); ++i) {
      if (i == 0 || i == I.size() - 1 ||
          I[i] == LIM) {  // ignore outside verts or tips
        vcids.push_back(LIM);
      } else {
        if (I[i - 1] == LIM || I[i + 1] == LIM)  // LIM-delimited tip
          vcids.push_back(LIM);
        else
          vcids.push_back(nvcid++);  // accepted location with existing neighbor
                                     // left and right
      }

      // for faces same as above but one additional skip for second to last vert
      // per yarn
      if (i == 0 || i >= I.size() - 2 || I[i] == LIM) {
        fcids.push_back(LIM);
      } else {
        if (I[i - 1] == LIM || I[i + 1] == LIM || I[i + 2] == LIM)
          fcids.push_back(LIM);
        else
          fcids.push_back(nfcid++);  // accepted location with existing neighbor
                                     // left and right
      }
    }
  }

  int NSEGS     = 8;  // number of cylinder sides
  float invSEG  = 1.0f / (NSEGS);
  int NVERTICES = NSEGS + 1;            // + 1 bc of radial seam
  float R       = m_model->getPYP().r;  // yarn radius
  float normalTwist =
      1.0f;  // TODO normal map params from app. as func params here?
  float normalNum = 4.0f;
  float normalLen = 6.0f;

  // allocate export data
  uint32_t num_total_verts = NVERTICES * nvcid;
  uint32_t num_total_faces = 2 * NSEGS * nfcid;
  FBXExportData data;
  data.vertices.resize(num_total_verts * 3);
  data.uv0.resize(num_total_verts * 2);
  data.uv1.resize(num_total_verts * 2);
  data.faces.resize(num_total_faces * 3);
  data.uvfaces.resize(num_total_faces * 3);

  // mediocre solution of copying back from gpu
  // TODO doesnt seem to work (always)? but in cpu mode export works?
  auto& Xws = m_soup.get_Xws().cpu();
  if (m_settings.gpu_compute) {
    auto data = m_soup.get_Xws().getGPUData();
    Xws.resize(data.size());
    threadutils::parallel_for(size_t(0), data.size(), [&](size_t i) {
      Xws[i] = data[i];
    });
  }
  // Xws = xyzt | arc | d1x d1y d1z | u v | rlocal
  threadutils::parallel_for(size_t(0), I.size(), [&](size_t i) {
    uint32_t vcid = vcids[i];
    if (vcid == LIM)
      return;

    // [vert_0] ---edge_A--- [vert_1] ---edge_B--- [vert_2]

    // indices of left, self and right yarn segment
    uint32_t left  = I[i - 1];
    uint32_t self  = I[i];
    uint32_t right = I[i + 1];
    auto& x0       = Xws[left];
    auto& x1       = Xws[self];
    auto& x2       = Xws[right];
    // DEBUG NOTE sometimes sort of randomly it failed here because it tries to access Xws[right] with right being LIM... this shouldnt happen?

    // edge rest length
    float lA = x1.a - x0.a;
    float lB = x2.a - x1.a;

    Vector3s tA = x1.mapX() - x0.mapX();
    Vector3s tB = x2.mapX() - x1.mapX();

    Vector3s tv = (lA * tA + lB * tB);                // vertex tangent
    Vector3s nv = (lA * x0.mapD() + lB * x1.mapD());  // vertex normal
    nv          = (nv - tv * tv.dot(nv) / tv.squaredNorm())
             .normalized();                   // orthonormalize
    Vector3s bv = tv.cross(nv).normalized();  // vertex binormal

    // TODO edge current length, and scale radius
    float rA = 1 * R;  // TODO vol.preserve
    float rB = 1 * R;  // TODO vol.preserve
    float r1 = (lA * rA + lB * rB) / (lA + lB);

    // if (i == 66) {
    //   // Debug::log(x0.mapX());
    //   // Debug::log(x1.mapX());
    //   // Debug::log(x2.mapX());
    //   // Debug::log(tv);
    //   // Debug::log(nv);
    //   // Debug::log(bv);
    // }

    for (int cs = 0; cs < NSEGS + 1;
         cs++) {  // iterate circular vertices, including duplicate end
      float alpha = cs * invSEG;  // 0 to 1, inclusive

      float a  = alpha * 2.0f * float(M_PI);
      float ca = std::cos(a);
      float sa = std::sin(a);
      // TODO edge twist theta ~> get average theta and shift a ?

      Vector3s n = ca * nv + sa * bv;
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

      uint32_t vertix               = vcid * NVERTICES + cs;
      data.vertices[3 * vertix + 0] = p[0];
      data.vertices[3 * vertix + 1] = p[1];
      data.vertices[3 * vertix + 2] = p[2];
      data.uv0[2 * vertix + 0]      = x1.u;
      data.uv0[2 * vertix + 1]      = x1.v;
      data.uv1[2 * vertix + 0]      = normalNum * (-alpha - normalTwist * x1.a / (r1 * normalLen)); // x1.a
      data.uv1[2 * vertix + 1]      = x1.a / (r1 * normalLen);
          // normalNum * (-alpha - 0.1591549f * normalTwist * x1.a / r1);
    }

    uint32_t fcid = fcids[i];
    if (fcid == LIM)
      return;

    uint32_t vcid_next = vcids[i + 1];  // assert == vcid+1?

    for (int cs = 0; cs < NSEGS; cs++) {  // iterate circular segments
      // in 2D pseudo indices, with i radial and f the circle/crosssection
      // index: { [f,i], [f+1,i+1], [f+1,i] } { [f,i], [f,i+1], [f+1,i+1] } note
      // that we model a radial seam explicitly, so cs+1 makes sense even for
      // the cs==NSEGS-1

      uint32_t facedofA = 3 * (fcid * 2 * NSEGS + 2 * cs + 0);
      uint32_t facedofB = 3 * (fcid * 2 * NSEGS + 2 * cs + 1);
      uint32_t vertix0  = vcid * NVERTICES + cs;       // "[f,cs]"
      uint32_t vertix1  = vcid_next * NVERTICES + cs;  // "[f+1,cs]"

      auto F0 = std::vector<uint32_t>{vertix0, vertix1 + 1, vertix1};
      // auto F0 = std::vector<uint32_t>{0,1,2};
      auto F1 = std::vector<uint32_t>{vertix0, vertix0 + 1, vertix1 + 1};

      // last poly index negated with ~
      data.faces[facedofA + 0] = F0[0];
      data.faces[facedofA + 1] = F0[1];
      data.faces[facedofA + 2] = ~int32_t(F0[2]);
      data.faces[facedofB + 0] = F1[0];
      data.faces[facedofB + 1] = F1[1];
      data.faces[facedofB + 2] = ~int32_t(F1[2]);

      data.uvfaces[facedofA + 0] = F0[0];
      data.uvfaces[facedofA + 1] = F0[1];
      data.uvfaces[facedofA + 2] = F0[2];
      data.uvfaces[facedofB + 0] = F1[0];
      data.uvfaces[facedofB + 1] = F1[1];
      data.uvfaces[facedofB + 2] = F1[2];
    }
  });

  return export_fbx(filename, data);
}

// TODO attempt to export X & U & U2, F (or other split) using cnpy
// and then load into blender python - using the fast setter funcs?
// it would be much better to export USD (but i cant even get the library
// running from vcpkg) and second best to export alembic (but i cant even find
// any code on that ...)


bool YarnMapper::export2fbx_cloth(const std::string& filename) {
  if (!m_initialized)
    return false;

  Mesh& mesh = m_meshProvider->getMesh();
  auto& X = mesh.X.cpu();
  auto& U = mesh.U.cpu();
  auto& F = mesh.F.cpu();
  auto& Fms = mesh.Fms.cpu();
  
  FBXExportData data;
  data.vertices.resize(X.size() * 3);
  data.uv0.resize(U.size() * 2);
  data.uv1.resize(U.size() * 2);
  data.faces.resize(F.size() * 3);
  data.uvfaces.resize(Fms.size() * 3);

  // Xws = xyzt | arc | d1x d1y d1z | u v | rlocal
  threadutils::parallel_for(size_t(0), X.size(), [&](size_t i) {
    data.vertices[3 * i + 0] = X[i].x;
    data.vertices[3 * i + 1] = X[i].y;
    data.vertices[3 * i + 2] = X[i].z;
  });
  threadutils::parallel_for(size_t(0), U.size(), [&](size_t i) {
    data.uv0[2 * i + 0]      = U[i].u;
    data.uv0[2 * i + 1]      = U[i].v;
    data.uv1[2 * i + 0]      = U[i].u;
    data.uv1[2 * i + 1]      = U[i].v;
  });
  threadutils::parallel_for(size_t(0), F.size(), [&](size_t i) {
    data.faces[3 * i + 0] = F[i].v0;
    data.faces[3 * i + 1] = F[i].v1;
    data.faces[3 * i + 2] = ~int32_t(F[i].v2);
    data.uvfaces[3 * i + 0] = Fms[i].v0;
    data.uvfaces[3 * i + 1] = Fms[i].v1;
    data.uvfaces[3 * i + 2] = Fms[i].v2;
  });

  return export_fbx(filename, data);
}

void YarnMapper::dbg_compare_bending() { // copy of deform_reference using both bending models to store new UVs
#ifdef MODEL4D
  Mesh& mesh = m_meshProvider->getMesh();
  int nverts = m_soup.num_vertices();
  auto& Xms  = m_soup.get_Xms().cpu();
  auto& Xws  = m_soup.get_Xws();
  auto& S    = mesh.strains.cpu();
  auto& Sv   = mesh.vertex_strains.cpu();
  
  auto& B0 = m_soup.get_B0().cpu();
  Xws.cpu().resize(Xms.size());
  threadutils::parallel_for(0, nverts, [&](int vix) {
    const auto& bary = B0[vix];
    int tri          = bary.tri;
    if (tri < 0)  // skip unassigned vertex
      return;

    Vector3s abc;
    abc << bary.a, bary.b, bary.c;

    Vector6s s;
    {
      auto ms_ixs = mesh.Fms.cpu()[tri].map();
      s = Sv[ms_ixs[0]].map() * abc[0] + Sv[ms_ixs[1]].map() * abc[1] +
          Sv[ms_ixs[2]].map() * abc[2];
    }

    /*const*/ auto& X = Xms[vix];

    float dbg0,dbg1;

    Vector4s g4D;
    {
      Vector6s s4D = s;
      if(m_settings.svdclamp > 0)
        hermitian_eig_clamp(s4D, m_settings.svdclamp);
      float sx = std::sqrt(s4D[0]) - 1;
      float sy = std::sqrt(s4D[2]) - 1;
      float sa = s4D[1]/((sx+1)*(sy+1));
      s4D.head<3>() << sx,sa,sy;
      g4D = m_model4D->deformation(s4D, X.pix);
    }

    Vector4s gLin;
    {
      Vector6s sLin = s;
      if (m_settings.linearized_bending > 0.001f) {
        sLin.head<3>() = sLin.head<3>() - m_settings.linearized_bending * 2 * X.h * sLin.tail<3>();
        sLin.tail<3>().setZero();
      }
      if(m_settings.svdclamp > 0)
        hermitian_eig_clamp(sLin, m_settings.svdclamp);
      float sx = std::sqrt(sLin[0]) - 1;
      float sy = std::sqrt(sLin[2]) - 1;
      float sa = sLin[1]/((sx+1)*(sy+1));
      sLin.head<3>() << sx,sa,sy;
      gLin = m_model->deformation(sLin, X.pix);
    }

    // pos error relative to radius
    // twist error relative to one full twist of 2pi
    float errx = (g4D.head<3>()-gLin.head<3>()).norm() / m_model->getPYP().r;
    float errt = std::abs(g4D(3)-gLin(3)) / (2*3.1415926f);


    Xws.cpu()[vix].u = errx;
    Xws.cpu()[vix].v = errt;
  });

  dbg_compare_keep_uv = true;
  #endif
}

void YarnMapper::dbg_compare_nsamples() { // copy of deform_reference current model and other fixed model

  // tmp load of comparison model
  std::string comparemodelfolder;
  if (m_settings.modelfolder.find("rib") != std::string::npos) {
    comparemodelfolder = "data/yarnmodels/num_samples/model_rib_31";
  } else {
    comparemodelfolder = "data/yarnmodels/num_samples/model_stock_31";
  }
  auto model2 = std::make_unique<Model>(comparemodelfolder);
  Debug::logf("Comparing\n  %s\n  %s\n",m_settings.modelfolder.c_str(), comparemodelfolder.c_str());

  Mesh& mesh = m_meshProvider->getMesh();
  int nverts = m_soup.num_vertices();
  auto& Xms  = m_soup.get_Xms().cpu();
  auto& Xws  = m_soup.get_Xws();
  auto& S    = mesh.strains.cpu();
  auto& Sv   = mesh.vertex_strains.cpu();
  
  auto& B0 = m_soup.get_B0().cpu();
  Xws.cpu().resize(Xms.size());
  threadutils::parallel_for(0, nverts, [&](int vix) {
    const auto& bary = B0[vix];
    int tri          = bary.tri;
    if (tri < 0)  // skip unassigned vertex
      return;

    Vector3s abc;
    abc << bary.a, bary.b, bary.c;

    Vector6s s;
    {
      auto ms_ixs = mesh.Fms.cpu()[tri].map();
      s = Sv[ms_ixs[0]].map() * abc[0] + Sv[ms_ixs[1]].map() * abc[1] +
          Sv[ms_ixs[2]].map() * abc[2];
    }

    /*const*/ auto& X = Xms[vix];

    if (m_settings.linearized_bending > 0.001f) {
      s.head<3>() = s.head<3>() - m_settings.linearized_bending * 2 * X.h * s.tail<3>();
      s.tail<3>().setZero();
    }
    if(m_settings.svdclamp > 0)
      hermitian_eig_clamp(s, m_settings.svdclamp);
    float sx = std::sqrt(s[0]) - 1;
    float sy = std::sqrt(s[2]) - 1;
    float sa = s[1]/((sx+1)*(sy+1));
    s.head<3>() << sx,sa,sy;

    Vector4s g1 = m_model->deformation(s, X.pix);
    Vector4s g2 = model2->deformation(s, X.pix);

    // pos error relative to radius
    // twist error relative to one full twist of 2pi
    float errx = (g1.head<3>()-g2.head<3>()).norm() / m_model->getPYP().r;
    float errt = std::abs(g1(3)-g2(3)) / (2*3.1415926f);


    Xws.cpu()[vix].u = errx;
    Xws.cpu()[vix].v = errt;
  });

  dbg_compare_keep_uv = true;
}