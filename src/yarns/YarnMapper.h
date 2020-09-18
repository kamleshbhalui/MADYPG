#ifndef __YARNMAPPER__H__
#define __YARNMAPPER__H__

#include "../utils/debug_includes.h"
#include "../utils/threadutils.h"
// #include "../mesh/AbstractMeshProvider.h"
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>  // iota

#include "../mesh/ObjSeqAnimation.h"
#include "Grid.h"
#include "PeriodicYarnPattern.h"
#include "YarnSoup.h"

class YarnMapper {
 public:
  struct Settings {
    std::string modelfolder = "";
    bool flat_normals       = false;
    bool shepard_weights    = true;
    // TODO objmeshprovider settings (and use those in its constructor)
    // TODO xpbdmeshprovider settings
  } m_settings;

  YarnMapper() : m_initialized(false) {}

  void initialize() { step(); }

  ~YarnMapper() {}

  void step() {
    Debug::Timer timer_2;
    Debug::Timer timer;

    if (m_initialized) {
      // step the mesh provider
      m_meshProvider->update();
      Debug::logf("Timer: meshProvider update %.3f ms\n",
                  timer.tock<microseconds>() * 0.001);
    } else {
      // set up mesh provider
      m_meshProvider = std::static_pointer_cast<AbstractMeshProvider>(
          std::make_shared<ObjSeqAnimation>("presim/test", true, true));

      // Load model / pyp
      m_pyp.deserialize(m_settings.modelfolder +
                        "/pyp");  // DEBUG hardcoded file
      m_pyp
          .rectangulize();  // TODO potentially assume that this is true for new
      // pyp, but it should take only a millisecond anyway

      Debug::logf("Timer: meshProvider & pyp init %.3f ms\n",
                  timer.tock<microseconds>() * 0.001);
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
      m_grid.fromTiling(mesh, m_pyp);
      m_grid.overlap_triangles(mesh);

      Debug::logf("Timer: grid setup %.3f ms\n",
                  timer.tock<microseconds>() * 0.001);

      if (!m_initialized) {
        // ... soup
        // NOTE about soup matrices: keep as 'good' matrix only the things
        // needed for gpu, the rest of the stuff could be joint into a vector of
        // structs of vertexdata ? dep on how its used ? . shader might want per
        // vertex:
        //    x y z (t)
        //    u_mesh v_mesh ; for meshspace texturing
        //    parametric_t(=acc.rest length) for yarnspace texturing/twist
        //    geom shader promotes and additionally produces parametric_c
        //    (around circle)

        m_soup.fill_from_grid(m_pyp, m_grid);

        // fancy: do some random uv displacement with multi-level 3D
        // displacement noise (sth like n levels with n strength values) ?
        // actually might be better to do in shader using a texture:
        // meshuv->noise, actually no. still want to do that in uv space (and
        // somehow account for yarns being pushed outside of uvmesh: tribary
        // choose closest tri)

        Debug::logf("Timer: soup tiling %.3f ms\n",
                    timer.tock<microseconds>() * 0.001);
      }

      mesh.compute_invDm();
      m_soup.assign_triangles(m_grid, mesh);
      Debug::logf("Timer: mesh invdm & soup triangles %.3f ms\n",
                  timer.tock<microseconds>() * 0.001);
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

        Debug::logf("Timer: soup cut & index list %.3f ms\n",
                    timer.tock<microseconds>() * 0.001);
      }

      mesh.compute_v2f_map(m_settings.shepard_weights);

      Debug::logf("Timer: ms change %.3f ms\n",
                  timer.tock<microseconds>() * 0.001);

    }  // end of MS change

    // @ WORLD SPACE MESH CHANGES

    mesh.compute_face_normals();  // TODO make optional in case obj file has
                                  // normals
    mesh.compute_vertex_normals();

    // DEBUG in lieu of masm for now just copy soup ms to ws
    // TODO load model (model class <--- pyp) / ws change / ... / masm / ...
    // mesh.strains_face
    // mesh.face2vertex(strains)
    // thisclasshere.masm
    // soup.assign_triangles again (at least redo bary)
    {
      int n     = m_soup.num_vertices();
      auto& Xms = m_soup.get_Xms();
      auto& Xws = m_soup.get_Xws();
      Xws.resize(Xms.rows(), Xms.cols());
      threadutils::parallel_for(0, n, [&](int i) { Xws.row(i) = Xms.row(i); });
    }

    shell_map(mesh, m_settings.flat_normals);

    Debug::logf("Timer: ws change %.3f ms\n",
                timer.tock<microseconds>() * 0.001);

    Debug::logf("Timer all: %.3f ms\n", timer_2.tock<microseconds>() * 0.001);
    m_initialized = true;
  }

  void shell_map(const Mesh& mesh, bool flat_normals = false) {
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

  const std::vector<uint32_t>& getIndices() const {
    return m_soup.getIndices();
  }
  const MatrixGLf& getVertexData() const { return m_soup.get_Xws(); }
  const std::shared_ptr<AbstractMeshProvider> getMeshSimulation() {
    return m_meshProvider;
  }
  float getRadius() { return m_pyp.r; }

  bool initialized() const { return m_initialized; }

 private:
  bool m_initialized;
  PeriodicYarnPattern m_pyp;
  Grid m_grid;
  YarnSoup m_soup;

  // simulation mesh stuff
  std::shared_ptr<AbstractMeshProvider> m_meshProvider;

  // yarn stuff
  std::vector<uint32_t> I;
};

#endif  // __YARNMAPPER__H__