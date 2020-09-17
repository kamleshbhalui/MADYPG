#ifndef __YARNMAPPER__H__
#define __YARNMAPPER__H__

#include "../utils/debug_includes.h"
#include "../utils/threadutils.h"
// #include "../mesh/AbstractMeshSimulation.h"
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>  // iota

#include "../mesh/TestMeshSimulation.h"
#include "Grid.h"
#include "PeriodicYarnPattern.h"
#include "YarnSoup.h"

class YarnMapper {
 public:
  YarnMapper(const std::string& modelfolder) {
    Debug::Timer timer;

    // set up mesh provider
    // TODO rename interface to provider?
    // DEBUG TODO ANIMATE
    m_meshSimulation = std::static_pointer_cast<AbstractMeshSimulation>(
        std::make_shared<TestMeshSimulation>());

    // Load model / pyp
    m_pyp.deserialize(modelfolder + "/pyp");  // DEBUG hardcoded file
    m_pyp.rectangulize();  // TODO potentially assume that this is true for new
    // pyp, but it should take only a millisecond anyway

    Debug::logf("Timer: startup %.3f ms\n", timer.tock<microseconds>() * 0.001);

    Mesh& mesh = m_meshSimulation->getMesh();
    m_grid.fromTiling(mesh, m_pyp);
    m_grid.overlap_triangles(mesh);

    Debug::logf("Timer: grid setup %.3f ms\n",
                timer.tock<microseconds>() * 0.001);

    // ... soup
    m_soup.fill_from_grid(m_pyp, m_grid);

    // fancy: do some random uv displacement with multi-level 3D displacement
    // noise (sth like n levels with n strength values) ?
    // actually might be better to do in shader using a texture: meshuv->noise,
    // actually no. still want to do that in uv space (and somehow account for
    // yarns being pushed outside of uvmesh: tribary choose closest tri)

    Debug::logf("Timer: soup tiling %.3f ms\n",
                timer.tock<microseconds>() * 0.001);

    // TODO cut outside / prune array (OR SHOULD I, NAIVE JUST KEEP MATRICES SAME SIZE AND SIMPLY NOT MAKE YARNLIST FOR OUTSIDE(tri=-1)) / ...
    mesh.compute_invDm();
    m_soup.assign_triangles(m_grid, mesh); // TODO consider split into subfunctions
    //. soup.cut_outside

    Debug::logf("Timer: mesh invdm & soup triangles %.3f ms\n",
                timer.tock<microseconds>() * 0.001);
    m_soup.generate_index_list(); // TODO prune by length while assembling! @ first parallel count: also sum up RL and prune (by not adding their indices and also deleting their edges), or just prune before assembly but that means I need to pass through yarns before.. (i guess if this happens once it doesnt matter, and its nicer however the code is more readable!!)

    Debug::logf("Timer: soup index list %.3f ms\n",
                timer.tock<microseconds>() * 0.001);

    // NOTE about soup matrices: keep as 'good' matrix only the things needed for gpu, the rest of the stuff could be joint into a vector of structs of vertexdata ? dep on how its used ?
    // . shader might want per vertex:
    //    x y z (t)
    //    u_mesh v_mesh ; for meshspace texturing
    //    parametric_t(=acc.rest length) for yarnspace texturing/twist
    //    geom shader promotes and additionally produces parametric_c (around circle)
    

    // TODO  ms change / ... / shell-map / ...
    // mesh.v2f (NEW!)
    // mesh.facenormals
    // mesh.face2vertex(normals)
    // thisclasshere.shellmap


    // TODO load model (model class <--- pyp) / ws change / ... / masm / ...
    // mesh.strains_face
    // mesh.face2vertex(strains)
    // thisclasshere.masm
    // soup.assign_triangles again (at least redo bary)

    // Debug::logf("Timer: yarnmapper setup %.3f ms\n",
    //             timer.tock<microseconds>() * 0.001);
  }

  ~YarnMapper() {}

  void step() { m_meshSimulation->update(); }

  const std::vector<uint32_t>& getIndices() const {
    return m_soup.getIndices();
  }
  const MatrixGLf& getVertexData() const { return m_soup.getVertexData(); }
  const std::shared_ptr<AbstractMeshSimulation> getMeshSimulation() {
    return m_meshSimulation;
  }
  float getRadius() { return m_pyp.r; }

 private:
  PeriodicYarnPattern m_pyp;
  Grid m_grid;
  YarnSoup m_soup;

  // simulation mesh stuff
  std::shared_ptr<AbstractMeshSimulation> m_meshSimulation;

  // yarn stuff
  std::vector<uint32_t> I;
};

#endif  // __YARNMAPPER__H__