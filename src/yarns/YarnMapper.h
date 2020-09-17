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
  YarnMapper() {
    Debug::Timer timer;

    // set up mesh provider
    // TODO rename interface to provider?
    // DEBUG
    m_meshSimulation = std::static_pointer_cast<AbstractMeshSimulation>(
        std::make_shared<TestMeshSimulation>());

    // Load model / pyp
    // DEBUG hardcoded file
    m_pyp.deserialize("models/model_stock/pyp");
    // m_pyp.rectangulize(); // TODO potentially assume that this is true for new
    // pyp, but it should take only a millisecond anyway

    Debug::logf("Timer: startup %.3f ms\n", timer.tock<microseconds>() * 0.001);

    const Mesh& mesh = m_meshSimulation->getMesh();
    m_grid.fromTiling(mesh, m_pyp);
    m_grid.overlap_triangles(mesh);
    // Grid <--- mesh, pyp
    //   get uv bounds + 0.1 * (px,py) robustness // OPT PARALLEL
    //   build grid from uv bounds
    //   tmp tri2cells: overlap v triangles and grid cells

    Debug::logf("Timer: grid setup %.3f ms\n",
                timer.tock<microseconds>() * 0.001);

    // ... soup
    m_soup.fill_from_grid(m_pyp, m_grid);

    // fancy: do some random uv displacement with multi-level 3D displacement
    // noise (sth like n levels with n strength values) ?
    // actually might be better to do in shader using a texture: meshuv->noise

    // .. cut prune etc ..

    m_soup.generate_index_list();

    Debug::logf("Timer: soup setup %.3f ms\n",
                timer.tock<microseconds>() * 0.001);

    // Debug::logf("Timer: yarnmapper setup %.3f ms\n",
    //             timer.tock<microseconds>() * 0.001);
  }

  ~YarnMapper() {}

  void step() {
    m_meshSimulation->update();
  }

  const std::vector<uint32_t>& getIndices() const { return m_soup.getIndices(); }
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