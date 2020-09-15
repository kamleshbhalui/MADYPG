#ifndef __YARNMAPPER__H__
#define __YARNMAPPER__H__

#include "../ThreadUtils.h"
// #include "../mesh/AbstractMeshSimulation.h"
#include "../mesh/TestMeshSimulation.h"
#include <numeric> // iota
#include <iostream>
#include <limits>
#include <memory>

class YarnMapper
{
public:
  int boost = 1;
  float rescale = 0.002f;//1.0f;//0.002f;
  float S = boost * 200.0f * rescale;
  float A = 0.75f * rescale;
  float f = 0.3333f / rescale;
  float DZ = 0.0f;
  int Ny = boost * boost * int(100 * 200 / S * rescale);
  int Nsub = 400;
  float wave = 1.0f;

  YarnMapper()
  {
    T = 0.0f;

    m_meshSimulation = std::static_pointer_cast<AbstractMeshSimulation>(std::make_shared<TestMeshSimulation>());//DEBUG

    I.reserve(Ny * (Nsub + 1 + 2)); // indices per yarn with Nsub verts, because of 1x primitiverestart and 2x adjacency
    uint32_t ix = 0;

    X.resize(2 * Ny * Nsub, 3);

    for (int j = 0; j < Ny; j++)
    {
      I.push_back(ix); // tip, self-adjacency
      for (int i = 0; i < Nsub; i++)
        I.push_back(ix++);
      I.push_back(ix - 1); // tip, self-adjacency
      I.push_back(std::numeric_limits<GLuint>::max());
    }
    for (int j = 0; j < Ny; j++)
    {
      I.push_back(ix); // tip, self-adjacency
      for (int i = 0; i < Nsub; i++)
        I.push_back(ix++);
      I.push_back(ix - 1); // tip, self-adjacency
      if (j < Ny - 1)
        I.push_back(std::numeric_limits<GLuint>::max());
    }

    printf("Test yarns with %d vertices\n", int(X.rows()));
    setPositions(T);
  }

  ~YarnMapper() {}

  void setPositions(float T)
  {
    float invS = 1.0f / S;

    // threadutils::parallel_for(0, Ny, [&](int j){});

    threadutils::parallel_for(0, Ny, [&](int j) {
      // for (int j = 0; j < Ny; j++) {
      float a1 = float(j) / (Ny - 1);
      for (int i = 0; i < Nsub; i++)
      {
        float a2 = float(i) / (Nsub - 1);
        Eigen::Vector3f v;
        v << (1 - a2) * (-0.5f * S) + a2 * 0.5f * S,
            (1 - a1) * (-0.5f * S) + a1 * 0.5f * S,
            DZ + A * std::sin(1.57f + j * 3.14f + a2 * f * S * 3.14f);
        v *= 0.75;
        Eigen::Vector3f dv;
        dv << 0.0f, 0.0f,
            5.0f * std::sin(v(0) * invS * 10.0f + T * 2 * 3.14f) + 2.5f * std::cos(0.5f + v(1) * invS * 10.0f + 0.333f * T * 2 * 3.14f) + 8.0f * std::cos((v(0) * invS + 2 * v(1) * invS + v(0) * v(1) * invS * invS) * 2.0f + T * 3.14f);
        X.row(j * Nsub + i) = v + rescale*dv*wave;
      }
      // }
    });
    threadutils::parallel_for(0, Ny, [&](int j) {
      // for (int j = 0; j < Ny; j++) {
      float a1 = float(j) / (Ny - 1);

      for (int i = 0; i < Nsub; i++)
      {
        float a2 = float(i) / (Nsub - 1);
        Eigen::Vector3f v;
        v << (1 - a1) * (-0.5f * S) + a1 * 0.5f * S,
            (1 - a2) * (-0.5f * S) + a2 * 0.5f * S,
            DZ + A * std::sin(-1.57f + j * 3.14f + a2 * f * S * 3.14f);
        v *= 0.75;
        Eigen::Vector3f dv;
        dv << 0.0f, 0.0f,
            5.0f * std::sin(v(0) * invS * 10.0f + T * 2 * 3.14f) + 2.5f * std::cos(0.5f + v(1) * invS * 10.0f + 0.333f * T * 2 * 3.14f) + 8.0f * std::cos((v(0) * invS + 2 * v(1) * invS + v(0) * v(1) * invS * invS) * 2.0f + T * 3.14f);
        X.row((Ny + j) * Nsub + i) = v + rescale*dv*wave;
      }
      // }
    });
  }

  void step()
  {
    // return;
    T += 0.01f;
    m_meshSimulation->update();
    setPositions(T);
  }

  const std::vector<uint32_t> &getIndices() { return I; }
  const MatrixGLf &getVertexData() { return X; }
  const std::shared_ptr<AbstractMeshSimulation> getMeshSimulation() { return m_meshSimulation; }

private:
  // yarn stuff
  MatrixGLf X;
  std::vector<uint32_t> I;

  // simulation mesh stuff
  std::shared_ptr<AbstractMeshSimulation> m_meshSimulation;

  // debug stuff
  float T;
};

#endif // __YARNMAPPER__H__