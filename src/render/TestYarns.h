#ifndef __TESTYARNS__H__
#define __TESTYARNS__H__

#include "YarnInterface.h"
#include <numeric> // iota
#include <iostream>
#include <limits>
#include "ThreadUtils.h"

class TestYarns : public YarnInterface
{
public:
  int boost = 1;
  float S = boost * 200.0f;
  float A = 0.75f;
  float f = 0.3333f;
  float DZ = 0.0f;
  int Ny = boost * boost * int(100 * 200 / S);
  int Nsub = 400;

  TestYarns()
  {
    T = 0.0f;

    I.reserve(Ny * (Nsub + 1 + 2)); // indices per yarn with Nsub verts, because of 1x primitiverestart and 2x adjacency
    uint32_t ix = 0;

    X.resize(2 * Ny * Nsub, 3);

    for (int j = 0; j < Ny; j++)
    {
      I.push_back(ix); // tip, self-adjacency
      for (int i = 0; i < Nsub; i++)
        I.push_back(ix++);
      I.push_back(ix-1); // tip, self-adjacency
      I.push_back(std::numeric_limits<GLuint>::max());
    }
    for (int j = 0; j < Ny; j++)
    {
      I.push_back(ix); // tip, self-adjacency
      for (int i = 0; i < Nsub; i++)
        I.push_back(ix++);
      I.push_back(ix-1); // tip, self-adjacency
      if (j < Ny - 1)
        I.push_back(std::numeric_limits<GLuint>::max());
    }

    printf("Test yarns with %d vertices\n", int(X.rows()));
    setPositions(T);
  }

  ~TestYarns() {}

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
        Eigen::Vector3f dv;
        dv << 0.0f, 0.0f,
            5.0f * std::sin(v(0) * invS * 10.0f + T * 2 * 3.14f) + 2.5f * std::cos(0.5f + v(1) * invS * 10.0f + 0.333f * T * 2 * 3.14f) + 8.0f * std::cos((v(0) * invS + 2 * v(1) * invS + v(0) * v(1) * invS * invS) * 2.0f + T * 3.14f);
        X.row(j * Nsub + i) = v + dv;
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
        Eigen::Vector3f dv;
        dv << 0.0f, 0.0f,
            5.0f * std::sin(v(0) * invS * 10.0f + T * 2 * 3.14f) + 2.5f * std::cos(0.5f + v(1) * invS * 10.0f + 0.333f * T * 2 * 3.14f) + 8.0f * std::cos((v(0) * invS + 2 * v(1) * invS + v(0) * v(1) * invS * invS) * 2.0f + T * 3.14f);
        X.row((Ny + j) * Nsub + i) = v + dv;
      }
      // }
    });
  }

  void step()
  {
    // return;
    T += 0.01f;
    setPositions(T);
  }

  const std::vector<uint32_t> &getIndices() { return I; }
  // const Eigen::Matrix<float,-1,-1> & getVertexData() { return X; }
  const Eigen::Matrix<float, -1, -1, Eigen::RowMajor> &getVertexData() { return X; }

  // const MatrixXXs & getVertexData() { return X; }
  // const Eigen::Matrix<scalar,Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & getVertexData() { return X; }

private:
  float T;
  // Eigen::Matrix<float,-1,-1> X;
  Eigen::Matrix<float, -1, -1, Eigen::RowMajor> X;
  // MatrixXXs  X;
  // Eigen::Matrix<scalar,Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  X;
  std::vector<unsigned int> I;
};

#endif // __TESTYARNS__H__

/*
_lines.clear();
        int N = 800;
        int Nv = 2400;
        float S = 350.0f;
        float amp = 2.0f;
        for (int j = 0; j < N; j++)
        {
          float A = float(j) / (N - 1);

          {
            float x = (1 - A) * -S + A * S;
            Vector3 d = Vector3(1.0, 0.0, 0.0);

            std::vector<Vector3> vertices;
            for (size_t i = 0; i < Nv; i++)
            {
              float a = float(i) / (Nv-1);

              vertices.push_back(((1 - a) * -S + a * S) * d + Vector3(0.0f, x, 40.0f + amp*std::sin(j*3.14f + a * 3.14f * 50.0f)));
            }
            _lines.emplace_back(_yarnGeometryShader);
            _lines.back().setVertices(vertices, true);
          }
          {
            float x = (1 - A) * -S + A * S;
            Vector3 d = Vector3(0.0, 1.0, 0.0);

            std::vector<Vector3> vertices;
            for (size_t i = 0; i < Nv; i++)
            {
              float a = float(i) / (Nv-1);

              vertices.push_back(((1 - a) * -S + a * S) * d + Vector3(x, 0.0f, 40.0f + amp*std::sin(3.14f * 0.33f + j*3.14f + a * 3.14f * 50.0f)));
            }
            _lines.emplace_back(_yarnGeometryShader);
            _lines.back().setVertices(vertices, true);
          }
        }
*/