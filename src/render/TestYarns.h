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
    float S = boost*200.0f;
    float A = 0.75f;
    float f = 0.3333f;
    float DZ = 40.0f;
    int Ny = boost*boost*int(100 * 200/S);
    int Nsub = 800;
    
    TestYarns()
    {
        T = 0.0f;

        I.reserve(Ny * (Nsub+1));
        uint32_t ix = 0;

        for (int j = 0; j < Ny; j++)
        {
            float a1 = float(j)/(Ny-1);

            int prev_size = X.rows();
            X.conservativeResize(X.rows() + Nsub, 3);
            for (int i = 0; i < Nsub; i++)
            {
              float a2 = float(i) / (Nsub-1);
              Eigen::Vector3f v;
            //   Vector3s v;
              v << (1-a2)*(-0.5f * S)+a2*0.5f*S,
                   (1-a1)*(-0.5f * S)+a1*0.5f*S,
                   DZ + A * std::sin(1.57f + j * 3.14f + a2 * f * S * 3.14f);
              X.row(prev_size + i) = v;
              I.push_back(ix++);
            }
            I.push_back(std::numeric_limits<GLuint>::max());
        }
        for (int j = 0; j < Ny; j++)
        {
            float a1 = float(j)/(Ny-1);

            int prev_size = X.rows();
            X.conservativeResize(X.rows() + Nsub, 3);
            for (int i = 0; i < Nsub; i++)
            {
              float a2 = float(i) / (Nsub-1);
              Eigen::Vector3f v;
            //   Vector3s v;
              v << (1-a1)*(-0.5f * S)+a1*0.5f*S,
                   (1-a2)*(-0.5f * S)+a2*0.5f*S,
                   DZ + A * std::sin(-1.57f + j * 3.14f + a2 * f * S * 3.14f);
              X.row(prev_size + i) = v;
              I.push_back(ix++);
            }
            if (j < Ny -1)
                I.push_back(std::numeric_limits<GLuint>::max());
        }

        // int N = 40000;
        // X.resize(N,3);
        // I.reserve(N);
        // uint32_t ix = 0;
        // for (int i = 0; i < N; i++)
        // {
        //     float a = float(i)/(N-1);
        //     X.row(i) << (1-a)*-100.0f+a*100.0f,0.0f,40.0f;
        //     I.push_back(ix++);
        // }
        // std::cout<<"Final "<<ix<<"\n";
        // std::cout<<"Final "<<X.row(N-1)<<"\n";

        printf("Test yarns with %d vertices\n", int(X.rows()));
        // printf("Test yarns with %d indices\n", int(I.size()));
        // std::cout<<I.back()<<"\n";
    }

    ~TestYarns() {}

    void step()
    {
        // return;
        T += 0.01f;


        // threadutils::parallel_for(0, Ny, [&](int j){});

        threadutils::parallel_for(0, Ny, [&](int j){
        // for (int j = 0; j < Ny; j++) {
            float a1 = float(j)/(Ny-1);
            for (int i = 0; i < Nsub; i++)
            {
              float a2 = float(i) / (Nsub-1);
              Eigen::Vector3f v;
              v << (1-a2)*(-0.5f * S)+a2*0.5f*S,
                   (1-a1)*(-0.5f * S)+a1*0.5f*S,
                   DZ + A * std::sin(1.57f + j * 3.14f + a2 * f * S * 3.14f);
              Eigen::Vector3f dv;
              dv << 0.0f, 0.0f, 5.0f*(std::sin(v(0)/S * 10.0f + T * 2 * 3.14f)+ 0.5f * std::cos(0.5f + v(1)/S * 10.0f + T * 2 * 3.14f));
              X.row(j*Nsub + i) = v + dv;
            }
        // }
        });
        threadutils::parallel_for(0, Ny, [&](int j){
        // for (int j = 0; j < Ny; j++) {
            float a1 = float(j)/(Ny-1);

            for (int i = 0; i < Nsub; i++)
            {
              float a2 = float(i) / (Nsub-1);
              Eigen::Vector3f v;
              v << (1-a1)*(-0.5f * S)+a1*0.5f*S,
                   (1-a2)*(-0.5f * S)+a2*0.5f*S,
                   DZ + A * std::sin(-1.57f + j * 3.14f + a2 * f * S * 3.14f);
              Eigen::Vector3f dv;
              dv << 0.0f, 0.0f, 5.0f*(std::sin(v(0)/S * 10.0f + T * 2 * 3.14f)+ 0.5f * std::cos(0.5f + v(1)/S * 10.0f + T * 2 * 3.14f));
              X.row((Ny+j)*Nsub + i) = v + dv;
            }
        // }
        });
        // * wave pattern dep on T and xy position ?
    }

    const std::vector<uint32_t> &getIndices() { return I; }
    // const Eigen::Matrix<float,-1,-1> & getVertexData() { return X; }
    const Eigen::Matrix<float,-1,-1, Eigen::RowMajor> & getVertexData() { return X; }

    // const MatrixXXs & getVertexData() { return X; }
    // const Eigen::Matrix<scalar,Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & getVertexData() { return X; }

private:
    float T;
    // Eigen::Matrix<float,-1,-1> X;
    Eigen::Matrix<float,-1,-1, Eigen::RowMajor> X;
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