#ifndef __LINES__H__
#define __LINES__H__

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/Magnum.h>
#include <numeric> // iota
#include <vector>
#include <limits>
#include "../EigenDefinitions.h"
#include <iostream>

namespace Magnum
{
  template <typename Shader>
  class YarnDrawable
  {
  public:
    struct Vertex
    {
      Vector3 position;
    };

    explicit YarnDrawable(
        Shader &shader,
        GL::BufferUsage vertexBufferUsage = GL::BufferUsage::StreamDraw, // probably change every frame
        GL::BufferUsage indexBufferUsage = GL::BufferUsage::StaticDraw)  // probably only ever change once
        : m_indexBufferUsage(indexBufferUsage),
          m_vertexBufferUsage(vertexBufferUsage), m_shader(shader)
    {
      m_mesh.setIndexBuffer(m_indexBuffer, 0, GL::MeshIndexType::UnsignedInt)
          // .setPrimitive(GL::MeshPrimitive::LineStrip);
          .setPrimitive(GL::MeshPrimitive::LineStripAdjacency);
      m_mesh.addVertexBuffer(
          this->m_vertexBuffer, 0, Shaders::Generic3D::Position{}); // something about memory layout of data in vertex buffer

      // enable breaking of linestrips within single index buffer by using the index GLuint::max
      glEnable(GL_PRIMITIVE_RESTART);
      glPrimitiveRestartIndex(std::numeric_limits<GLuint>::max());
    }

    void setIndices(const std::vector<UnsignedInt> &indices)
    {
      // this->m_indices = indices;
      m_mesh.setCount(indices.size());
      m_indexBuffer.setData(indices, m_indexBufferUsage);
    }

    template< typename Matrix>
    void setVertices(const Matrix &M)
    {
      // assume row major ? or col major? 

      // Eigen::Map<Eigen::RowVectorXf> v({M.data()}, M.size());
// (size_t) sizeof(typename Matrix::Scalar) *
// std::cout<< (uint32_t(sizeof(typename Matrix::Scalar)) * uint32_t(M.size()))<<"\n";
// std::cout<< uint32_t(sizeof(typename Matrix::Scalar))<<" * "<<uint32_t(M.size())<<"\n";
//       if (M.innerSize() == M.outerStride()) {
//         std::cout<< "YES\n";
//       }else{
//         std::cout<< "NO\n";

//       }

//       std::cout<< M.rows()<<" "<<M.cols()<<"\n";
//       int N = 40000;
//       for (int i = 0; i < N; i++)
//       {
//         std::cout<<M(i,0)<<" "<<M(i,1)<<" "<<M(i,2)<<"\n";
//         std::cout<<M.data()[3*i+0]<<" "<<M.data()[3*i+1]<<" "<<M.data()[3*i+2]<<"\n";
//       }
//       std::cout<< "\n";
      

      this->m_vertexBuffer.setData({M.data(),
        uint32_t(M.size())}, this->m_vertexBufferUsage);
    }

    void setVertices(const std::vector<Vector3> &vertices,
                     bool recalculateIndices = false)
    {
      this->m_vertexBuffer.setData(vertices, this->m_vertexBufferUsage);

      if (recalculateIndices)
        this->recalculateIndices(vertices.size());
    }

    // void draw(const Matrix4 &transformationMatrix, Camera3D &camera) {
    //   (*m_shader).setTransformationProjectionMatrix(camera.projectionMatrix() *
    //                                                 transformationMatrix);
    //   GL::Renderer::setLineWidth(m_size);
    //   m_mesh.draw(*m_shader);
    // }
    void draw(const Matrix4 &V)
    {
      // assumed that camera proj already set
      // Matrix4 M = Matrix4::scaling(Vector3(1.0f));
      // Matrix4 MV = V * M;
      const Matrix4& MV = V;
      m_shader.setTransformation(MV)
          // .setNormalMatrix(MV.normalMatrix())
          // .setProjection(camera.projectionMatrix())
          // .setProjection(_projection)
          .setDiffuseColor(Color4(1.0)) // more like tint
          // .setDiffuseColor(Color4(1.0,0.5,0.5,1.0))
          .setRadius(m_radius)
          .draw(m_mesh);
    }

    float m_radius = 0.001f;

  protected:
    GL::BufferUsage m_indexBufferUsage, m_vertexBufferUsage;
    GL::Buffer m_indexBuffer, m_vertexBuffer;
    Shader &m_shader;
    GL::Mesh m_mesh;

    void recalculateIndices(int n)
    {
      std::vector<UnsignedInt> indices(n);
      std::iota(std::begin(indices), std::end(indices),
                0); // increasing integer values
      m_mesh.setCount(indices.size());
      m_indexBuffer.setData(indices, m_indexBufferUsage);
    }

  private:
    std::vector<UnsignedInt> m_indices;
    std::vector<Vertex> m_vertices;
      // _yarnGeometryShader.bindTexture(_matcap);

  };

} // namespace Magnum

#endif // __LINES__H__