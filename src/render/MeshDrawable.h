#ifndef __MESHDRAWABLE__H_
#define __MESHDRAWABLE__H_

#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>
#include <vector>
#include <iostream>
#include "../EigenDefinitions.h"
// #include "../mesh/Mesh.h"

namespace Magnum
{
  template <typename Shader>
  class MeshDrawable
  {
  public:
    struct Vertex
    {
      Vector3 position;
    };

    explicit MeshDrawable(
        Shader &shader,
        GL::BufferUsage vertexBufferUsage = GL::BufferUsage::StreamDraw, // probably change every frame
        GL::BufferUsage indexBufferUsage = GL::BufferUsage::StaticDraw)  // probably only ever change once
        : m_indexBufferUsage(indexBufferUsage),
          m_vertexBufferUsage(vertexBufferUsage), m_shader(shader)
    {
      m_glmesh.setIndexBuffer(m_indexBuffer, 0, GL::MeshIndexType::UnsignedInt)
          .setPrimitive(GL::MeshPrimitive::Triangles);
      m_glmesh.addVertexBuffer(
          this->m_vertexBuffer, 0, Shaders::Generic3D::Position{}); // something about memory layout of data in vertex buffer
    }

    void setIndices(const MatrixGLi &M)
    {
      // this->m_indices = indices;
      m_glmesh.setCount(M.rows() * M.cols());
      m_indexBuffer.setData({M.data(),
                             uint32_t(M.size())},
                            m_indexBufferUsage);
    }

    void setVertices(const MatrixGLf &M)
    {
      this->m_vertexBuffer.setData({M.data(),
                                    uint32_t(M.size())},
                                   this->m_vertexBufferUsage);
    }

    void draw(const Matrix4 &V)
    {
      GL::Renderer::setLineWidth(m_linewidth);

      // assumed that camera proj already set
      // Matrix4 M = Matrix4::scaling(Vector3(1.0f));
      // Matrix4 MV = V * M;
      const Matrix4 &MV = V;
      m_shader.setTransformation(MV)
          // .setNormalMatrix(MV.normalMatrix())
          // .setProjection(camera.projectionMatrix())
          // .setProjection(_projection)
          .setDiffuseColor(m_color) // more like tint
          // .setDiffuseColor(Color4(1.0,0.5,0.5,1.0))
          // .setRadius(m_radius)
          .draw(m_glmesh);
    }

  private:
    GL::BufferUsage m_indexBufferUsage, m_vertexBufferUsage;
    GL::Buffer m_indexBuffer, m_vertexBuffer;
    Shader &m_shader;
    GL::Mesh m_glmesh;
    Color4 m_color = Color4(Color3(0.9f),1.0f);
    float m_linewidth = 2.0f;
  };

} // namespace Magnum

#endif // __MESHDRAWABLE__H_