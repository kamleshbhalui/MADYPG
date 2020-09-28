#ifndef __OBSMESHDRAWABLE__H_
#define __OBSMESHDRAWABLE__H_

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
  class ObsMeshDrawable
  {
  public:
    struct Vertex
    {
      Vector3 position;
    };

    explicit ObsMeshDrawable(
        Shader &shader,
        GL::BufferUsage vertexBufferUsage = GL::BufferUsage::StaticDraw,
        GL::BufferUsage indexBufferUsage = GL::BufferUsage::StaticDraw)
        : m_indexBufferUsage(indexBufferUsage),
          m_vertexBufferUsage(vertexBufferUsage), m_shader(shader)
    {
      m_glmesh.setIndexBuffer(m_indexBuffer, 0, GL::MeshIndexType::UnsignedInt)
          .setPrimitive(GL::MeshPrimitive::Triangles);
      m_glmesh.addVertexBuffer(
          this->m_vertexBuffer, 0, Shaders::Generic3D::Position{});
    }

    void setIndices(const MatrixGLi &M)
    {
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
      const Matrix4 &MV = V;
      m_shader.setTransformation(MV)
          .draw(m_glmesh);
    }

  private:
    GL::BufferUsage m_indexBufferUsage, m_vertexBufferUsage;
    GL::Buffer m_indexBuffer, m_vertexBuffer;
    Shader &m_shader;
    GL::Mesh m_glmesh;
  };

} // namespace Magnum

#endif // __OBSMESHDRAWABLE__H_