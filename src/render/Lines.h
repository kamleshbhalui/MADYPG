#ifndef __LINES__H__
#define __LINES__H__

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/Magnum.h>
#include <numeric> // iota
#include <vector>

namespace Magnum
{
  template <typename Shader>
  class Lines
  {
  public:
    struct Vertex
    {
      Vector3 position;
    };

    explicit Lines(
        Shader &shader,
        GL::BufferUsage vertexBufferUsage = GL::BufferUsage::StreamDraw, // probably change every frame
        GL::BufferUsage indexBufferUsage = GL::BufferUsage::StaticDraw)  // probably only ever change once
        : m_indexBufferUsage(indexBufferUsage),
          m_vertexBufferUsage(vertexBufferUsage), m_shader(shader)
    {
      m_mesh.setIndexBuffer(m_indexBuffer, 0, GL::MeshIndexType::UnsignedInt)
          .setPrimitive(GL::MeshPrimitive::LineStrip);
      m_mesh.addVertexBuffer(
          this->m_vertexBuffer, 0, Shaders::Generic3D::Position{}); // something about memory layout of data in vertex buffer
    }

    void setIndices(const std::vector<UnsignedInt> &indices)
    {
      this->m_indices = indices;
      m_mesh.setCount(this->m_indices.size());
      m_indexBuffer.setData(this->m_indices, m_indexBufferUsage);
    }

    void setVertices(const std::vector<Vector3> &vertices,
                     bool recalculateIndices = false)
    {
      printf("VERTS, TODO avoid storing local vertices and indices ...\n");
      this->m_vertices.resize(vertices.size());
      for (size_t i = 0; i < vertices.size(); ++i)
      {
        this->m_vertices[i].position = vertices[i];
      }
      this->m_vertexBuffer.setData(this->m_vertices, this->m_vertexBufferUsage);

      if (recalculateIndices)
        this->recalculateIndices(this->m_vertices.size());
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
          .setDiffuseColor(Color4(0.9))
          .setRadius(m_radius)
          .draw(m_mesh);
    }

    float m_radius = 1.0f;

  protected:
    GL::BufferUsage m_indexBufferUsage, m_vertexBufferUsage;
    GL::Buffer m_indexBuffer, m_vertexBuffer;
    Shader &m_shader;
    GL::Mesh m_mesh;

    void recalculateIndices(int n)
    {
      m_indices.resize(n);
      std::iota(std::begin(m_indices), std::end(m_indices),
                0); // increasing integer values
      m_mesh.setCount(this->m_indices.size());
      m_indexBuffer.setData(this->m_indices, m_indexBufferUsage);
    }

  private:
    std::vector<UnsignedInt> m_indices;
    std::vector<Vertex> m_vertices;
  };

} // namespace Magnum

#endif // __LINES__H__