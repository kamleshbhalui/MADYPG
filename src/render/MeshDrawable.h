#ifndef __MESHDRAWABLE__H_
#define __MESHDRAWABLE__H_

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/Generic.h>

#include <iostream>
#include <vector>

#include "../EigenDefinitions.h"
#include "../mesh/Mesh.h"

namespace Magnum {
template <typename Shader>
// Class for drawing cloth meshes
// using references to vertex/index(=face) buffer from YarnMapper->Mesh 
class MeshDrawable {
 public:
  explicit MeshDrawable(Shader &shader, VectorBuffer<Mesh::Face> &indexBuffer,
                        VectorBuffer<Mesh::WSVertex> &vertexBuffer)
      : m_shader(shader) {
    m_glmesh
        .setIndexBuffer(indexBuffer.gpu(), 0, GL::MeshIndexType::UnsignedInt)
        .setPrimitive(GL::MeshPrimitive::Triangles);
    m_glmesh.setCount(indexBuffer.getGPUSize() *
                      3);  // *3 because of triangle faces
    m_glmesh.addVertexBuffer(
        vertexBuffer.gpu(), 0,
        Shaders::Generic3D::Position{});  // something about memory layout of
                                          // data in vertex buffer
  }

  void updateIndexCount(int32_t Nindices) {
    if (m_glmesh.count() != Nindices) {
      m_glmesh.setCount(Nindices);
    }
  }

  void draw(const Matrix4 &V) {
    GL::Renderer::setLineWidth(m_linewidth);

    // assumed that camera proj already set
    // Matrix4 M = Matrix4::scaling(Vector3(1.0f));
    // Matrix4 MV = V * M;
    const Matrix4 &MV = V;
    m_shader
        .setTransformation(MV)
        .draw(m_glmesh);
  }

 private:
  GL::BufferUsage m_indexBufferUsage, m_vertexBufferUsage;
  GL::Buffer m_indexBuffer, m_vertexBuffer;
  Shader &m_shader;
  GL::Mesh m_glmesh;
  float m_linewidth = 2.0f;
};

}  // namespace Magnum

#endif  // __MESHDRAWABLE__H_