#ifndef __OBSMESHDRAWABLE__H_
#define __OBSMESHDRAWABLE__H_

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/Generic.h>

#include <iostream>
#include <vector>

#include "../EigenDefinitions.h"
#include "../mesh/AbstractMeshProvider.h"
// #include "../mesh/Mesh.h"
#include "../utils/debug_logging.h"

namespace Magnum {
template <typename Shader>
class ObsMeshDrawable {
 public:
  explicit ObsMeshDrawable(
      Shader &shader,
      GL::BufferUsage vertexBufferUsage = GL::BufferUsage::StaticDraw,
      GL::BufferUsage indexBufferUsage  = GL::BufferUsage::StaticDraw)
      : m_indexBufferUsage(indexBufferUsage),
        m_vertexBufferUsage(vertexBufferUsage),
        m_shader(shader) {
    m_glmesh.setIndexBuffer(m_indexBuffer, 0, GL::MeshIndexType::UnsignedInt)
        .setPrimitive(GL::MeshPrimitive::Triangles);
    m_glmesh.addVertexBuffer(this->m_vertexBuffer, 0,
                             Shaders::Generic3D::Position{});                             
  }

  // void setIndices(const MatrixGLi &M) {
  //   m_glmesh.setCount(M.rows() * M.cols());
  //   m_indexBuffer.setData({M.data(), uint32_t(M.size())}, m_indexBufferUsage);
  // }

  // void setVertices(const MatrixGLf &M) {
  //   this->m_vertexBuffer.setData({M.data(), uint32_t(M.size())},
  //                                this->m_vertexBufferUsage);
  // }
  // TODO POTENTIALLY SET OUTSIDE / IN YARNMAPPER
  void setIndices(const VectorBuffer<Mesh::Face> &buf)
  { 
    m_glmesh.setCount(buf.cpu().size()*3); // total number of indices!
    m_indexBuffer.setData({&buf.cpu()[0],
                            uint32_t(buf.cpu().size())},
                          m_indexBufferUsage);
  }

  void setVertices(const VectorBuffer<Mesh::WSVertex> &buf)
  { 
    this->m_vertexBuffer.setData({&buf.cpu()[0],
                            uint32_t(buf.cpu().size())},
                                  this->m_vertexBufferUsage);
  }

  void draw(const Matrix4 &V, const Trafo &m) {
    Matrix4 MV =
        V *
        Matrix4::translation(Vector3{m.translation[0],
                                     m.translation[1],
                                     m.translation[2]}) *
        Matrix4::scaling(Vector3(m.scale)) *
        Matrix4::rotation(Deg(m.axis_angle[0]),
                          Vector3{m.axis_angle[1], m.axis_angle[2],
                                  m.axis_angle[3]});
    m_shader.setTransformation(MV).draw(m_glmesh);
  }

 private:
  GL::BufferUsage m_indexBufferUsage, m_vertexBufferUsage;
  GL::Buffer m_indexBuffer, m_vertexBuffer;
  Shader &m_shader;
  GL::Mesh m_glmesh;
};

}  // namespace Magnum

#endif  // __OBSMESHDRAWABLE__H_