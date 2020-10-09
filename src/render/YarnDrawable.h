#ifndef __LINES__H__
#define __LINES__H__

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Magnum.h>
#include <Magnum/Shaders/Generic.h>

#include <iostream>
#include <limits>
#include <vector>

#include "../EigenDefinitions.h"
#include "../yarns/YarnSoup.h"
#include "shaders/YarnShader.h"

namespace Magnum {
class YarnDrawable {
 public:

  explicit YarnDrawable(
      YarnShader &shader, VectorBuffer<VertexWSData> &vertexBufferRef, VectorBuffer<uint32_t> &indexBuffer)
      : m_shader(shader) {
    m_mesh
        .setIndexBuffer(indexBuffer.gpu(), 0, GL::MeshIndexType::UnsignedInt)
        .setCount(indexBuffer.getGPUSize())
        // .setPrimitive(GL::MeshPrimitive::LineStrip);
        .setPrimitive(GL::MeshPrimitive::LineStripAdjacency);
    m_mesh.addVertexBuffer(
        vertexBufferRef.gpu(), 0, YarnShader::Position{},
        // this->m_vertexBuffer, 0, YarnShader::Position{},
        4,  // skip twist variable for now
        YarnShader::TextureCoordinates{},
        YarnShader::Radius{});  // something about memory layout of data in
                                // vertex buffer

    // enable breaking of linestrips within single index buffer by using the
    // index GLuint::max
    glEnable(GL_PRIMITIVE_RESTART); // TODO DO THIS SOMEWHERE ELSE AND EXPOSE GLUINT
    glPrimitiveRestartIndex(std::numeric_limits<GLuint>::max());
  }

  void draw(const Matrix4 &V) {
    // assumed that non-obj specific things like camera proj already set
    // Matrix4 MV = V * M;
    const Matrix4 &MV = V;
    m_shader
        .setTransformation(MV)
        // .setNormalMatrix(MV.normalMatrix())
        .setDiffuseColor(Color4(1.0))  // more like tint
        .setRadius(m_radius)
        .draw(m_mesh);
  }

  float m_radius = 0.001f;

 protected:
  YarnShader &m_shader;
  GL::Mesh m_mesh;

 private:
  std::vector<UnsignedInt> m_indices;
};

}  // namespace Magnum

#endif  // __LINES__H__