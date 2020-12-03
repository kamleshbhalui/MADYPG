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
        YarnShader::Arc{},
        YarnShader::Director{},
        YarnShader::TextureCoordinates{},
        YarnShader::Radius{}
        // ,YarnShader::Pad{}
        );  // something about memory layout of data in
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
        .setNormalMatrix(MV.normalMatrix())
        .setDiffuseColor(Color4(1.0))  // more like tint
        .setRadius(m_radius)
        .setPlyTwist(m_nmtwist)
        .setPlyNum(m_nmnum)
        .setPlyHeight(m_nmheight)
        .setPlyLength(m_nmlen)
        .draw(m_mesh);
  }

  float m_radius = 0.001f;
  float m_nmtwist = 1.0f;
  float m_nmnum = 1.0f;
  float m_nmheight = 0.1f;
  float m_nmlen = 1.0f;

 protected:
  YarnShader &m_shader;
  GL::Mesh m_mesh;

 private:
  std::vector<UnsignedInt> m_indices;
};

}  // namespace Magnum

#endif  // __LINES__H__