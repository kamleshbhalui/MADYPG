#ifndef __LINES__H__
#define __LINES__H__

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Shaders/Generic.h>
#include "shaders/YarnShader.h"
#include <Magnum/Magnum.h>
#include <numeric> // iota
#include <vector>
#include <limits>
#include "../EigenDefinitions.h"
#include <iostream>

namespace Magnum
{
  class YarnDrawable
  {
  public:
    struct Vertex
    {
      Vector3 position;
      float radius;
    };

    explicit YarnDrawable(
        YarnShader &shader,
        GL::BufferUsage vertexBufferUsage = GL::BufferUsage::StreamDraw, // probably change every frame
        GL::BufferUsage indexBufferUsage = GL::BufferUsage::StaticDraw)  // probably only ever change once
        : m_indexBufferUsage(indexBufferUsage),
          m_vertexBufferUsage(vertexBufferUsage), m_shader(shader)
    {
      m_mesh.setIndexBuffer(m_indexBuffer, 0, GL::MeshIndexType::UnsignedInt)
          // .setPrimitive(GL::MeshPrimitive::LineStrip);
          .setPrimitive(GL::MeshPrimitive::LineStripAdjacency);
      m_mesh.addVertexBuffer(
          this->m_vertexBuffer, 0, YarnShader::Position{},
          4, // skip twist variable for now
          YarnShader::TextureCoordinates{},
          YarnShader::Radius{}
          ); // something about memory layout of data in vertex buffer

      // enable breaking of linestrips within single index buffer by using the index GLuint::max
      glEnable(GL_PRIMITIVE_RESTART);
      glPrimitiveRestartIndex(std::numeric_limits<GLuint>::max());
    /* @ compute shader 
      // only works if vertexbufferdata is correctly padded !!! (maybe not for std430 layout, when making struct just be floats, eg float val_x, val_y, val_z...)
      m_vertexBuffer.bind(GL::Buffer::Target::ShaderStorage, inputV_index_in_shader?)
      m_indexBuffer.bind(GL::Buffer::Target::ShaderStorage, inputI_index_in_shader?)
      outputBuffer.bind(GL::Buffer::Target::ShaderStorage, outputV_index_in_shader?)
      outputBuffer.bind(GL::Buffer::Target::ShaderStorage, outputF_index_in_shader?)
      compshader.dispatchCompute({Nverts,1,1}) //local size 1,1,1 ie iterate each vertex separately

      shader will check each vertex, if it is start (ix==0) or (ix-1 == RESTART) fill left cap, if (ix+1=last) or (ix+1 == RESTART) fill right, if ix==RESTART skip computation (and set output such that there is no face). slightly more complicated with more neighbor lookups maybe

      cast outputbuffer.data() to vector<struct> and serialize

      // https://github.com/mosra/magnum-examples/pull/91/files
    */
    }

    void setIndices(const std::vector<UnsignedInt> &indices)
    {
      // this->m_indices = indices;
      m_mesh.setCount(indices.size());
      m_indexBuffer.setData(indices, m_indexBufferUsage);
    }

    template< typename Matrix>
    void setVertices(const Matrix &M) // assume row major M
    {
      // for (int i = 0; i < std::min(int(M.rows()),100); i++)
      // {
      //   std::cout<<M.row(i)<<"\n";
      // }
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
    YarnShader &m_shader;
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