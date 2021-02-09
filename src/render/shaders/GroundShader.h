#ifndef __GROUNDSHADER__H__
#define __GROUNDSHADER__H__

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/GL.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/Shaders/Generic.h>

namespace Magnum {

// class to specifically draw a textured ground plane
// NOTE: shader and hardcoded quad mesh drawable combined
class GroundShader : public Magnum::GL::AbstractShaderProgram {
 public:
  using Position = Magnum::Shaders::Generic3D::Position;

  enum {
    AlbedoOutput    = 0,
    PositionsOutput = 1,
    NormalsOutput   = 2,
  };

  explicit GroundShader(Magnum::NoCreateT)
      : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

  explicit GroundShader();

  // GroundShader& setTransformation(const Matrix4& transformation);
  // GroundShader& setProjection(const Matrix4& projection);
  void renderQuad(const Matrix4& transformation, const Matrix4& projection);
  GroundShader& bindTexture(GL::Texture2D& texture);

  float dY    = -0.01f;
  float scale = 1.0f;

 private:
  enum : Int { TextureUnit = 0 };
  Int _transformationUniform, _projectionUniform;
  Int _dYUniform, _scaleUniform;
  GL::Mesh m_glmesh{Magnum::NoCreate};
};

}  // namespace Magnum

#endif