#ifndef __MESHSHADER__H__
#define __MESHSHADER__H__

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/GL/GL.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/Shaders/Generic.h>

namespace Magnum {

class MeshShader : public Magnum::GL::AbstractShaderProgram {
 public:
  using Position = Magnum::Shaders::Generic3D::Position;
  // using Normal = Magnum::Shaders::Generic3D::Normal;
  // using TextureCoordinates = Magnum::Shaders::Generic3D::TextureCoordinates;

  enum {
    AlbedoOutput    = 0,
    PositionsOutput = 1,
    NormalsOutput   = 2,
  };

  explicit MeshShader(Magnum::NoCreateT)
      : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

  explicit MeshShader();

  MeshShader& setTransformation(const Matrix4& transformation);
  MeshShader& setProjection(const Matrix4& projection);
  MeshShader& bindMatCap(GL::Texture2D& texture);

 private:
  enum : Int { TextureUnit_Matcap = 0 };

  Int _transformationUniform, _normalMatrixUniform, _projectionUniform;
};

}  // namespace Magnum

#endif