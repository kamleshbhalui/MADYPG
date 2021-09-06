#ifndef __MESHSHADER__H__
#define __MESHSHADER__H__

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/GL.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/Shaders/Generic.h>

namespace Magnum {

class MeshShader : public Magnum::GL::AbstractShaderProgram {
 public:
  using Position = Magnum::Shaders::Generic3D::Position;
  // using Normal = Magnum::Shaders::Generic3D::Normal;
  using TextureCoordinates = Magnum::Shaders::Generic3D::TextureCoordinates;

  enum {
    AlbedoOutput    = 0,
    PositionsOutput = 1,
    NormalsOutput   = 2,
  };

  enum SSBO {
    FmsBuffer = 0,
    UBuffer   = 1,
  };

  explicit MeshShader(Magnum::NoCreateT)
      : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

  explicit MeshShader();

  MeshShader& setTransformation(const Matrix4& transformation);
  MeshShader& setProjection(const Matrix4& projection);
  MeshShader& bindMatCap(GL::Texture2D& texture);
  MeshShader& bindClothTexture(GL::Texture2D& texture);
  MeshShader& setTextureScale(float scale);
  MeshShader& setTextureOffset(const Vector2& offset);
  MeshShader& bindUBuffer(Magnum::GL::Buffer& U);
  MeshShader& bindFmsBuffer(Magnum::GL::Buffer& Fms);

 private:
  enum : Int { TextureUnit_Matcap = 0, TextureUnit_ClothTexture = 1 };

  Int _transformationUniform, _normalMatrixUniform, _projectionUniform,
      _uvscaleUniform, _uvoffsetUniform;
};

}  // namespace Magnum

#endif