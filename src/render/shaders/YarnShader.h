#ifndef _YARNSHADER_H_
#define _YARNSHADER_H_

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/GL/GL.h>
#include <Magnum/GL/Texture.h>

namespace Magnum
{

  class YarnShader : public Magnum::GL::AbstractShaderProgram
  {
  public:
    using Position = GL::Attribute<0, Vector4>;
    using Arc = GL::Attribute<1, Float>;
    using Director = GL::Attribute<2, Vector3>;
    using TextureCoordinates = GL::Attribute<3, Vector2>;
    using Radius = GL::Attribute<4, Float>;

    enum
    {
      AlbedoOutput = 0,
      PositionsOutput = 1,
      NormalsOutput = 2,
    };

    explicit YarnShader(Magnum::NoCreateT) : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

    explicit YarnShader();

    YarnShader &setTransformation(const Matrix4 &transformation);
    YarnShader& setNormalMatrix(const Matrix3x3& normalMatrix);
    YarnShader &setProjection(const Matrix4 &projection);
    YarnShader &setDiffuseColor(const Color4 &color);
    YarnShader &setRadius(float radius);
    YarnShader &setNormalTwist(float speed);
    YarnShader &setNormalNum(float num);
    YarnShader &setNormalHeight(float height);
    YarnShader &bindMatCap(GL::Texture2D &texture);
    YarnShader &bindClothTexture(GL::Texture2D &texture);
    YarnShader &bindNormalMap(GL::Texture1D &texture);
    YarnShader &setTextureScale(float scale);

  private:
    enum : Int
    {
      TextureUnit_Matcap = 0,
      TextureUnit_ClothTexture = 1,
      TextureUnit_NormalMap = 2
    };

    Int _transformationUniform,
        _normalMatrixUniform,
        _projectionUniform,
        _diffuseColorUniform,
        _radiusUniform,
        _normalTwistUniform,
        _normalNumUniform,
        _normalHeightUniform,
        _texscaleUniform;
  };

} // namespace Magnum

#endif