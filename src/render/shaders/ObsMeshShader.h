#ifndef __OBSMESHSHADER__H__
#define __OBSMESHSHADER__H__

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/GL/GL.h>
#include <Magnum/GL/Texture.h>

namespace Magnum {

class ObsMeshShader : public Magnum::GL::AbstractShaderProgram {
public:

    using Position = Magnum::Shaders::Generic3D::Position;
    // using Normal = Magnum::Shaders::Generic3D::Normal;
    // using TextureCoordinates = Magnum::Shaders::Generic3D::TextureCoordinates;

    enum {
        AlbedoOutput = 0,
        PositionsOutput = 1,
        NormalsOutput = 2,
    };

    explicit ObsMeshShader(Magnum::NoCreateT) : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

    explicit ObsMeshShader();

    ObsMeshShader& setTransformation(const Matrix4& transformation);
    // MeshShader& setNormalMatrix(const Matrix3x3& normalMatrix);
    ObsMeshShader& setProjection(const Matrix4& projection);
    ObsMeshShader &bindMatCap(GL::Texture2D &texture);

private:
    enum : Int
    {
      TextureUnit_Matcap = 0
    };

    Int _transformationUniform,
        _normalMatrixUniform,
        _projectionUniform;

};

}

#endif