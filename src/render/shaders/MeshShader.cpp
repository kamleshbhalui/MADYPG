#include "MeshShader.h"

#include <Corrade/Utility/Resource.h>
#include <Corrade/Containers/Reference.h>

#include <Magnum/GL/Shader.h>
#include <Magnum/GL/Version.h>
#include <Magnum/GL/Context.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Matrix3.h>
#include <Magnum/Math/Color.h>

namespace Magnum
{
  MeshShader::MeshShader()
  {
    MAGNUM_ASSERT_GL_VERSION_SUPPORTED(GL::Version::GL330);

    const Utility::Resource rs{"ssao-data"};

    GL::Shader vert{GL::Version::GL330, GL::Shader::Type::Vertex};
    GL::Shader geom{GL::Version::GL330, GL::Shader::Type::Geometry};
    GL::Shader frag{GL::Version::GL330, GL::Shader::Type::Fragment};

    vert.addSource(rs.get("MeshShader.vert"));

    geom.addSource(rs.get("MeshShader.geom"));

    frag.addSource(rs.get("MeshShader.frag"));

    CORRADE_INTERNAL_ASSERT_OUTPUT(GL::Shader::compile({vert, geom, frag}));

    attachShaders({vert, geom, frag});

    CORRADE_INTERNAL_ASSERT_OUTPUT(link());

    bindAttributeLocation(Position::Location, "position");
    // bindAttributeLocation(Normal::Location, "normal");

    bindFragmentDataLocation(AlbedoOutput, "color");
    bindFragmentDataLocation(PositionsOutput, "position");
    bindFragmentDataLocation(NormalsOutput, "normal");

    _transformationUniform = uniformLocation("transformation");
    // _normalMatrixUniform = uniformLocation("normalMatrix");
    _projectionUniform = uniformLocation("projection");
    _diffuseColorUniform = uniformLocation("diffuseColor");
    // _radiusUniform = uniformLocation("radius");
    // setUniform(uniformLocation("matcap"), TextureUnit);
  }

  MeshShader &MeshShader::setTransformation(const Matrix4 &transformation)
  {
    setUniform(_transformationUniform, transformation);
    return *this;
  }

  // MeshShader& MeshShader::setNormalMatrix(const Matrix3x3& normalMatrix){
  //    setUniform(_normalMatrixUniform, normalMatrix);
  //    return *this;
  // }

  MeshShader &MeshShader::setProjection(const Matrix4 &projection)
  {
    setUniform(_projectionUniform, projection);
    return *this;
  }

  MeshShader &MeshShader::setDiffuseColor(const Color4 &color)
  {
    setUniform(_diffuseColorUniform, color);
    return *this;
  }

  // MeshShader &MeshShader::setRadius(float radius)
  // {
  //     setUniform(_radiusUniform, radius);
  //     return *this;
  // }

  // MeshShader &MeshShader::bindTexture(GL::Texture2D &texture)
  // {
  //     texture.bind(TextureUnit);
  //     return *this;
  // }

} // namespace Magnum