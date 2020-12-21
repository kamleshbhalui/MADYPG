#include "GroundShader.h"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/Resource.h>
#include <Magnum/GL/Context.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/GL/Version.h>
#include <Magnum/Math/Matrix4.h>

using namespace Magnum;

GroundShader::GroundShader() {
  MAGNUM_ASSERT_GL_VERSION_SUPPORTED(GL::Version::GL330);

  const Utility::Resource rs{"ssao-data"};

  GL::Shader vert{GL::Version::GL330, GL::Shader::Type::Vertex};
  GL::Shader geom{GL::Version::GL330, GL::Shader::Type::Geometry};
  GL::Shader frag{GL::Version::GL330, GL::Shader::Type::Fragment};

  vert.addSource(rs.get("GroundShader.vert"));
  geom.addSource(rs.get("GroundShader.geom"));
  frag.addSource(rs.get("GroundShader.frag"));

  CORRADE_INTERNAL_ASSERT_OUTPUT(GL::Shader::compile({vert, geom, frag}));

  attachShaders({vert, geom, frag});

  CORRADE_INTERNAL_ASSERT_OUTPUT(link());

  bindAttributeLocation(Position::Location, "position");
  // bindAttributeLocation(Normal::Location, "normal");

  bindFragmentDataLocation(AlbedoOutput, "color");
  bindFragmentDataLocation(PositionsOutput, "position");
  bindFragmentDataLocation(NormalsOutput, "normal");

  _transformationUniform = uniformLocation("transformation");
  _projectionUniform     = uniformLocation("projection");
  _dYUniform             = uniformLocation("dY");
  _scaleUniform          = uniformLocation("scale");
  setUniform(uniformLocation("tex"), TextureUnit);

  Vector3 X[4]     = {Vector3(-1, 0, 1), Vector3(1, 0, 1), Vector3(1, 0, -1),
                  Vector3(-1, 0, -1)};
  UnsignedInt I[6] = {0, 1, 2, 0, 2, 3};
  GL::Buffer vBuf, iBuf;
  vBuf.setData(X);
  iBuf.setData(I);

  m_glmesh = GL::Mesh();
  m_glmesh.setPrimitive(GL::MeshPrimitive::Triangles)
      .setCount(4 * 3)  // total number of indices!
      .addVertexBuffer(vBuf, 0, Shaders::Generic3D::Position{})
      .setIndexBuffer(iBuf, 0, GL::MeshIndexType::UnsignedInt);
}

GroundShader &GroundShader::bindTexture(GL::Texture2D &texture) {
  texture.bind(TextureUnit);
  return *this;
}

void GroundShader::renderQuad(const Matrix4 &transformation,
                              const Matrix4 &projection) {
  setUniform(_transformationUniform, transformation);
  setUniform(_projectionUniform, projection);
  setUniform(_dYUniform, dY);
  setUniform(_scaleUniform, scale);
  this->draw(m_glmesh);
}
