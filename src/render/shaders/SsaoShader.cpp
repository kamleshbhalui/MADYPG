#include "SsaoShader.h"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/FormatStl.h>
#include <Corrade/Utility/Resource.h>
#include <Magnum/GL/Context.h>
#include <Magnum/GL/ImageFormat.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/GL/Version.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/Math/Matrix4.h>

#include <random>

using namespace Magnum;

enum {
  PositionUnit  = 0,
  NormalUnit    = 1,
  NoiseUnit     = 2,
  OcclusionUnit = 3,
};

SsaoShader::SsaoShader(UnsignedInt sampleCount) {
  Containers::Array<Vector3> randomSamples(Containers::NoInit, sampleCount);
  std::default_random_engine engine;
  std::normal_distribution<float> ndistr(0.0f, 0.4f);
  std::uniform_real_distribution<float> distr(0.01f, 1.0f);
  constexpr float c = 0.1f;  // minimum z value of random sample, to avoid self
                             // shadowing of inclined planes

  auto lerp = [](float a, float b, float f) { return a + f * (b - a); };

  for (UnsignedInt i = 0; i < sampleCount; ++i) {
    Vector3 p{ndistr(engine), ndistr(engine), ndistr(engine)};
    if (p.z() < 0)
      p.z() *= -1.f;
    // p = Math::fmod(p, Vector3{1});
    // randomSamples[i] = p;
    // continue;
    p = p.normalized();
    if (p.z() < c) {  // minimum cosine / z value
      Vector2 p2 =
          Vector2{ndistr(engine), ndistr(engine)}.normalized() * (1 - c * c);
      p = Vector3{p2.x(), p2.y(), c};
    }
    p *= distr(engine);  // hemisphere volume
    float scale = float(i) / (sampleCount - 1);
    p *= lerp(
        0.1f, 1.0f,
        scale * scale);  // sample closer to origin // artistically, this
                         // concentrates AO to obj edges compared to spreading
                         // darkness along the entire radius distance
    randomSamples[i] = p;
  }

  Utility::Resource rs{"ssao-data"};

  MAGNUM_ASSERT_GL_VERSION_SUPPORTED(GL::Version::GL330);
  GL::Shader vert{GL::Version::GL330, GL::Shader::Type::Vertex};
  GL::Shader frag{GL::Version::GL330, GL::Shader::Type::Fragment};

  vert.addSource(rs.get("FullScreenTriangle.vert"));
#ifdef MSAA
  frag.addSource("#define MSAA " + std::to_string(MSAA) + "\n");
#endif
  frag.addSource(
          Utility::formatString("#define SAMPLE_COUNT {}\n", sampleCount))
      .addSource(rs.get("Ssao.frag"));

  CORRADE_INTERNAL_ASSERT_OUTPUT(GL::Shader::compile({vert, frag}));
  attachShaders({vert, frag});

  CORRADE_INTERNAL_ASSERT_OUTPUT(link());

  bindFragmentDataLocation(AmbientOcclusionOutput, "ambientOcclusion");

  setUniform(uniformLocation("positions"), PositionUnit);
  setUniform(uniformLocation("normals"), NormalUnit);
  setUniform(uniformLocation("noise"), NoiseUnit);

  _projectionUniform = uniformLocation("projection");
  _biasUniform       = uniformLocation("bias");
  _radiusUniform     = uniformLocation("radius");
  _samplesUniform    = uniformLocation("samples");

  setUniform(_samplesUniform, randomSamples);
}

#ifdef MSAA
SsaoShader& SsaoShader::bindPositionTexture(GL::MultisampleTexture2D& texture) {
  texture.bind(PositionUnit);
  return *this;
}

SsaoShader& SsaoShader::bindNormalTexture(GL::MultisampleTexture2D& texture) {
  texture.bind(NormalUnit);
  return *this;
}

#else
SsaoShader& SsaoShader::bindPositionTexture(GL::Texture2D& texture) {
  texture.bind(PositionUnit);
  return *this;
}

SsaoShader& SsaoShader::bindNormalTexture(GL::Texture2D& texture) {
  texture.bind(NormalUnit);
  return *this;
}

#endif

SsaoShader& SsaoShader::bindOcclusionTexture(GL::Texture2D& texture) {
  texture.bindImage(OcclusionUnit, 0, GL::ImageAccess::WriteOnly,
                    GL::ImageFormat::R32F);
  return *this;
}

SsaoShader& SsaoShader::bindNoiseTexture(Magnum::GL::Texture2D& texture) {
  texture.bind(NoiseUnit);
  return *this;
}

SsaoShader& SsaoShader::setProjectionMatrix(const Matrix4& projection) {
  setUniform(_projectionUniform, projection);
  return *this;
}

SsaoShader& SsaoShader::setSampleRadius(Float radius) {
  setUniform(_radiusUniform, radius);
  return *this;
}

SsaoShader& SsaoShader::setBias(Float bias) {
  setUniform(_biasUniform, bias);
  return *this;
}
