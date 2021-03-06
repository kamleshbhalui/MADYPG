/*
    This file is part of Magnum.

    Original authors — credit is appreciated but not required:

        2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 —
            Vladimír Vondruš <mosra@centrum.cz>
        2020 — Janos Meny <janos.meny@gmail.com>

    This is free and unencumbered software released into the public domain.

    Anyone is free to copy, modify, publish, use, compile, sell, or distribute
    this software, either in source code form or as a compiled binary, for any
    purpose, commercial or non-commercial, and by any means.

    In jurisdictions that recognize copyright laws, the author or authors of
    this software dedicate any and all copyright interest in the software to
    the public domain. We make this dedication for the benefit of the public
    at large and to the detriment of our heirs and successors. We intend this
    dedication to be an overt act of relinquishment in perpetuity of all
    present and future rights to this software under copyright law.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "SsaoApplyShader.h"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/Resource.h>
#include <Magnum/GL/Context.h>
#include <Magnum/GL/ImageFormat.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/GL/Version.h>
#include <Magnum/Math/Matrix4.h>

using namespace Magnum;

enum {
  PositionUnit         = 0,
  NormalUnit           = 1,
  AlbedoUnit           = 2,
  AmbientOcclusionUnit = 3
};

SsaoApplyShader::SsaoApplyShader(Flag flag) {
  MAGNUM_ASSERT_GL_VERSION_SUPPORTED(GL::Version::GL330);

  Utility::Resource rs{"ssao-data"};

  GL::Shader vert{GL::Version::GL330, GL::Shader::Type::Vertex};
  GL::Shader frag{GL::Version::GL330, GL::Shader::Type::Fragment};

  vert.addSource(rs.get("FullScreenTriangle.vert"));
#ifdef MSAA
  frag.addSource("#define MSAA " + std::to_string(MSAA) + "\n");
#endif
  frag.addSource(flag == Flag::DrawAmbientOcclusion
                     ? "#define DRAW_OCCLUSION_FACTOR\n"
                     : "")
      .addSource(rs.get("SsaoApply.frag"));

  CORRADE_INTERNAL_ASSERT_OUTPUT(GL::Shader::compile({vert, frag}));

  attachShaders({vert, frag});

  CORRADE_INTERNAL_ASSERT_OUTPUT(link());

  // _aoBlurRadiusUniform  = uniformLocation("ao_blur_radius");
  // _aoBlurFeatureUniform = uniformLocation("ao_blur_feature");
  // setUniform(uniformLocation("positionTexture"), PositionUnit);

  _aoPowUniform = uniformLocation("ao_pow");

  if (flag != Flag::DrawAmbientOcclusion) {
    // setUniform(uniformLocation("normalTexture"), NormalUnit);
    setUniform(uniformLocation("albedoTexture"), AlbedoUnit);
  }

  setUniform(uniformLocation("ambientOcclusionTexture"), AmbientOcclusionUnit);
}

#ifdef MSAA
SsaoApplyShader& SsaoApplyShader::bindAlbedoTexture(
    GL::MultisampleTexture2D& texture) {
  texture.bind(AlbedoUnit);
  return *this;
}
SsaoApplyShader& SsaoApplyShader::bindPositionTexture(
    GL::MultisampleTexture2D& texture) {
  texture.bind(PositionUnit);
  return *this;
}
// SsaoApplyShader& SsaoApplyShader::bindNormalTexture(GL::MultisampleTexture2D&
// texture) {
//   texture.bind(NormalUnit);
//   return *this;
// }
#else
SsaoApplyShader& SsaoApplyShader::bindAlbedoTexture(GL::Texture2D& texture) {
  texture.bind(AlbedoUnit);
  return *this;
}
SsaoApplyShader& SsaoApplyShader::bindPositionTexture(GL::Texture2D& texture) {
  texture.bind(PositionUnit);
  return *this;
}
// SsaoApplyShader& SsaoApplyShader::bindNormalTexture(GL::Texture2D& texture) {
//   texture.bind(NormalUnit);
//   return *this;
// }
#endif

SsaoApplyShader& SsaoApplyShader::bindOcclusionTexture(GL::Texture2D& texture) {
  texture.bind(AmbientOcclusionUnit);
  return *this;
}

SsaoApplyShader& SsaoApplyShader::setLightPosition(const Vector3& position) {
  setUniform(_lightPositionUniform, position);
  return *this;
}

SsaoApplyShader& SsaoApplyShader::setLightColor(const Color3& color) {
  setUniform(_lightColorUniform, color);
  return *this;
}

SsaoApplyShader& SsaoApplyShader::setShininess(Float shininess) {
  setUniform(_shininessUniform, shininess);
  return *this;
}

SsaoApplyShader& SsaoApplyShader::setSpecularColor(const Color3& color) {
  setUniform(_specularColorUniform, color);
  return *this;
}

SsaoApplyShader& SsaoApplyShader::setAOBlurRadius(Int radius) {
  setUniform(_aoBlurRadiusUniform, radius);
  return *this;
}

SsaoApplyShader& SsaoApplyShader::setAOBlurFeature(Float feature) {
  setUniform(_aoBlurFeatureUniform, feature);
  return *this;
}

SsaoApplyShader& SsaoApplyShader::setAOPow(Float pow) {
  setUniform(_aoPowUniform, pow);
  return *this;
}
