#include "DeformShader.h"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/FormatStl.h>
#include <Corrade/Utility/Resource.h>
#include <Magnum/GL/Context.h>
#include <Magnum/GL/GL.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/GL/Version.h>
#include <Magnum/Math/Vector3.h>

#include <cmath>
#include <string>

#include "ShaderSettings.h"

using namespace Magnum;

DeformShader::DeformShader() {
  const Utility::Resource rs{"compute-shaders"};

  MAGNUM_ASSERT_GL_VERSION_SUPPORTED(GL::Version::GL430);
  GL::Shader comp{GL::Version::GL430, GL::Shader::Type::Compute};
  comp.addSource("#define DEFORMSHADER_WRKGRPSIZE " +
                 std::to_string(DEFORMSHADER_WRKGRPSIZE) + "\n");
  comp.addSource(rs.get("DeformShader.comp"));

  CORRADE_INTERNAL_ASSERT_OUTPUT(comp.compile());
  attachShaders({comp});

  CORRADE_INTERNAL_ASSERT_OUTPUT(link());

  _deformUniform   = uniformLocation("deform_reference");
  _numVertsUniform = uniformLocation("num_vertices");
}

void DeformShader::compute(size_t N, Magnum::GL::Buffer &Xws,
                           Magnum::GL::Buffer &Xms, Magnum::GL::Buffer &B0,
                           Magnum::GL::Buffer &mS, Magnum::GL::Buffer &mFms,
                           Magnum::GL::Buffer &texHeader,
                           Magnum::GL::Buffer &texData,
                           float deform_reference) {
  Xws.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::XwsBuffer);
  Xms.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::XmsBuffer);
  B0.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::Bary0Buffer);
  mS.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::MeshStrainBuffer);
  mFms.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::MeshFmsBuffer);
  texHeader.bind(Magnum::GL::Buffer::Target::ShaderStorage,
                 SSBO::TexAxesBuffer);
  texData.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::TexDataBuffer);
  setUniform(_deformUniform, deform_reference);
  setUniform(_numVertsUniform, uint32_t(N));

  dispatchCompute(Magnum::Vector3ui{
      uint32_t(std::ceil(N * 1.0f / DEFORMSHADER_WRKGRPSIZE)), 1, 1});
}