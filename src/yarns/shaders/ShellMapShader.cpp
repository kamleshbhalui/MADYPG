#include "ShellMapShader.h"

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

ShellMapShader::ShellMapShader() {
  const Utility::Resource rs{"compute-shaders"};

  MAGNUM_ASSERT_GL_VERSION_SUPPORTED(GL::Version::GL430);
  GL::Shader comp{GL::Version::GL430, GL::Shader::Type::Compute};
  comp.addSource("#define SHELLMAPSHADER_WRKGRPSIZE " +
                 std::to_string(SHELLMAPSHADER_WRKGRPSIZE) + "\n");
  comp.addSource(rs.get("ShellMapShader.comp"));

  CORRADE_INTERNAL_ASSERT_OUTPUT(comp.compile());
  attachShaders({comp});

  CORRADE_INTERNAL_ASSERT_OUTPUT(link());

  // _applyUniform   = uniformLocation("apply");
  _flatNormalsUniform = uniformLocation("flat_normals");
  _numVertsUniform    = uniformLocation("num_vertices");
  _phongUniform       = uniformLocation("phong_deformation");

  // bindFragmentDataLocation(VertexIOBuffer, "ambientOcclusion");
  // setUniform(_samplesUniform, randomSamples);
}

void ShellMapShader::compute(size_t N, Magnum::GL::Buffer &Xws,
                             Magnum::GL::Buffer &B0, Magnum::GL::Buffer &DinvU,
                             Magnum::GL::Buffer &NNv, Magnum::GL::Buffer &mX,
                             Magnum::GL::Buffer &mF, Magnum::GL::Buffer &mFms,
                             Magnum::GL::Buffer &mdefF,
                             Magnum::GL::Buffer &mdefFv, Magnum::GL::Buffer &U,
                             bool flat_normals, float phong) {
  Xws.bind(Magnum::GL::Buffer::Target::ShaderStorage, VertexBuffer);
  B0.bind(Magnum::GL::Buffer::Target::ShaderStorage, BaryBuffer);
  DinvU.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshDiUBuffer);
  NNv.bind(Magnum::GL::Buffer::Target::ShaderStorage, NormalBuffer);
  mX.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshXBuffer);
  mF.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshFBuffer);
  mFms.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshFmsBuffer);
  mdefF.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshdefFBuffer);
  U.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshUBuffer);
  if (phong > 0) {
    mdefFv.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshdefFvBuffer);
  }
  // setUniform(_applyUniform, apply);
  setUniform(_flatNormalsUniform, flat_normals);
  setUniform(_numVertsUniform, uint32_t(N));
  setUniform(_phongUniform, phong);

  dispatchCompute(Magnum::Vector3ui{
      uint32_t(std::ceil(N * 1.0f / SHELLMAPSHADER_WRKGRPSIZE)), 1, 1});
}