#include "ShellMapShader.h"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/FormatStl.h>
#include <Corrade/Utility/Resource.h>
#include <Magnum/GL/Context.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/GL/Version.h>

using namespace Magnum;

ShellMapShader::ShellMapShader()
    {
  const Utility::Resource rs{"compute-shaders"};

  MAGNUM_ASSERT_GL_VERSION_SUPPORTED(GL::Version::GL430);
  GL::Shader comp{GL::Version::GL430, GL::Shader::Type::Compute};
  comp.addSource(rs.get("ShellMapShader.comp"));

  CORRADE_INTERNAL_ASSERT_OUTPUT(comp.compile());
  attachShaders({comp});

  CORRADE_INTERNAL_ASSERT_OUTPUT(link());

  
  _applyUniform   = uniformLocation("apply");
  _flatNormalsUniform = uniformLocation("flat_normals");

  // bindFragmentDataLocation(VertexIOBuffer, "ambientOcclusion");
  // setUniform(_samplesUniform, randomSamples);
}