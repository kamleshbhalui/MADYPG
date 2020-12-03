#include "DeformShader.h"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/FormatStl.h>
#include <Corrade/Utility/Resource.h>
#include <Magnum/GL/Context.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/GL/Version.h>

#include <string>

using namespace Magnum;

DeformShader::DeformShader()
    {
  const Utility::Resource rs{"compute-shaders"};

  MAGNUM_ASSERT_GL_VERSION_SUPPORTED(GL::Version::GL430);
  GL::Shader comp{GL::Version::GL430, GL::Shader::Type::Compute};
  comp.addSource("#define DEFORMSHADER_WRKGRPSIZE " +  std::to_string(DEFORMSHADER_WRKGRPSIZE) + "\n");
  comp.addSource(rs.get("DeformShader.comp"));

  CORRADE_INTERNAL_ASSERT_OUTPUT(comp.compile());
  attachShaders({comp});

  CORRADE_INTERNAL_ASSERT_OUTPUT(link());

  _deformUniform   = uniformLocation("deform_reference");
  _numVertsUniform = uniformLocation("num_vertices");
}