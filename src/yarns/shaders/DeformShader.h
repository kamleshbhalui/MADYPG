#ifndef __DeformShader__H__
#define __DeformShader__H__

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/GL.h>
#include <Magnum/GL/Renderer.h>
#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Vector3.h>
#include "../../EigenDefinitions.h"
#include "../../utils/threadutils.h"
#include <iostream>
#include "ShaderSettings.h"
#include <cmath>

class DeformShader : public Magnum::GL::AbstractShaderProgram {
 public:
  enum SSBO {
    XwsBuffer = 0,
    XmsBuffer = 1,
    Bary0Buffer = 2,
    MeshStrainBuffer = 3,
    MeshFmsBuffer = 4,
    TexAxesBuffer = 5,
    TexDataBuffer = 6
  };
  // enum UBO {
  //   TexAxesBuffer = 0 // UBO not possible because layout std430 not supported in UBO, and std140 is stupid padding each entry in a float[] to 16 bytes..
  // };

  explicit DeformShader(Magnum::NoCreateT)
      : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

  explicit DeformShader();

  // TODO into cpp
  void compute(size_t N, Magnum::GL::Buffer &Xws, Magnum::GL::Buffer &Xms, Magnum::GL::Buffer &B0, Magnum::GL::Buffer &mS, Magnum::GL::Buffer &mFms, Magnum::GL::Buffer &texHeader, Magnum::GL::Buffer &texData, float deform_reference) {
    Xws.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::XwsBuffer);
    Xms.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::XmsBuffer);
    B0.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::Bary0Buffer);
    mS.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::MeshStrainBuffer);
    mFms.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::MeshFmsBuffer);
    texHeader.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::TexAxesBuffer);
    texData.bind(Magnum::GL::Buffer::Target::ShaderStorage, SSBO::TexDataBuffer);
    setUniform(_deformUniform, deform_reference);
    setUniform(_numVertsUniform, uint32_t(N));

    dispatchCompute(Magnum::Vector3ui{uint32_t(std::ceil(N * 1.0f / DEFORMSHADER_WRKGRPSIZE)), 1, 1}); // TODO TEST SIZE

    // Magnum::GL::Renderer::setMemoryBarrier(Magnum::GL::Renderer::MemoryBarrier::ShaderStorage); // the version of Magnum used in this project defines the wrong constant here, defaulting to direct opengl instead
    glMemoryBarrier(GLbitfield(GL_SHADER_STORAGE_BARRIER_BIT));
  }

 private:
    Magnum::Int _deformUniform, _numVertsUniform;
};

#endif  // __DeformShader__H__