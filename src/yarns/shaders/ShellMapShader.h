#ifndef __SHELLMAPSHADER__H__
#define __SHELLMAPSHADER__H__

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/GL.h>
#include <Magnum/GL/Renderer.h>
#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Vector3.h>
#include "../../EigenDefinitions.h"
#include "../../utils/threadutils.h"
#include <iostream>
#define SHELLMAPSHADER_WRKGRPSIZE 32 // for nvidia?
#include <cmath>

class ShellMapShader : public Magnum::GL::AbstractShaderProgram {
 public:
  enum {
    VertexBuffer = 0,
    BaryBuffer = 1,
    NormalBuffer = 2,
    MeshXBuffer = 3,
    MeshFBuffer = 4,
    MeshFmsBuffer = 5,
    MeshDiUBuffer = 6,
    MeshdefFBuffer = 7,
    MeshdefFvBuffer = 8,
    MeshUBuffer = 9
  };

  explicit ShellMapShader(Magnum::NoCreateT)
      : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

  explicit ShellMapShader();

  // TODO into cpp
  void compute(size_t N, Magnum::GL::Buffer &Xws, Magnum::GL::Buffer &B0, Magnum::GL::Buffer &DinvU, Magnum::GL::Buffer &NNv, Magnum::GL::Buffer &mX, Magnum::GL::Buffer &mF, Magnum::GL::Buffer &mFms, Magnum::GL::Buffer &mdefF, Magnum::GL::Buffer &mdefFv, Magnum::GL::Buffer &U, bool flat_normals, float phong) {
    Xws.bind(Magnum::GL::Buffer::Target::ShaderStorage, VertexBuffer);
    B0.bind(Magnum::GL::Buffer::Target::ShaderStorage, BaryBuffer);
    DinvU.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshDiUBuffer);
    NNv.bind(Magnum::GL::Buffer::Target::ShaderStorage, NormalBuffer);
    mX.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshXBuffer);
    mF.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshFBuffer);
    mFms.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshFmsBuffer);
    mdefF.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshdefFBuffer);
    U.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshUBuffer);
    if(phong > 0) {
      mdefFv.bind(Magnum::GL::Buffer::Target::ShaderStorage, MeshdefFvBuffer);
    }
    // setUniform(_applyUniform, apply);
    setUniform(_flatNormalsUniform, flat_normals);
    setUniform(_numVertsUniform, uint32_t(N));
    setUniform(_phongUniform, phong);

    dispatchCompute(Magnum::Vector3ui{uint32_t(std::ceil(N * 1.0f / SHELLMAPSHADER_WRKGRPSIZE)), 1, 1}); // TODO TEST SIZE

    Magnum::GL::Renderer::setMemoryBarrier(Magnum::GL::Renderer::MemoryBarrier::VertexAttributeArray);
    // Magnum::GL::Renderer::setMemoryBarrier(Magnum::GL::Renderer::MemoryBarrier::ShaderStorage);
  }

 private:
    Magnum::Int _applyUniform, _flatNormalsUniform, _numVertsUniform, _phongUniform;
};

#endif  // __SHELLMAPSHADER__H__