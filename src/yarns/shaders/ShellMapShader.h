#ifndef __SHELLMAPSHADER__H__
#define __SHELLMAPSHADER__H__

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/GL/Buffer.h>

class ShellMapShader : public Magnum::GL::AbstractShaderProgram {
 public:
  enum {
    VertexBuffer    = 0,
    BaryBuffer      = 1,
    NormalBuffer    = 2,
    MeshXBuffer     = 3,
    MeshFBuffer     = 4,
    MeshFmsBuffer   = 5,
    MeshDiUBuffer   = 6,
    MeshdefFBuffer  = 7,
    MeshdefFvBuffer = 8,
    MeshUBuffer     = 9
  };

  explicit ShellMapShader(Magnum::NoCreateT)
      : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

  explicit ShellMapShader();

  void compute(size_t N, Magnum::GL::Buffer &Xws, Magnum::GL::Buffer &B0,
               Magnum::GL::Buffer &DinvU, Magnum::GL::Buffer &NNv,
               Magnum::GL::Buffer &mX, Magnum::GL::Buffer &mF,
               Magnum::GL::Buffer &mFms, Magnum::GL::Buffer &mdefF,
               Magnum::GL::Buffer &mdefFv, Magnum::GL::Buffer &U,
               bool flat_normals, float phong);

 private:
  Magnum::Int _applyUniform, _flatNormalsUniform, _numVertsUniform,
      _phongUniform;
};

#endif  // __SHELLMAPSHADER__H__