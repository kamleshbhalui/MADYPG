#ifndef __DeformShader__H__
#define __DeformShader__H__

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/GL/Buffer.h>

class DeformShader : public Magnum::GL::AbstractShaderProgram {
 public:
  enum SSBO {
    XwsBuffer        = 0,
    XmsBuffer        = 1,
    Bary0Buffer      = 2,
    MeshStrainBuffer = 3,
    MeshFmsBuffer    = 4,
    TexAxesBuffer    = 5,
    TexDataBuffer    = 6
  };
  // enum UBO {
  //   TexAxesBuffer = 0 // UBO not possible because layout std430 not supported
  //   in UBO, and std140 is stupid padding each entry in a float[] to 16
  //   bytes..
  // };

  explicit DeformShader(Magnum::NoCreateT)
      : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

  explicit DeformShader();

  void compute(size_t N, Magnum::GL::Buffer &Xws, Magnum::GL::Buffer &Xms,
               Magnum::GL::Buffer &B0, Magnum::GL::Buffer &mS,
               Magnum::GL::Buffer &mFms, Magnum::GL::Buffer &texHeader,
               Magnum::GL::Buffer &texData, float deform_reference,
               float linearized_bending, float min_eigval);

 private:
  Magnum::Int _deformUniform, _linearizedBending, _minEig, _numVertsUniform;
};

#endif  // __DeformShader__H__