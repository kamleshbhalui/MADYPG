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

class ShellMapShader : public Magnum::GL::AbstractShaderProgram {
 public:
  enum {
    VertexBuffer = 0,
    TriBuffer = 1,
    BaryBuffer = 2
  };

  explicit ShellMapShader(Magnum::NoCreateT)
      : Magnum::GL::AbstractShaderProgram{Magnum::NoCreate} {};

  explicit ShellMapShader();


  void compute(size_t N, Magnum::GL::Buffer &vertexBuffer, Magnum::GL::Buffer &triBuffer, Magnum::GL::Buffer &baryBuffer) {
    vertexBuffer.bind(Magnum::GL::Buffer::Target::ShaderStorage, VertexBuffer);
    triBuffer.bind(Magnum::GL::Buffer::Target::ShaderStorage, TriBuffer);
    baryBuffer.bind(Magnum::GL::Buffer::Target::ShaderStorage, BaryBuffer);

    dispatchCompute(Magnum::Vector3ui{N, 1, 1});


    Magnum::GL::Renderer::setMemoryBarrier(Magnum::GL::Renderer::MemoryBarrier::VertexAttributeArray);
    // Magnum::GL::Renderer::setMemoryBarrier(Magnum::GL::Renderer::MemoryBarrier::ShaderStorage);
  }

  // void compute(MatrixGLf& M) {
  //   // std::cout<<"Msize "<<M.size()<<" "<<M.rows()<<" "<<M.cols()<<"\n";

  //   // TODO USE SUBDATA? BC CHRISTIAN SAYS ELSE SLOW REALLOC

  //   #warning todo subdata
  //   m_vertexbuffer.setData({M.data(), uint32_t(M.size())},
  //                          Magnum::GL::BufferUsage::StreamDraw);

  //   m_vertexbuffer.bind(Magnum::GL::Buffer::Target::ShaderStorage, VertexIOBuffer);

  //   dispatchCompute(Magnum::Vector3ui{Magnum::UnsignedInt(M.rows()), 1, 1});
    
  //   // Magnum::GL::Renderer::setMemoryBarrier(Magnum::GL::Renderer::MemoryBarrier::ShaderStorage); // TODO ?
  //   // typename RowMajorMatrix::Scalar* matdata = M.data();
  //   // matdata = reinterpret_cast<typename RowMajorMatrix::Scalar*>(m_vertexbuffer.data().data());
  //   // TODO NEED TO GET DATA BACK FROM GPU BUFFER? HOW TO HACK INTO EIGEN?

  //   // float* data = M.data(); // get the pointer
  //   // new (&M) Eigen::Map<MatrixGLf>(NULL); // swap the mapped array with anything else
  //   // do something with data

  //   // Eigen::Map<MatrixGLf> tst((float*)m_vertexbuffer.data().data(),M.rows(),M.cols());
  //   // for (size_t i = 0; i < 10; i++)
  //   // {
  //   //   std::cout<<M.row(i)<<"\n";
  //   //   std::cout<<tst.row(i)<<"\n";
  //   // }

  //   // TODO SLOW
                           
  //   /* 
  //   it has to allocate a CPU-side memory and then copy there from the GPU (plus stalling the GPU pipeline if inevitable)
  //    instead of this you could map() the buffer memory, which might avoid the other copy, but only under certain conditions (ideally you'll need to call setStorage() on the buffer with correct flags, the classic setData() don't guarantee a buffer is allocated in CPU-accessible memory
  //   */
  //   Magnum::Containers::Array<char> data = m_vertexbuffer.data(); // THIS MAKES A COPY, MIGHT BE AVOIDABLE
  //   Magnum::Containers::ArrayView<float> floats = Magnum::Containers::arrayCast<float>(data); // the arrayCast is just a reinterpret_cast underneath, so no copy
  //   // std::cout<<floats[0]<<"\n";
  //   // std::cout<<"Fsize "<<floats.size()<<"\n";

  //   //when you copy from the floats, you'll have a second copy
  //   // what you could do instead of this is creating an Eigen matrix using the memory of floats, and this is something i'll try to merge today to make it easier: https://github.com/mosra/magnum-integration/pull/74
  //   // for (size_t i = 0; i < M.size(); i++) // TODO BETTER? PARALLEL?
  //     // for (size_t j = 0; j < M.cols(); j++) {
  //       // std::cout<<"ij"<<i<<" "<<j<<"\n";
  //       // std::cout<<"M(ij)="<<M(i,j)<<"; f[ij]="<<floats[i*M.cols()+j]<<"\n";
  //       // M(i,j) = floats[i*M.cols()+j];
  //     // }
  //     // M.data()[i] = floats[i];
  //   threadutils::parallel_for(size_t(0),size_t(M.size()),[&](size_t i){
  //     M.data()[i] = floats[i];
  //   });
  // }

 private:
  Magnum::GL::Buffer m_vertexbuffer;
};

#endif  // __SHELLMAPSHADER__H__