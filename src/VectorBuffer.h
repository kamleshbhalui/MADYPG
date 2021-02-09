#ifndef __VECTORBUFFER__H__
#define __VECTORBUFFER__H__

#include <Magnum/GL/Buffer.h>
#include <Corrade/Containers/Array.h>
#include <assert.h>

// #include <Eigen/Dense>
// #include <Eigen/Core>
// #include <type_traits>

// Class that interfaces cpu and gpu data, where the cpu side is a vector
// of some struct T, and the gpu side is a Magnum GL::Buffer.
template <typename T>
class VectorBuffer {
 public:
  // reference to the cpu-data vector
  std::vector<T>& cpu() { return m_data; }
  const std::vector<T>& cpu() const { return m_data; }
  // reference to the gpu-data GL::Buffer
  Magnum::GL::Buffer& gpu() { return m_buffer; }
  const Magnum::GL::Buffer& gpu() const { return m_buffer; }

  // upload the data from the cpu to the gpu buffer
  void bufferData(Magnum::GL::BufferUsage bufferUsage =
                      Magnum::GL::BufferUsage::StreamDraw) {
    bool realloc = (m_allocatedGPUsize != m_data.size());
    bufferData(realloc, bufferUsage);
  }
  void bufferData(bool realloc, Magnum::GL::BufferUsage bufferUsage =
                                    Magnum::GL::BufferUsage::StreamDraw) {
    if (realloc) {
      m_buffer.setData({m_data.data(), uint32_t(m_data.size())}, bufferUsage);
      m_allocatedGPUsize = m_data.size();
    } else
      m_buffer.setSubData(0, {m_data.data(), uint32_t(m_data.size())});
  }

  // allocate memory on the GPU
  void allocateGPU(size_t size, Magnum::GL::BufferUsage bufferUsage = Magnum::GL::BufferUsage::StreamDraw) {
    m_buffer.setData({nullptr, sizeof(T)*size}, bufferUsage);
    m_allocatedGPUsize = size;
  }

  // return (a copy of the) data on the GPU 
  auto getGPUData() {
    return Corrade::Containers::arrayCast<T>(m_buffer.data());
  }

  size_t getGPUSize() const { return m_allocatedGPUsize; }
  size_t getCPUSize() const { return m_data.size(); }
  // void allocateCPU(size_t n) { m_data.resize(n); }
  // void freeCPU() { m_data.clear(); }

  // NOTE: deprecate matrix-mapping of rows / entire arrays
  // replaced through map-functions in respective structs T

  // // access (type-casted part of) a data row in the cpu array
  // template <typename _scalar, int _nelems, int _offset = 0>
  // auto row(size_t row) {
  //   assert(row < m_data.size() && "trying to access row outside of vector");
  //   return Eigen::Map<Eigen::Matrix<_scalar, _nelems, 1, Eigen::ColMajor>,
  //                     Eigen::Unaligned>(
  //       reinterpret_cast<_scalar*>(&m_data[row]) + _offset, _nelems);
  // }
  // // access (type-casted part of) a data row in the cpu array
  // template <typename _scalar, int _nelems, int _offset = 0>
  // auto row(size_t row) const {
  //   assert(row < m_data.size() && "trying to access row outside of vector");
  //   return Eigen::Map<const Eigen::Matrix<_scalar, _nelems, 1, Eigen::ColMajor>,
  //                     Eigen::Unaligned>(
  //       reinterpret_cast<const _scalar*>(&m_data[row]) + _offset, _nelems);
  // }

  // template <typename _scalar, int _ncols, int _offset = 0>
  // auto matrixView() {
  //   // assert vec not empty?
  //   // ...

  //   // estimate stride based on struct memory size compared to scalar size
  //   constexpr int _stride = sizeof(T) / sizeof(_scalar);

  //   // in case of single column have to use colmajor vector, so this expression
  //   // checks if _cols == 1 and defines either a rowmajor matrix type or a
  //   // colmajor vector
  //   typedef typename std::conditional<
  //       _ncols == 1,
  //       Eigen::Map<Eigen::Matrix<_scalar, Eigen::Dynamic, 1, Eigen::ColMajor>,
  //                  Eigen::Unaligned, Eigen::Stride<1, _stride>>,
  //       Eigen::Map<
  //           Eigen::Matrix<_scalar, Eigen::Dynamic, _ncols, Eigen::RowMajor>,
  //           Eigen::Unaligned, Eigen::Stride<_stride, 1>>>::type MapType;
  //   return MapType(reinterpret_cast<_scalar*>(&m_data[0]) + _offset, m_data.size(),
  //                  _ncols);
  // }
  // template <typename _scalar, int _ncols, int _offset = 0>
  // auto matrixView() const {
  //   constexpr int _stride = sizeof(T) / sizeof(_scalar);
  //   typedef typename std::conditional<
  //       _ncols == 1,
  //       Eigen::Map<const Eigen::Matrix<_scalar, Eigen::Dynamic, 1, Eigen::ColMajor>,
  //                  Eigen::Unaligned, Eigen::Stride<1, _stride>>,
  //       Eigen::Map<
  //           const Eigen::Matrix<_scalar, Eigen::Dynamic, _ncols, Eigen::RowMajor>,
  //           Eigen::Unaligned, Eigen::Stride<_stride, 1>>>::type MapType;
  //   return MapType(reinterpret_cast<const _scalar*>(&m_data[0]) + _offset, m_data.size(),
  //                  _ncols);
  // }

  VectorBuffer() {}
  VectorBuffer(VectorBuffer<T> && rhs) = default;
  // VectorBuffer(const VectorBuffer<T> & rhs) = default;

 private:
  Magnum::GL::Buffer m_buffer;
  size_t m_allocatedGPUsize = 0u;
  std::vector<T> m_data;
};

#endif  // __VECTORBUFFER__H__