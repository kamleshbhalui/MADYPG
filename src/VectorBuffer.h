#ifndef __VECTORBUFFER__H__
#define __VECTORBUFFER__H__

#include <Magnum/GL/Buffer.h>
#include <Corrade/Containers/Array.h>
#include <assert.h>

// #include <Eigen/Dense>
// #include <Eigen/Core>
#include <type_traits>
// template <typename _scalar, int _ncols, int _offset>
// static auto map(std::vector<MyStruct>& data) {
//   //assert vec not empty

//   //estimate stride based on struct memory size compared to scalar size
//   constexpr int _stride = sizeof(MyStruct)/sizeof(_scalar);

//   // in case of single column have to use colmajor vector, so this expression
//   // checks if _cols == 1 and defines either a rowmajor matrix type or a
//   colmajor vector typedef typename std::conditional<_ncols == 1,
//   Eigen::Map<Eigen::Matrix<_scalar, Eigen::Dynamic, 1, Eigen::ColMajor>,
//   Eigen::Unaligned, Eigen::Stride<1, _stride>>,
//   Eigen::Map<Eigen::Matrix<_scalar, Eigen::Dynamic, _ncols, Eigen::RowMajor>,
//   Eigen::Unaligned, Eigen::Stride<_stride, 1>>>::type MapType; return
//   MapType(reinterpret_cast<_scalar*>(&data[0])+_offset, data.size(), _ncols);
// }
// template <typename _scalar, int _ncols, int _offset, typename T>
// auto structvec2map(std::vector<T>& data) {
//   //assert vec not empty
//   constexpr int _stride = sizeof(T)/sizeof(_scalar);
//   return Eigen::Map<Eigen::Matrix<_scalar, Eigen::Dynamic, _ncols,
//   Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<_stride,
//   1>>(reinterpret_cast<_scalar*>(&data[0])+_offset, data.size(), _ncols);
// }
// auto map = MyStruct::map<float,4,0>(data);
// Map<Vector3f>(&data.back().x).normalize();

// T is expected to be a struct type defining the per-element data
template <typename T>
class VectorBuffer {
 public:
  std::vector<T>& cpu() { return m_data; }
  const std::vector<T>& cpu() const { return m_data; }
  Magnum::GL::Buffer& gpu() { return m_buffer; }
  const Magnum::GL::Buffer& gpu() const { return m_buffer; }

  // upload the data from the cpu to the gpu buffer
  void bufferData(bool realloc, Magnum::GL::BufferUsage bufferUsage =
                                    Magnum::GL::BufferUsage::StreamDraw) {
    if (realloc) {
      m_buffer.setData({m_data.data(), uint32_t(m_data.size())}, bufferUsage);
      m_allocatedGPUsize = m_data.size();
    } else
      m_buffer.setSubData(0, {m_data.data(), uint32_t(m_data.size())});
  }
  void bufferData(Magnum::GL::BufferUsage bufferUsage =
                      Magnum::GL::BufferUsage::StreamDraw) {
    bool realloc = (m_allocatedGPUsize != m_data.size());
    bufferData(realloc, bufferUsage);
  }
  void allocateGPU(size_t size, Magnum::GL::BufferUsage bufferUsage = Magnum::GL::BufferUsage::StreamDraw) {
    m_buffer.setData({nullptr, sizeof(T)*size}, bufferUsage);
    m_allocatedGPUsize = size;
  }

  auto getGPUData() {
    return Corrade::Containers::arrayCast<T>(m_buffer.data());
  }

  size_t getGPUSize() const { return m_allocatedGPUsize; }
  size_t getCPUSize() const { return m_data.size(); }
  // void allocateCPU(size_t n) { m_data.resize(n); }
  // void freeCPU() { m_data.clear(); }

  // TODO DEPRECATE row & matrixView BELOW, REPLACE WITH EXPLICIT GETTER AND STRUCT MAPS

  // access (type-casted part of) a data row in the cpu array
  template <typename _scalar, int _nelems, int _offset = 0>
  auto row(size_t row) {
    assert(row < m_data.size() && "trying to access row outside of vector");
    return Eigen::Map<Eigen::Matrix<_scalar, _nelems, 1, Eigen::ColMajor>,
                      Eigen::Unaligned>(
        reinterpret_cast<_scalar*>(&m_data[row]) + _offset, _nelems);
  }
  // access (type-casted part of) a data row in the cpu array
  template <typename _scalar, int _nelems, int _offset = 0>
  auto row(size_t row) const {
    assert(row < m_data.size() && "trying to access row outside of vector");
    return Eigen::Map<const Eigen::Matrix<_scalar, _nelems, 1, Eigen::ColMajor>,
                      Eigen::Unaligned>(
        reinterpret_cast<const _scalar*>(&m_data[row]) + _offset, _nelems);
  }

  template <typename _scalar, int _ncols, int _offset = 0>
  auto matrixView() {
    // assert vec not empty?
    // ...

    // estimate stride based on struct memory size compared to scalar size
    constexpr int _stride = sizeof(T) / sizeof(_scalar);

    // in case of single column have to use colmajor vector, so this expression
    // checks if _cols == 1 and defines either a rowmajor matrix type or a
    // colmajor vector
    typedef typename std::conditional<
        _ncols == 1,
        Eigen::Map<Eigen::Matrix<_scalar, Eigen::Dynamic, 1, Eigen::ColMajor>,
                   Eigen::Unaligned, Eigen::Stride<1, _stride>>,
        Eigen::Map<
            Eigen::Matrix<_scalar, Eigen::Dynamic, _ncols, Eigen::RowMajor>,
            Eigen::Unaligned, Eigen::Stride<_stride, 1>>>::type MapType;
    return MapType(reinterpret_cast<_scalar*>(&m_data[0]) + _offset, m_data.size(),
                   _ncols);
  }
  template <typename _scalar, int _ncols, int _offset = 0>
  auto matrixView() const {
    constexpr int _stride = sizeof(T) / sizeof(_scalar);
    typedef typename std::conditional<
        _ncols == 1,
        Eigen::Map<const Eigen::Matrix<_scalar, Eigen::Dynamic, 1, Eigen::ColMajor>,
                   Eigen::Unaligned, Eigen::Stride<1, _stride>>,
        Eigen::Map<
            const Eigen::Matrix<_scalar, Eigen::Dynamic, _ncols, Eigen::RowMajor>,
            Eigen::Unaligned, Eigen::Stride<_stride, 1>>>::type MapType;
    return MapType(reinterpret_cast<const _scalar*>(&m_data[0]) + _offset, m_data.size(),
                   _ncols);
  }

  VectorBuffer() {}
  VectorBuffer(VectorBuffer<T> && rhs) = default;
  // VectorBuffer(const VectorBuffer<T> & rhs) = default;

 private:
  Magnum::GL::Buffer m_buffer;
  size_t m_allocatedGPUsize = 0u;
  std::vector<T> m_data;
};

#endif  // __VECTORBUFFER__H__