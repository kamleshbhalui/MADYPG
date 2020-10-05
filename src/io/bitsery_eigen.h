#ifndef BITSERY_EXT_EIGEN_H
#define BITSERY_EXT_EIGEN_H

#include <bitsery/details/adapter_common.h>
#include <bitsery/details/serialization_common.h>
#include <bitsery/traits/core/traits.h>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace bitsery {
namespace ext {
namespace Eigen {
class Matrix {
 public:
  template <typename Ser, typename Scalar, int _Rows, int _Cols, int _Options,
            int _MaxRows, int _MaxCols, typename Fnc>
  void serialize(Ser &ser,
                 const ::Eigen::Matrix<Scalar, _Rows, _Cols, _Options, _MaxRows,
                                       _MaxCols> &matrix,
                 Fnc &&fnc) const {
    uint32_t rows  = matrix.rows();
    uint32_t cols  = matrix.cols();
    uint32_t elems = rows * cols;

    auto &writer = ser.adapter();
    writer.template writeBytes<4>(static_cast<uint32_t>(rows));
    writer.template writeBytes<4>(static_cast<uint32_t>(cols));

    static_assert(details::IsFundamentalType<Scalar>::value,
                  "Value must be integral, float or enum type.");
    using TValue = typename details::IntegralFromFundamental<Scalar>::TValue;
    writer.template writeBuffer<sizeof(TValue)>(
        reinterpret_cast<const TValue *>(matrix.data()), elems);
  }

  template <typename Des, typename Scalar, int _Rows, int _Cols, int _Options,
            int _MaxRows, int _MaxCols, typename Fnc>
  void deserialize(Des &des,
                   ::Eigen::Matrix<Scalar, _Rows, _Cols, _Options, _MaxRows,
                                   _MaxCols> &matrix,
                   Fnc &&fnc) const {
    auto &reader  = des.adapter();
    uint32_t rows = 0u, cols = 0u;
    reader.template readBytes<4>(rows);
    reader.template readBytes<4>(cols);
    uint32_t elems = rows * cols;

    matrix.resize(rows, cols);
    static_assert(details::IsFundamentalType<Scalar>::value,
                  "Value must be integral, float or enum type.");
    using TValue = typename details::IntegralFromFundamental<Scalar>::TValue;
    reader.template readBuffer<sizeof(TValue)>(
        reinterpret_cast<TValue *>(matrix.data()), elems);
  }

 private:
};

}  // namespace Eigen
}  // namespace ext

namespace traits {
template <typename T>
struct ExtensionTraits<ext::Eigen::Matrix, T> {
  using TValue                                = void;
  static constexpr bool SupportValueOverload  = false;
  static constexpr bool SupportObjectOverload = true;
  static constexpr bool SupportLambdaOverload = false;
};
}  // namespace traits

}  // namespace bitsery

#endif  // BITSERY_EXT_EIGEN_H
