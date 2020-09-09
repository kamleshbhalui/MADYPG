#ifndef MASTERGIT_H
#define MASTERGIT_H

#include <Magnum/Math/TypeTraits.h>
#include <type_traits>
#include <Magnum/Tags.h>

namespace Magnum {
namespace Math {

template <class T>
inline typename std::enable_if<IsScalar<T>::value, T>::type fmod(T a, T b) {
  return T(std::fmod(UnderlyingTypeOf<T>(a), UnderlyingTypeOf<T>(b)));
}

template <std::size_t size, class T>
inline Vector<size, T> fmod(const Vector<size, T>& a,
                            const Vector<size, T>& b) {
  Vector<size, T> out{Magnum::NoInit};
  for (std::size_t i = 0; i != size; ++i) out[i] = Math::fmod(a[i], b[i]);
  return out;
}
}
}
#endif  // MASTERGIT_H
