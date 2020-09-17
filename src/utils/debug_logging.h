#ifndef _DEBUG_LOGGING_H_
#define _DEBUG_LOGGING_H_

#include <assert.h>
#include <iostream>
#include <vector>
#include <deque>
#include <stdarg.h>
#include <Eigen/Sparse>
#include "debug_macros.h"

#ifndef DEBUG_VERBOSE_LEVEL
// #ifdef NDEBUG // release mode
#define DEBUG_VERBOSE_LEVEL 3 // define from compiler
#endif

namespace Debug {

// convenience type extensions
// std::ostream& operator<< (std::ostream& stream, const Magnum::Vector3& val);
// std::ostream& operator<< (std::ostream& stream, const Magnum::Vector2& val);

// printf wrapper
void logf(const char *format, ...);
void warningf(const char *format, ...);
void errorf(const char *format, ...);

// space-separated outstream
template <typename... Args> void log(Args &&... args) {
#if DEBUG_VERBOSE_LEVEL >= 3
  std::cout << "Log: ";
  int dummy[] = {0, (std::cout << std::forward<Args>(args) << " ", 0)...};
  DECLARE_UNUSED(dummy)
  std::cout << std::endl;
#endif
}

template <typename... Args> void warning(Args &&... args) {
#if DEBUG_VERBOSE_LEVEL >= 2
  std::cout << "Warning: ";
  int dummy[] = {0, (std::cout << std::forward<Args>(args) << " ", 0)...};
  DECLARE_UNUSED(dummy)
  std::cout << std::endl;
#endif
}

template <typename... Args> void error(Args &&... args) {
#if DEBUG_VERBOSE_LEVEL >= 1
  std::cerr << "Error: ";
  int dummy[] = {0, (std::cerr << std::forward<Args>(args) << " ", 0)...};
  DECLARE_UNUSED(dummy)
  std::cerr << std::endl;
#endif
}

template <typename T>
std::ostream& operator<< (std::ostream& stream, const std::vector<T>& vec) {
  stream << "[";
  if(vec.size() < 1)
    return stream << "]";
  for (int i = 0; i < int(vec.size()) - 1; i++) {
    stream << vec[i] << ", ";
  }
  stream << vec.back() << "]";
  return stream;
}
template <typename T>
std::ostream& operator<< (std::ostream& stream, const std::deque<T>& vec) {
  stream << "[";
  if(vec.size() < 1)
    return stream << "]";
  for (int i = 0; i < int(vec.size()) - 1; i++) {
    stream << vec[i] << ", ";
  }
  stream << vec.back() << "]";
  return stream;
}

template <typename Scalar, typename StorageIndex>
std::ostream& operator<< (std::ostream& stream, const Eigen::Triplet<Scalar, StorageIndex>& triplet) {
  stream << "(" << triplet.row() << ", " << triplet.col() << ", " << triplet.value() << ")";
  return stream;
}

template <typename A, typename B>
std::ostream& operator<< (std::ostream& stream, const std::pair<A, B>& pair) {
  stream << "(" << pair.first << ", " << pair.second << ")";
  return stream;
}

void msgassert(const std::string& msg, bool cond);

} // namespace Debug

#endif /* end of include guard: _DEBUG_LOGGING_H_ */
