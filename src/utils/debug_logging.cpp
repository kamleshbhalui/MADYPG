#include "debug_logging.h"
#include <alloca.h>
// convenience type extensions
// std::ostream& Debug::operator<< (std::ostream& stream, const Magnum::Vector3&
// val) {
//   stream << "("<< val[0] << " " << val[1] << " " << val[2] << ")";
//   return stream;
// }
// std::ostream& Debug::operator<< (std::ostream& stream, const Magnum::Vector2&
// val) {
//   stream << "("<< val[0] << " " << val[1] << ")";
//   return stream;
// }

void Debug::logf(const char *format, ...) {
#if DEBUG_VERBOSE_LEVEL >= 3
  std::cout << "Log: ";
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
#endif
}

void Debug::warningf(const char *format, ...) {
#if DEBUG_VERBOSE_LEVEL >= 2
  std::cout << "Warning: ";
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
#endif
}

void Debug::errorf(const char *format, ...) {
// NOTE not just using "vfprintf(stderr,format, args);"
// because pybind redirects std::cerr, but not stderr
#if DEBUG_VERBOSE_LEVEL >= 1
  std::cerr << "Error: ";
  va_list args1, args2;
  va_start(args1, format);
  va_copy(args2, args1);
  // char buf[1 + vsnprintf(nullptr, 0, format, args1)];
  int buf_len = vsnprintf(nullptr, 0, format, args1);
  char *buf = (char *)alloca(buf_len); // allocate dynamically on stack, gets freed automatically
  va_end(args1);
  // vsnprintf(buf, sizeof(buf), format, args2);
  vsnprintf(buf, buf_len, format, args2);
  va_end(args2);
  std::cerr << buf;
#endif
}


bool Debug::msgassert(const std::string& msg, bool cond) {
  if(!cond)
    std::cerr << "Assert failed: " << msg << "\n";
  #ifndef NDEBUG // debug mode
  assert(cond);
  #endif
  return cond;
}