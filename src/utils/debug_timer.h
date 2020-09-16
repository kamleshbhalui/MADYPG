#ifndef __DEBUG_TIMER__H__
#define __DEBUG_TIMER__H__

#include <chrono>

typedef std::chrono::microseconds microseconds;
typedef std::chrono::milliseconds milliseconds;
typedef std::chrono::seconds seconds;

namespace Debug {
class Timer {
  using TimePoint = std::chrono::high_resolution_clock::time_point;

 public:
  Timer() {
    tick();
    _tock();
  }
  ~Timer() {}

  // set t_prev to now
  void tick() { t_prev = std::chrono::high_resolution_clock::now(); }

  // advance time and return duration
  template <typename T = std::chrono::seconds>
  auto tock() {
    _tock();
    return std::chrono::duration_cast<T>(dt).count();
  }

 private:
  TimePoint t_prev;
  std::chrono::duration<double> dt;  // defaults to duration in seconds

  // compute duration and set t_prev to now
  void _tock() {
    TimePoint t = std::chrono::high_resolution_clock::now();
    dt = std::chrono::duration_cast<std::chrono::duration<double>>(t - t_prev);
    t_prev = t;
  }
};
}  // namespace Debug

#endif  // __DEBUG_TIMER__H__