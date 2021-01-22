#ifndef __DEBUG_TIMER__H__
#define __DEBUG_TIMER__H__

#include <chrono>
#include <string>
#include <unordered_map>
#include <deque>
#include <vector>

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

template <int N, typename duration_type>
class MovingAverageTimer {
 public:
  using TimePoint = std::chrono::high_resolution_clock::time_point;

  MovingAverageTimer() {
    tick();
    _tock();
  }
  ~MovingAverageTimer() {}

  // // explicitly add labeled timer, to enforce order
  // void declare(const std::string & label) {
  //   // m_map.try_emplace(label, {});
  //   auto& entry   = m_map[label];
  // }

  void tick() { t_prev = std::chrono::high_resolution_clock::now(); }

  void tock(const std::string& label) {
    int duration = _tock();
    tockDuration(label, duration);
  }

  void tockDuration(const std::string& label, int duration) { 
    auto& entry   = m_map[label];
    auto& A       = entry.value;
    auto& samples = entry.samples;

    int n = int(samples.size());
    if (n < N) {
      A = (A * n + duration) / (n + 1);  // incrementing average
    } else {
      A = A + (duration - samples.front()) * invN;  // moving average
      samples.pop_front();
    }
    samples.push_back(duration);
    tick();  // don't count tock execution in next timer usage
  }

  double getAverage(const std::string& label) const {
    // double a = 0.0;
    // for (auto& v : m_map[label].samples) {
    //   a += v;
    // }

    auto search = m_map.find(label);
    if (search != m_map.end()) {
        return search->second.value;
    } else {
        return -1;
    }
    // return m_map[label].value; // non-const lookup
  }

  std::vector<std::pair<std::string, double>> getAverageList() {
    std::vector<std::pair<std::string, double>> pairs;
    pairs.reserve(m_map.size());
    for (auto kv : m_map) {
      pairs.push_back(std::make_pair(kv.first,kv.second.value));
    }
    return pairs;
  }

 private:
  struct Average {
    double value = 0.0;  // ensure 0 init
    std::deque<int> samples;
  };

  std::unordered_map<std::string, Average> m_map;
  constexpr static double invN = 1.0 / N;
  TimePoint t_prev;
  std::chrono::duration<double> dt;

  auto _tock() {
    TimePoint t = std::chrono::high_resolution_clock::now();
    // auto dt = std::chrono::duration_cast<double>(t - t_prev);
    return std::chrono::duration_cast<duration_type>(t - t_prev).count();
  }
};
}  // namespace Debug

#endif  // __DEBUG_TIMER__H__