#ifndef _DEBUG_SINGLETON_H_
#define _DEBUG_SINGLETON_H_

#include <boost/format.hpp>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

// TODO maybe add second template parameter with default value
// e.g. Map<typename T, int ID = 0> such that we can have multiple maps
// -> Debug::get<bool> ("valuename")
// -> Debug::get<bool, 1> ("valuename") in different map

// https://stackoverflow.com/questions/1008019/c-singleton-design-pattern
namespace Debug {
template <typename A>
class Map {
 public:
  static Map &getInstance() {
    static Map instance;  // Guaranteed to be destroyed.
                          // Instantiated on first use.
    return instance;
  }

  // only private constructor
  Map(Map const &) = delete;
  void operator=(Map const &) = delete;

  static A &get(const std::string &id) { return getInstance().m_map[id]; }
  static A &get(const std::string &id, const A &defaultvalue) {
    auto &map = getInstance().m_map;

    auto got = map.find(id);
    if (got == map.end()) {
      auto &val = map[id];
      val       = defaultvalue;
      return val;
    } else {
      return got->second;
    }
  }
  static auto &getMap() { return getInstance().m_map; }

  // A &operator()(const std::string &id) { return m_map[id]; }

  static void print() {
    std::cout << "----------Map Entries---------------\n";
    for (const auto elem : getMap()) {
      std::cout << elem.first << " " << elem.second << "\n";
    }
    std::cout << "------------------------------------\n";
  }

 private:
  std::unordered_map<std::string, A> m_map;

  Map() {}
};

// Debug::get<bool>("valuename")
template <typename T>
static T &get(const std::string &id) {
  return Map<T>::get(id);
}
template <typename T>
static T &get(const std::string &id, const T &defaultvalue) {
  return Map<T>::get(id, defaultvalue);
}
//
// class FormattedOutput {
// public:
//   static FormattedOutput &getInstance() {
//     static FormattedOutput instance;
//     return instance;
//   }
//
//   // only private constructor
//   FormattedOutput(FormattedOutput const &) = delete;
//   void operator=(FormattedOutput const &) = delete;
//
//   static std::string &get(int id) { return getInstance().m_map[id]; }
//   static auto &getMap() { return getInstance().m_map; }
//
// private:
//   std::map<int, std::string> m_map;
//
//   FormattedOutput() {}
// };
//
// template <typename T>
// static void output(int id, const std::string &fmt, const T &value) {
//   FormattedOutput::get(id) = str(boost::format(fmt) % value);
// }

}  // namespace Debug

#endif /* end of include guard: _DEBUG_SINGLETON_H_ */
