#ifndef __PERIODICYARNPATTERN__H__
#define __PERIODICYARNPATTERN__H__

#include <string>
#include "../EigenDefinitions.h"

struct PeriodicYarnPattern
{
  void deserialize(const std::string &filename);
  void rectangulize();

  bool isPeriodicEdge(int eix);

  // debugging functionality to generate uint32_t::max separated yarnlists of yarns without periodic edges
  std::vector<uint32_t> compute_simple_yarns();


  // compute_periodic_yarns() USED IN MODEL FOR DEFO LOOKUP

  float px = 0;
  float py = 0;
  float r = 0;
  MatrixGLf Q;  // vertex data [x y z t]
  MatrixXXRMi E;  // periodic edges [v0, v1, di, dj]
  MatrixGLf RL; // restlengths
  Vector2s Qmin;
};
using PYP = PeriodicYarnPattern; // short alias

#endif // __PERIODICYARNPATTERN__H__