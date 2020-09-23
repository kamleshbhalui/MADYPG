#ifndef __PERIODICYARNPATTERN__H__
#define __PERIODICYARNPATTERN__H__

#include <string>
#include "../EigenDefinitions.h"

struct PeriodicYarnPattern
{
  void deserialize(const std::string &filename);
  void rectangulize();

  bool isPeriodicEdge(int eix);

  void recompute_VE_table();

  // debugging functionality to generate uint32_t::max separated yarnlists of yarns without periodic edges
  std::vector<uint32_t> compute_simple_yarns();

  void compute_parametric(); // TODO maybe do this in python model creation

  float px = 0;
  float py = 0;
  float r = 0;
  MatrixGLf Q;  // vertex data [x y z t]
  MatrixXXRMi E;  // periodic edges [v0, v1, di, dj]
  std::vector<scalar> RL; // restlengths
  Vector2s Qmin;
  std::vector<int> param_v2y;
  std::vector<scalar> param_v2t;
  std::vector<std::deque<int>> param_y2v;
  std::vector<std::vector<scalar>> param_y2t;
  MatrixXXRMi VE; // vertex edge table [eix_prev, eix_next]
};
using PYP = PeriodicYarnPattern; // short alias

#endif // __PERIODICYARNPATTERN__H__