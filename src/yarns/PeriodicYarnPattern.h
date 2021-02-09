#ifndef __PERIODICYARNPATTERN__H__
#define __PERIODICYARNPATTERN__H__

#include <string>

#include "../EigenDefinitions.h"

// Periodic yarn pattern data, with vertices, periodic edges, etc.
struct PeriodicYarnPattern {
  // load pattern data from file
  void deserialize(const std::string &filename);
  // 'choose' vertex location such that the pattern is contained in a rectangle
  void rectangulize();
  // recompute a vertex-edge table (for each vertex, index of incident edges)
  void recompute_VE_table();

  // debugging functionality to generate yarnlists of yarns without periodic
  // edges (concatenated and separated by uint32_t::max)
  std::vector<uint32_t> compute_simple_yarns();

  bool isPeriodicEdge(int eix);

  float px = 0;
  float py = 0;
  float r  = 0;
  MatrixGLf Q;             // vertex data [x y z t]
  MatrixXXRMi E;           // periodic edges [v0, v1, di, dj]
  std::vector<scalar> RL;  // restlengths (NOTE: stored per vertex outgoing
                           // edge, not per edge index!)
  MatrixGLf RefD1;  // ref directors (NOTE: stored per vertex outgoing edge, not
                    // per edge index!)
  Vector2s Qmin;
  MatrixXXRMi VE;  // vertex edge table [eix_prev, eix_next]
};
using PYP = PeriodicYarnPattern;  // Periodic Yarn Pattern (alias)

#endif  // __PERIODICYARNPATTERN__H__