#ifndef __FRAME__H__
#define __FRAME__H__

#include "../EigenDefinitions.h"
#include "Trafo.h"

// single frame of animation
// including cloth mesh, obstacle meshes, obstacle transformations
// empty matrices correspond to no change from previous frames
struct Frame {
  MatrixGLf cloth_V, cloth_U;
  MatrixGLi cloth_F, cloth_Fms;
  std::vector<MatrixGLf> obs_V;
  std::vector<MatrixGLi> obs_F;
  std::vector<Trafo> obs_trafo;
};

#endif  // __FRAME__H__