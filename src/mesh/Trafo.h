#ifndef __TRAFO__H__
#define __TRAFO__H__

#include "../EigenDefinitions.h"

// arcsim-like obstacle transformation struct
struct Trafo {
  Vector4s axis_angle;  // angle, axis
  scalar scale;
  Vector3s translation;
  Trafo() {
    axis_angle << 0, 0, 0, 1;
    translation.setZero();
    scale = 1;
  }
};

#endif  // __TRAFO__H__
