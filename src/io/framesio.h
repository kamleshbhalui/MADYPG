#ifndef __FRAMESIO__H__
#define __FRAMESIO__H__

#define MAX_OBS_COUNT 10000
#define MAX_FRAMES 10000
#define MAX_VERTICES_OR_FACES 10000000  // arbitrary 10mio, alt: std limit uint?

#include <bitsery/bitsery.h>
#include <bitsery/traits/vector.h>

#include "../mesh/Frame.h"
#include "bitsery_eigen.h"

namespace bitsery {
template <typename S>
void serialize(S& ser, Trafo& o) {
  ser.template value<sizeof(scalar)>(o.scale);
  ser.ext(o.translation, ext::Eigen::Matrix{});
  ser.ext(o.axis_angle, ext::Eigen::Matrix{});
}

template <typename S>
void serialize(S& ser, Frame& o) {
  ser.ext(o.cloth_V, ext::Eigen::Matrix{});
  ser.ext(o.cloth_F, ext::Eigen::Matrix{});
  ser.ext(o.cloth_U, ext::Eigen::Matrix{});
  ser.ext(o.cloth_Fms, ext::Eigen::Matrix{});

  ser.container(o.obs_trafo, MAX_OBS_COUNT);
  ser.container(o.obs_V, MAX_OBS_COUNT, [](S& s, MatrixGLf& matrix) {
    s.ext(matrix, ext::Eigen::Matrix{});
  });
  ser.container(o.obs_F, MAX_OBS_COUNT, [](S& s, MatrixGLi& matrix) {
    s.ext(matrix, ext::Eigen::Matrix{});
  });
}

template <typename S>
void serialize(S& ser, std::vector<Frame>& o) {
  ser.container(o, MAX_FRAMES, [](S& s, Frame& frame) { s.object(frame); });
}
}  // namespace bitsery

bool deserialize_frames(const std::string& filepath,
                        std::vector<Frame>& frames);

#endif  // __FRAMESIO__H__