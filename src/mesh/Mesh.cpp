#include "Mesh.h"

#include "../utils/threadutils.h"

bool barycentric_inside(const Vector3s& abc) {
  return abc[0] > 0 && abc[1] > 0 && abc[2] > 0;
}

void Mesh::compute_invDm() {
  invDm.resize(Fms.rows());
  threadutils::parallel_for(size_t(0), invDm.size(), [&](size_t i) {
    auto ixs = Fms.row(i);
    Matrix2s Dm;
    Dm.col(0) = U.row(ixs[1]) - U.row(ixs[0]);
    Dm.col(1) = U.row(ixs[2]) - U.row(ixs[0]);
    invDm[i]  = Dm.inverse();
  });
}

Vector3s Mesh::barycentric_ms(int tri, const Vector2s& p) const {
  Vector3s abc;
  abc.tail<2>() = invDm[tri] * (p - U.block<1,2>(Fms(tri, 0),0).transpose());
  abc[0]        = scalar(1) - abc[1] - abc[2];
  return abc;
}
