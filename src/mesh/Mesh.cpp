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
  abc.tail<2>() = invDm[tri] * (p - U.block<1, 2>(Fms(tri, 0), 0).transpose());
  abc[0]        = scalar(1) - abc[1] - abc[2];
  return abc;
}

void Mesh::compute_face_normals() {
  normals.resize(F.rows());
  threadutils::parallel_for(size_t(0), normals.size(), [&](size_t i) {
    auto ixs    = F.row(i);
    Vector3s e0 = X.row(ixs[1]) - X.row(ixs[0]);
    Vector3s e1 = X.row(ixs[2]) - X.row(ixs[0]);
    normals[i]  = (e0.cross(e1)).normalized();
  });
}

void Mesh::compute_v2f_map(bool shepard_weights) {
  // doing this in materials space, i.e. using U and Fms

  // in serial: push all face indices to their respective verts
  // NOTE maybe parallelizable with concurrent deques
  v2f.resize(U.rows());
  for (int f = 0; f < int(Fms.rows()); ++f) {
    for (size_t i = 0; i < 3; ++i) {
      int v = int(Fms(f, i));
      v2f[v].push_back(std::make_pair(f, scalar(1)));
    }
  }

  if (shepard_weights) {
    // Weighted Shepard interpolation from 'Phong Deformation' [James 2020]:
    // for vertex i, with incident faces k
    // with vertex-to-centroid vector r_ki = rhat_ki * |r_ki| // |rhat_ki| = 1
    //   (sum_k(rhat_ki rhat_ki^T) + eps I) lam = -sum_k(rhat_ki) // eps = 1
    //   w_ki' = (1 + lam . rhat_ki) / |r_ki|

    AlignedVector<Vector2s> centroids(Fms.rows());
    scalar _third = scalar(1) / scalar(3);
    threadutils::parallel_for(size_t(0), centroids.size(), [&](size_t i) {
      auto ixs     = Fms.row(i);
      centroids[i] = _third * (U.row(ixs[0]) + U.row(ixs[1]) + U.row(ixs[2]));
    });

    // compute weights
    scalar eps = 1;
    threadutils::parallel_for(size_t(0), v2f.size(), [&](size_t i) {
      auto& v2f_i = v2f[i];

      // compute vertex-to-centroid vector data and assemble system
      Matrix2s A = Matrix2s::Identity() * eps;
      Vector2s b = Vector2s::Zero();
      AlignedVector<Vector2s> rhats_i(v2f_i.size());
      std::vector<scalar> rnorms_i(v2f_i.size());
      for (size_t k = 0; k < v2f_i.size(); ++k) {
        Vector2s r_ki = centroids[v2f_i[k].first] - U.row(i).transpose();
        rnorms_i[k]   = r_ki.norm();
        r_ki /= rnorms_i[k];
        rhats_i[k] = r_ki;

        A += r_ki * r_ki.transpose();
        b -= r_ki;
      }

      Vector2s lam = A.inverse() * b;

      for (size_t k = 0; k < v2f_i.size(); ++k) {
        scalar w_ki = (1 + lam.dot(rhats_i[k])) / rnorms_i[k];
        v2f_i[k].second *= w_ki;
      }
    });
  }

  // normalize weights
  threadutils::parallel_for(size_t(0), v2f.size(), [&](size_t i) {
    scalar Wtotal = 0.0;
    for (const auto& fw : v2f[i]) {
      Wtotal += fw.second;
    }
    scalar wnorm = scalar(1) / Wtotal;
    assert(v2f[i].size() > 0);  // vertex without any incident face
    for (auto& fw : v2f[i]) {
      fw.second *= wnorm;
    }
  });
}

void Mesh::compute_vertex_normals() {
  const auto& fdata = normals;
  auto& vdata       = vertex_normals;

  vdata.resize(v2f.size());

  threadutils::parallel_for(size_t(0), v2f.size(), [&](size_t i) {
    Vector3s n = Vector3s::Zero();
    for (const auto& fw : v2f[i]) {
      n += fdata[fw.first] * fw.second;
    }
    vdata[i] = n.normalized();  // NOTE: could consider not normalizing until
                                // barycentric interpolation?
  });
}
