#include "Mesh.h"

#include <algorithm>

#include "../utils/threadutils.h"
#include "../utils/debug_logging.h" // DEBUG REMOVE

bool barycentric_inside(const Vector3s &abc) {
  return abc[0] > 0 && abc[1] > 0 && abc[2] > 0;
}


scalar signed_angle(const Vector3s& u, const Vector3s& v, const Vector3s& e) {
  return std::atan2(e.dot(u.cross(v)), u.dot(v));
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

Vector3s Mesh::barycentric_ms(int tri, const Vector2s &p) const {
  Vector3s abc;
  abc.tail<2>() = invDm[tri] * (p - U.block<1, 2>(Fms(tri, 0), 0).transpose());
  abc[0]        = scalar(1) - abc[1] - abc[2];
  return abc;
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
      auto &v2f_i = v2f[i];

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
    for (const auto &fw : v2f[i]) {
      Wtotal += fw.second;
    }
    scalar wnorm = scalar(1) / Wtotal;
    assert(v2f[i].size() > 0);  // vertex without any incident face
    for (auto &fw : v2f[i]) {
      fw.second *= wnorm;
    }
  });
}

void Mesh::compute_face_adjacency() {
  // Following "A note on a linear time algorithm for constructing adjacency
  // graphs of 3D FEA data"

  Vector3i none;
  none << -1, -1, -1;
  Fmsadj.resize(Fms.rows(),
                none);  // NOTE as with many other inits, non parallel resize...

  // edge list: (ov1, ov2, f, fe)
  //  where ov1,ov2 are ordered vertex indices, f is the face index,
  //  and fe < 3 is the index of the edge within the face
  AlignedVector<Vector4i> edge_list;
  edge_list.resize(Fms.rows() * 3);

  threadutils::parallel_for(0, int(Fms.rows()), [&](int f) {
    // store edge with sorted vertex indices
    std::vector<int> vs = {int(Fms(f, 0)), int(Fms(f, 1)), int(Fms(f, 2))};
    edge_list[3 * f + 0] << std::min(vs[0], vs[1]), std::max(vs[0], vs[1]), f,
        0;
    edge_list[3 * f + 1] << std::min(vs[1], vs[2]), std::max(vs[1], vs[2]), f,
        1;
    edge_list[3 * f + 2] << std::min(vs[2], vs[0]), std::max(vs[2], vs[0]), f,
        2;
  });

  // sort edge list, s.t. same edges of neighboring faces are adjacent in list
  threadutils::parallel_sort(edge_list,
                             [](const Vector4i &a, const Vector4i &b) {
                               return a[0] == b[0] ? a[1] < b[1] : a[0] < b[0];
                             });

  // compare list-adjacent edges to check for adjacent faces
  threadutils::parallel_for(0, int(edge_list.size()) - 1, [&](int i) {
    const auto &ea = edge_list[i];
    const auto &eb = edge_list[i + 1];
    bool adj       = ea[0] == eb[0] && ea[1] == eb[1];  // matching vertices
    if (adj) {  // faces ea[2] and eb[2] are adjacent
      // store opposite vertex
      // for a face with vertices 0,1,2, we have opposite vertices n0,n1,n2
      // where ni is opposite the edge(i,(i+1)%3), or -1 if it doesnt exist
      Fmsadj[ea[2]][ea[3]] = int(Fms(eb[2], (eb[3] + 2) % 3));
      Fmsadj[eb[2]][eb[3]] = int(Fms(ea[2], (ea[3] + 2) % 3));
    }
  });
}

void Mesh::compute_face_data() {
  normals.resize(F.rows());
  strains.resize(F.rows());
  threadutils::parallel_for(size_t(0), normals.size(), [&](size_t f) {
    auto ixs    = F.row(f);
    Vector3s e01 = X.row(ixs[1]) - X.row(ixs[0]);
    Vector3s e02 = X.row(ixs[2]) - X.row(ixs[0]);

    Vector3s n  = (e01.cross(e02));
    scalar invA = 1 / n.norm();
    n *= invA;
    normals[f] = n;

    MatrixNMs<3, 2> F_;
    F_.col(0) = e01;
    F_.col(1) = e02;
    F_ *= invDm[f];
    Vector6s &s = strains[f];

    // in plane
    scalar C00 = F_.col(0).squaredNorm();
    scalar C01 = F_.col(0).dot(F_.col(1));
    scalar C11 = F_.col(1).squaredNorm();
    s[0]       = std::sqrt(C00) - 1;          // s_x
    s[2]       = std::sqrt(C11) - 1;          // s_y
    s[1]       = C01 / std::sqrt(C00 * C11);  // s_a

    // bending
    // II = F.T Lam F = sum_i thetai / (2 A li) F^T ti o F^t ti  , with ti of
    // length li
    const auto &adj = Fmsadj[f];
    s.tail<3>().setZero();
    for (int i = 0; i < 3; ++i) {
      if (adj[i] < 0)
        continue;
      Vector3s ei = X.row(ixs[(i+1)%3])-X.row(ixs[i]);
      // NOTE alternatively could precompute normals and then just look up neighbor normal, but then might also want to precompute invA such as to not have to compute normal again 
      Vector3s ni = Vector3s(X.row(adj[i])-X.row(ixs[i])).cross(ei);
      Vector2s FTti = F_.transpose() * ei.cross(n);
      scalar invli = 1 / ei.norm();
      scalar theta = -signed_angle(n,ni,ei * invli);
      scalar c = theta * scalar(0.5) * invA * invli;
      #warning check II sign against mathematica / arcsim. 
      // NOTE @ warning: seems like II should be positive for bowl shape
      // it might be that its just about the order of stuff passed to signed_angle
      // if i change the sign here, also change it in python!
      s[3] += c * FTti[0] * FTti[0];
      s[4] += c * FTti[0] * FTti[1];
      s[5] += c * FTti[1] * FTti[1];
    }

    // TODO consider eigendecomp of II already here
  });

  // Vector6s min = strains[0];
  // Vector6s max = strains[0];
  // for (size_t i = 0; i < strains.size(); i++) {
  //   for (size_t j = 0; j < 6; j++)
  //   {
  //     if (strains[i](j) < min(j))
  //       min(j) = strains[i](j);
  //     if (strains[i](j) > max(j))
  //       max(j) = strains[i](j);
  //   }
    
  // }
  // Debug::log("STRAIN MIN", min.transpose());
  // Debug::log("STRAIN MAX", max.transpose());
}

void Mesh::compute_vertex_normals() {
  const auto &fdata = normals;
  auto &vdata       = vertex_normals;

  vdata.resize(v2f.size());

  threadutils::parallel_for(size_t(0), v2f.size(), [&](size_t i) {
    vdata[i].setZero();
    for (const auto &fw : v2f[i]) {
      vdata[i] += fdata[fw.first] * fw.second;
    }
    vdata[i].normalize();  // NOTE: could consider not normalizing until
                           // barycentric interpolation?
  });
}

void Mesh::compute_vertex_strains() {
  const auto &fdata = strains;
  auto &vdata       = vertex_strains;

  vdata.resize(v2f.size());

  threadutils::parallel_for(size_t(0), v2f.size(), [&](size_t i) {
    vdata[i].setZero();
    for (const auto &fw : v2f[i]) {
      vdata[i] += fdata[fw.first] * fw.second;
    }
  });
}
