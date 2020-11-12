#include "Mesh.h"

#include <algorithm>

#include "../utils/debug_logging.h"  // DEBUG REMOVE
#include "../utils/threadutils.h"

bool barycentric_inside(const Vector3s &abc) {
  return abc[0] > 0 && abc[1] > 0 && abc[2] > 0;
}

scalar signed_angle(const Vector3s &u, const Vector3s &v, const Vector3s &e) {
  return std::atan2(e.dot(u.cross(v)), u.dot(v));
}

void Mesh::compute_invDm() {
  auto& invDmU_cpu = invDmU.cpu();
  invDmU_cpu.resize(Fms.getCPUSize());
  threadutils::parallel_for(size_t(0), invDmU_cpu.size(), [&](size_t i) {
    auto &ixs = Fms.cpu()[i];
    Matrix2s Dm;
    Vector2s U0 = U.row<float, 2>(ixs.v0);
    Dm.col(0) = U.row<float, 2>(ixs.v1) - U0;
    Dm.col(1) = U.row<float, 2>(ixs.v2) - U0;
    invDmU_cpu[i].mapDinv() = Dm.inverse();
    invDmU_cpu[i].mapU0() = U0;
  });
}

Vector3s Mesh::barycentric_ms(int tri, const Vector2s &p) const {
  const auto& invDmU_cpu = invDmU.cpu();
  Vector3s abc;
  // auto u0       = U.row<float, 2>(Fms.cpu()[tri].v0);
  const auto& iDmU = invDmU_cpu[tri];
  abc.tail<2>() = iDmU.mapDinv() * (p - iDmU.mapU0());
  abc[0]        = scalar(1) - abc[1] - abc[2];
  return abc;
}

void Mesh::compute_v2f_map(bool shepard_weights) {
  // doing this in materials space, i.e. using U and Fms
  auto &Uc   = U.cpu();
  auto &Fmsc = Fms.cpu();

  // in serial: push all face indices to their respective verts
  // NOTE maybe parallelizable with concurrent deques
  v2f.resize(Uc.size());
  int nfaces = int(Fmsc.size());
  for (int f = 0; f < nfaces; ++f) {
    auto face = Fmsc[f].map();
    for (size_t i = 0; i < 3; ++i) {
      int v = int(face(i));
      v2f[v].push_back(std::make_pair(f, scalar(1)));
    }
  }

  if (shepard_weights) {
    // Weighted Shepard interpolation from 'Phong Deformation' [James 2020]:
    // for vertex i, with incident faces k
    // with vertex-to-centroid vector r_ki = rhat_ki * |r_ki| // |rhat_ki| = 1
    //   (sum_k(rhat_ki rhat_ki^T) + eps I) lam = -sum_k(rhat_ki) // eps = 1
    //   w_ki' = (1 + lam . rhat_ki) / |r_ki|

    AlignedVector<Vector2s> centroids(Fmsc.size());
    scalar _third = scalar(1) / scalar(3);
    threadutils::parallel_for(size_t(0), centroids.size(), [&](size_t i) {
      auto ixs = Fmsc[i].map();
      centroids[i] =
          _third * (Uc[ixs[0]].map() + Uc[ixs[1]].map() + Uc[ixs[2]].map());
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
        Vector2s r_ki = centroids[v2f_i[k].first] - Uc[i].map();
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

  auto &Fmsc = Fms.cpu();

  Vector3i none;
  none << -1, -1, -1;
  Fmsadj.resize(
      Fmsc.size());  // NOTE as with many other inits, non parallel resize...
  threadutils::parallel_for(size_t(0), Fmsadj.size(), [&](size_t f) {
    Fmsadj[f].faces << -1, -1, -1;  // init as having no neighbor
  });

  // edge list: (ov1, ov2, f, fe)
  //  where ov1,ov2 are ordered vertex indices, f is the face index,
  //  and fe < 3 is the index of the edge within the face
  AlignedVector<Vector4i> edge_list;
  edge_list.resize(Fmsc.size() * 3);

  threadutils::parallel_for(size_t(0), Fmsc.size(), [&](size_t f) {
    auto &face = Fmsc[f];
    // store edge with sorted vertex indices
    edge_list[3 * f + 0] << std::min(face.v0, face.v1),
        std::max(face.v0, face.v1), f, 0;
    edge_list[3 * f + 1] << std::min(face.v1, face.v2),
        std::max(face.v1, face.v2), f, 1;
    edge_list[3 * f + 2] << std::min(face.v2, face.v0),
        std::max(face.v2, face.v0), f, 2;
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
    if (adj) {
      // faces ea[2] and eb[2] are adjacent
      // -> store face and (local) opposite vertex
      // for a face with vertices 0,1,2, we have opposite vertices n0,n1,n2
      // where ni is opposite the edge(i,(i+1)%3), or -1 if it doesnt exist
      Fmsadj[ea[2]].faces[ea[3]] = eb[2];            // adjacent face
      Fmsadj[ea[2]].opp[ea[3]]   = (eb[3] + 2) % 3;  // opp vert ix
      Fmsadj[eb[2]].faces[eb[3]] = ea[2];
      Fmsadj[eb[2]].opp[eb[3]]   = (ea[3] + 2) % 3;
      // NOTE this assumes that the vertex order within Fms and F is the same!
      // but given that we load faces from obj files as x0/u0 x1/u1 x2/u2
      // this should be fine
    }
  });
}

void Mesh::compute_face_data() {
  auto &Fc = F.cpu();
  auto Xc  = X.matrixView<float, 3>();
  auto& invDmU_cpu = invDmU.cpu();
  normals.cpu().resize(Fc.size());
  auto& defgrads = defF.cpu();
  defgrads.resize(Fc.size());
  // normals.resize(Fc.size());
  strains.resize(Fc.size());
  threadutils::parallel_for(size_t(0), normals.cpu().size(), [&](size_t f) {
    auto ixs     = Fc[f].map();
    Vector3s e01 = Xc.row(ixs[1]) - Xc.row(ixs[0]);
    Vector3s e02 = Xc.row(ixs[2]) - Xc.row(ixs[0]);

    Vector3s n  = (e01.cross(e02));
    scalar invA = 1 / n.norm();
    n *= invA;
    normals.cpu()[f].map() = n;

    MatrixNMs<3, 2> defoF;
    defoF.col(0) = e01;
    defoF.col(1) = e02;
    defoF *= invDmU_cpu[f].mapDinv();
    Vector6s &s = strains[f];

    defgrads[f].map() = defoF;

    // in plane
    scalar C00 = defoF.col(0).squaredNorm();
    scalar C01 = defoF.col(0).dot(defoF.col(1));
    scalar C11 = defoF.col(1).squaredNorm();
    s[0]       = std::sqrt(C00) - 1;          // s_x
    s[2]       = std::sqrt(C11) - 1;          // s_y
    s[1]       = C01 / std::sqrt(C00 * C11);  // s_a

    // bending
    // II = F.T Lam F = sum_i thetai / (2 A li) F^T ti o F^t ti  , with ti of
    // length li
    const auto &adj = Fmsadj[f];
    s.tail<3>().setZero();
    for (int i = 0; i < 3; ++i) {
      if (adj.faces[i] < 0)
        continue;
      Vector3s ei = Xc.row(ixs[(i + 1) % 3]) - Xc.row(ixs[i]);
      // NOTE alternatively could precompute normals and then just look up
      // neighbor normal, but then might also want to precompute invA such as to
      // not have to compute normal again
      Vector3s ni =
          Vector3s(Xc.row(Fc[adj.faces[i]].map()[adj.opp[i]]) - Xc.row(ixs[i]))
              .cross(ei);
      Vector2s FTti = defoF.transpose() * ei.cross(n);
      scalar invli  = 1 / ei.norm();
      scalar theta  = -signed_angle(n, ni, ei * invli);
      scalar c      = theta * scalar(0.5) * invA * invli;
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
  const auto &fdata = normals.cpu();
  auto &vdata       = vertex_normals.cpu();

  vdata.resize(v2f.size());

  threadutils::parallel_for(size_t(0), v2f.size(), [&](size_t i) {
    auto v = vdata[i].map();
    v.setZero();
    for (const auto &fw : v2f[i]) {
      v += fdata[fw.first].map() * fw.second;
    }
    v.normalize();  // NOTE: could consider not normalizing until
                           // barycentric interpolation?
  });
}

void Mesh::compute_vertex_defF() {
  const auto &fdata = defF.cpu();
  auto &vdata       = vertex_defF.cpu();

  vdata.resize(v2f.size());

  threadutils::parallel_for(size_t(0), v2f.size(), [&](size_t i) {
    auto v = vdata[i].map();
    v.setZero();
    for (const auto &fw : v2f[i]) {
      v += fdata[fw.first].map() * fw.second;
    }
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
