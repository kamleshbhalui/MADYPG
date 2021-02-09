#ifndef __MESH__H__
#define __MESH__H__

#include "../EigenDefinitions.h"
#include "../VectorBuffer.h"

bool barycentric_inside(const Vector3s& abc);
// compute signed angle from u to v over e, e has unit length and is orthogonal
// to u,v
scalar signed_angle(const Vector3s& u, const Vector3s& v, const Vector3s& e);

// Mesh class with VectorBuffers for material-space and world-space vertices
// and faces, per-face / per-vertex normals / fundamental forms / deformation
// gradients, etc.. as well as methods to compute those
class Mesh {
 public:
  // uv coords
  struct MSVertex {
    float u, v;
    auto map() {
      return Eigen::Map<Eigen::Matrix<float, 2, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<float*>(this), 2);
    }
    auto map() const {
      return Eigen::Map<const Eigen::Matrix<float, 2, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<const float*>(this),
                                          2);
    }
  };
  // position xyz
  struct WSVertex {
    float x, y, z;
    auto map() {
      return Eigen::Map<Eigen::Matrix<float, 3, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<float*>(this), 3);
    }
    auto map() const {
      return Eigen::Map<const Eigen::Matrix<float, 3, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<const float*>(this),
                                          3);
    }
  };
  // vertex indices: v0, v1, v2
  struct Face {
    uint32_t v0, v1, v2;
    auto map() {
      return Eigen::Map<Eigen::Matrix<uint32_t, 3, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<uint32_t*>(this), 3);
    }
    auto map() const {
      return Eigen::Map<const Eigen::Matrix<uint32_t, 3, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(
          reinterpret_cast<const uint32_t*>(this), 3);
    }
  };
  // FEM-matrix (inverse of (e1, e2))
  // and first mat-space vertex U0 in a triangle
  struct DinvU {
    float Dinv11, Dinv21, Dinv12, Dinv22, U0x, U0y;
    auto mapDinv() {
      return Eigen::Map<Eigen::Matrix<float, 2, 2, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<float*>(this), 2, 2);
    }
    auto mapDinv() const {
      return Eigen::Map<const Eigen::Matrix<float, 2, 2, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<const float*>(this),
                                          2, 2);
    }
    auto mapU0() {
      return Eigen::Map<Eigen::Matrix<float, 2, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<float*>(this) + 4, 2,
                                          1);
    }
    auto mapU0() const {
      return Eigen::Map<const Eigen::Matrix<float, 2, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(
          reinterpret_cast<const float*>(this) + 4, 2, 1);
    }
  };
  // Deformation Gradient F
  struct FDefo {
    float F11, F21, F31, F12, F22, F32;
    auto map() {
      return Eigen::Map<Eigen::Matrix<float, 3, 2, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<float*>(this), 3, 2);
    }
    auto map() const {
      return Eigen::Map<const Eigen::Matrix<float, 3, 2, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<const float*>(this),
                                          3, 2);
    }
  };
  // Ixx, Ixy, Iyy, IIxx, IIxy, IIyy;
  struct Strain {
    float Ixx, Ixy, Iyy, IIxx, IIxy, IIyy;
    auto map() {
      return Eigen::Map<Eigen::Matrix<float, 6, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<float*>(this), 6);
    }
    auto map() const {
      return Eigen::Map<const Eigen::Matrix<float, 6, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<const float*>(this),
                                          6);
    }
  };

  Mesh() {}
  ~Mesh() {}
  Mesh(Mesh&& rhs) = default;

  // compute barycentric coordinates of point p in material-space triangle with
  // index tri
  Vector3s barycentric_ms(int tri, const Vector2s& p) const;
  // compute FEM-matrix
  void compute_invDm();
  void compute_face_adjacency();
  // compute vertex-to-face weights
  void compute_v2f_map(bool shepard_weights = true);
  // compute normals, deformation gradient, I, II
  void compute_face_data();
  void compute_vertex_normals();
  void compute_vertex_defF();
  void compute_vertex_strains();

  bool empty() const {
    return (U.cpu().size() == 0) || (Fms.cpu().size() == 0);
  }

  VectorBuffer<WSVertex> X;
  VectorBuffer<MSVertex> U;
  VectorBuffer<Face> F, Fms;  // world- and material-space faces
  // NOTE: Fms is generally different from F due to mesh seams
  //  where Fms & U encode the individual flat pieces of a cloth
  //  F & X encode the current configuration _after gluing the pieces
  //  together_, ie in F some vertices are merged compared to Fms

  VectorBuffer<DinvU> invDmU;             // FEM-matrix (inverse of (e1, e2))
  VectorBuffer<WSVertex> normals;         // world-space face normals
  VectorBuffer<WSVertex> vertex_normals;  // world-space vertex normals
  VectorBuffer<FDefo> defF;               // per-face deformation gradient
  VectorBuffer<FDefo> vertex_defF;        // per-vertex deformation gradient
  VectorBuffer<Strain> strains;           // per-face I,II (flattened)
  VectorBuffer<Strain> vertex_strains;    // per-vertex I,II
  std::vector<std::deque<std::pair<int, scalar>>> v2f; // vertex-to-face weights

 private:
  struct FaceAdjacency {
    Vector3i faces;  // adjacent faces
    Vector3i opp;    // local opposite vertex index in resp. adjacent face
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };
  AlignedVector<FaceAdjacency> Fmsadj;  // (ordered) vertices opposite of the
                                        // face's edges (or -1 if non-existent)
};

#endif  // __MESH__H__