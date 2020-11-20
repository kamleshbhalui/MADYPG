#ifndef __MESH__H__
#define __MESH__H__

#include "../EigenDefinitions.h"
#include "../VectorBuffer.h"

bool barycentric_inside(const Vector3s& abc);
// compute signed angle from u to v over e, e has unit length and is orthogonal
// to u,v
scalar signed_angle(const Vector3s& u, const Vector3s& v, const Vector3s& e);


class Mesh {
 public:
struct MSVertex {
  float u,v;
  auto map() {
    return Eigen::Map<Eigen::Matrix<float, 2, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this), 2);
  }
  auto map() const {
    return Eigen::Map<const Eigen::Matrix<float, 2, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<const float*>(this), 2);
  }
};
struct WSVertex {
  float x,y,z;
  auto map() {
    return Eigen::Map<Eigen::Matrix<float, 3, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this), 3);
  }
  auto map() const {
    return Eigen::Map<const Eigen::Matrix<float, 3, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<const float*>(this), 3);
  }
};
struct Face {
  uint32_t v0,v1,v2;
  auto map() {
    return Eigen::Map<Eigen::Matrix<uint32_t, 3, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<uint32_t*>(this), 3);
  }
  auto map() const {
    return Eigen::Map<const Eigen::Matrix<uint32_t, 3, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<const uint32_t*>(this), 3);
  }
};
struct DinvU {
  float Dinv11,Dinv21,Dinv12,Dinv22,U0x,U0y;
  auto mapDinv() {
    return Eigen::Map<Eigen::Matrix<float, 2, 2, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this), 2, 2);
  }
  auto mapDinv() const {
    return Eigen::Map<const Eigen::Matrix<float, 2, 2, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<const float*>(this), 2, 2);
  }
  auto mapU0() {
    return Eigen::Map<Eigen::Matrix<float, 2, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this)+4, 2, 1);
  }
  auto mapU0() const {
    return Eigen::Map<const Eigen::Matrix<float, 2, 1, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<const float*>(this)+4, 2, 1);
  }
};
struct FDefo {
  float F11,F12,F21,F22,F31,F32;
  auto map() {
    return Eigen::Map<Eigen::Matrix<float, 3, 2, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<float*>(this), 3, 2);
  }
  auto map() const {
    return Eigen::Map<const Eigen::Matrix<float, 3, 2, Eigen::ColMajor>, Eigen::Unaligned>(reinterpret_cast<const float*>(this), 3, 2);
  }
};

  Mesh() {}
  ~Mesh() {}
  Mesh(Mesh && rhs) = default;


  // MatrixGLi &getMSFaces() { return (Fms.rows() == 0 && F.rows() > 0) ? F :
  // Fms; } const MatrixGLi &getMSFaces() const { return (Fms.rows() == 0 &&
  // F.rows() > 0) ? F : Fms; } MatrixGLi &getWSFaces() { return F; } const
  // MatrixGLi &getWSFaces() const { return F; }

  Vector3s barycentric_ms(int tri, const Vector2s& p) const;
  void compute_invDm();
  void compute_v2f_map(bool shepard_weights = true);
  void compute_face_adjacency();
  void compute_face_data(float svdclamp);  // normals and strains
  void compute_vertex_normals();
  void compute_vertex_defF();
  void compute_vertex_strains();

  bool empty() const { return (U.cpu().size() == 0) || (Fms.cpu().size() == 0); }

  VectorBuffer<WSVertex> X;
  VectorBuffer<MSVertex> U;
  VectorBuffer<Face> F, Fms;

  // // private:
  // MatrixGLf X, U;    // vertex positions and uv coordinates
  // MatrixGLi F, Fms;  // tri or quad faces {[v0, v1, v2(, v3)],...} for
                     // worldspace faces and material space faces
// AlignedVector<Matrix2s> invDm; 
  VectorBuffer<DinvU> invDmU; // FEM-matrix (inverse of (e1, e2))
  // AlignedVector<Vector3s> normals;         // world-space face normals
  // AlignedVector<Vector3s> vertex_normals;  // world-space vertex normals
  VectorBuffer<WSVertex> normals;
  VectorBuffer<WSVertex> vertex_normals;
  VectorBuffer<FDefo> defF;
  VectorBuffer<FDefo> vertex_defF;
  AlignedVector<Vector6s> strains;
  AlignedVector<Vector6s> vertex_strains;
  std::vector<std::deque<std::pair<int, scalar>>> v2f;

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