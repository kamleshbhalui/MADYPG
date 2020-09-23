#ifndef __MESH__H__
#define __MESH__H__

#include "../EigenDefinitions.h"

bool barycentric_inside(const Vector3s& abc);
// compute signed angle from u to v over e, e has unit length and is orthogonal
// to u,v
scalar signed_angle(const Vector3s& u, const Vector3s& v, const Vector3s& e);

class Mesh {
 public:
  Mesh() {}
  ~Mesh() {}

  // MatrixGLi &getMSFaces() { return (Fms.rows() == 0 && F.rows() > 0) ? F :
  // Fms; } const MatrixGLi &getMSFaces() const { return (Fms.rows() == 0 &&
  // F.rows() > 0) ? F : Fms; } MatrixGLi &getWSFaces() { return F; } const
  // MatrixGLi &getWSFaces() const { return F; }

  Vector3s barycentric_ms(int tri, const Vector2s& p) const;
  void compute_invDm();
  void compute_v2f_map(bool shepard_weights = true);
  void compute_face_adjacency();
  void compute_face_data();  // normals and strains
  void compute_vertex_normals();
  void compute_vertex_strains();

  // private:
  MatrixGLf X, U;    // vertex positions and uv coordinates
  MatrixGLi F, Fms;  // tri or quad faces {[v0, v1, v2(, v3)],...} for
                     // worldspace faces and material space faces

  AlignedVector<Matrix2s> invDm;           // FEM-matrix (inverse of (e1, e2))
  AlignedVector<Vector3s> normals;         // world-space face normals
  AlignedVector<Vector3s> vertex_normals;  // world-space vertex normals
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