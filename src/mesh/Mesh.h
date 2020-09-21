#ifndef __MESH__H__
#define __MESH__H__

#include "../EigenDefinitions.h"

bool barycentric_inside(const Vector3s& abc);

class Mesh
{
public:
  Mesh() {}
  ~Mesh() {}

  // MatrixGLi &getMSFaces() { return (Fms.rows() == 0 && F.rows() > 0) ? F : Fms; }
  // const MatrixGLi &getMSFaces() const { return (Fms.rows() == 0 && F.rows() > 0) ? F : Fms; }
  // MatrixGLi &getWSFaces() { return F; }
  // const MatrixGLi &getWSFaces() const { return F; }

  void compute_invDm();
  Vector3s barycentric_ms(int tri, const Vector2s& p) const;
  void compute_face_normals();
  void compute_v2f_map(bool shepard_weights = true);
  void compute_vertex_normals();

  // private:
  MatrixGLf X, U;   // vertex positions and uv coordinates
  MatrixGLi F, Fms; // tri or quad faces {[v0, v1, v2(, v3)],...} for worldspace faces and material space faces

  AlignedVector<Matrix2s> invDm; // FEM-matrix (inverse of (e1, e2))
  AlignedVector<Vector3s> normals; // world-space face normals
  AlignedVector<Vector3s> vertex_normals; // world-space vertex normals
  AlignedVector<Vector6s> strains;
  AlignedVector<Vector6s> vertex_strains;
  std::vector<std::deque<std::pair<int, scalar>>> v2f;

private:


};

#endif // __MESH__H__