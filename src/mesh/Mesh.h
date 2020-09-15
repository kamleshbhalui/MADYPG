#ifndef __MESH__H__
#define __MESH__H__

#include "../EigenDefinitions.h"

class Mesh
{
public:
  Mesh() {}
  ~Mesh() {}

  // MatrixGLi &getMSFaces() { return (Fms.rows() == 0 && F.rows() > 0) ? F : Fms; }
  // const MatrixGLi &getMSFaces() const { return (Fms.rows() == 0 && F.rows() > 0) ? F : Fms; }
  // MatrixGLi &getWSFaces() { return F; }
  // const MatrixGLi &getWSFaces() const { return F; }

  // private:
  MatrixGLf X, U;   // vertex positions and uv coordinates
  MatrixGLi F, Fms; // tri or quad faces {[v0, v1, v2(, v3)],...} for worldspace faces and material space faces
};

#endif // __MESH__H__