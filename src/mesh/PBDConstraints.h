#ifndef __PBDCONSTRAINTS__H__
#define __PBDCONSTRAINTS__H__

#include <Simulation/Constraints.h>

namespace PBD {

class UVFEMTriangleConstraint : public Constraint
{
public:
  static int TYPE_ID;
  Real m_area;
  Matrix2r m_invRestMat;

  UVFEMTriangleConstraint() : Constraint(3) {}
  virtual int &getTypeId() const { return TYPE_ID; }

  virtual bool initConstraint(SimulationModel &model, const Vector2r& u0, const Vector2r& u1, const Vector2r& u2);
  virtual bool solvePositionConstraint(SimulationModel &model, const unsigned int iter);
};

bool solve_UVFEMTriangleConstraint(
  const Vector3r &p0, Real invMass0, 
  const Vector3r &p1, Real invMass1,
  const Vector3r &p2, Real invMass2,
  const Real &area,
  const Matrix2r &invRestMat,
  const Real youngsModulusX,
  const Real youngsModulusY,
  const Real youngsModulusShear,
  const Real poissonRatioXY,
  const Real poissonRatioYX,
  Vector3r &corr0, Vector3r &corr1, Vector3r &corr2);

class UVStrainTriangleConstraint : public Constraint
{
public:
  static int TYPE_ID;
  // Real m_area;
  Matrix2r m_invRestMat;

  UVStrainTriangleConstraint() : Constraint(3) {}
  virtual int &getTypeId() const { return TYPE_ID; }

  virtual bool initConstraint(SimulationModel &model, const Vector2r& u0, const Vector2r& u1, const Vector2r& u2);
  virtual bool solvePositionConstraint(SimulationModel &model, const unsigned int iter);
};

bool solve_UVStrainTriangleConstraint(
  const Vector3r &p0, Real invMass0, 
  const Vector3r &p1, Real invMass1,
  const Vector3r &p2, Real invMass2,
  const Matrix2r &invRestMat,
  const Real xxStiffness, 
  const Real yyStiffness, 
  const Real xyStiffness,
  const bool normalizeStretch,
  const bool normalizeShear,
  Vector3r &corr0, Vector3r &corr1, Vector3r &corr2);
}

#endif // __PBDCONSTRAINTS__H__	
