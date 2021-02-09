#include "PBDConstraints.h"

#include <PositionBasedDynamics/PositionBasedDynamics.h>
#include <Simulation/IDFactory.h>
#include <Simulation/SimulationModel.h>

#include "../utils/debug_logging.h"

using namespace PBD;

const Real eps = static_cast<Real>(1e-6);

/*
int UVFEMTriangleConstraint::TYPE_ID = IDFactory::getId();
bool UVFEMTriangleConstraint::initConstraint(SimulationModel &model, const
Vector2r& u0, const Vector2r& u1, const Vector2r& u2)
{
  // m_bodies[0] = particle1;
  // m_bodies[1] = particle2;
  // m_bodies[2] = particle3;

  // ParticleData &pd = model.getParticles();

  // Vector3r &x1 = pd.getPosition0(particle1);
  // Vector3r &x2 = pd.getPosition0(particle2);
  // Vector3r &x3 = pd.getPosition0(particle3);
  // Vector3r x1,x2,x3;
  // x1 << u0[0], u0[1], 0;
  // x2 << u1[0], u1[1], 0;
  // x3 << u2[0], u2[1], 0;

  // Vector3r &p0 = x1;
  // Vector3r &p1 = x2;
  // Vector3r &p2 = x3;
  // Real& area = m_area;
  // Matrix2r& invRestMat = m_invRestMat;

  // { // adapted from: PositionBasedDynamics::init_FEMTriangleConstraint
  //   // TODO REPLACE WITH AREA AND DETERMINANT FROM UVS SPACE
  //   // probably feed into constructor and init: uv_tri id ...or idk


  //   Vector3r normal0 = (p1 - p0).cross(p2 - p0);
  //   area = normal0.norm() * static_cast<Real>(0.5);

  //   Vector3r axis0_1 = p1 - p0;
  //   axis0_1.normalize();
  //   Vector3r axis0_2 = normal0.cross(axis0_1);
  //   axis0_2.normalize();

  //   Vector2r p[3];
  //   p[0] = Vector2r(p0.dot(axis0_2), p0.dot(axis0_1));
  //   p[1] = Vector2r(p1.dot(axis0_2), p1.dot(axis0_1));
  //   p[2] = Vector2r(p2.dot(axis0_2), p2.dot(axis0_1));

  //   Matrix2r P;
  //   P(0, 0) = p[0][0] - p[2][0];
  //   P(1, 0) = p[0][1] - p[2][1];
  //   P(0, 1) = p[1][0] - p[2][0];
  //   P(1, 1) = p[1][1] - p[2][1];

  //   const Real det = P.determinant();
  //   const Real eps = static_cast<Real>(1e-6);
  //   if (fabs(det) > eps)
  //   {
  //     invRestMat = P.inverse();
  //     return true;
  //   }
  //   return false;
  // }

  Matrix2r D;
  D.col(0) = u1 - u0;
  D.col(1) = u2 - u0;
  // D.col(0) *= 0; // TODO DEBUG why does this resolve crash
  // D.col(0) = u2 - u0; D.col(1) = u1 - u0;
  // D.col(0) = u0 - u2; D.col(1) = u1 - u2;
  // D.col(0) = u1 - u2; D.col(1) = u0 - u2;
  // D.col(0) = u0 - u1; D.col(1) = u2 - u1;
  // D.col(0) = u2 - u1; D.col(1) = u0 - u1;
  const Real det = D.determinant();
  if (fabs(det) > eps)
  {
    m_invRestMat = D.inverse();
    return true;
  }
  return false;
}

bool UVFEMTriangleConstraint::solvePositionConstraint(SimulationModel &model,
const unsigned int iter)
{
  ParticleData &pd = model.getParticles();

  const unsigned i1 = m_bodies[0];
  const unsigned i2 = m_bodies[1];
  const unsigned i3 = m_bodies[2];

  Vector3r &x1 = pd.getPosition(i1);
  Vector3r &x2 = pd.getPosition(i2);
  Vector3r &x3 = pd.getPosition(i3);

  const Real invMass1 = pd.getInvMass(i1);
  const Real invMass2 = pd.getInvMass(i2);
  const Real invMass3 = pd.getInvMass(i3);

  Vector3r corr1, corr2, corr3;
  // const bool res = PositionBasedDynamics::solve_FEMTriangleConstraint(
  const bool res = solve_UVFEMTriangleConstraint(
    x1, invMass1,
    x2, invMass2,
    x3, invMass3,
    m_area,
    m_invRestMat,
    model.getValue<Real>(SimulationModel::CLOTH_STIFFNESS_XX),
    model.getValue<Real>(SimulationModel::CLOTH_STIFFNESS_YY),
    model.getValue<Real>(SimulationModel::CLOTH_STIFFNESS_XY),
    model.getValue<Real>(SimulationModel::CLOTH_POISSON_RATIO_XY),
    model.getValue<Real>(SimulationModel::CLOTH_POISSON_RATIO_YX),
    corr1, corr2, corr3);

  if (res)
  {
    if (invMass1 != 0.0f)
      x1 += corr1;
    if (invMass2 != 0.0f)
      x2 += corr2;
    if (invMass3 != 0.0f)
      x3 += corr3;
  }
  return res;
}

bool PBD::solve_UVFEMTriangleConstraint(
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
        Vector3r &corr0, Vector3r &corr1, Vector3r &corr2)
{
        // Orthotropic elasticity tensor
        Matrix3r C;
        C.setZero();
        C(0, 0) = youngsModulusX / (static_cast<Real>(1.0) -
poissonRatioXY*poissonRatioYX); C(0, 1) = youngsModulusX*poissonRatioYX /
(static_cast<Real>(1.0) - poissonRatioXY*poissonRatioYX); C(1, 1) =
youngsModulusY / (static_cast<Real>(1.0) - poissonRatioXY*poissonRatioYX); C(1,
0) = youngsModulusY*poissonRatioXY / (static_cast<Real>(1.0) -
poissonRatioXY*poissonRatioYX); C(2, 2) = youngsModulusShear;

        // Determine \partial x/\partial m_i
        Eigen::Matrix<Real, 3, 2> F;
  F.col(0) = p1 - p0;
  F.col(1) = p2 - p0;
  F *= invRestMat;

        // const Vector3r p13 = p0 - p2;
        // const Vector3r p23 = p1 - p2;
        // F(0,0) = p13[0] * invRestMat(0,0) + p23[0] * invRestMat(1,0);
        // F(0,1) = p13[0] * invRestMat(0,1) + p23[0] * invRestMat(1,1);
        // F(1,0) = p13[1] * invRestMat(0,0) + p23[1] * invRestMat(1,0);
        // F(1,1) = p13[1] * invRestMat(0,1) + p23[1] * invRestMat(1,1);
        // F(2,0) = p13[2] * invRestMat(0,0) + p23[2] * invRestMat(1,0);
        // F(2,1) = p13[2] * invRestMat(0,1) + p23[2] * invRestMat(1,1);

        // epsilon = 0.5(F^T * F - I)
        Matrix2r epsilon;
        epsilon(0,0) = static_cast<Real>(0.5)*(F(0,0) * F(0,0) + F(1,0) * F(1,0)
+ F(2,0) * F(2,0) - static_cast<Real>(1.0));		// xx epsilon(1,1) =
static_cast<Real>(0.5)*(F(0,1) * F(0,1) + F(1,1) * F(1,1) + F(2,1) * F(2,1) -
static_cast<Real>(1.0));		// yy epsilon(0,1) =
static_cast<Real>(0.5)*(F(0,0) * F(0,1) + F(1,0) * F(1,1) + F(2,0) * F(2,1));
// xy epsilon(1,0) = epsilon(0,1);

        // P(F) = det(F) * C*E * F^-T => E = green strain
        Matrix2r stress;
        stress(0,0) = C(0,0) * epsilon(0,0) + C(0,1) * epsilon(1,1) + C(0,2) *
epsilon(0,1); stress(1,1) = C(1,0) * epsilon(0,0) + C(1,1) * epsilon(1,1) +
C(1,2) * epsilon(0,1); stress(0,1) = C(2,0) * epsilon(0,0) + C(2,1) *
epsilon(1,1) + C(2,2) * epsilon(0,1); stress(1,0) = stress(0,1);

        const Eigen::Matrix<Real, 3, 2> piolaKirchhoffStres = F * stress;

        Real psi = 0.0;
        for (unsigned char j = 0; j < 2; j++)
                for (unsigned char k = 0; k < 2; k++)
                        psi += epsilon(j,k) * stress(j,k);
        psi = static_cast<Real>(0.5)*psi;
        Real energy = area*psi;

        // compute gradient
        Eigen::Matrix<Real, 3, 2> H = area * piolaKirchhoffStres *
invRestMat.transpose(); Vector3r gradC[3]; for (unsigned char j = 0; j < 3; ++j)
        {
                // gradC[0][j] = H(j,0);
                // gradC[1][j] = H(j,1);
                gradC[1][j] = H(j,0);
                gradC[2][j] = H(j,1);
        }
        // gradC[2] = -gradC[0] - gradC[1];
        gradC[0] = -gradC[1] - gradC[2];

  // gradC[0] = Vector3r::Zero();
  // gradC[1] = Vector3r::Zero();
  // gradC[2] = Vector3r::Zero();


        Real sum_normGradC = invMass0 * gradC[0].squaredNorm();
        sum_normGradC += invMass1 * gradC[1].squaredNorm();
        sum_normGradC += invMass2 * gradC[2].squaredNorm();

        // exit early if required
        if (fabs(sum_normGradC) > eps)
        {
                // compute scaling factor
                const Real s = energy / sum_normGradC;

                // update positions
                corr0 = -(s*invMass0) * gradC[0];
                corr1 = -(s*invMass1) * gradC[1];
                corr2 = -(s*invMass2) * gradC[2];


    // static bool dot = true;
    // if(dot && invMass0 > 0.0001){
    //   Debug::log("THING",invMass0, corr0.transpose(), gradC[0].transpose());
    //   dot = false;
    // }

                return true;
        }

        return false;
}
*/

int UVStrainTriangleConstraint::TYPE_ID = IDFactory::getId();
bool UVStrainTriangleConstraint::initConstraint(SimulationModel &model,
                                                const Vector2r &u0,
                                                const Vector2r &u1,
                                                const Vector2r &u2) {
  // Matrix2r D;
  // D.col(0) = u1 - u0;
  // D.col(1) = u2 - u0;

  // const Real det = D.determinant();
  // if (fabs(det) > eps)
  // {
  //   m_invRestMat = D.inverse();
  //   return true;
  // }
  // return false;

  Real a = u1[0] - u0[0];
  Real b = u2[0] - u0[0];
  Real c = u1[1] - u0[1];
  Real d = u2[1] - u0[1];

  // inverse
  Real det = a * d - b * c;
  if (fabs(det) < eps)
    return false;

  Real s             = static_cast<Real>(1.0) / det;
  m_invRestMat(0, 0) = d * s;
  m_invRestMat(0, 1) = -b * s;
  m_invRestMat(1, 0) = -c * s;
  m_invRestMat(1, 1) = a * s;

  return true;
}

bool UVStrainTriangleConstraint::solvePositionConstraint(
    SimulationModel &model, const unsigned int iter) {
  ParticleData &pd = model.getParticles();

  const unsigned i1 = m_bodies[0];
  const unsigned i2 = m_bodies[1];
  const unsigned i3 = m_bodies[2];

  Vector3r &x1 = pd.getPosition(i1);
  Vector3r &x2 = pd.getPosition(i2);
  Vector3r &x3 = pd.getPosition(i3);

  const Real invMass1 = pd.getInvMass(i1);
  const Real invMass2 = pd.getInvMass(i2);
  const Real invMass3 = pd.getInvMass(i3);

  Vector3r corr1, corr2, corr3;
  // const bool res = PositionBasedDynamics::solve_FEMTriangleConstraint(
  const bool res = solve_UVStrainTriangleConstraint(
      x1, invMass1, x2, invMass2, x3, invMass3, m_invRestMat,
      model.getValue<Real>(SimulationModel::CLOTH_STIFFNESS_XX),
      model.getValue<Real>(SimulationModel::CLOTH_STIFFNESS_YY),
      model.getValue<Real>(SimulationModel::CLOTH_STIFFNESS_XY),
      // model.getValue<bool>(SimulationModel::CLOTH_NORMALIZE_STRETCH),
      // model.getValue<bool>(SimulationModel::CLOTH_NORMALIZE_SHEAR),
      false, false, corr1, corr2, corr3);

  if (res) {
    if (invMass1 != 0.0f)
      x1 += corr1;  //*0.001*invMass1;
    if (invMass2 != 0.0f)
      x2 += corr2;  //*0.001*invMass2;
    if (invMass3 != 0.0f)
      x3 += corr3;  //*0.001*invMass3;

    // if(invMass1 > 0)
    // 	std::cout<<"   "<< invMass1<<"\n";
  }
  return res;
}

bool PBD::solve_UVStrainTriangleConstraint(
    const Vector3r &p0, Real invMass0, const Vector3r &p1, Real invMass1,
    const Vector3r &p2, Real invMass2, const Matrix2r &invRestMat,
    const Real xxStiffness, const Real yyStiffness, const Real xyStiffness,
    const bool normalizeStretch, const bool normalizeShear, Vector3r &corr0,
    Vector3r &corr1, Vector3r &corr2) {
  Vector2r c[2];
  c[0] = Vector2r(invRestMat(0, 0), invRestMat(1, 0));
  c[1] = Vector2r(invRestMat(0, 1), invRestMat(1, 1));

  Vector2r r[3];

  corr0.setZero();
  corr1.setZero();
  corr2.setZero();

  // constraint seems to be:
  // C(x) = 1/k ||S - I||^2 or similar? S = (X12+dX12 - (X0+dX0)) D^-1

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j <= i; j++) {
      // 			r[0] = Vector3r(p1[0] - p0[0], p2[0] - p0[0],
      // 0.0);  // Jacobi 			r[1] = Vector3r(p1[1] - p0[1], p2[1] - p0[1], 0.0);
      // 			r[2] = Vector3r(p1[2] - p0[2], p2[2] - p0[2],
      // 0.0);

      r[0] =
          Vector2r((p1[0] + corr1[0]) - (p0[0] + corr0[0]),
                   (p2[0] + corr2[0]) - (p0[0] + corr0[0]));  // Gauss - Seidel
      r[1] = Vector2r((p1[1] + corr1[1]) - (p0[1] + corr0[1]),
                      (p2[1] + corr2[1]) - (p0[1] + corr0[1]));
      r[2] = Vector2r((p1[2] + corr1[2]) - (p0[2] + corr0[2]),
                      (p2[2] + corr2[2]) - (p0[2] + corr0[2]));

      Real Sij = 0.0;
      for (int k = 0; k < 3; k++) Sij += r[k].dot(c[i]) * r[k].dot(c[j]);

      Vector3r d[3];  // grad_xi C
      d[0] = Vector3r(0.0, 0.0, 0.0);

      for (int k = 0; k < 2; k++) {
        d[k + 1] = Vector3r(r[0].dot(c[j]), r[1].dot(c[j]), r[2].dot(c[j])) *
                   invRestMat(k, i);
        d[k + 1] += Vector3r(r[0].dot(c[i]), r[1].dot(c[i]), r[2].dot(c[i])) *
                    invRestMat(k, j);
        d[0] -= d[k + 1];
      }

      if (i != j && normalizeShear) {
        Real fi2 = 0.0;
        Real fj2 = 0.0;
        for (int k = 0; k < 3; k++) {
          fi2 += r[k].dot(c[i]) * r[k].dot(c[i]);
          fj2 += r[k].dot(c[j]) * r[k].dot(c[j]);
        }
        Real fi = sqrt(fi2);
        Real fj = sqrt(fj2);

        d[0]   = Vector3r(0.0, 0.0, 0.0);
        Real s = Sij / (fi2 * fi * fj2 * fj);
        for (int k = 0; k < 2; k++) {
          d[k + 1] /= fi * fj;
          d[k + 1] -= fj * fj *
                      Vector3r(r[0].dot(c[i]), r[1].dot(c[i]), r[2].dot(c[i])) *
                      invRestMat(k, i) * s;
          d[k + 1] -= fi * fi *
                      Vector3r(r[0].dot(c[j]), r[1].dot(c[j]), r[2].dot(c[j])) *
                      invRestMat(k, j) * s;
          d[0] -= d[k + 1];
        }
        Sij = Sij / (fi * fj);
      }

      Real lambda = invMass0 * d[0].squaredNorm() +
                    invMass1 * d[1].squaredNorm() +
                    invMass2 * d[2].squaredNorm();

      if (lambda == 0.0f)
        continue;

      if (i == 0 && j == 0) {
        if (normalizeStretch) {
          Real s = sqrt(Sij);
          lambda = static_cast<Real>(2.0) * s * (s - static_cast<Real>(1.0)) /
                   lambda * xxStiffness;
        } else {
          lambda = (Sij - static_cast<Real>(1.0)) / lambda * xxStiffness;
        }
      } else if (i == 1 && j == 1) {
        if (normalizeStretch) {
          Real s = sqrt(Sij);
          lambda = static_cast<Real>(2.0) * s * (s - static_cast<Real>(1.0)) /
                   lambda * yyStiffness;
        } else {
          lambda = (Sij - static_cast<Real>(1.0)) / lambda * yyStiffness;
        }
      } else {
        lambda = Sij / lambda * xyStiffness;
      }

      corr0 -= lambda * invMass0 * d[0];
      corr1 -= lambda * invMass1 * d[1];
      corr2 -= lambda * invMass2 * d[2];

      // if (i+j<1 && invMass0 > 0.001f){
      // 	std::cout<<"  "<<corr0[0]<<"  "<<corr0[1]<<"  "<<corr0[2]<<"
      // "<<invMass0<<"\n";
      // }
    }
  }
  return true;
}

/*
int UVIsometricBendingConstraint::TYPE_ID = IDFactory::getId();
bool UVIsometricBendingConstraint::initConstraint(SimulationModel &model, const
unsigned int particle1, const unsigned int particle2, const unsigned int
particle3, const unsigned int particle4)
{
        m_bodies[0] = particle1;
        m_bodies[1] = particle2;
        m_bodies[2] = particle3;
        m_bodies[3] = particle4;

        ParticleData &pd = model.getParticles();

        const Vector3r &x1 = pd.getPosition0(particle1);
        const Vector3r &x2 = pd.getPosition0(particle2);
        const Vector3r &x3 = pd.getPosition0(particle3);
        const Vector3r &x4 = pd.getPosition0(particle4);

        return PositionBasedDynamics::init_IsometricBendingConstraint(x1, x2,
x3, x4, m_Q);
}

bool UVIsometricBendingConstraint::solvePositionConstraint(SimulationModel
&model, const unsigned int iter)
{
        ParticleData &pd = model.getParticles();

        const unsigned i1 = m_bodies[0];
        const unsigned i2 = m_bodies[1];
        const unsigned i3 = m_bodies[2];
        const unsigned i4 = m_bodies[3];

        Vector3r &x1 = pd.getPosition(i1);
        Vector3r &x2 = pd.getPosition(i2);
        Vector3r &x3 = pd.getPosition(i3);
        Vector3r &x4 = pd.getPosition(i4);

        const Real invMass1 = pd.getInvMass(i1);
        const Real invMass2 = pd.getInvMass(i2);
        const Real invMass3 = pd.getInvMass(i3);
        const Real invMass4 = pd.getInvMass(i4);

        Vector3r corr1, corr2, corr3, corr4;
        const bool res =
PositionBasedDynamics::solve_IsometricBendingConstraint( x1, invMass1, x2,
invMass2, x3, invMass3, x4, invMass4, m_Q,
                model.getValue<Real>(SimulationModel::CLOTH_BENDING_STIFFNESS),
                corr1, corr2, corr3, corr4);

        if (res)
        {
                if (invMass1 != 0.0)
                        x1 += corr1;
                if (invMass2 != 0.0)
                        x2 += corr2;
                if (invMass3 != 0.0)
                        x3 += corr3;
                if (invMass4 != 0.0)
                        x4 += corr4;
        }
        return res;
}
*/