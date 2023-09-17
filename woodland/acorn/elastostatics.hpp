// Derived from https://github.com/ambrad/dc3dm; see
//     https://github.com/ambrad/dc3dm/blob/main/LICENSE
// for details of the Eclipse Public License - v 1.0.

#ifndef INCLUDE_WOODLAND_ACORN_ELASTOSTATICS
#define INCLUDE_WOODLAND_ACORN_ELASTOSTATICS

#include "woodland/acorn/acorn.hpp"
#include "woodland/acorn/dc3d.h"

namespace woodland {
namespace acorn {
namespace okada {

inline int dc3d0(
  // in
  double alpha,
  // Obs point relative to the source point. z <= 0.
  double x, double y, double z,
  // Src fault. depth >= 0.
  double depth, double dipdeg,
  double pot1, double pot2, double pot3, double pot4,
  // out
  double* u,  // [ux uy uz]
  double* du, // [uxx uyx uzx uxy uyy uzy uxz uyz uzz]
  bool want_fullspace = false)
{
#ifdef WOODLAND_ACORN_HAVE_DC3D
  const char space = want_fullspace ? 'f' : 'h';
  int iret;
  dc3d0_(&space, &alpha, &x, &y, &z, &depth, &dipdeg, &pot1, &pot2, &pot3, &pot4,
         u, u+1, u+2, du, du+1, du+2, du+3, du+4, du+5, du+6, du+7, du+8, &iret);
  return iret;
#else
  return -1;
#endif
}

inline int dc3d(
  // in
  double alpha,
  // Obs point relative to the source point. z <= 0.
  double x, double y, double z,
  // Src fault. depth >= 0.
  double depth, double dipdeg,
  double al1, double al2, double aw1, double aw2,
  double disl1, double disl2, double disl3,
  // out
  double* u,  // [ux uy uz]
  double* du, // [uxx uyx uzx uxy uyy uzy uxz uyz uzz]
  bool want_fullspace = false)
{
#ifdef WOODLAND_ACORN_HAVE_DC3D
  const char space = want_fullspace ? 'f' : 'h';
  int iret;
  dc3d_(&space, &alpha, &x, &y, &z,
        &depth, &dipdeg, &al1, &al2, &aw1, &aw2, &disl1, &disl2, &disl3,
        u, u+1, u+2, du, du+1, du+2, du+3, du+4, du+5, du+6, du+7, du+8, &iret);
  return iret;
#else
  return -1;
#endif
}

} // namespace okada

inline Real mu_nu_to_lambda (const Real mu, const Real nu)
{ return 2*mu*nu/(1 - 2*nu); }

inline Real lambda_mu_to_alpha (const Real lambda, const Real mu)
{ return (lambda + mu)/(lambda + 2*mu); }

inline void
du_to_tau (const Real lambda, const Real mu, const Real* du, Real* tau) {
  const Real theta = du[0] + du[4] + du[8];
  tau[0] = lambda*theta + 2.0*mu*du[0];
  tau[1] = mu*(du[1] + du[3]);
  tau[2] = mu*(du[2] + du[6]);
  tau[3] = lambda*theta + 2.0*mu*du[4];
  tau[4] = mu*(du[5] + du[7]);
  tau[5] = lambda*theta + 2.0*mu*du[8];
}

inline void sigma_matvec (
  // Sigma in 6-vector form.
  const Real s[6],
  // y = sigma x.
  const Real x[3], Real y[3])
{
  y[0] = s[0]*x[0] + s[1]*x[1] + s[2]*x[2];
  y[1] = s[1]*x[0] + s[3]*x[1] + s[4]*x[2];
  y[2] = s[2]*x[0] + s[4]*x[1] + s[5]*x[2];
}

inline void sigma_transform (
  // Sigma in 6-vector form.
  const Real sigma[6],
  // Coordinate basis vectors.
  const Real x[3], const Real y[3], const Real z[3],
  // Transformed sigma.
  Real sigma_t[6])
{
  const auto dot = [&] (const Real a[3], const Real b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  };
  Real v[3];
  sigma_matvec(sigma, x, v);
  sigma_t[0] = dot(x, v);
  sigma_t[1] = dot(y, v);
  sigma_t[2] = dot(z, v);
  sigma_matvec(sigma, y, v);
  sigma_t[3] = dot(y, v);
  sigma_t[4] = dot(z, v);
  sigma_matvec(sigma, z, v);
  sigma_t[5] = dot(z, v);
}

void call_okada(
  const bool fullspace, const bool point,
  const Real lam, const Real mu,
  const int n,
  const Real* srcs, const Real* rcvs,
  const Real* dislocs, const Real* nmls,
  const Real* strike_dip_dims,
  Real* sigmas);

} // namespace acorn
} // namespace woodland

#endif
