#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace acorn {
namespace fs3d {

typedef Matvec<3,Real> mv3;

/* Solution for the case
   z = x - y = (x1,0,0)
   nml = v = (v1,0,v3)
*/
template <typename RealT>
void calc_sigma_point_znml (
  const Real lam, const Real mu, const RealT& v1, const RealT& v3,
  const RealT d[3], const RealT& x1, RealT sigma[6])
{
  const RealT
    d1 = d[0], d2 = d[1], d3 = d[2];
  const Real
    lm = lam + mu,
    l2m = lam + 2*mu;
  const RealT
    R3 = std::abs(x1*x1*x1),
    fac1 = 1.0/(M_PI*l2m*mu*R3),
    fac2 = 1.0/(M_PI*l2m*R3);
  const auto T_x_0  = -1.0/2.0*(2*l2m*mu + lam*(l2m - lm))*fac1;
  const auto T_x_8  = (1.0/4.0)*(l2m + 2*lm)*fac2;
  const auto T_x_10 = T_x_8 ;
  const auto T_x_20 = T_x_8 ;
  const auto T_x_24 = (1.0/2.0)*(-lam*(l2m - lm) + lm*mu)*fac1;
  const auto T_x_28 = (1.0/4.0)*(lam*(l2m - lm) + 2*lm*mu)*fac1;
  const auto T_x_36 = (1.0/2.0)*(-l2m + lm)*fac2;
  const auto T_x_44 = (1.0/4.0)*(l2m - lm)*fac2;
  const auto T_x_52 = (1.0/4.0)*(lam*(l2m - lm) - lm*mu)*fac1;
  const auto T_x_56 = T_x_28;
  const auto T_x_60 = T_x_36;
  const auto T_x_70 = T_x_44;
  const auto T_x_72 = T_x_36;
  const auto T_x_80 = (1.0/4.0)*(lam*(l2m - lm) + mu*(2*l2m - 3*lm))*fac1;
  RealT u_x[9];
  u_x[0] = T_x_0 *d1*v1 + T_x_24*d3*v3;
  u_x[1] = T_x_10*d2*v1;
  u_x[2] = T_x_8 *d1*v3 + T_x_20*d3*v1;
  u_x[3] = T_x_36*d2*v1;
  u_x[4] = T_x_28*d1*v1 + T_x_52*d3*v3;
  u_x[5] = T_x_44*d2*v3;
  u_x[6] = T_x_60*d1*v3 + T_x_72*d3*v1;
  u_x[7] = T_x_70*d2*v3;
  u_x[8] = T_x_56*d1*v1 + T_x_80*d3*v3;
  const RealT udiv = u_x[0] + u_x[4] + u_x[8];
  sigma[0] = lam*udiv + mu*(u_x[0] + u_x[0]);
  sigma[1] = mu*(u_x[1] + u_x[3]);
  sigma[2] = mu*(u_x[2] + u_x[6]);
  sigma[3] = lam*udiv + mu*(u_x[4] + u_x[4]);
  sigma[4] = mu*(u_x[5] + u_x[7]);
  sigma[5] = lam*udiv + mu*(u_x[8] + u_x[8]);}

template <typename RealT>
void calc_sigma_point (const Real lam, const Real mu, const RealT y[3],
                       const RealT v[3], const RealT d[3], const RealT x[3],
                       RealT sigma[6]) {
  Real z[3], R[9], dlcl[3];
  mv3::subtract(x, y, z);
  form_rotation_given_x_then_z(z, v, R);
  const auto x1 = mv3::dot(R, z);
  const auto v1 = mv3::dot(R, v);
  const auto v3 = mv3::dot(R+6, v);
  mv3::matvec(R, d, dlcl);
  calc_sigma_point_znml(lam, mu, v1, v3, dlcl, x1, sigma);
  rotate_sym_tensor_3x3_RtAR(R, sigma);
}

template void
calc_sigma_point<Real>(const Real lam, const Real mu, const Real y[3], const Real v[3],
                       const Real d[3], const Real x[3], Real sigma[6]);

#ifdef WOODLAND_ACORN_VECTORIZE
# message "ERROR: WOODLAND_ACORN_VECTORIZE is not supported yet."
error
#endif

} // namespace fs3d
} // namespace acorn
} // namespace woodland
