#include "woodland/acorn/hs3d.hpp"

#include <cmath>

namespace woodland {
namespace acorn {
namespace hs3d {

template <typename RealT>
void calc_halfspace_term_A (
  const Real lam, const Real mu, const RealT v[3], const RealT d[3],
  const RealT& y3, const RealT& x1, const RealT& x3, RealT sigma[6])
{
  const Real
    lm = lam + mu,
    l2m = lam + 2*mu;
  const RealT
    d1 = d[0], d2 = d[1], d3 = d[2],
    v1 = v[0], v2 = v[1], v3 = v[2],
    z1 = x1, z3 = -(x3 + y3),
    z11 = z1*z1, z33 = z3*z3, z13 = z1*z3,
    R2 = z11 + z33,
    R = std::sqrt(R2),
    R3 = R2*R;
  const RealT
    R4 = R2*R2, R5 = R4*R, R7 = R5*R2,
    z333 = z33*z3,
    z1111 = z11*z11, z3333 = z33*z33,
    den1 = M_PI*l2m*R7,
    den2 = M_PI*l2m*mu*R7,
    den3 = M_PI*l2m*R5,
    den4 = M_PI*l2m*mu*R5;
  const auto T_x_0  = (1.0/4.0)*(2*l2m*lam*z1111 + l2m*lam*z11*z33 - l2m*lam*z3333 + 4*l2m*mu*z1111 + 2*l2m*mu*z11*z33 - 2*l2m*mu*z3333 - 2*lam*lm*z1111 - lam*lm*z11*z33 + lam*lm*z3333 - 12*lm*mu*z11*z33 + 3*lm*mu*z3333)/den2;
  const auto T_x_2  = -3.0/8.0*z13*(lam*(5*lm*z11 + 5*lm*z33 - 4*lm*R2 + (2*l2m - 3*lm)*R2) + 2*mu*(5*lm*z11 + (2*l2m - 3*lm)*R2))/den2;
  const auto T_x_4  = (1.0/8.0)*(-6*lm*z11 + lm*R2 - (2*l2m - lm)*R2)/den3;
  const auto T_x_6  = (3.0/8.0)*z13*(10*lm*z11 - 3*lm*R2 + (2*l2m - 3*lm)*R2)/den1;
  const auto T_x_8  = (1.0/8.0)*(-30*lm*z11*z33 + 2*lm*R4 + (2*l2m - lm)*R4 + 3*R2*(lm*z11 - z33*(2*l2m - lm)))/den1;
  const auto T_x_10 = T_x_4 ;
  const auto T_x_12 = (1.0/4.0)*(2*l2m*lam*z11 - l2m*lam*z33 - 2*lam*lm*z11 + lam*lm*z33 - 2*lm*mu*z11 + lm*mu*z33)/den4;
  const auto T_x_14 = (3.0/4.0)*z13*(-l2m*lam + lam*lm + lm*mu)/den4;
  const auto T_x_16 = -3.0/4.0*lm*z13/den3;
  const auto T_x_18 = T_x_6 ;
  const auto T_x_20 = T_x_8 ;
  const auto T_x_22 = T_x_16;
  const auto T_x_24 = (1.0/4.0)*(2*l2m*lam*z1111 + l2m*lam*z11*z33 - l2m*lam*z3333 - 2*lam*lm*z1111 - lam*lm*z11*z33 + lam*lm*z3333 - 2*lm*mu*z1111 + 11*lm*mu*z11*z33 - 2*lm*mu*z3333)/den2;
  const auto T_x_26 = (3.0/8.0)*z13*(-lam*(5*lm*z11 + 5*lm*z33 - 4*lm*R2 + (2*l2m - 3*lm)*R2) + 2*lm*mu*(3*z11 - 2*z33))/den2;
  const auto T_x_28 = (1.0/4.0)*(-lam*(l2m - lm)*R2 + lm*mu*(-2*z11 + z33))/den4;
  const auto T_x_30 = (1.0/8.0)*(-3*lm*z11 + lm*R2 + 3*z11*(2*l2m - lm) + (-2*l2m + lm)*R2)/den3;
  const auto T_x_32 = (3.0/4.0)*z13*(-l2m + lm)/den3;
  const auto T_x_34 = T_x_16;
  const auto T_x_36 = T_x_30;
  const auto T_x_38 = T_x_32;
  const auto T_x_40 = (1.0/4.0)*(-lam*(l2m - lm) + mu*(-2*l2m + 3*lm))/(M_PI*l2m*mu*R3);
  const auto T_x_42 = (3.0/4.0)*z13*(l2m - lm)/den3;
  const auto T_x_44 = (1.0/8.0)*(3*lm*z33 - lm*R2 - 3*z33*(2*l2m - lm) - (-2*l2m + lm)*R2)/den3;
  const auto T_x_46 = T_x_16;
  const auto T_x_48 = T_x_42;
  const auto T_x_50 = T_x_44;
  const auto T_x_52 = (1.0/4.0)*(-lam*(l2m - lm)*R2 + lm*mu*(z11 - 2*z33))/den4;
  const auto T_x_54 = (3.0/8.0)*z1*(lam*(5*lm*z11*z3 + 5*lm*z333 - 4*lm*z3*R2 + R2*(-2*lm*z3 + z3*(2*l2m - lm))) + 2*lm*mu*z3*(2*z11 - 3*z33))/den2;
  const auto T_x_56 = (1.0/8.0)*(-lam*(15*lm*z11*z33 + 15*lm*z3333 + 2*lm*R4 - 3*lm*R2*(z11 + 2*z33) + 3*z3*R2*(-5*lm*z3 + z3*(2*l2m - lm)) + (-2*l2m + 3*lm)*R4) + 2*lm*mu*(-15*z11*z33 + 2*R4))/den2;
  const auto T_x_58 = T_x_16;
  const auto T_x_60 = (1.0/4.0)*(2*l2m*z1111 + l2m*z11*z33 - l2m*z3333 - 2*lm*z1111 + 11*lm*z11*z33 - 2*lm*z3333)/den1;
  const auto T_x_62 = (3.0/8.0)*z1*(-10*lm*z333 + 3*lm*z3*R2 - R2*(-2*lm*z3 + z3*(2*l2m - lm)))/den1;
  const auto T_x_64 = T_x_16;
  const auto T_x_66 = (3.0/4.0)*z1*(l2m*lam*z3 + lam*lm*x3 + lam*lm*y3 - lm*mu*z3)/den4;
  const auto T_x_68 = -1.0/4.0*(-l2m*lam*z11 + 2*l2m*lam*z33 + 6*lam*lm*x3*z3 + 6*lam*lm*y3*z3 + lam*lm*z11 + 4*lam*lm*z33 + lm*mu*z11 - 2*lm*mu*z33)/den4;
  const auto T_x_70 = (1.0/8.0)*(-6*lm*z33 + lm*R2 - (2*l2m - lm)*R2)/den3;
  const auto T_x_72 = T_x_60;
  const auto T_x_74 = T_x_62;
  const auto T_x_76 = T_x_70;
  const auto T_x_78 = (3.0/8.0)*z1*(lam*(5*lm*z11*z3 + 5*lm*z333 - 4*lm*z3*R2 + R2*(-2*lm*z3 + z3*(2*l2m - lm))) + 2*mu*(5*lm*z333 + R2*(-2*lm*z3 + z3*(2*l2m - lm))))/den2;
  const auto T_x_80 = (1.0/8.0)*(-lam*(15*lm*z11*z33 + 15*lm*z3333 + 2*lm*R4 - 3*lm*R2*(z11 + 2*z33) + 3*z3*R2*(-5*lm*z3 + z3*(2*l2m - lm)) + (-2*l2m + 3*lm)*R4) + 2*mu*(-15*lm*z3333 + 3*z3*R2*(5*lm*z3 + z3*(-2*l2m + lm)) + (2*l2m - 3*lm)*R4))/den2;
  RealT u_x[9];
  u_x[0] = T_x_0 *d1*v1 + T_x_6 *d1*v3 + T_x_12*d2*v2 + T_x_18*d3*v1 + T_x_24*d3*v3;
  u_x[1] = T_x_4 *d1*v2 + T_x_10*d2*v1 + T_x_16*d2*v3 + T_x_22*d3*v2;
  u_x[2] = T_x_2 *d1*v1 + T_x_8 *d1*v3 + T_x_14*d2*v2 + T_x_20*d3*v1 + T_x_26*d3*v3;
  u_x[3] = T_x_30*d1*v2 + T_x_36*d2*v1 + T_x_42*d2*v3 + T_x_48*d3*v2;
  u_x[4] = T_x_28*d1*v1 + T_x_34*d1*v3 + T_x_40*d2*v2 + T_x_46*d3*v1 + T_x_52*d3*v3;
  u_x[5] = T_x_32*d1*v2 + T_x_38*d2*v1 + T_x_44*d2*v3 + T_x_50*d3*v2;
  u_x[6] = T_x_54*d1*v1 + T_x_60*d1*v3 + T_x_66*d2*v2 + T_x_72*d3*v1 + T_x_78*d3*v3;
  u_x[7] = T_x_58*d1*v2 + T_x_64*d2*v1 + T_x_70*d2*v3 + T_x_76*d3*v2;
  u_x[8] = T_x_56*d1*v1 + T_x_62*d1*v3 + T_x_68*d2*v2 + T_x_74*d3*v1 + T_x_80*d3*v3;
  const RealT udiv = u_x[0] + u_x[4] + u_x[8];
  sigma[0] = lam*udiv + mu*(u_x[0] + u_x[0]);
  sigma[1] = mu*(u_x[1] + u_x[3]);
  sigma[2] = mu*(u_x[2] + u_x[6]);
  sigma[3] = lam*udiv + mu*(u_x[4] + u_x[4]);
  sigma[4] = mu*(u_x[5] + u_x[7]);
  sigma[5] = lam*udiv + mu*(u_x[8] + u_x[8]);
}

template void calc_halfspace_term_A<Real> (
  const Real lam, const Real mu, const Real v[3], const Real d[3],
  const Real& y3, const Real& x1, const Real& x3, Real sigma[6]);

#ifdef WOODLAND_ACORN_VECTORIZE
# message "ERROR: WOODLAND_ACORN_VECTORIZE is not supported yet."
error
#endif

} // namespace hs3d
} // namespace acorn
} // namespace woodland
