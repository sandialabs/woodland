#include "woodland/acorn/hs3d.hpp"

#include <cmath>

namespace woodland {
namespace acorn {
namespace hs3d {

template <typename RealT>
void calc_halfspace_term_C (
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
    z111 = z11*z1,
    R2 = z11 + z33,
    R = std::sqrt(R2),
    R3 = R2*R;
  const RealT
    R4 = R2*R2, R5 = R4*R, R7 = R5*R2, R9 = R7*R2,
    z1111 = z11*z11,
    z333 = z33*z3,
    z3333 = z33*z33,
    c2 = M_PI*l2m,
    c1 = c2*mu;
  const auto T_x_0  = lam*((1.0/4.0)*x3*(-15*z11*z3*(2 - lm/l2m)/R7 + 3*z3*(2 - lm/l2m)/R5 + 105*lm*y3*z11*z33/(l2m*R9) - 15*lm*y3*z11/(l2m*R7) - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) + 15*lm*z11*z3/(l2m*R7) - 3*lm*z3/(l2m*R5))/(M_PI*mu) - 15.0/4.0*lm*x3*y3*z11/(c1*R7) + (1.0/4.0)*lm*x3*y3*(105*z1111/R9 - 90*z11/R7 + 9/R5)/(c1) + (3.0/4.0)*lm*x3*y3/(c1*R5)) + (1.0/2.0)*lm*x3*y3*(105*z1111/R9 - 90*z11/R7 + 9/R5)/(c2);
  const auto T_x_2  = lam*((1.0/4.0)*x3*(15*z1*z33*(2 - lm/l2m)/R7 - 3*z1*(2 - lm/l2m)/R5 - 105*lm*y3*z1*z333/(l2m*R9) + 45*lm*y3*z13/(l2m*R7) - 15*lm*z1*z33/(l2m*R7) + 3*lm*z1/(l2m*R5))/(M_PI*mu) + (1.0/4.0)*(3*z13*(2 - lm/l2m)/R5 - 15*lm*y3*z1*z33/(l2m*R7) + 3*lm*y3*z1/(l2m*R5) - 3*lm*z13/(l2m*R5))/(M_PI*mu) + (15.0/4.0)*lm*x3*y3*z13/(c1*R7) + (1.0/4.0)*lm*x3*y3*(-105*z111*z3/R9 + 45*z13/R7)/(c1) + (3.0/4.0)*lm*y3*z1/(c1*R5) + (1.0/4.0)*lm*y3*(-15*z111/R7 + 9*z1/R5)/(c1)) + (1.0/2.0)*lm*x3*y3*(-105*z111*z3/R9 + 45*z13/R7)/(c2) + (1.0/2.0)*lm*y3*(-15*z111/R7 + 9*z1/R5)/(c2);
  const auto T_x_4  = mu*(-15.0/4.0*lm*x3*y3*z11/(c1*R7) + (1.0/4.0)*lm*x3*y3*(-15*z11/R7 + 3/R5)/(c1) + (3.0/4.0)*lm*x3*y3/(c1*R5));
  const auto T_x_6  = mu*((1.0/4.0)*x3*(-15*z111*(2 - lm/l2m)/R7 + 9*z1*(2 - lm/l2m)/R5 + 105*lm*y3*z111*z3/(l2m*R9) - 45*lm*y3*z13/(l2m*R7))/(M_PI*mu) + (1.0/4.0)*lm*x3*y3*(105*z111*z3/R9 - 45*z13/R7)/(c1) + (1.0/4.0)*lm*x3*(15*z111/R7 - 9*z1/R5)/(c1));
  const auto T_x_8  = mu*((1.0/4.0)*x3*(15*z11*z3*(2 - lm/l2m)/R7 - 3*z3*(2 - lm/l2m)/R5 - 105*lm*y3*z11*z33/(l2m*R9) + 15*lm*y3*z11/(l2m*R7) + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/(M_PI*mu) + (1.0/4.0)*(3*z11*(2 - lm/l2m)/R5 - (2 - lm/l2m)/R3 - 15*lm*y3*z11*z3/(l2m*R7) + 3*lm*y3*z3/(l2m*R5))/(M_PI*mu) + (1.0/4.0)*lm*x3*y3*(-105*z11*z33/R9 + 15*z11/R7 + 15*z33/R7 - 3/R5)/(c1) + (1.0/4.0)*lm*x3*(-15*z11*z3/R7 + 3*z3/R5)/(c1) + (1.0/4.0)*lm*y3*(-15*z11*z3/R7 + 3*z3/R5)/(c1) + (1.0/4.0)*lm*(-3*z11/R5 + (1.0/R3))/(c1));
  const auto T_x_10 = T_x_4 ;
  const auto T_x_12 = lam*((1.0/4.0)*x3*(-15*z11*z3*(2 - lm/l2m)/R7 + 3*z3*(2 - lm/l2m)/R5 + 105*lm*y3*z11*z33/(l2m*R9) - 15*lm*y3*z11/(l2m*R7) - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) + 15*lm*z11*z3/(l2m*R7) - 3*lm*z3/(l2m*R5))/(M_PI*mu) - 15.0/4.0*lm*x3*y3*z11/(c1*R7) + (1.0/4.0)*lm*x3*y3*(105*z1111/R9 - 90*z11/R7 + 9/R5)/(c1) + (3.0/4.0)*lm*x3*y3/(c1*R5)) + mu*(-15.0/2.0*lm*x3*y3*z11/(c1*R7) + (3.0/2.0)*lm*x3*y3/(c1*R5));
  const auto T_x_14 = lam*((1.0/4.0)*x3*(15*z1*z33*(2 - lm/l2m)/R7 - 3*z1*(2 - lm/l2m)/R5 - 105*lm*y3*z1*z333/(l2m*R9) + 45*lm*y3*z13/(l2m*R7) - 15*lm*z1*z33/(l2m*R7) + 3*lm*z1/(l2m*R5))/(M_PI*mu) + (1.0/4.0)*(3*z13*(2 - lm/l2m)/R5 - 15*lm*y3*z1*z33/(l2m*R7) + 3*lm*y3*z1/(l2m*R5) - 3*lm*z13/(l2m*R5))/(M_PI*mu) + (15.0/4.0)*lm*x3*y3*z13/(c1*R7) + (1.0/4.0)*lm*x3*y3*(-105*z111*z3/R9 + 45*z13/R7)/(c1) + (3.0/4.0)*lm*y3*z1/(c1*R5) + (1.0/4.0)*lm*y3*(-15*z111/R7 + 9*z1/R5)/(c1)) + mu*((15.0/2.0)*lm*x3*y3*z13/(c1*R7) + (3.0/2.0)*lm*y3*z1/(c1*R5));
  const auto T_x_16 = mu*((1.0/4.0)*x3*(3*z1*(2 - lm/l2m)/R5 - 15*lm*y3*z13/(l2m*R7))/(M_PI*mu) - 15.0/4.0*lm*x3*y3*z13/(c1*R7) - 3.0/4.0*lm*x3*z1/(c1*R5));
  const auto T_x_18 = T_x_6 ;
  const auto T_x_20 = T_x_8 ;
  const auto T_x_22 = T_x_16;
  const auto T_x_24 = lam*((1.0/4.0)*x3*(-15*z11*z3*(2 - lm/l2m)/R7 + 3*z3*(2 - lm/l2m)/R5 + 105*lm*y3*z11*z33/(l2m*R9) - 15*lm*y3*z11/(l2m*R7) - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) + 15*lm*z11*z3/(l2m*R7) - 3*lm*z3/(l2m*R5))/(M_PI*mu) - 15.0/4.0*lm*x3*y3*z11/(c1*R7) + (1.0/4.0)*lm*x3*y3*(105*z1111/R9 - 90*z11/R7 + 9/R5)/(c1) + (3.0/4.0)*lm*x3*y3/(c1*R5)) + (1.0/2.0)*x3*(-15*z11*z3*(2 - lm/l2m)/R7 + 3*z3*(2 - lm/l2m)/R5 + 105*lm*y3*z11*z33/(l2m*R9) - 15*lm*y3*z11/(l2m*R7) - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) + 15*lm*z11*z3/(l2m*R7) - 3*lm*z3/(l2m*R5))/M_PI;
  const auto T_x_26 = lam*((1.0/4.0)*x3*(15*z1*z33*(2 - lm/l2m)/R7 - 3*z1*(2 - lm/l2m)/R5 - 105*lm*y3*z1*z333/(l2m*R9) + 45*lm*y3*z13/(l2m*R7) - 15*lm*z1*z33/(l2m*R7) + 3*lm*z1/(l2m*R5))/(M_PI*mu) + (1.0/4.0)*(3*z13*(2 - lm/l2m)/R5 - 15*lm*y3*z1*z33/(l2m*R7) + 3*lm*y3*z1/(l2m*R5) - 3*lm*z13/(l2m*R5))/(M_PI*mu) + (15.0/4.0)*lm*x3*y3*z13/(c1*R7) + (1.0/4.0)*lm*x3*y3*(-105*z111*z3/R9 + 45*z13/R7)/(c1) + (3.0/4.0)*lm*y3*z1/(c1*R5) + (1.0/4.0)*lm*y3*(-15*z111/R7 + 9*z1/R5)/(c1)) + (1.0/2.0)*x3*(15*z1*z33*(2 - lm/l2m)/R7 - 3*z1*(2 - lm/l2m)/R5 - 105*lm*y3*z1*z333/(l2m*R9) + 45*lm*y3*z13/(l2m*R7) - 15*lm*z1*z33/(l2m*R7) + 3*lm*z1/(l2m*R5))/M_PI + (1.0/2.0)*(3*z13*(2 - lm/l2m)/R5 - 15*lm*y3*z1*z33/(l2m*R7) + 3*lm*y3*z1/(l2m*R5) - 3*lm*z13/(l2m*R5))/M_PI;
  const auto T_x_28 = lam*((1.0/4.0)*x3*(3*z3*(2 - lm/l2m)/R5 - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) - 3*lm*z3/(l2m*R5))/(M_PI*mu) - 15.0/4.0*lm*x3*y3*z11/(c1*R7) + 3*lm*x3*y3/(c1*R5)) + mu*(-15.0/2.0*lm*x3*y3*z11/(c1*R7) + (3.0/2.0)*lm*x3*y3/(c1*R5));
  const auto T_x_30 = T_x_4 ;
  const auto T_x_32 = mu*((15.0/2.0)*lm*x3*y3*z13/(c1*R7) + (3.0/2.0)*lm*y3*z1/(c1*R5));
  const auto T_x_34 = T_x_16;
  const auto T_x_36 = T_x_4 ;
  const auto T_x_38 = T_x_32;
  const auto T_x_40 = lam*((1.0/4.0)*x3*(3*z3*(2 - lm/l2m)/R5 - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) - 3*lm*z3/(l2m*R5))/(M_PI*mu) - 15.0/4.0*lm*x3*y3*z11/(c1*R7) + 3*lm*x3*y3/(c1*R5)) + (9.0/2.0)*lm*x3*y3/(c2*R5);
  const auto T_x_42 = T_x_16;
  const auto T_x_44 = mu*((1.0/4.0)*x3*(-3*z3*(2 - lm/l2m)/R5 + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/(M_PI*mu) + (1.0/4.0)*(-(2 - lm/l2m)/R3 + 3*lm*y3*z3/(l2m*R5))/(M_PI*mu) + (1.0/4.0)*lm*x3*y3*(15*z33/R7 - 3/R5)/(c1) + (3.0/4.0)*lm*x3*z3/(c1*R5) + (3.0/4.0)*lm*y3*z3/(c1*R5) + (1.0/4.0)*lm/(c1*R3));
  const auto T_x_46 = T_x_16;
  const auto T_x_48 = T_x_16;
  const auto T_x_50 = T_x_44;
  const auto T_x_52 = lam*((1.0/4.0)*x3*(3*z3*(2 - lm/l2m)/R5 - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) - 3*lm*z3/(l2m*R5))/(M_PI*mu) - 15.0/4.0*lm*x3*y3*z11/(c1*R7) + 3*lm*x3*y3/(c1*R5)) + (1.0/2.0)*x3*(3*z3*(2 - lm/l2m)/R5 - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) - 3*lm*z3/(l2m*R5))/M_PI;
  const auto T_x_54 = lam*(-1.0/4.0*x3*(-3*z1*(2 - lm/l2m)/R5 - 15*lm*y3*z13/(l2m*R7))/(M_PI*mu) - 1.0/4.0*x3*(15*z111*(2 - lm/l2m)/R7 - 9*z1*(2 - lm/l2m)/R5 + 105*lm*y3*z111*z3/(l2m*R9) - 45*lm*y3*z13/(l2m*R7))/(M_PI*mu) - 1.0/4.0*lm*x3*y3*(105*z1*z333/R9 - 45*z13/R7)/(c1) - 1.0/4.0*lm*x3*(15*z1*z33/R7 - 3*z1/R5)/(c1)) - 1.0/2.0*x3*(15*z111*(2 - lm/l2m)/R7 - 9*z1*(2 - lm/l2m)/R5 + 105*lm*y3*z111*z3/(l2m*R9) - 45*lm*y3*z13/(l2m*R7))/M_PI;
  const auto T_x_56 = lam*(-1.0/4.0*x3*(3*z3*(2 - lm/l2m)/R5 + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*x3*(-15*z11*z3*(2 - lm/l2m)/R7 + 3*z3*(2 - lm/l2m)/R5 - 105*lm*y3*z11*z33/(l2m*R9) + 15*lm*y3*z11/(l2m*R7) + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*((2 - lm/l2m)/R3 + 3*lm*y3*z3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*(-3*z11*(2 - lm/l2m)/R5 + (2 - lm/l2m)/R3 - 15*lm*y3*z11*z3/(l2m*R7) + 3*lm*y3*z3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*lm*x3*y3*(-105*z3333/R9 + 90*z33/R7 - 9/R5)/(c1) - 1.0/4.0*lm*x3*(-15*z333/R7 + 9*z3/R5)/(c1) - 1.0/4.0*lm*y3*(-15*z333/R7 + 9*z3/R5)/(c1) - 1.0/4.0*lm*(-3*z33/R5 + (1.0/R3))/(c1)) - 1.0/2.0*x3*(-15*z11*z3*(2 - lm/l2m)/R7 + 3*z3*(2 - lm/l2m)/R5 - 105*lm*y3*z11*z33/(l2m*R9) + 15*lm*y3*z11/(l2m*R7) + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/M_PI - 1.0/2.0*(-3*z11*(2 - lm/l2m)/R5 + (2 - lm/l2m)/R3 - 15*lm*y3*z11*z3/(l2m*R7) + 3*lm*y3*z3/(l2m*R5))/M_PI;
  const auto T_x_58 = -1.0/2.0*x3*(-3*z1*(2 - lm/l2m)/R5 - 15*lm*y3*z13/(l2m*R7))/M_PI;
  const auto T_x_60 = mu*(-1.0/4.0*x3*(15*z11*z3*(2 - lm/l2m)/R7 - 3*z3*(2 - lm/l2m)/R5 + 105*lm*y3*z11*z33/(l2m*R9) - 15*lm*y3*z11/(l2m*R7) - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) + 15*lm*z11*z3/(l2m*R7) - 3*lm*z3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*lm*x3*y3*(105*z11*z33/R9 - 15*z11/R7 - 15*z33/R7 + 3/R5)/(c1));
  const auto T_x_62 = mu*(-1.0/4.0*x3*(-15*z1*z33*(2 - lm/l2m)/R7 + 3*z1*(2 - lm/l2m)/R5 - 105*lm*y3*z1*z333/(l2m*R9) + 45*lm*y3*z13/(l2m*R7) - 15*lm*z1*z33/(l2m*R7) + 3*lm*z1/(l2m*R5))/(M_PI*mu) - 1.0/4.0*(-3*z13*(2 - lm/l2m)/R5 - 15*lm*y3*z1*z33/(l2m*R7) + 3*lm*y3*z1/(l2m*R5) - 3*lm*z13/(l2m*R5))/(M_PI*mu) - 1.0/4.0*lm*x3*y3*(-105*z1*z333/R9 + 45*z13/R7)/(c1) - 1.0/4.0*lm*y3*(-15*z1*z33/R7 + 3*z1/R5)/(c1));
  const auto T_x_64 = T_x_58;
  const auto T_x_66 = lam*(-1.0/4.0*x3*(-3*z1*(2 - lm/l2m)/R5 - 15*lm*y3*z13/(l2m*R7))/(M_PI*mu) - 1.0/4.0*x3*(15*z111*(2 - lm/l2m)/R7 - 9*z1*(2 - lm/l2m)/R5 + 105*lm*y3*z111*z3/(l2m*R9) - 45*lm*y3*z13/(l2m*R7))/(M_PI*mu) - 1.0/4.0*lm*x3*y3*(105*z1*z333/R9 - 45*z13/R7)/(c1) - 1.0/4.0*lm*x3*(15*z1*z33/R7 - 3*z1/R5)/(c1)) - 1.0/2.0*x3*(-3*z1*(2 - lm/l2m)/R5 - 15*lm*y3*z13/(l2m*R7))/M_PI;
  const auto T_x_68 = lam*(-1.0/4.0*x3*(3*z3*(2 - lm/l2m)/R5 + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*x3*(-15*z11*z3*(2 - lm/l2m)/R7 + 3*z3*(2 - lm/l2m)/R5 - 105*lm*y3*z11*z33/(l2m*R9) + 15*lm*y3*z11/(l2m*R7) + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*((2 - lm/l2m)/R3 + 3*lm*y3*z3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*(-3*z11*(2 - lm/l2m)/R5 + (2 - lm/l2m)/R3 - 15*lm*y3*z11*z3/(l2m*R7) + 3*lm*y3*z3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*lm*x3*y3*(-105*z3333/R9 + 90*z33/R7 - 9/R5)/(c1) - 1.0/4.0*lm*x3*(-15*z333/R7 + 9*z3/R5)/(c1) - 1.0/4.0*lm*y3*(-15*z333/R7 + 9*z3/R5)/(c1) - 1.0/4.0*lm*(-3*z33/R5 + (1.0/R3))/(c1)) - 1.0/2.0*x3*(3*z3*(2 - lm/l2m)/R5 + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/M_PI - 1.0/2.0*((2 - lm/l2m)/R3 + 3*lm*y3*z3/(l2m*R5))/M_PI;
  const auto T_x_70 = mu*(-1.0/4.0*x3*(-3*z3*(2 - lm/l2m)/R5 - 15*lm*y3*z33/(l2m*R7) + 3*lm*y3/(l2m*R5) - 3*lm*z3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*lm*x3*y3*(-15*z33/R7 + 3/R5)/(c1));
  const auto T_x_72 = T_x_60;
  const auto T_x_74 = T_x_62;
  const auto T_x_76 = T_x_70;
  const auto T_x_78 = lam*(-1.0/4.0*x3*(-3*z1*(2 - lm/l2m)/R5 - 15*lm*y3*z13/(l2m*R7))/(M_PI*mu) - 1.0/4.0*x3*(15*z111*(2 - lm/l2m)/R7 - 9*z1*(2 - lm/l2m)/R5 + 105*lm*y3*z111*z3/(l2m*R9) - 45*lm*y3*z13/(l2m*R7))/(M_PI*mu) - 1.0/4.0*lm*x3*y3*(105*z1*z333/R9 - 45*z13/R7)/(c1) - 1.0/4.0*lm*x3*(15*z1*z33/R7 - 3*z1/R5)/(c1)) + mu*(-1.0/2.0*lm*x3*y3*(105*z1*z333/R9 - 45*z13/R7)/(c1) - 1.0/2.0*lm*x3*(15*z1*z33/R7 - 3*z1/R5)/(c1));
  const auto T_x_80 = lam*(-1.0/4.0*x3*(3*z3*(2 - lm/l2m)/R5 + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*x3*(-15*z11*z3*(2 - lm/l2m)/R7 + 3*z3*(2 - lm/l2m)/R5 - 105*lm*y3*z11*z33/(l2m*R9) + 15*lm*y3*z11/(l2m*R7) + 15*lm*y3*z33/(l2m*R7) - 3*lm*y3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*((2 - lm/l2m)/R3 + 3*lm*y3*z3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*(-3*z11*(2 - lm/l2m)/R5 + (2 - lm/l2m)/R3 - 15*lm*y3*z11*z3/(l2m*R7) + 3*lm*y3*z3/(l2m*R5))/(M_PI*mu) - 1.0/4.0*lm*x3*y3*(-105*z3333/R9 + 90*z33/R7 - 9/R5)/(c1) - 1.0/4.0*lm*x3*(-15*z333/R7 + 9*z3/R5)/(c1) - 1.0/4.0*lm*y3*(-15*z333/R7 + 9*z3/R5)/(c1) - 1.0/4.0*lm*(-3*z33/R5 + (1.0/R3))/(c1)) + mu*(-1.0/2.0*lm*x3*y3*(-105*z3333/R9 + 90*z33/R7 - 9/R5)/(c1) - 1.0/2.0*lm*x3*(-15*z333/R7 + 9*z3/R5)/(c1) - 1.0/2.0*lm*y3*(-15*z333/R7 + 9*z3/R5)/(c1) - 1.0/2.0*lm*(-3*z33/R5 + (1.0/R3))/(c1));
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

template void calc_halfspace_term_C<Real> (
  const Real lam, const Real mu, const Real v[3], const Real d[3],
  const Real& y3, const Real& x1, const Real& x3, Real sigma[6]);

#ifdef WOODLAND_ACORN_VECTORIZE
# message "ERROR: WOODLAND_ACORN_VECTORIZE is not supported yet."
error
#endif

} // namespace hs3d
} // namespace acorn
} // namespace woodland
