#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/vectorization.hpp"

namespace woodland {
namespace acorn {
namespace fs3d {

template <typename RealT>
void calc_sigma_point (const Real lam, const Real mu, const RealT y[3],
                       const RealT v[3], const RealT d[3], const RealT x[3],
                       RealT sigma[6]) {
  const Real
    lm = lam + mu,
    l2m = lam + 2*mu,
    l3m = lam + 3*mu;
  const RealT
    x1 = x[0], x2 = x[1], x3 = x[2],
    y1 = y[0], y2 = y[1], y3 = y[2],
    z1 = x1 - y1,
    z2 = x2 - y2,
    z3 = x3 - y3,
    z11 = square(z1),
    z111 = z1*z11,
    z12 = z1*z2,
    z13 = z1*z3,
    z22 = square(z2),
    z222 = z2*z22,
    z23 = z2*z3,
    z33 = square(z3),
    z333 = z3*z33,
    z11pz22 = z11 + z22,
    z11pz33 = z11 + z33,
    z22pz33 = z22 + z33,
    R2 = z11 + z22 + z33,
    R4 = square(R2),
    R7 = acorn::pow(R2, 7.0/2.0),
    den1 = M_PI*R7*l2m*mu,
    den2 = M_PI*R7*l2m;
  RealT T_x[81];
T_x[0] = (1.0/8.0)*(-lam*(2*R4*(2*lam + lm + 4*mu) + 3*R2*(-2*R2*l2m - 8*l2m*z11 - 2*lm*z11 - lm*z22 - lm*z33 + lm*z22pz33) + 15*z1*(lm*z1*z22 + lm*z1*z33 + z1*(2*R2*l2m - lm*z22pz33))) + 2*mu*(-4*R4*l2m + 3*R2*(2*R2*l2m + 8*l2m*z11 - lm*z22pz33) - 15*z11*(2*R2*l2m - lm*z22pz33)))/den1;
T_x[1] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z12 + 2*l3m*z12 + 4*lm*z12) + 5*z2*(lm*z1*z22 + lm*z1*z33 + z1*(2*R2*l2m - lm*z22pz33))) + 2*mu*z12*(-10*R2*l2m + 2*R2*(3*lam + 7*mu) + 5*lm*z22pz33))/den1;
T_x[2] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z13 + 2*l3m*z13 + 4*lm*z13) + 5*z3*(lm*z1*z22 + lm*z1*z33 + z1*(2*R2*l2m - lm*z22pz33))) + 2*mu*z13*(-10*R2*l2m + 2*R2*(3*lam + 7*mu) + 5*lm*z22pz33))/den1;
T_x[3] = (3.0/8.0)*(R2*(4*l2m*z12 + 2*l3m*z12 + 3*lm*z12) - 5*z1*(lm*z11*z2 + z2*(2*R2*l2m - lm*z22pz33)))/den2;
T_x[4] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z22 + lm*z11 + lm*z22 - lm*z22pz33) - 15*z2*(lm*z11*z2 + z2*(2*R2*l2m - lm*z22pz33)))/den2;
T_x[5] = (3.0/8.0)*z3*(R2*(4*l3m*z2 + lm*z2) - 5*lm*z11*z2 - 5*z2*(2*R2*l2m - lm*z22pz33))/den2;
T_x[6] = (3.0/8.0)*(R2*(4*l2m*z13 + 2*l3m*z13 + 3*lm*z13) - 5*z1*(lm*z11*z3 + z3*(2*R2*l2m - lm*z22pz33)))/den2;
T_x[7] = (3.0/8.0)*z2*(R2*(4*l3m*z3 + lm*z3) - 5*lm*z11*z3 - 5*z3*(2*R2*l2m - lm*z22pz33))/den2;
T_x[8] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z33 + lm*z11 + lm*z33 - lm*z22pz33) - 15*z3*(lm*z11*z3 + z3*(2*R2*l2m - lm*z22pz33)))/den2;
T_x[9] = (3.0/8.0)*(R2*(4*l2m*z12 + 2*l3m*z12 + 3*lm*z12) - 5*z1*(lm*z11*z2 + z2*(2*R2*l2m - lm*z22pz33)))/den2;
T_x[10] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z22 + lm*z11 + lm*z22 - lm*z22pz33) - 15*z2*(lm*z11*z2 + z2*(2*R2*l2m - lm*z22pz33)))/den2;
T_x[11] = (3.0/8.0)*z3*(R2*(4*l3m*z2 + lm*z2) - 5*lm*z11*z2 - 5*z2*(2*R2*l2m - lm*z22pz33))/den2;
T_x[12] = (1.0/8.0)*(-lam*(2*R4*(2*lam + lm + 4*mu) + 3*R2*(-2*R2*l2m - 8*l2m*z11 - 2*lm*z11 - lm*z22 - lm*z33 + lm*z22pz33) + 15*z1*(lm*z1*z22 + lm*z1*z33 + z1*(2*R2*l2m - lm*z22pz33))) + 2*lm*mu*(-R4 + 3*R2*z11pz22 - 15*z11*z22))/den1;
T_x[13] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z12 + 2*l3m*z12 + 4*lm*z12) + 5*z2*(lm*z1*z22 + lm*z1*z33 + z1*(2*R2*l2m - lm*z22pz33))) + 2*lm*mu*z1*(3*R2*z2 - 5*z222))/den1;
T_x[14] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z13 + 2*l3m*z13 + 4*lm*z13) + 5*z3*(lm*z1*z22 + lm*z1*z33 + z1*(2*R2*l2m - lm*z22pz33))) + 2*lm*mu*z13*(R2 - 5*z22))/den1;
T_x[15] = (3.0/4.0)*lm*z23*(R2 - 5*z11)/den2;
T_x[16] = (3.0/8.0)*lm*z1*(2*R2*z3 - 10*z22*z3)/den2;
T_x[17] = (3.0/8.0)*lm*z1*(2*R2*z2 - 10*z2*z33)/den2;
T_x[18] = (3.0/8.0)*(R2*(4*l2m*z13 + 2*l3m*z13 + 3*lm*z13) - 5*z1*(lm*z11*z3 + z3*(2*R2*l2m - lm*z22pz33)))/den2;
T_x[19] = (3.0/8.0)*z2*(R2*(4*l3m*z3 + lm*z3) - 5*lm*z11*z3 - 5*z3*(2*R2*l2m - lm*z22pz33))/den2;
T_x[20] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z33 + lm*z11 + lm*z33 - lm*z22pz33) - 15*z3*(lm*z11*z3 + z3*(2*R2*l2m - lm*z22pz33)))/den2;
T_x[21] = (3.0/4.0)*lm*z23*(R2 - 5*z11)/den2;
T_x[22] = (3.0/8.0)*lm*z1*(2*R2*z3 - 10*z22*z3)/den2;
T_x[23] = (3.0/8.0)*lm*z1*(2*R2*z2 - 10*z2*z33)/den2;
T_x[24] = (1.0/8.0)*(-lam*(2*R4*(2*lam + lm + 4*mu) + 3*R2*(-2*R2*l2m - 8*l2m*z11 - 2*lm*z11 - lm*z22 - lm*z33 + lm*z22pz33) + 15*z1*(lm*z1*z22 + lm*z1*z33 + z1*(2*R2*l2m - lm*z22pz33))) + 2*lm*mu*(-R4 + 3*R2*z11pz33 - 15*z11*z33))/den1;
T_x[25] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z12 + 2*l3m*z12 + 4*lm*z12) + 5*z2*(lm*z1*z22 + lm*z1*z33 + z1*(2*R2*l2m - lm*z22pz33))) + 2*lm*mu*z12*(R2 - 5*z33))/den1;
T_x[26] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z13 + 2*l3m*z13 + 4*lm*z13) + 5*z3*(lm*z1*z22 + lm*z1*z33 + z1*(2*R2*l2m - lm*z22pz33))) + 2*lm*mu*z1*(3*R2*z3 - 5*z333))/den1;
T_x[27] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z12 + 2*l3m*z12 + 4*lm*z12) + 5*z1*(lm*z11*z2 + lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz33))) + 2*lm*mu*z2*(3*R2*z1 - 5*z111))/den1;
T_x[28] = (1.0/8.0)*(-lam*(2*R4*(2*lam + lm + 4*mu) + 3*R2*(-2*R2*l2m - 8*l2m*z22 - lm*z11 - 2*lm*z22 - lm*z33 + lm*z11pz33) + 15*z2*(lm*z11*z2 + lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz33))) + 2*lm*mu*(-R4 + 3*R2*z11pz22 - 15*z11*z22))/den1;
T_x[29] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z23 + 2*l3m*z23 + 4*lm*z23) + 5*z3*(lm*z11*z2 + lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz33))) + 2*lm*mu*z23*(R2 - 5*z11))/den1;
T_x[30] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z11 + lm*z11 + lm*z22 - lm*z11pz33) - 15*z1*(lm*z1*z22 + z1*(2*R2*l2m - lm*z11pz33)))/den2;
T_x[31] = (3.0/8.0)*(R2*(4*l2m*z12 + 2*l3m*z12 + 3*lm*z12) - 5*z2*(lm*z1*z22 + z1*(2*R2*l2m - lm*z11pz33)))/den2;
T_x[32] = (3.0/8.0)*z3*(R2*(4*l3m*z1 + lm*z1) - 5*lm*z1*z22 - 5*z1*(2*R2*l2m - lm*z11pz33))/den2;
T_x[33] = (3.0/8.0)*lm*z2*(2*R2*z3 - 10*z11*z3)/den2;
T_x[34] = (3.0/4.0)*lm*z13*(R2 - 5*z22)/den2;
T_x[35] = (3.0/8.0)*lm*z2*(2*R2*z1 - 10*z1*z33)/den2;
T_x[36] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z11 + lm*z11 + lm*z22 - lm*z11pz33) - 15*z1*(lm*z1*z22 + z1*(2*R2*l2m - lm*z11pz33)))/den2;
T_x[37] = (3.0/8.0)*(R2*(4*l2m*z12 + 2*l3m*z12 + 3*lm*z12) - 5*z2*(lm*z1*z22 + z1*(2*R2*l2m - lm*z11pz33)))/den2;
T_x[38] = (3.0/8.0)*z3*(R2*(4*l3m*z1 + lm*z1) - 5*lm*z1*z22 - 5*z1*(2*R2*l2m - lm*z11pz33))/den2;
T_x[39] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z12 + 2*l3m*z12 + 4*lm*z12) + 5*z1*(lm*z11*z2 + lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz33))) + 2*mu*z12*(-10*R2*l2m + 2*R2*(3*lam + 7*mu) + 5*lm*z11pz33))/den1;
T_x[40] = (1.0/8.0)*(-lam*(2*R4*(2*lam + lm + 4*mu) + 3*R2*(-2*R2*l2m - 8*l2m*z22 - lm*z11 - 2*lm*z22 - lm*z33 + lm*z11pz33) + 15*z2*(lm*z11*z2 + lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz33))) + 2*mu*(-4*R4*l2m + 3*R2*(2*R2*l2m + 8*l2m*z22 - lm*z11pz33) - 15*z22*(2*R2*l2m - lm*z11pz33)))/den1;
T_x[41] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z23 + 2*l3m*z23 + 4*lm*z23) + 5*z3*(lm*z11*z2 + lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz33))) + 2*mu*z23*(-10*R2*l2m + 2*R2*(3*lam + 7*mu) + 5*lm*z11pz33))/den1;
T_x[42] = (3.0/8.0)*z1*(R2*(4*l3m*z3 + lm*z3) - 5*lm*z22*z3 - 5*z3*(2*R2*l2m - lm*z11pz33))/den2;
T_x[43] = (3.0/8.0)*(R2*(4*l2m*z23 + 2*l3m*z23 + 3*lm*z23) - 5*z2*(lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz33)))/den2;
T_x[44] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z33 + lm*z22 + lm*z33 - lm*z11pz33) - 15*z3*(lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz33)))/den2;
T_x[45] = (3.0/8.0)*lm*z2*(2*R2*z3 - 10*z11*z3)/den2;
T_x[46] = (3.0/4.0)*lm*z13*(R2 - 5*z22)/den2;
T_x[47] = (3.0/8.0)*lm*z2*(2*R2*z1 - 10*z1*z33)/den2;
T_x[48] = (3.0/8.0)*z1*(R2*(4*l3m*z3 + lm*z3) - 5*lm*z22*z3 - 5*z3*(2*R2*l2m - lm*z11pz33))/den2;
T_x[49] = (3.0/8.0)*(R2*(4*l2m*z23 + 2*l3m*z23 + 3*lm*z23) - 5*z2*(lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz33)))/den2;
T_x[50] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z33 + lm*z22 + lm*z33 - lm*z11pz33) - 15*z3*(lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz33)))/den2;
T_x[51] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z12 + 2*l3m*z12 + 4*lm*z12) + 5*z1*(lm*z11*z2 + lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz33))) + 2*lm*mu*z12*(R2 - 5*z33))/den1;
T_x[52] = (1.0/8.0)*(-lam*(2*R4*(2*lam + lm + 4*mu) + 3*R2*(-2*R2*l2m - 8*l2m*z22 - lm*z11 - 2*lm*z22 - lm*z33 + lm*z11pz33) + 15*z2*(lm*z11*z2 + lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz33))) + 2*lm*mu*(-R4 + 3*R2*z22pz33 - 15*z22*z33))/den1;
T_x[53] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z23 + 2*l3m*z23 + 4*lm*z23) + 5*z3*(lm*z11*z2 + lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz33))) + 2*lm*mu*z2*(3*R2*z3 - 5*z333))/den1;
T_x[54] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z13 + 2*l3m*z13 + 4*lm*z13) + 5*z1*(lm*z11*z3 + lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz22))) + 2*lm*mu*z3*(3*R2*z1 - 5*z111))/den1;
T_x[55] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z23 + 2*l3m*z23 + 4*lm*z23) + 5*z2*(lm*z11*z3 + lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz22))) + 2*lm*mu*z23*(R2 - 5*z11))/den1;
T_x[56] = (1.0/8.0)*(-lam*(2*R4*(2*lam + lm + 4*mu) + 3*R2*(-2*R2*l2m - 8*l2m*z33 - lm*z11 - lm*z22 - 2*lm*z33 + lm*z11pz22) + 15*z3*(lm*z11*z3 + lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz22))) + 2*lm*mu*(-R4 + 3*R2*z11pz33 - 15*z11*z33))/den1;
T_x[57] = (3.0/8.0)*lm*z3*(2*R2*z2 - 10*z11*z2)/den2;
T_x[58] = (3.0/8.0)*lm*z3*(2*R2*z1 - 10*z1*z22)/den2;
T_x[59] = (3.0/4.0)*lm*z12*(R2 - 5*z33)/den2;
T_x[60] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z11 + lm*z11 + lm*z33 - lm*z11pz22) - 15*z1*(lm*z1*z33 + z1*(2*R2*l2m - lm*z11pz22)))/den2;
T_x[61] = (3.0/8.0)*z2*(R2*(4*l3m*z1 + lm*z1) - 5*lm*z1*z33 - 5*z1*(2*R2*l2m - lm*z11pz22))/den2;
T_x[62] = (3.0/8.0)*(R2*(4*l2m*z13 + 2*l3m*z13 + 3*lm*z13) - 5*z3*(lm*z1*z33 + z1*(2*R2*l2m - lm*z11pz22)))/den2;
T_x[63] = (3.0/8.0)*lm*z3*(2*R2*z2 - 10*z11*z2)/den2;
T_x[64] = (3.0/8.0)*lm*z3*(2*R2*z1 - 10*z1*z22)/den2;
T_x[65] = (3.0/4.0)*lm*z12*(R2 - 5*z33)/den2;
T_x[66] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z13 + 2*l3m*z13 + 4*lm*z13) + 5*z1*(lm*z11*z3 + lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz22))) + 2*lm*mu*z13*(R2 - 5*z22))/den1;
T_x[67] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z23 + 2*l3m*z23 + 4*lm*z23) + 5*z2*(lm*z11*z3 + lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz22))) + 2*lm*mu*z3*(3*R2*z2 - 5*z222))/den1;
T_x[68] = (1.0/8.0)*(-lam*(2*R4*(2*lam + lm + 4*mu) + 3*R2*(-2*R2*l2m - 8*l2m*z33 - lm*z11 - lm*z22 - 2*lm*z33 + lm*z11pz22) + 15*z3*(lm*z11*z3 + lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz22))) + 2*lm*mu*(-R4 + 3*R2*z22pz33 - 15*z22*z33))/den1;
T_x[69] = (3.0/8.0)*z1*(R2*(4*l3m*z2 + lm*z2) - 5*lm*z2*z33 - 5*z2*(2*R2*l2m - lm*z11pz22))/den2;
T_x[70] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z22 + lm*z22 + lm*z33 - lm*z11pz22) - 15*z2*(lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz22)))/den2;
T_x[71] = (3.0/8.0)*(R2*(4*l2m*z23 + 2*l3m*z23 + 3*lm*z23) - 5*z3*(lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz22)))/den2;
T_x[72] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z11 + lm*z11 + lm*z33 - lm*z11pz22) - 15*z1*(lm*z1*z33 + z1*(2*R2*l2m - lm*z11pz22)))/den2;
T_x[73] = (3.0/8.0)*z2*(R2*(4*l3m*z1 + lm*z1) - 5*lm*z1*z33 - 5*z1*(2*R2*l2m - lm*z11pz22))/den2;
T_x[74] = (3.0/8.0)*(R2*(4*l2m*z13 + 2*l3m*z13 + 3*lm*z13) - 5*z3*(lm*z1*z33 + z1*(2*R2*l2m - lm*z11pz22)))/den2;
T_x[75] = (3.0/8.0)*z1*(R2*(4*l3m*z2 + lm*z2) - 5*lm*z2*z33 - 5*z2*(2*R2*l2m - lm*z11pz22))/den2;
T_x[76] = (1.0/8.0)*(-R4*(2*lam + lm + 6*mu) + 3*R2*(2*R2*l2m + 4*l3m*z22 + lm*z22 + lm*z33 - lm*z11pz22) - 15*z2*(lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz22)))/den2;
T_x[77] = (3.0/8.0)*(R2*(4*l2m*z23 + 2*l3m*z23 + 3*lm*z23) - 5*z3*(lm*z2*z33 + z2*(2*R2*l2m - lm*z11pz22)))/den2;
T_x[78] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z13 + 2*l3m*z13 + 4*lm*z13) + 5*z1*(lm*z11*z3 + lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz22))) + 2*mu*z13*(-10*R2*l2m + 2*R2*(3*lam + 7*mu) + 5*lm*z11pz22))/den1;
T_x[79] = (3.0/8.0)*(-lam*(-R2*(4*l2m*z23 + 2*l3m*z23 + 4*lm*z23) + 5*z2*(lm*z11*z3 + lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz22))) + 2*mu*z23*(-10*R2*l2m + 2*R2*(3*lam + 7*mu) + 5*lm*z11pz22))/den1;
T_x[80] = (1.0/8.0)*(-lam*(2*R4*(2*lam + lm + 4*mu) + 3*R2*(-2*R2*l2m - 8*l2m*z33 - lm*z11 - lm*z22 - 2*lm*z33 + lm*z11pz22) + 15*z3*(lm*z11*z3 + lm*z22*z3 + z3*(2*R2*l2m - lm*z11pz22))) + 2*mu*(-4*R4*l2m + 3*R2*(2*R2*l2m + 8*l2m*z33 - lm*z11pz22) - 15*z33*(2*R2*l2m - lm*z11pz22)))/den1;
  RealT u_x[9];
  u_x[0] = T_x[0]*d[0]*v[0] + T_x[3]*d[0]*v[1] + T_x[6]*d[0]*v[2] + T_x[9]*d[1]*v[0] + T_x[12]*d[1]*v[1] + T_x[15]*d[1]*v[2] + T_x[18]*d[2]*v[0] + T_x[21]*d[2]*v[1] + T_x[24]*d[2]*v[2];
  u_x[1] = T_x[1]*d[0]*v[0] + T_x[4]*d[0]*v[1] + T_x[7]*d[0]*v[2] + T_x[10]*d[1]*v[0] + T_x[13]*d[1]*v[1] + T_x[16]*d[1]*v[2] + T_x[19]*d[2]*v[0] + T_x[22]*d[2]*v[1] + T_x[25]*d[2]*v[2];
  u_x[2] = T_x[2]*d[0]*v[0] + T_x[5]*d[0]*v[1] + T_x[8]*d[0]*v[2] + T_x[11]*d[1]*v[0] + T_x[14]*d[1]*v[1] + T_x[17]*d[1]*v[2] + T_x[20]*d[2]*v[0] + T_x[23]*d[2]*v[1] + T_x[26]*d[2]*v[2];
  u_x[3] = T_x[27]*d[0]*v[0] + T_x[30]*d[0]*v[1] + T_x[33]*d[0]*v[2] + T_x[36]*d[1]*v[0] + T_x[39]*d[1]*v[1] + T_x[42]*d[1]*v[2] + T_x[45]*d[2]*v[0] + T_x[48]*d[2]*v[1] + T_x[51]*d[2]*v[2];
  u_x[4] = T_x[28]*d[0]*v[0] + T_x[31]*d[0]*v[1] + T_x[34]*d[0]*v[2] + T_x[37]*d[1]*v[0] + T_x[40]*d[1]*v[1] + T_x[43]*d[1]*v[2] + T_x[46]*d[2]*v[0] + T_x[49]*d[2]*v[1] + T_x[52]*d[2]*v[2];
  u_x[5] = T_x[29]*d[0]*v[0] + T_x[32]*d[0]*v[1] + T_x[35]*d[0]*v[2] + T_x[38]*d[1]*v[0] + T_x[41]*d[1]*v[1] + T_x[44]*d[1]*v[2] + T_x[47]*d[2]*v[0] + T_x[50]*d[2]*v[1] + T_x[53]*d[2]*v[2];
  u_x[6] = T_x[54]*d[0]*v[0] + T_x[57]*d[0]*v[1] + T_x[60]*d[0]*v[2] + T_x[63]*d[1]*v[0] + T_x[66]*d[1]*v[1] + T_x[69]*d[1]*v[2] + T_x[72]*d[2]*v[0] + T_x[75]*d[2]*v[1] + T_x[78]*d[2]*v[2];
  u_x[7] = T_x[55]*d[0]*v[0] + T_x[58]*d[0]*v[1] + T_x[61]*d[0]*v[2] + T_x[64]*d[1]*v[0] + T_x[67]*d[1]*v[1] + T_x[70]*d[1]*v[2] + T_x[73]*d[2]*v[0] + T_x[76]*d[2]*v[1] + T_x[79]*d[2]*v[2];
  u_x[8] = T_x[56]*d[0]*v[0] + T_x[59]*d[0]*v[1] + T_x[62]*d[0]*v[2] + T_x[65]*d[1]*v[0] + T_x[68]*d[1]*v[1] + T_x[71]*d[1]*v[2] + T_x[74]*d[2]*v[0] + T_x[77]*d[2]*v[1] + T_x[80]*d[2]*v[2];
  const RealT udiv = u_x[0] + u_x[4] + u_x[8];
  sigma[0] = lam*udiv + mu*(u_x[0] + u_x[0]);
  sigma[1] = mu*(u_x[1] + u_x[3]);
  sigma[2] = mu*(u_x[2] + u_x[6]);
  sigma[3] = lam*udiv + mu*(u_x[4] + u_x[4]);
  sigma[4] = mu*(u_x[5] + u_x[7]);
  sigma[5] = lam*udiv + mu*(u_x[8] + u_x[8]);
}

template void
calc_sigma_point<Real>(const Real lam, const Real mu, const Real y[3], const Real v[3],
                       const Real d[3], const Real x[3], Real sigma[6]);

#ifdef WOODLAND_ACORN_VECTORIZE
template void
calc_sigma_point<RealPack>(const Real lam, const Real mu, const RealPack y[3],
                           const RealPack v[3], const RealPack d[3],
                           const RealPack x[3], RealPack sigma[6]);
#endif

} // namespace fs3d
} // namespace acorn
} // namespace woodland
