#ifndef INCLUDE_WOODLAND_ACORN_TRIANGLE
#define INCLUDE_WOODLAND_ACORN_TRIANGLE

#include "woodland/acorn/acorn.hpp"
#include "woodland/acorn/matvec.hpp"

namespace woodland {
namespace acorn {

/* Operations on a single triangle in 2D space.

   A triangle is a list of three (x,y) 2-vectors and is usually called 't'. A
   point in x-y space is usually called 'p'.

   The barycentric coordinate is usually called 'lam', for 'lambda'. The
   barycentric matrix that transforms from (x,y) to lam is 'b'.

   The Clough-Tocher code is adapted from my package dc3dm. calc_coefs follows
   section 2 of ref [1].

   [1] P. Alfeld, A trivariate Clough-Tocher scheme, 1984.
 */

struct Triangle2D {
  typedef Matvec2d<Real> mv2;

  // The triangle is (x1,x2,x3).
  static void calc_barycentric_matrix (
    const Real x0[2], const Real x1[2], const Real x2[2], Real b[4])
  {
    mv2::subtract(x0, x2, b  );
    mv2::subtract(x1, x2, b+2);
    const auto det = b[0]*b[3] - b[1]*b[2];
    const auto tmp = b[0];
    b[0] =  b[3] / det;
    b[1] = -b[1] / det;
    b[2] = -b[2] / det;
    b[3] =   tmp / det;    
  }

  static void calc_barycentric_matrix (const Real t[6], Real b[4]) {
    calc_barycentric_matrix(t, t+2, t+4, b);
  }

  // x2 is the third vertex in the original t.
  static void xy_to_barycentric (const Real x2[2], const Real b[4],
                                 const Real p[2], Real lam[3]) {
    Real dx[2];
    mv2::subtract(p, x2, dx);
    lam[0] = b[0] * dx[0] + b[2] * dx[1];
    lam[1] = b[1] * dx[0] + b[3] * dx[1];
    lam[2] = 1.0 - lam[0] - lam[1];
  }

  static void barycentric_to_xy (
    const Real x0[2], const Real x1[2], const Real x2[2], const Real lam[3],
    Real p[2])
  { mv2::sum3(lam[0], x0, lam[1], x1, lam[2], x2, p); }

  static void barycentric_to_xy (const Real t[6], const Real lam[3], Real p[2]) {
    barycentric_to_xy(t, t+2, t+4, lam, p);
  }

  static void calc_grad_lam (const Real b[4], const Real lam[3],
                             // (lam[0]_x, lam[0]_y, ..., lam_2[]_y)
                             Real grad[6]) {
    grad[0] = b[0]; grad[1] = b[2];
    grad[2] = b[1]; grad[3] = b[3];
    grad[4] = -(grad[0] + grad[2]);
    grad[5] = -(grad[1] + grad[3]);
  }

  // Check if point p is in the triangle t. If label is non-null and p is not,
  // then print information.
  static bool is_p_in_t(const Real t[6], const Real b[4], const Real p[2],
                        const char* const label = nullptr,
                        Real tol = -1);

  // >= 0 for CCW list of points.
  static Real calc_signed_area(const Real t[6]);
  static Real calc_signed_area(const Real a[2], const Real b[2],
                               const Real c[2]);

  static int unittest();
};

} // namespace acorn
} // namespace woodland

#endif
