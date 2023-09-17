#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/macros.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/dbg.hpp"

#include <cassert>
#include <cstdio>
#include <algorithm>

namespace woodland {
namespace acorn {

bool Triangle2D::is_p_in_t (const Real t[6], const Real b[4], const Real p[2],
                            const char* const label, Real tol) {
  if (tol < 0) tol = 10*Triangle2D::mv2::eps;
  Real lam[3];
  Triangle2D::xy_to_barycentric(t+4, b, p, lam);
  bool in = true;
  for (int i = 0; i < 3; ++i)
    if (lam[i] > 1 + tol || lam[i] < -tol)
      in = false;
  if ( ! in && label)
    fprintf(stderr, "p not in t for '%s'\n"
            "  tri: %22.15e %22.15e\n %28.15e %22.15e\n %28.15e %22.15e\n"
            "    p: %22.15e %22.15e\n"
            "  lam: %22.15e %22.15e %22.15e\n",
            label, t[0], t[1], t[2], t[3], t[4], t[5], p[0], p[1],
            lam[0], lam[1], lam[2]);
  return in;
}

Real Triangle2D::calc_signed_area (const Real t[6]) {
  return calc_signed_area(t, t+2, t+4);
}

Real Triangle2D::calc_signed_area (const Real u[2], const Real v[2],
                                   const Real w[2]) {
  Real a[2], b[2];
  mv2::subtract(v, u, a);
  mv2::subtract(w, u, b);
  return (a[0]*b[1] - a[1]*b[0])/2;
}

int Triangle2D::unittest ()  {
  int nerr = 0;

  {
    const Real t[] = {-1, -0.3, 2, 0.1, 0, 1}, p[] = {-0.5, 0};
    Real b[4], lam[3], p1[2];

    calc_barycentric_matrix(t, b);
    xy_to_barycentric(t+4, b, p, lam);
    barycentric_to_xy(t, lam, p1);

    if (lam[0] < 0 || lam[1] < 0 || lam[2] < 0) ++nerr;
    if (std::max(std::abs(p1[0] - p[0]), std::abs(p1[1] - p[1])) > 5*mv2::eps)
      ++nerr;
    if ( ! is_p_in_t(t, b, p1, "Triangle2D::unittest")) ++nerr;

    Real pout[] = {-1.1, -0.2};
    if (is_p_in_t(t, b, pout)) ++nerr;
  }

  {
    const Real t[] = {0,0, 1,0, 1,1};
    if (std::abs(calc_signed_area(t) - 0.5) > 10*mv2::eps) ++nerr;
  }
  {
    const Real t[] = {0,0, 0.5,1, 1,0};
    if (std::abs(calc_signed_area(t) + 0.5) > 10*mv2::eps) ++nerr;
  }

  return nerr;
}

} // namespace acorn
} // namespace woodland
