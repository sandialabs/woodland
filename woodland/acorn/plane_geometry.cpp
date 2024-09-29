#include "woodland/acorn/plane_geometry.hpp"

#include <algorithm>

#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace acorn {
namespace plane {

typedef Matvec2d<Real> mv2;

bool is_ccw (const Pt a, const Pt b, const Pt c) {
  Pt u, v;
  mv2::subtract(b, a, u);
  mv2::subtract(c, a, v);
  return mv2::cross_given_z0(u, v) >= 0;
}

bool is_ccw (const Polygon& p, const bool assume_noncollinear) {
  Pt a, b;
  mv2::subtract(p.xys+2, p.xys, a);
  mv2::subtract(p.xys+4, p.xys, b);
  auto z = mv2::cross_given_z0(a, b);
  if (assume_noncollinear) return z > 0;
  for (int iv = 6; iv < 2*p.n; ++iv) {
    mv2::subtract(p.xys+iv, p.xys, b);
    const auto z1 = mv2::cross_given_z0(a, b);
    if (std::abs(z1) > std::abs(z)) z = z1;
  }
  return z > 0;
}

bool is_point_to_right_of_ray (const Pt p, const Pt v1, const Pt v2) {
  return is_ccw(v1, p, v2);
}

bool is_point_inside_convex_polygon (const Polygon& p, const Pt c) {
  const bool p_ccw = is_ccw(p);
  for (int i = 0; i < p.n; ++i) {
    const auto v1 = &p.xys[2*i];
    const auto v2 = &p.xys[2*((i+1) % p.n)];
    const bool c_inside =
      p_ccw != plane::is_point_to_right_of_ray(c, v1, v2);
    if ( ! c_inside) return false;
  }
  return true;
}

void project_p_onto_line (const Pt p, const Pt v1, const Pt v2,
                          Real& a, Pt projection) {
  // Project p onto the line by finding a such that
  //     dot((v1 + v2)/2 + a (v2 - v1) - p, v2 - v1) = 0.
  Pt e;
  mv2::subtract(v2, v1, e);
  // Make the solution finite-precision invariant to order of v1 and v2.
  Pt vmid;
  mv2::add(v1, v2, vmid);
  mv2::scale(0.5, vmid);
  const auto pe = mv2::dot(p, e), ee = mv2::dot(e, e), vme = mv2::dot(vmid, e);
  a = (pe - vme)/ee;
  mv2::axpbyz(1, vmid, a, e, projection);
}

Real calc_shortest_distance2_between_point_and_line (
  const Pt p, const Pt v1, const Pt v2)
{
  Real a;
  Pt projection;
  project_p_onto_line(p, v1, v2, a, projection);
  return mv2::distance2(projection, p);
}

Real inscribe_largest_circle (const Polygon& p, const Pt s,
                              int* support_sides) {
  Real radius2 = calc_shortest_distance2_between_point_and_line(
    s, &p.xys[0], &p.xys[2]);
  const int nd = 2*p.n;
  if (support_sides) *support_sides = 0;
  for (int i = 1; i < p.n; ++i) {
    const auto r2 = calc_shortest_distance2_between_point_and_line(
      s, &p.xys[2*i], &p.xys[(2*(i+1)) % nd]);
    if (r2 < radius2) {
      radius2 = r2;
      if (support_sides) *support_sides = i;
    }
  }
  return radius2;
}

static void solve_trig_eq (const Real d, const Real e, const Real f,
                           Real soln[2][2]) {
  // Algebra.
  Real x[2]; {
    const auto d2 = square(d);
    solve_quadratic_eq_for_real_roots(square(e) + d2, -2*e*f, square(f) - d2, x,
                                      true /* assume_real */);
  }
  // Enumerate the four candidates.
  Real cand[4][2];
  cand[0][1] = x[0]; cand[0][0] = std::sqrt(1 - square(x[0]));
  cand[1][1] = x[0]; cand[1][0] = -cand[0][0];
  cand[2][1] = x[1]; cand[2][0] = std::sqrt(1 - square(x[1]));
  cand[3][1] = x[1]; cand[3][0] = -cand[2][0];
  // Select the two candidates that satisfy the trig eq.
  Real residual[4];
  for (int i = 0; i < 4; ++i)
    residual[i] = std::abs(d*cand[i][0] + e*cand[i][1] - f);
  int idxs[2];
  idxs[0] = 0;
  Real rmin = residual[0], rmax = rmin;
  for (int i = 1; i < 4; ++i) {
    if (residual[i] < rmin) {
      idxs[0] = i;
      rmin = residual[i];
    }
    if (residual[i] > rmax) rmax = residual[i];
  }
  rmin = rmax;
  idxs[1] = -1;
  for (int i = 0; i < 4; ++i) {
    if (i == idxs[0]) continue;
    if (residual[i] < rmin) {
      idxs[1] = i;
      rmin = residual[i];
    }
  }
  if (idxs[1] == -1) {
    // All remaining residuals are the same. Choose the candidate that is
    // farthest from the first solution.
    Real max_dist2 = 0;
    for (int i = 0; i < 4; ++i) {
      if (i == idxs[0]) continue;
      const auto d2 = mv2::distance2(cand[idxs[0]], cand[i]);
      if (d2 > max_dist2) {
        max_dist2 = d2;
        idxs[1] = i;
      }
    }
    if (idxs[1] == -1) {
      // All solutions are the same.
      idxs[1] = idxs[0];
    }
  }
  assert(residual[idxs[0]] < 1e8*mv2::eps);
  assert(residual[idxs[1]] < 1e8*mv2::eps);
  for (int i = 0; i < 2; ++i)
    mv2::copy(cand[idxs[i]], soln[i]);
}

bool calc_tangent_points (const Pt a, const Pt cc, const Real radius2,
                          Pt pts[2]) {
  Pt c;
  mv2::subtract(cc, a, c);
  const auto s2 = mv2::norm22(c), L2 = s2 - radius2;
  if (L2 < 0) return false;
  const auto L = std::sqrt(L2);
  solve_trig_eq(L*c[0], L*c[1], L2, pts);
  for (int i = 0; i < 2; ++i) mv2::axpy(1, a, mv2::scale(L, pts[i]));
  return true;
}

bool intersect_line_line (const Pt xa, const Pt xb, const Pt ya, const Pt yb,
                          Pt p, Real* const p_xb) {
  /* Let
         p(a) := xa + a (xb - xa).
     Solve
         z := perpcw2(yb - ya)
         z'p(a) = z'ya
     for a:
         z'(xa + a (xb - xa) - ya) = 0
         a = z'(ya - xa) / z'(xb - xa).
   */
  Pt tmp, z;
  mv2::subtract(yb, ya, tmp);
  mv2::perpcw2(tmp, z);
  mv2::subtract(ya, xa, tmp);
  const auto num = mv2::dot(z, tmp);
  mv2::subtract(xb, xa, tmp);
  const auto den = mv2::dot(z, tmp);
  if (den == 0) return false;
  const auto a = num/den;
  mv2::axpbyz(1, xa, a, tmp, p);
  if (p_xb) {
    const auto f = -a/den;
    // a I + (xb - xa) a_xb
    p_xb[0] = a + tmp[0]*f*z[0];
    p_xb[1] =     tmp[0]*f*z[1];
    p_xb[2] =     tmp[1]*f*z[0];
    p_xb[3] = a + tmp[1]*f*z[1];
  }
  return true;
}

// Intersection of xa + a (xb - xa) with the nearer of the two points to xa on
// the circle, with radius^2 as input. Return false if the line does not
// intersect the circle, otherwise true.
bool intersect_line_circle_nearer (
  const Pt xa, const Pt xb, const Pt ctr, const Real radius2, Pt p,
  const bool assume_intersect, Real* const p_xb)
{
  /* Let
         p(a) := xa + a (xb - xa).
     Solve
         norm(p(a) - ctr, 2)^2 = r^2
     for a and take the nearer of the two solutions to xa:
         c := xa - ctr
         d := xb - xa
         0 = (c + a d)'(c + a d) - r^2
           = (c'c - r^2) + 2 a c'd + a^2 d'd
     => a = [-2 c'd +/- sqrt(4 (c'd)^2 - 4 (c'c - r^2) d'd)]/(2 d'd).
   */
  assert(radius2 >= 0);
  Pt c, d;
  mv2::subtract(xa, ctr, c);
  mv2::subtract(xb, xa , d);
  Real as[2];
  const auto dd = mv2::dot(d, d), dc = mv2::dot(d, c), cc = mv2::dot(c, c);
  if ( ! solve_quadratic_eq_for_real_roots(dd, 2*dc, cc - radius2, as,
                                           assume_intersect))
    return false;
  const auto a = std::abs(as[0]) <= std::abs(as[1]) ? as[0] : as[1];
  mv2::axpbyz(1, xa, a, d, p);
  if (p_xb) {
    // If norm(d) > 0 and c'd /= 0 (both of which are needed to make a solvable
    // problem), then den /= 0.
    const auto den = a*dd + dc;
    assert(den != 0);
    /* 0 = a^2 d'd + 2 a d'c + (c'c - r^2)
         = 2 a a_d (d'd) + 2 a^2 d + 2 a_d d'c + 2 a c
         = [a (d'd) + d'c] a_d + a^2 d + a c
       => a_y = -a (a d + c)/(a (d'd) + d'c)
     */
    const auto f = -a/den;
    Pt a_y;
    for (int i = 0; i < 2; ++i)
      a_y[i] = f*(a*d[i] + c[i]);
    // a I + (xb - xa) a_xb
    p_xb[0] = a + d[0]*a_y[0];
    p_xb[1] =     d[0]*a_y[1];
    p_xb[2] =     d[1]*a_y[0];
    p_xb[3] = a + d[1]*a_y[1];
  }
  return true;
}

void calc_nearest_point_on_polygon_boundary_to_point (
  const Polygon& p, const Pt pt, Pt pt_nearest, int& on_vertex, int& on_edge)
{
  Real dist2 = 0;
  for (int i = 0; i < p.n; ++i) {
    const auto v1 = &p.xys[2*i], v2 = &p.xys[2*((i+1) % p.n)];
    Real a;
    Pt pt_proj;
    project_p_onto_line(pt, v1, v2, a, pt_proj);
    int ovi = -1, oei = -1;
    if (a <= -0.5) {
      ovi = i;
      mv2::copy(v1, pt_proj);
    } else if (a >= 0.5) {
      ovi = i+1;
      mv2::copy(v2, pt_proj);
    } else {
      oei = i;
    }
    const Real idist2 = mv2::distance2(pt, pt_proj);
    if (i == 0 || idist2 < dist2) {
      dist2 = idist2;
      mv2::copy(pt_proj, pt_nearest);
      on_vertex = ovi;
      on_edge = oei;
    }
  }
}

static int test_inscribe_largest_circle () {
  int ne = 0;
  {
    // Convex polygon.
    const Real xys[] = {1.1,-2, 1.1,1, -0.9,1, -0.9,-1.5};
    static const int nv = sizeof(xys)/sizeof(*xys)/2;
    Polygon p(xys, nv);
    // Iterate over some points.
    for (const Real y : {0.0, 0.5, 1.0}) {
      const Real s[] = {0.1, y};
      for (const Real f : {0.1, 0.8, 1.3}) {
        // Rotate the problem a few times to check for rotational invariance.
        const auto cth = std::cos(f*M_PI), sth = std::sin(f*M_PI);
        Real rxys[2*nv], rs[2];
        Polygon rp(rxys, nv);
        std::copy(xys, xys + 2*nv, rxys);
        std::copy(s, s+2, rs);
        rotate_vector_2(cth, sth, rs);
        for (int i = 0; i < nv; ++i) rotate_vector_2(cth, sth, &rxys[2*i]);
        // Calculate and check.
        const auto radius2 = inscribe_largest_circle(rp, rs);
        if (std::abs(std::sqrt(radius2) - (1-y)) > 10*mv2::eps) ++ne;
      }
    }
  }
  if (ne) printf("test_inscribe_largest_circle failed\n");
  return ne;
}

static int test_intersect_line_line () {
  int ne = 0;
  const Pt a = {1, 1}, b = {-1, -0.5}, c = {0.5, 0.7};
  {
    Pt p;
    if ( ! intersect_line_line(a, b, a, c, p)) ++ne;
    if (reldif(2, p, a) > mv2::eps) ++ne;
    intersect_line_line(b, c, a, c, p);
    if (reldif(2, p, c) > mv2::eps) ++ne;
    if (intersect_line_line(a, b, a, b, p)) ++ne;
  }
  for (int i = 0, n = 100; i < n; ++i) {
    const auto th = 2*M_PI*Real(i)/n;
    const Pt xa = {std::cos(th), std::sin(th)};
    const Pt xb = {-(i+1)*xa[0], -(i+1)*xa[1]};
    const Pt ya = {0.1, -1}, yb = {-0.1, 1};
    Pt p;
    Real p_xb[4];
    if ( ! intersect_line_line(xa, xb, ya, yb, p, p_xb)) ++ne;
    { // test Jacobian
      Real p_xb_fd[4];
      for (int d = 0; d < 2; ++d) {
        static const Real e = 1e-4;
        Pt p_pe, p_me, xb_e, delta;
        mv2::copy(xb, xb_e);
        xb_e[d] += e;
        intersect_line_line(xa, xb_e, ya, yb, p_pe);
        xb_e[d] = xb[d] - e;
        intersect_line_line(xa, xb_e, ya, yb, p_me);
        mv2::subtract(p_pe, p_me, delta);
        p_xb_fd[  d] = 0.5*delta[0]/e;
        p_xb_fd[2+d] = 0.5*delta[1]/e;
      }
      if (reldif(4, p_xb_fd, p_xb) > 1e-7) ++ne;
    }
    if (absmax(2, p) > 5*mv2::eps) ++ne;
  }
  if (ne) printf("test_intersect_line_line failed\n");
  return ne;
}

static int test_intersect_line_circle_nearer () {
  int ne = 0;
  const Pt a = {-4, 3}, cc = {-2.5, 1.3};
  const Real r2 = square(0.9), da2cc2 = mv2::distance2(a, cc);
  for (int i = 0, n = 100; i < n; ++i) {
    const auto th = 2*M_PI*Real(i)/n;
    const Pt d = {std::cos(th), std::sin(th)};
    for (const Real scale : {0.1, 4.0}) {
      Pt b, p;
      Real p_b[4];
      mv2::axpbyz(1, a, scale, d, b);
      const auto dpline2 =
        calc_shortest_distance2_between_point_and_line(cc, a, b);
      const auto dpline2order =
        calc_shortest_distance2_between_point_and_line(cc, b, a);
      if (dpline2order != dpline2) ++ne;
      if (intersect_line_circle_nearer(a, b, cc, r2, p, false, p_b)) {
        if (std::abs(mv2::distance2(cc, p) - r2) > 20*mv2::eps) ++ne;
        if (mv2::distance2(a, p) > da2cc2) ++ne;
        if (dpline2 > r2) ++ne;
        { // test Jacobian
          Real p_b_fd[4];
          const Real e = 1e-7*mv2::distance(a, b);
          for (int d = 0; d < 2; ++d) {
            Pt p_pe, p_me, b_e, delta;
            mv2::copy(b, b_e);
            b_e[d] += e;
            intersect_line_circle_nearer(a, b_e, cc, r2, p_pe);
            b_e[d] = b[d] - e;
            intersect_line_circle_nearer(a, b_e, cc, r2, p_me);
            mv2::subtract(p_pe, p_me, delta);
            p_b_fd[  d] = 0.5*delta[0]/e;
            p_b_fd[2+d] = 0.5*delta[1]/e;
          }
          if (reldif(4, p_b_fd, p_b) > 1e-7) ++ne;
        }
      } else {
        if (dpline2 < r2) ++ne;
      }
    }
  }
  if (ne) printf("test_intersect_line_circle_nearer failed\n");
  return ne;
}

static int test_calc_nearest_point_on_polygon_boundary_to_point () {
  int ne = 0;
  // The polygon doesn't have to be convex.
  const int n = 7;
  Real vtxs[2*n];
  for (int i = 0; i < n; ++i) {
    const Real theta = Real(i)*2*M_PI/n;
    const Real r = 1 + (i % 3);
    vtxs[2*i+0] = r*std::cos(theta);
    vtxs[2*i+1] = r*std::sin(theta);
  }
  const Polygon poly(vtxs, n);
  for (int i = 0; i < n; ++i) {
    int on_vertex, on_edge;
    Pt pn;
    calc_nearest_point_on_polygon_boundary_to_point(poly, &vtxs[2*i],
                                                    pn, on_vertex, on_edge);
    if (on_vertex != i) ++ne;
    if (mv2::distance2(&vtxs[2*i], pn) > 0) ++ne;
  }
  const int ns = 20;
  for (int i = 0; i < ns; ++i) {
    const Real a = Real(i)/(ns-1);
    const Real y = -4*(1-a) + 4*a;
    for (int j = 0; j < ns; ++j) {
      const Real b = Real(i)/(ns-1);
      const Real x = -4*(1-b) + 4*b;
      const Pt pt{x,y};
      Real vdist = 10;
      for (int k = 0; k < n; ++k)
        vdist = std::min(vdist, mv2::distance(pt, &vtxs[2*k]));
      int on_vertex, on_edge;
      Pt pn;
      calc_nearest_point_on_polygon_boundary_to_point(poly, pt,
                                                      pn, on_vertex, on_edge);
      const Real dist = mv2::distance(pt, pn);
      if (dist > vdist) ++ne;
      if (on_edge >= 0) {
        const Real d1 = mv2::distance(pt, &vtxs[2*on_edge]);
        const Real d2 = mv2::distance(pt, &vtxs[2*((on_edge+1) % n)]);
        if (dist > std::min(d1, d2)) ++ne;
      }
      if (on_vertex >= 0 && mv2::distance(pn, &vtxs[2*on_vertex]) > 0) ++ne;
      Pt pnt;
      mv2::copy(pn, pnt);
      calc_nearest_point_on_polygon_boundary_to_point(poly, pnt,
                                                      pn, on_vertex, on_edge);
      if (mv2::distance2(pnt, pn) > 1e-15) ++ne;
    }
  }
  if (ne)
    printf("test_calc_nearest_point_on_polygon_boundary_to_point failed\n");
  return ne;
}

static int test_misc () {
  int ne = 0;
  {
    const Real xys[] = {0,0, 1,0, 1,1};
    const Polygon p(xys, 3);
    if ( ! is_ccw(p)) ++ne;
    if ( ! is_ccw(p), false) ++ne;
    if ( ! is_ccw(xys, xys+2, xys+4)) ++ne;
    if (is_ccw(xys+2, xys, xys+4)) ++ne;
  }
  {
    const Real xys[] = {0,0, 1,1, 1,0};
    const Polygon p(xys, 3);
    if (is_ccw(p)) ++ne;
    if (is_ccw(p), false) ++ne;
  }
  return ne;
}

int unittest () {
  int ne = 0;
  ne += test_inscribe_largest_circle();
  ne += test_intersect_line_line();
  ne += test_intersect_line_circle_nearer();
  ne += test_calc_nearest_point_on_polygon_boundary_to_point();
  ne += test_misc();
  return ne;
}

} // namespace plane
} // namespace acorn
} // namespace woodland
