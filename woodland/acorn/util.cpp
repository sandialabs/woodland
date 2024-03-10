#include "woodland/acorn/util.hpp"
#include "woodland/acorn/macros.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/dbg.hpp"

#include <cassert>
#include <cstdio>
#include <cstdarg>
#include <limits>

namespace woodland {
namespace acorn {

typedef Matvec<3,Real> mv3;

double urand () { return rand() / double(RAND_MAX); }

void rotate_vector_3 (const Real u[3], const Real a, Real v[3]) {
  Real u_cross_v[3], ucv_cross_u[3];
  mv3::cross(u, v, u_cross_v);
  mv3::cross(u_cross_v, u, ucv_cross_u);
  const Real u_dot_v = mv3::dot(u, v);
  const Real c = std::cos(a), s = std::sin(a);
  mv3::sum3(u_dot_v, u, c, ucv_cross_u, s, u_cross_v, v);
}

void form_rotation_3 (const Real u[3], const Real a, Real R[9]) {
  const Real c = std::cos(a), omc = 1-c, s = std::sin(a);
  const Real x = u[0], y = u[1], z = u[2];
  R[0] = c + x*x*omc;
  R[1] = x*y*omc - z*s;
  R[2] = x*z*omc + y*s;
  R[3] = y*x*omc + z*s;
  R[4] = c + y*y*omc;
  R[5] = y*z*omc - x*s;
  R[6] = z*x*omc - y*s;
  R[7] = z*y*omc + x*s;
  R[8] = c + z*z*omc;
}

static int test_rotate () {
  int ne = 0;
  {
    Real axis[] = {-0.1, 0.7, 0.3};
    mv3::normalize(axis);
    const Real theta = -0.3*M_PI;
    const Real v[] = {0.2, -1.3, 0.7};
    Real vr[3];
    mv3::copy(v, vr);
    rotate_vector_3(axis, theta, vr);
    const auto nv2 = mv3::norm22(v), nvr2 = mv3::norm22(vr);
    if (reldif(nv2, nvr2) > 2*mv3::eps) ++ne;
    const auto d = mv3::dot(axis, v), dr = mv3::dot(axis, vr);
    if (reldif(d, dr) > 2*mv3::eps) ++ne;
    Real R[9];
    form_rotation_3(axis, theta, R);
    Real vr1[3];
    mv3::matvec(R, v, vr1);
    if (reldif(3, vr, vr1) > 2*mv3::eps) ++ne;
  }
  {
    const Real axis[] = {0, 1, 0};
    Real v[] = {0, 0, 1};
    rotate_vector_3(axis, M_PI/2, v);
    const Real xhat[] = {1, 0, 0}, zhat[] = {0, 0, 1};
    if (std::abs(mv3::dot(xhat, v) - 1) > 2*mv3::eps) ++ne;
    rotate_vector_3(axis, -M_PI/2, v);
    if (std::abs(mv3::dot(zhat, v) - 1) > 2*mv3::eps) ++ne;
  }
  return ne;
}

void xyb_to_FlatSegments (const int ne, const Real* const xb, const Real* const yb,
                          std::vector<FlatSegment>& els) {
  els.resize(ne);
  /*ompparfor*/ for (int i = 0; i < ne; ++i) {
    auto& e = els[i];
    e.xc = (xb[i] + xb[i+1])/2;
    e.dx = (xb[i+1] - xb[i])/2;
    e.yc = (yb[i] + yb[i+1])/2;
    e.dy = (yb[i+1] - yb[i])/2;
    e.reset_hl();
  }
}

void calc_parabola_coefs (const Real dx[2], const Real y[3], Real c[3]) {
  assert(dx[0] < 0 && dx[1] > 0);
  const Real
    d1 = dx[0], d12 = square(d1),
    d2 = dx[1], d22 = square(d2),
    det = d1*d22 - d12*d2,
    dy1 = y[0] - y[1], dy2 = y[2] - y[1];
  assert(det < 0);
  c[0] = y[1];
  c[1] = ( d22*dy1 - d12*dy2)/det;
  c[2] = (-d2 *dy1 + d1 *dy2)/det;
}

static int test_calc_parabola_coefs () {
  int nerr = 0;
  const Real eps = std::numeric_limits<Real>::epsilon(), tol = 10*eps;
  {
    int ne = 0;
    const Real x[] = {-1, 0, 2}, y[] = {1, 1, 1};
    const Real dx[] = {x[0]-x[1], x[2]-x[1]};
    Real c[3];
    calc_parabola_coefs(dx, y, c);
    if (reldif(c[0], 1.0) > tol) ++ne;
    if (std::abs(c[1]) > tol) ++ne;
    if (std::abs(c[2]) > tol) ++ne;
    if (ne) printf("test_calc_parabola_coefs: const test failed\n");
    nerr += ne;
  }
  {
    int ne = 0;
    const Real x[] = {-1, 0, 2.5}, y[] = {-3.3, -1, 3};
    const Real dx[] = {x[0]-x[1], x[2]-x[1]};
    const Real ctrue[] = {-1, 2.1, -0.2};
    Real c[3];
    calc_parabola_coefs(dx, y, c);
    for (int i = 0; i < 3; ++i)
      if (reldif(c[i], ctrue[i]) > tol) ++ne;
    if (ne) printf("test_calc_parabola_coefs: general test failed\n");
    nerr += ne;
  }
  return nerr;
}

Real cosxm1 (Real x) {
  // The Taylor series to fourth order of y(x) is
  //     1/2 x^2 - 1/24 x^4.
  // Switch between the Taylor series to second order and the trig when
  //     1/2 x^2 - [1 - cos(x)] = eps
  //     => 1/2 x^2 - [1/2 x^2 - 1/24 x^4] = eps
  //     => x = (24 eps)^(1/4).
  static const auto
    eps = std::numeric_limits<Real>::epsilon(),
    x_co = std::pow(24*eps, 0.25);
  x = std::abs(x);
  return (x <= x_co ?
          -0.5*square(x) :
          std::cos(x) - 1);
}

// See, e.g., https://people.csail.mit.edu/bkph/articles/Quadratics.pdf.
bool solve_quadratic_eq_for_real_roots (const Real a, const Real b, const Real c,
                                        Real x[2], const bool assume_real) {
  const auto d2 = square(b) - 4*a*c;
  if (d2 < 0 and ! assume_real) return false;
  const auto d = d2 <= 0 ? 0 : std::sqrt(d2);
  if (b >= 0) {
    const auto t = -(b + d);
    x[0] = t/(2*a);
    x[1] = t == 0 ? x[0] : (2*c)/t;
  } else {
    const auto t = d - b;
    x[1] = t/(2*a);
    x[0] = t == 0 ? x[1] : (2*c)/t;
  }
  return true;
}

static int test_quadratic_eq () {
  static const auto eps = std::numeric_limits<Real>::epsilon();
  int ne = 0;
  for (const Real a : {-1.0, 0.5})
    for (const Real b : {-2, 1}) {
      const Real y = 0.3, c = -(a*sqrt(y) + b*y);
      Real x[2];
      const auto stat = solve_quadratic_eq_for_real_roots(a, b, c, x);
      if ( ! stat) ++ne;
      for (int i = 0; i < 2; ++i) {
        const auto res = std::abs(c + x[i]*(b + x[i]*a));
        if (res > 10*eps) {
          printf("test_quadratic_eq: %d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",
                 int(stat), a, b, c, x[0], x[1], res);
          ++ne;
        }
      }
    }
  return ne;
}

static int test_cosxm1 () {
  static const auto eps = std::numeric_limits<Real>::epsilon();
  int ne = 0;
  for (const Real xbase : {2, 3, 7})
    for (const int p : {0, 3, 7, 10, 14}) {
      const Real
        x = xbase*std::pow(10, -p),
        y_sf = std::cos(x) - 1,
        y = cosxm1(x);
      if (std::abs(y - y_sf) > eps) ++ne;
    }
  if (ne) printf("test_cosxm1: ne %d\n", ne);
  return ne;
}

Real sqrt1pxm1 (const Real x) {
  static const auto
    eps = std::numeric_limits<Real>::epsilon(),
    x_co = std::sqrt(eps);
  return (std::abs(x) < x_co ?
          (-1.0/8.0*x + 1.0/2.0)*x :
          std::sqrt(1 + x) - 1);
}

static int test_sqrt1pxm1 () {
#if 0
  // Decided cutoff based on this:
  for (int p = 0; p <= 20; ++p) {
    const Real
      x = std::pow(10, -p),
      y_sf = std::sqrt(1 + x) - 1,
      y = sqrt1pxm1(x);
    printf("%2d %22.15e %22.15e %9.2e %9.2e\n", p, y_sf, y, y - y_sf, y-0.5*x);
  }
#endif
  static const auto eps = std::numeric_limits<Real>::epsilon();
  int ne = 0;
  for (const Real xbase : {2, 3, 7})
    for (const int p : {0, 3, 7, 10, 14}) {
      const Real
        x = xbase*std::pow(10, -p),
        y_sf = std::sqrt(1 + x) - 1,
        y = sqrt1pxm1(x);
      if (std::abs(y - y_sf) > eps) ++ne;
    }
  if (ne) printf("test_sqrt1pxm1: ne %d\n", ne);
  return ne;
}

void matvec (const Real x[3], const Real y[3], const Real z[3],
             const Real u[3], Real v[3]) {
  v[0] = x[0]*u[0] + x[1]*u[1] + x[2]*u[2];
  v[1] = y[0]*u[0] + y[1]*u[1] + y[2]*u[2];
  v[2] = z[0]*u[0] + z[1]*u[1] + z[2]*u[2];
}

void tmatvec (const Real x[3], const Real y[3], const Real z[3],
              const Real u[3], Real v[3]) {
  v[0] = x[0]*u[0] + y[0]*u[1] + z[0]*u[2];
  v[1] = x[1]*u[0] + y[1]*u[1] + z[1]*u[2];
  v[2] = x[2]*u[0] + y[2]*u[1] + z[2]*u[2];
}

std::vector<std::string> split (std::string s, const std::string& delim) {
  std::vector<std::string> toks;
  const auto dsz = delim.size();
  for (;;) {
    const auto p = s.find(delim);
    if (p == std::string::npos) break;
    toks.push_back(s.substr(0, p));
    s.erase(0, p + dsz);
  }
  toks.push_back(s);
  return toks;
}

static int test_split () {
  const auto toks = split("  hi how 3.14?  7", " ");
  int ne = 0;
  if (toks.size() != 7) ++ne;
  const int idxs[] = {2,3,4,6};
  const char* strs[] = {"hi", "how", "3.14?", "7"};
  for (size_t i = 0; i < sizeof(idxs)/sizeof(*idxs); ++i)
    if (toks[idxs[i]] != strs[i]) ++ne;
  return ne;
}

void Sprinter::add (const char* format, ...) {
  va_list args;
  va_start(args, format);
  for (;;) {
    const int ncur = vsnprintf(&buf[n], nbuf-n, format, args);
    if (ncur < nbuf-n) {
      n += ncur;
      break;
    }
    nbuf *= 2;
    buf.resize(nbuf);
  }
  va_end(args);
  assert(n < nbuf);
}

void Sprinter::out (FILE* stream, const bool newline) {
  const char* format = newline ? "%s\n" : "%s";
  fprintf(stream, format, buf.data());
}

int util_test () {
  int nerr = 0;
  nerr += test_calc_parabola_coefs();
  nerr += test_quadratic_eq();
  nerr += test_cosxm1();
  nerr += test_sqrt1pxm1();
  nerr += test_rotate();
  nerr += test_split();
  return nerr;
}

} // namespace acorn
} // namespace woodland
