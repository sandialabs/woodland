#ifndef INCLUDE_WOODLAND_ACORN_MATVEC
#define INCLUDE_WOODLAND_ACORN_MATVEC

#include <cassert>
#include <cmath>
#include <limits>

namespace woodland {
namespace acorn {

template <int dim, typename Real>
struct Matvec {
  typedef const Real* const CRPtr;
  typedef Real* const RPtr;

  static constexpr Real eps = std::numeric_limits<Real>::epsilon();

  static void assign (const Real value, RPtr d) {
    for (int i = 0; i < dim; ++i) d[i] = value;
  }

  static void copy (CRPtr s, RPtr d) {
    for (int i = 0; i < dim; ++i) d[i] = s[i];
  }

  static Real dot (CRPtr a, CRPtr b) {
    Real d = 0;
    for (int i = 0; i < dim; ++i) d += a[i]*b[i];
    return d;
  }

  // Row-major dxd A.
  static void matvec (CRPtr A, CRPtr x, RPtr y) {
    for (int d = 0; d < dim; ++d) y[d] = dot(A+d*dim, x);
  }
  static void tmatvec (CRPtr A, CRPtr x, RPtr y) {
    for (int d = 0; d < dim; ++d) y[d] = 0;
    for (int c = 0; c < dim; ++c)
      for (int r = 0; r < dim; ++r)
        y[r] += A[c*dim + r]*x[c];
  }

  static Real norm22 (CRPtr v) { return dot(v, v); }
  static Real norm2  (CRPtr v) { return std::sqrt(norm22(v)); }

  static void normalize (RPtr v) {
    const Real nrm = std::sqrt(norm22(v));
    for (int i = 0; i < dim; ++i) v[i] /= nrm;
  };

  static RPtr scale (const Real scale, RPtr v) {
    for (int i = 0; i < dim; ++i) v[i] *= scale;
    return v;
  }

  // c = a + b
  static void add (CRPtr a, CRPtr b, RPtr c) {
    for (int i = 0; i < dim; ++i) c[i] = a[i] + b[i];
  }

  // c = a - b
  static void subtract (CRPtr a, CRPtr b, RPtr c) {
    for (int i = 0; i < dim; ++i) c[i] = a[i] - b[i];
  }

  // y += a x
  static void axpy (const Real a, CRPtr x, RPtr y) {
    for (int i = 0; i < dim; ++i) y[i] += a*x[i];
  }

  // z = a x + b y
  static void axpbyz (const Real a, CRPtr x, const Real b, CRPtr y, RPtr z) {
    for (int i = 0; i < dim; ++i) z[i] = a*x[i] + b*y[i];
  }

  // z = a u + b v
  static void sum2 (const Real a, CRPtr u, const Real b, CRPtr v, RPtr z) {
    for (int i = 0; i < dim; ++i) z[i] = a*u[i] + b*v[i];
  }

  // z = a u + b v + c w
  static void sum3 (const Real a, CRPtr u, const Real b, CRPtr v, const Real c,
                    CRPtr w, RPtr z) {
    for (int i = 0; i < dim; ++i) z[i] = a*u[i] + b*v[i] + c*w[i];
  }

  static void perpcw2 (CRPtr v, RPtr n) {
    assert(dim == 2);
    n[0] =  v[1];
    n[1] = -v[0];
  }

  static void perpccw2 (CRPtr v, RPtr n) {
    assert(dim == 2);
    n[0] = -v[1];
    n[1] =  v[0];
  }

  // c = a x b
  static RPtr cross (CRPtr a, CRPtr b, RPtr c) {
    assert(dim == 3);
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
  }

  // Return the z component of [a; 0] x [b; 0].
  static Real cross_given_z0 (CRPtr a, CRPtr b) {
    assert(dim == 2);
    return a[0]*b[1] - a[1]*b[0];
  }

  // norm(a - b, 2)^2
  static Real distance2 (CRPtr a, CRPtr b) {
    Real c[dim];
    subtract(a, b, c);
    return norm22(c);
  }

  // norm(a - b, 2)
  static Real distance (CRPtr a, CRPtr b) {
    return std::sqrt(distance2(a, b));
  }

  static int unittest();
};

template <typename Real> using Matvec2d = Matvec<2,Real>;

} // namespace acorn
} // namespace woodland

#endif
