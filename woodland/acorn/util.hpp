#ifndef INCLUDE_WOODLAND_ACORN_UTIL
#define INCLUDE_WOODLAND_ACORN_UTIL

#include "woodland/acorn/acorn.hpp"
#include "woodland/acorn/macros.hpp"

#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <memory>

namespace woodland {
namespace acorn {

template <typename T> inline T square (const T& x) { return x*x; }
template <typename T> inline T cube (const T& x) { return x*x*x; }
template <typename T> inline T zero (const T&) { return T(0); }

double urand();

template <typename T> void copy (const int n, const T* const src, T* const dst) {
  assert(n >= 0);
  for (int i = 0; i < n; ++i) dst[i] = src[i];
}

template <typename Real> inline Real atanx0 (const Real y, const Real x) {
  return (x == 0) ? 0 : std::atan(y/x);
}

template <typename Real> inline Real reldif (const Real a, const Real b) {
  return std::abs(b - a)/std::max(std::abs(a), std::abs(b));
}

template <typename Real>
inline Real absmax (const int n, const Real* a) {
  Real num = 0;
  for (int i = 0; i < n; ++i)
    num = std::abs(a[i]);
  return num;
}

template <typename Real>
inline Real absdif (const int n, const Real* a, const Real* b) {
  Real num = 0;
  for (int i = 0; i < n; ++i)
    num = std::max(num, std::abs(b[i] - a[i]));
  return num;
}

template <typename Real>
inline Real reldif (const int n, const Real* a, const Real* b) {
  Real num = 0, den = 0;
  for (int i = 0; i < n; ++i) {
    num = std::max(num, std::abs(b[i] - a[i]));
    den = std::max(den, std::max(std::abs(a[i]), std::abs(b[i])));
  }
  return num/den;
}

// Compute y(x) = sum_{i=0}^{n-1} c[i] x^i.
template <typename Real> inline
Real evalpoly (const int n, const Real* const c, const Real& x) {
  Real y = c[n-1];
  for (int i = 2; i <= n; ++i)
    y = y*x + c[n-i];
  return y;
}

template <int N, typename Real> inline
Real evalpoly (const Real* const c, const Real& x) { return evalpoly(N, c, x); }

// For R = [c -s; s c] and A = [a11 a12; a12 a22], return R A R'.
inline void rotate_sym_tensor_2x2 (const Real c, const Real s,
                                   Real& a11, Real& a12, Real& a22) {
  const Real
    f1 = a11*c - a12*s,
    f2 = a12*c - a22*s,
    f3 = a12*s + a22*c,
    f4 = a11*s + a12*c;
  a11 = f1*c - f2*s;
  a12 = f2*c + f1*s;
  a22 = f3*c + f4*s;
}

// [v1; v2] = [c -s; s c] [v1; v2]
template <typename T>
inline void rotate_vector_2 (const T c, const T s, T& v1, T& v2) {
  const T t  = c*v1 - s*v2;
  v2 = s*v1 + c*v2;
  v1 = t;
}

template <typename T>
inline void rotate_vector_2 (const T c, const T s, T v[2])
{ rotate_vector_2(c, s, v[0], v[1]); }

template <typename T>
inline void rotate_transp_vector_2 (const T c, const T s, T v[2])
{ rotate_vector_2(c, -s, v); }

// axis must have unit length.
void rotate_vector_3(const Real axis[3], const Real angle, Real v[3]);

// For row-major R, R v = rotate_vector_3(axis, angle, v).
void form_rotation_3(const Real axis[3], const Real angle, Real R[9]);

// For row-major R = [x; y; z] and T = [t11 t12 t13; t12 t22 t23; t13 t23 t33],
// return R'T R.
template <typename T> inline void
rotate_sym_tensor_3x3_RtAR (const T* const x, const T* const y, const T* const z,
                            T& t11, T& t12, T& t13,
                            T& t22, T& t23, T& t33) {
  const T
    f1 = t11*x[0] + t12*y[0] + t13*z[0],
    f2 = t12*x[0] + t22*y[0] + t23*z[0],
    f3 = t13*x[0] + t23*y[0] + t33*z[0],
    f4 = t11*x[1] + t12*y[1] + t13*z[1],
    f5 = t12*x[1] + t22*y[1] + t23*z[1],
    f6 = t13*x[1] + t23*y[1] + t33*z[1],
    f7 = t11*x[2] + t12*y[2] + t13*z[2],
    f8 = t12*x[2] + t22*y[2] + t23*z[2],
    f9 = t13*x[2] + t23*y[2] + t33*z[2];
  t11 = x[0]*f1 + y[0]*f2 + z[0]*f3;
  t12 = x[1]*f1 + y[1]*f2 + z[1]*f3;
  t13 = x[2]*f1 + y[2]*f2 + z[2]*f3;
  t22 = x[1]*f4 + y[1]*f5 + z[1]*f6;
  t23 = x[2]*f4 + y[2]*f5 + z[2]*f6;
  t33 = x[2]*f7 + y[2]*f8 + z[2]*f9;
}

template <typename T> inline void
rotate_sym_tensor_3x3_RtAR (const T* const x, const T* const y, const T* const z,
                            T* const t) {
  rotate_sym_tensor_3x3_RtAR(x, y, z, t[0], t[1], t[2], t[3], t[4], t[5]);
}

// Row-major R.
template <typename T> inline void
rotate_sym_tensor_3x3_RtAR (const T* const R, T* const t) {
  rotate_sym_tensor_3x3_RtAR(R, R+3, R+6, t[0], t[1], t[2], t[3], t[4], t[5]);
}

// For row-major R = [x; y; z] and T = [t11 t12 t13; t12 t22 t23; t13 t23 t33],
// return R T R'.
template <typename T> inline void
rotate_sym_tensor_3x3_RARt (const T* const x, const T* const y, const T* const z,
                            T& t11, T& t12, T& t13,
                            T& t22, T& t23, T& t33) {
  const T
    f1 = t11*x[0] + t12*x[1] + t13*x[2],
    f2 = t12*x[0] + t22*x[1] + t23*x[2],
    f3 = t13*x[0] + t23*x[1] + t33*x[2],
    f4 = t11*y[0] + t12*y[1] + t13*y[2],
    f5 = t12*y[0] + t22*y[1] + t23*y[2],
    f6 = t13*y[0] + t23*y[1] + t33*y[2],
    f7 = t11*z[0] + t12*z[1] + t13*z[2],
    f8 = t12*z[0] + t22*z[1] + t23*z[2],
    f9 = t13*z[0] + t23*z[1] + t33*z[2];
  t11 = x[0]*f1 + x[1]*f2 + x[2]*f3;
  t12 = y[0]*f1 + y[1]*f2 + y[2]*f3;
  t13 = z[0]*f1 + z[1]*f2 + z[2]*f3;
  t22 = y[0]*f4 + y[1]*f5 + y[2]*f6;
  t23 = z[0]*f4 + z[1]*f5 + z[2]*f6;
  t33 = z[0]*f7 + z[1]*f8 + z[2]*f9;
}

template <typename T> inline void
rotate_sym_tensor_3x3_RARt (const T* const x, const T* const y, const T* const z,
                            T* const t) {
  rotate_sym_tensor_3x3_RARt(x, y, z, t[0], t[1], t[2], t[3], t[4], t[5]);
}

template <typename T> inline void
rotate_sym_tensor_3x3_RARt (const T* const R, T* const t) {
  rotate_sym_tensor_3x3_RARt(R, R+3, R+6, t[0], t[1], t[2], t[3], t[4], t[5]);
}

// dx[2] = {x0-x1, x2-x1} for x0 < x1 = 0 < x2.
// y[] = {f(x0), f(x1=0), f(x2)}.
void calc_parabola_coefs(const Real dx[2], const Real y[3], Real c[3]);

// cos(x) - 1
Real cosxm1(const Real x);
// sqrt(1 + x) - 1
Real sqrt1pxm1(const Real x);

// Solve a x^2 + b x + c = 0. Return false if the roots are complex, true
// otherwise. If assume_real, then if the discriminant is < 0, set it to 0. This
// permits roundoff-level negative discriminants if the caller knows the
// discriminant is in fact 0 in that case.
bool solve_quadratic_eq_for_real_roots(const Real a, const Real b, const Real c,
                                       Real x[2],
                                       const bool assume_real = false);

// v = [x; y; z] u.
void matvec (const Real x[3], const Real y[3], const Real z[3],
             const Real u[3], Real v[3]);
// v = [x; y; z]' u.
void tmatvec (const Real x[3], const Real y[3], const Real z[3],
              const Real u[3], Real v[3]);

// Form v so that |dot(u,v)|/(norm(u) norm(v)) <= sqrt(1/3).
void form_separated_vector(const Real u[3], Real v[3]);

// Form row-major R such that R(1,:) = normalize(x).
void form_rotation_given_x(const Real x[3], Real R[9]);

// Form row-major R such that R(1,:) = normalize(x) and R(2,:)'z = 0.
void form_rotation_given_x_then_z(const Real x[3], const Real z[3], Real R[9]);

// A = [a, c; c b].
void sym_2x2_eig(const Real a, const Real b, const Real c,
                 Real lams[2], Real v1[2], Real v2[2]);

std::vector<std::string> split(std::string s, const std::string& delim);

template <typename T> int nbytes () { return 0; }
template <> inline int nbytes<float > () { return 4; }
template <> inline int nbytes<double> () { return 8; }

struct File {
  typedef std::shared_ptr<File> Ptr;
  File (FILE* fid_ = nullptr) : fid(fid_) {}
  ~File () { if (fid) std::fclose(fid); }
  FILE* id () { return fid; }
  // fid == nullptr is permitted.
  static File::Ptr create(FILE* fid);
  // If fopen fails, nullptr is returned. On error and !quiet, print a message.
  static File::Ptr create(const std::string& filename, const std::string& mode,
                          const bool quiet=false);
private:
  FILE* fid;
};

bool file_exists(const std::string& filename);

struct Sprinter {
  typedef std::shared_ptr<Sprinter> Ptr;
  
  // Write.
  void wr(const char* format, ...);

  // Use this class in one of two ways. Mixing them between reset calls is an
  // error. Reset to switch APIs. At rest, the buffer is cleared.

  // API 1. Gather one large string, then output it or return it.
  Sprinter(const int nbuf_ = 128);
  void reset(const int nbuf_ = 128);
  void out(FILE* stream = stdout, const bool newline = true) const;
  std::string str() const;
  void clear();

  // API 2. Buffer to a string if bufsz > 0; periodically write to the stream.
  Sprinter(const File::Ptr& stream, const int bufsz = -1);
  void reset(const File::Ptr& stream, const int bufsz = -1);
  // The file is not assured to be complete until flush is called. reset and
  // ~Sprinter call flush if needed.
  void flush();

  ~Sprinter();

private:
  bool str_api;
  int n, nbuf;
  std::vector<char> buf;
  int bufsz;
  File::Ptr fid;

  void clear_();
};

struct SilenceThrowPrint {
  SilenceThrowPrint();
  ~SilenceThrowPrint();
};

int util_test();

} // namespace acorn
} // namespace woodland

#endif
