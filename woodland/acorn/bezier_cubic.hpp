#ifndef INCLUDE_WOODLAND_ACORN_BEZIER_CUBIC
#define INCLUDE_WOODLAND_ACORN_BEZIER_CUBIC

#include "woodland/acorn/compose_quadrature.hpp"

#include "woodland/acorn/macros.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/solver1d.hpp"

#include <cmath>
#include <limits>

namespace woodland {
namespace acorn {

/* Cubic Bezier curve.
     Arrays of vectors are 2xn, with the 2D vectors the fast index. Array names
   end with 's' for plural.
     Arrays of coefficients are 4xn.
     As in all other docs, () indexing is 1-based and [] indexing is 0-based.
     Each curve p(t) is parameterized by t in [0,1], with p(0) = p1 and p(1) = p2.
     Functions over n values use ompfor. Functions over just m do not.
 */
template <typename Real>
struct BezierCubic {
  typedef const Real* const CRPtr;
  typedef Real* const RPtr;
  typedef Real Vec2[2];
  static constexpr Real eps = std::numeric_limits<Real>::epsilon();

private:
  using v2 = Matvec2d<Real>;

public:
  // Initialize the coefficients cs(1:8) from p1,n1,p2,n2s(1:2) and optional d.
  //   p1,2 are the curve end points. The curve unit normals at p1,2 are n1,2.
  //   d = norm(p2 - p1, 2).
  //   The coefficients are chosen so that the arc covers <(=, an edge case)
  // half a circle.
  static void init_from_nml(CRPtr p1, CRPtr n1, CRPtr p2, CRPtr n2, RPtr c,
                            Real d = -1);

  // p,ns(1:2,1:n+1) define n adjacent curves.
  static void init_from_nml(const int n, CRPtr ps, CRPtr ns, RPtr cs,
                            CRPtr ds = nullptr);

  // Evaluate p(t).
  static void eval_p  (CRPtr c, const Real& t, RPtr p);
  static void eval_p  (CRPtr c, const int m, CRPtr ts, RPtr ps);
  // Evaluate p_t(t).
  static void eval_pt (CRPtr c, const Real& t, Vec2 pt);
  static void eval_pt (CRPtr c, const int m, CRPtr ts, RPtr ps);
  // Evaluate p_tt(t).
  static void eval_ptt(CRPtr c, const Real& t, Vec2 ptt);
  static void eval_ptt(CRPtr c, const int m, CRPtr ts, RPtr ps);
  // Evaluate ps(1:2,blk) = p(ts(1:m)) for n m-length blocks. Note the same
  // ts(1:m) is used for each of the n blocks.
  static void eval_p_one_ts  (const int n, CRPtr cs, const int m, CRPtr ts, RPtr ps);
  static void eval_pt_one_ts (const int n, CRPtr cs, const int m, CRPtr ts, RPtr ps);
  static void eval_ptt_one_ts(const int n, CRPtr cs, const int m, CRPtr ts, RPtr ps);

  // Compute R, the signed radius of curvature. The curvature kappa = 1/R. See
  //     https://en.wikipedia.org/wiki/Curvature#In_terms_of_a_general_parametrization
  static Real calc_radius_of_curvature(CRPtr c, const Real& t);
  
  // Compute arclength(t1,t2) and optionally arclength_t2(t1,t2).
  static void calc_arclength(const Quadrature& q, CRPtr c,
                             const Real& t1, const Real& t2,
                             Real& al, Real* al_d = nullptr);

  // Find tc in [0, 1] such that arclength(0,tc) = arclength(0,1)/2.
  static Real solve_for_ctr_t(const Quadrature& q, CRPtr c,
                              const Real tol = 1e4*eps);

  // Find t in [0, 1] such that for p = eval_p(c, t, p), p[0] = x. x must be in
  // [c[0], c[6]].
  static Real calc_t_for_x(CRPtr c, const Real x, const Real tol = 1e4*eps);

  static int unittest();
};

} // namespace acorn
} // namespace woodland

#endif
