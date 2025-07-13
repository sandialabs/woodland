#include "woodland/squirrel/bezier_cubic.hpp"

#include "woodland/acorn/util.hpp"

#include "woodland/squirrel/solver1d.hpp"

namespace woodland {
namespace squirrel {

using acorn::square;
using acorn::cube;
using acorn::Quadrature;

template <typename Real>
void BezierCubic<Real>
::init_from_tan (CRPtr p1, CRPtr m1, CRPtr p2, CRPtr m2, const Real d, RPtr c) {
  /* Derivation:
       Given end points p1, p2 and end-point tangent vectors m1, m2.
       Let c1, c2 be the two internal control points.
       Let r = 1-t.
       For t in [0,1],
         p(t) = p1 r^3 + c1 r^2 t + c2 r t^2 + p2 t^3.
         p_t(t) = -3 p1 r^2 + r (r - 2 t) c1 + t (2 r - t) c2 + 3 p2 t^2
         p_t(0) = d m1 = -3 p1 + c1 => c1 = 3 p1 + d m1
         p_t(1) = d m2 =  3 p2 - c2 => c2 = 3 p2 - d m2.
   */
  v2::copy(p1, c  );
  v2::copy(p2, c+6);
  v2::axpbyz(3, p1,  d, m1, c+2);
  v2::axpbyz(3, p2, -d, m2, c+4);
}

template <typename Real>
void BezierCubic<Real>
::init_from_tan (const int n, CRPtr ps, CRPtr ms, CRPtr ds, RPtr cs) {
  ompfor for (int i = 0; i < n; ++i)
    init_from_tan(&ps[2*i], &ms[2*i], &ps[2*(i+1)], &ms[2*(i+1)], ds[i],
                  &cs[8*i]);
}

template <typename Real>
void BezierCubic<Real>
::eval_p (CRPtr c, const Real& t, RPtr p) {
  const auto t2 = t*t, t3 = t2*t, r = 1-t, r2 = r*r, r3 = r2*r;
  for (int d = 0; d < 2; ++d)
    p[d] = c[d]*r3 + c[2+d]*r2*t + c[4+d]*r*t2 + c[6+d]*t3;
}

template <typename Real>
void BezierCubic<Real>
::eval_pt (CRPtr c, const Real& t, Vec2 pt) {
  const auto t2 = t*t, t2p = 2*t, t3p = 3*t2;
  const auto r = 1-t, r2 = r*r, r2p = -2*r, r3p = -3*r2;
  for (int d = 0; d < 2; ++d)
    pt[d] = c[d]*r3p + c[2+d]*(r2p*t + r2) + c[4+d]*(r*t2p - t2) + c[6+d]*t3p;
}

template <typename Real>
void BezierCubic<Real>
::eval_ptt (CRPtr c, const Real& t, Vec2 ptt) {
  const auto t2p = 2*t, t2pp = 2, t3pp = 6*t;
  const auto r = 1-t, r2p = -2*r, r2pp = 2, r3pp = 6*r;
  for (int d = 0; d < 2; ++d)
    ptt[d] = (c[d]*r3pp + c[2+d]*(r2pp*t + 2*r2p) + c[4+d]*(r*t2pp - 2*t2p) +
              c[6+d]*t3pp);
}

template <typename Real>
void BezierCubic<Real>
::eval_p (CRPtr c, const int m, CRPtr ts, RPtr ps) {
  for (int j = 0; j < m; ++j)
    eval_p(c, ts[j], &ps[2*j]);
}

template <typename Real>
void BezierCubic<Real>
::eval_pt (CRPtr c, const int m, CRPtr ts, RPtr ps) {
  for (int j = 0; j < m; ++j)
    eval_pt(c, ts[j], &ps[2*j]);
}

template <typename Real>
void BezierCubic<Real>
::eval_ptt (CRPtr c, const int m, CRPtr ts, RPtr ps) {
  for (int j = 0; j < m; ++j)
    eval_ptt(c, ts[j], &ps[2*j]);
}

template <typename Real>
void BezierCubic<Real>
::eval_p_one_ts (const int n, CRPtr cs, const int m, CRPtr ts, RPtr ps) {
  ompfor for (int i = 0; i < n; ++i)
    eval_p(&cs[8*i], m, ts, &ps[2*m*i]);
}

template <typename Real>
void BezierCubic<Real>
::eval_pt_one_ts (const int n, CRPtr cs, const int m, CRPtr ts, RPtr ps) {
  ompfor for (int i = 0; i < n; ++i)
    eval_pt(&cs[8*i], m, ts, &ps[2*m*i]);
}

template <typename Real>
void BezierCubic<Real>
::eval_ptt_one_ts (const int n, CRPtr cs, const int m, CRPtr ts, RPtr ps) {
  ompfor for (int i = 0; i < n; ++i)
    eval_ptt(&cs[8*i], m, ts, &ps[2*m*i]);
}

template <typename Real>
Real BezierCubic<Real>
::calc_radius_of_curvature (CRPtr c, const Real& t) {
  Vec2 pt, ptt;
  eval_pt (c, t, pt );
  eval_ptt(c, t, ptt);
  return (std::sqrt(cube(square(pt[0]) + square(pt[1])))/
          (pt[0]*ptt[1] - pt[1]*ptt[0]));
}

template <typename Real>
void BezierCubic<Real>
::calc_arclength (const Quadrature& q, CRPtr c, const Real& t1, const Real& t2,
                  Real& al_, Real* al_d_) {
  const Real* gx, * gw;
  q.get_xw(gx, gw);
  Real al = 0, al_d = 0;
  const bool deriv = al_d_;
  for (int iq = 0; iq < q.nq; ++iq) {
    const auto a = (1 + gx[iq])/2;
    assert(a > -1e-6 and a < 1+1e-6);
    const auto t = (1-a)*t1 + a*t2;
    Vec2 pt;
    eval_pt(c, t, pt);
    const auto f = std::sqrt(square(pt[0]) + square(pt[1]));
    al += f*gw[iq];
    if (deriv) {
      Vec2 ptt;
      eval_ptt(c, t, ptt);
      al_d += (f + (pt[0]*ptt[0] + pt[1]*ptt[1])*(t2-t1)*a/f)*gw[iq]/2;
    }
  }
  al *= (t2 - t1)/2;
  al_ = al;
  if (deriv) *al_d_ = al_d;
}

template <typename Real>
struct CtrFn : Solver1d<Real>::Fn {
  typedef typename BezierCubic<Real>::CRPtr CRPtr;

  Quadrature q;
  CRPtr c;
  Real al;

  CtrFn (const Quadrature& q_, CRPtr c_)
    : q(q_), c(c_)
  { BezierCubic<Real>::calc_arclength(q, c, 0, 1, al); }

  void eval (const Real& t, Real& f, Real& fp) const override {
    BezierCubic<Real>::calc_arclength(q, c, 0, t, f, &fp);
    f -= al/2;
  }
};

template <typename Real>
Real BezierCubic<Real>
::solve_for_ctr_t (const Quadrature& q, CRPtr c, const Real tol) {
  typedef Solver1d<Real> S;
  CtrFn<Real> fn(q, c);
  typename S::Info info;
  typename S::Tols tols(tol, tol);
  S::fzero(tols, fn, 0, 1, -fn.al/2, fn.al - fn.al/2, info);
  assert(info.success());
  return info.x;
}

template <typename Real>
struct XFn : Solver1d<Real>::Fn {
  typedef typename BezierCubic<Real>::CRPtr CRPtr;

  CRPtr c;
  Real x;

  XFn (CRPtr c_, const Real x_)
    : c(c_), x(x_)
  {}

  void eval (const Real& t, Real& f, Real& fp) const override {
    Real p[2];
    BezierCubic<Real>::eval_p(c, t, p);
    f = p[0] - x;
    BezierCubic<Real>::eval_pt(c, t, p);
    fp = p[0];
  }
};

template <typename Real>
Real BezierCubic<Real>::calc_t_for_x (CRPtr c, const Real x, const Real tol) {
  assert(x >= c[0] && x <= c[6]);
  if (x == c[0]) return 0;
  if (x == c[6]) return 1;
  typedef Solver1d<Real> S;
  XFn<Real> fn(c, x);
  typename S::Info info;
  typename S::Tols tols(tol, tol);
  S::fzero(tols, fn, 0, 1, c[0] - x, c[6] - x, info);
  assert(info.success());
  return info.x;  
}

template <typename Real>
int BezierCubic<Real>::unittest () {
  static constexpr Real tol = 1e4*std::numeric_limits<Real>::epsilon();
  int nerr = 0;
  BezierCubic<Real> bc;
  Quadrature q;
  for (int tno = 0; tno < 1; ++tno){
    const Real p1[] = {0,0}, p2[] = {1,0}, m1[] = {1,0}, m2 [] = {1,0};
    Real c[8];
    bc.init_from_tan(p1, m1, p2, m2, p2[0] - p1[0], c);
    Real t = 0.5, p[2], pt[2], ptt[2];
    bc.eval_p  (c, t, p  );
    bc.eval_pt (c, t, pt );
    bc.eval_ptt(c, t, ptt);
    if (std::abs(p  [0] - 0.5) > tol) ++nerr;
    if (std::abs(p  [1])       > tol) ++nerr;
    if (std::abs(pt [0] - 1)   > tol) ++nerr;
    if (std::abs(pt [1])       > tol) ++nerr;
    if (std::abs(ptt[0])       > tol) ++nerr;
    if (std::abs(ptt[1])       > tol) ++nerr;
    const auto tc = bc.solve_for_ctr_t(q, c);
    if (std::abs(tc - 0.5) > tol) ++nerr;
    {
      const Real x = 0.7;
      const auto tx = bc.calc_t_for_x(c, x);
      Vec2 px;
      bc.eval_p(c, tx, px);
      if (std::abs(px[0] - x) > tol) ++nerr;
    }
  }
  for (const Real n10 : {0.0, -0.1, 0.1, -0.3, -1.5}) {
    Real p1[] = {-0.1,-0.5}, p2[] = {1,0.3}, m1[] = {1,-n10}, m2[] = {0.1,5}, c[8];
    bc.init_from_tan(p1, m1, p2, m2, p2[0] - p1[0], c);
    const auto tc = bc.solve_for_ctr_t(q, c, tol);
    Real al, alh;
    bc.calc_arclength(q, c, 0, 1, al);
    bc.calc_arclength(q, c, 0, tc, alh);
    if (std::abs(alh - al/2) > tol) ++nerr;
    {
      const Real x = 0.3*p1[0] + 0.7*p2[0];
      const auto tx = bc.calc_t_for_x(c, x);
      Vec2 px;
      bc.eval_p(c, tx, px);
      if (std::abs(px[0] - x) > tol) ++nerr;
    }
  }
  {
    typedef acorn::Matvec2d<Real> v2;
    const Real ths[] = {0, 0.01*M_PI}, fac = 3.2;
    Real ps[2][2], ns[2][2], ms[2][2];
    for (int i = 0; i < 2; ++i) {
      ps[i][0] = fac*std::cos(ths[i]); ps[i][1] = fac*std::sin(ths[i]);
      v2::copy(ps[i], ns[i]);
      v2::normalize(ns[i]);
      ms[i][0] = -fac*std::sin(ths[i]); ms[i][1] = fac*std::cos(ths[i]);
    }
    for (int tno = 0; tno < 1; ++tno) {
      Real c[8];
      bc.init_from_tan(ps[0], ms[0], ps[1], ms[1], ths[1]-ths[0], c);
      for (const Real t : {0.1, 0.5, 0.8}) {
        const auto R = bc.calc_radius_of_curvature(c, t);
        if (std::abs(R - fac) > 1e-3) ++nerr;
      }
    }
  }
  return nerr;
}

template struct BezierCubic<double>;

} // namespace squirrel
} // namespace woodland
