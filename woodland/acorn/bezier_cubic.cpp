#include "woodland/acorn/bezier_cubic.hpp"
#include "woodland/acorn/util.hpp"

namespace woodland {
namespace acorn {

template <typename Real>
void BezierCubic<Real>
::init_from_nml (CRPtr p1, CRPtr n1, CRPtr p2, CRPtr n2, RPtr c, Real d) {
  /* Derivation:
     given end points p1, p2 and end-point normal unit vectors n1, n2.
     let c1, c2 be the two internal control points.
     let r = 1-t.
     let perpcw(n) = [n(2); -n(1)].
     for t in [0,1],
     p(t) = p1 r^3 + c1 r^2 t + c2 r t^2 + p2 t^3.
     p_t(t) = -3 p1 r^2 + r (r - 2 t) c1 + t (2 r - t) c2 + 3 p2 t^2
     p_t(0) = d perpcw(n1) = -3 p1 + c1 => c1 = 3 p1 + d perpcw(n1)
     p_t(1) = d perpcw(n2) =  3 p2 - c2 => c2 = 3 p2 - d perpcw(n2)
     rule: choose path from p1 to p2 that is the shorter of the two
     possible. this leads to inserting the following sign factor:
     sgn = sign(dot(p_t(0), p2-p1))
     c1 = 3 p1 + sgn d perpcw(n1)
     c2 = 3 p2 - sgn d perpcw(n2)
  */
  assert(std::abs(v2::norm22(n1) - 1) < 1e2*eps);
  assert(std::abs(v2::norm22(n2) - 1) < 1e2*eps);
  Vec2 pd1; v2::perpcw2(n1, pd1);
  Vec2 pd2; v2::perpcw2(n2, pd2);
  Vec2 p2m1; v2::subtract(p2, p1, p2m1);
  const auto sgn = v2::dot(pd1, p2m1) >= 0 ? 1 : -1;
  if (d < 0) d = std::sqrt(v2::norm22(p2m1));
  else assert(std::abs(d - std::sqrt(v2::norm22(p2m1))) < 1e2*d*eps);
  v2::copy(p1, c  );
  v2::copy(p2, c+6);
  v2::axpbyz(3, p1,  sgn*d, pd1, c+2);
  v2::axpbyz(3, p2, -sgn*d, pd2, c+4);
}

template <typename Real>
void BezierCubic<Real>
::init_from_nml (const int n, CRPtr ps, CRPtr ns, RPtr cs, CRPtr ds) {
  ompfor for (int i = 0; i < n; ++i)
    init_from_nml(&ps[2*i], &ns[2*i], &ps[2*(i+1)], &ns[2*(i+1)], &cs[8*i],
                  ds ? ds[i] : -1);
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
  {
    Real p1[] = {0,0}, p2[] = {1,0}, n1[] = {0,1}, n2[] = {0,1}, c[8];
    bc.init_from_nml(p1, n1, p2, n2, c);
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
    Real p1[] = {-0.1,-0.5}, p2[] = {1,0.3}, n1[] = {n10,1}, n2[] = {0.5,1}, c[8];
    Matvec2d<Real>::normalize(n1);
    Matvec2d<Real>::normalize(n2);
    bc.init_from_nml(p1, n1, p2, n2, c);
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
    typedef Matvec2d<Real> v2;
    const Real ths[] = {0, 0.01*M_PI}, fac = 3.2;
    Real ps[2][2], ns[2][2], c[8];
    for (int i = 0; i < 2; ++i) {
      ps[i][0] = fac*std::cos(ths[i]); ps[i][1] = fac*std::sin(ths[i]);
      v2::copy(ps[i], ns[i]);
      v2::normalize(ns[i]);
    }
    bc.init_from_nml(ps[0], ns[0], ps[1], ns[1], c);
    for (const Real t : {0.1, 0.5, 0.8}) {
      const auto R = bc.calc_radius_of_curvature(c, t);
      if (std::abs(R - fac) > 1e-3) ++nerr;
    }
  }
  return nerr;
}

template struct BezierCubic<double>;

} // namespace acorn
} // namespace woodland
