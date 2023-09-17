#include "woodland/acorn/hfp.hpp"

#include <cassert>
#include <cstdio>
#include <limits>

#include "woodland/acorn/compose_quadrature.hpp"
#include "woodland/acorn/compose_lagrange_polynomial.hpp"

#include "woodland/acorn/util.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace acorn {
namespace hfp {

CallerPolyP::CallerPolyP () : n(0), xs(nullptr), ys(nullptr) {}

CallerPolyP::CallerPolyP (const int n_, CRPtr xs_, CRPtr ys_)
  : n(n_), xs(xs_), ys(ys_)
{
  assert(n > 0);
}

Real CallerPolyP::eval (const Real x) const {
  assert(n > 0);
  return eval_lagrange_poly(n, xs, ys, x);
}

void CallerPolyP::eval_p_c (const Real c, Real& pc, Real& pdc) const {
  pc  = eval_lagrange_poly           (n, xs, ys, c);
  pdc = eval_lagrange_poly_derivative(n, xs, ys, c);  
}

bool CallerPolyP::in_support (const Real x) const {
  return n > 1 && x >= xs[0] && x <= xs[n-1];
}

static Real
integrate_third_term (const Options& o, const CallerP& p,
                      const Real pc, const Real pdc,
                      const Real a, const Real b, const Real c) {
  assert(a <= b);
  assert(is_gl_supported(o.gl_np));

  if (a == b) return 0;

  if (a < c && c < b)
    return (integrate_third_term(o, p, a, c, c) +
            integrate_third_term(o, p, c, b, c));

  // At this point, one of a or b is the singular point or the singular point is
  // outside of [a,b].
  const auto qnp = o.gl_np;
  const auto* const qxs = get_x_gl(qnp);
  const auto* const qws = get_w_gl(qnp);
  Real integral = 0;
  for (int i = 0; i < qnp; ++i) {
    const auto xepc = a + (b - a)*((1 + qxs[i])/2);
    const auto xe = xepc - c;
    assert(xe != 0);
    assert(p.in_support(xepc));
    integral += qws[i]*((p.eval(xepc) - pc - pdc*xe) / square(xe));
  }
  return integral*((b - a)/2);
}

Real calc_hfp (const Options& o, const CallerP& p,
               const Real a, const Real b, const Real c) {
  assert(a <= b);
  Real pc, pdc;
  p.eval_p_c(c, pc, pdc);
  const auto third_term = integrate_third_term(o, p, pc, pdc, a, b, c);
  if (c > b || c < a) {                   // proper integral
    return (pc*(1/(a - c) - 1/(b - c)) +
            pdc*std::log((b - c)/(a - c)) +
            third_term);
  } else if (a == 0 && c == 0 && b > a) { // case 2
    return (-pc/b +
            pdc*std::log(b) +
            third_term);
  } else {                                // case 1
    return (pc*(1/(a - c) - 1/(b - c)) +
            pdc*std::log((b - c)/(c - a)) +
            third_term);
  }
}

Real integrate_third_term (const Options& o, const CallerP& p,
                           const Real a, const Real b, const Real c) {
  assert(a <= b);
  Real pc, pdc;
  p.eval_p_c(c, pc, pdc);
  return integrate_third_term(o, p, pc, pdc, a, b, c);
}

static void setup_problem (const int nd, CRPtr ds,
                           const Real a, const Real b, const Real c,
                           const int n, RPtr xs, RPtr ys) {
  assert(a <= b);
  const Real lo = std::min(a, c), hi = std::max(b, c);
  const Real* qxs = get_x_gll(n);
  for (int i = 0; i < n; ++i) {
    const auto x = lo + (hi - lo)*(1 + qxs[i])/2;
    xs[i] = x;
    Real p = 0;
    for (int i = nd-1; i >= 0; --i) p = ds[i] + (x-c)*p;
    ys[i] = p;
  }  
}

static Real calc_hfp (const Options& o, CRPtr xs, CRPtr ys, const int n,
                      const Real a, const Real b, const Real c) {
  CallerPolyP p(n, xs, ys);
  return calc_hfp(o, p, a, b, c);
}

/* Test problem:
       p(x) = sum_{i=0}^{n-1} di (x-c)^i.
   Then
       p (c) = d0
       p'(c) = d1
       p(x+c) = sum_{i=0}^{n-1} di x^i
       (p(x+c) - p(c) - p'(c) x)/x^2 = sum_{i=2}^{n-1} di x^(i-2)
   and
       int_(a-c)^(b-c) (p(x+c) - p(c) - p'(c) x)/x^2 dx
         = int_(a-c)^(b-c) sum_{i=2}^{n-1} di x^(i-2) dx
         = sum_{i=2}^{n-1} 1/(i-1) di x^{i-1} |_(a-c)^(b-c).
 */
static int
test_integrate_third_term (const int nd, CRPtr ds,
                           const Real a, const Real b, const Real c) {
  int ne = 0;
  Real int_num = 0;
  {
    static const int n = 20;
    Real xs[n], ys[n];
    setup_problem(nd, ds, a, b, c, n, xs, ys);
    Options o;
    CallerPolyP p(n, xs, ys);
    int_num = integrate_third_term(o, p, a, b, c); 
  }
  Real int_ana_a = 0, int_ana_b = 0;
  for (int i = nd-1; i >= 2; --i) {
    const auto f = ds[i]/(i-1);
    int_ana_a = (a - c)*(f + int_ana_a);
    int_ana_b = (b - c)*(f + int_ana_b);
  }
  const auto int_ana = int_ana_b - int_ana_a;
  const auto re = (std::abs(int_ana - int_num)/
                   std::max(std::abs(int_ana), std::abs(int_num)));
  if (re > 1e4*std::numeric_limits<Real>::epsilon()) {
    ++ne;
    printf("test_integrate_third_term:\n  ds (%d):", nd);
    for (int i = 0; i < nd; ++i) printf(" %10.3e", ds[i]);
    printf("\n  int_ana %22.15e\n  int_num %22.15e\n  re      %10.3e\n",
           int_ana, int_num, re);
  }
  return ne;
}

int test_calc_hfp () {
  const auto eps = std::numeric_limits<Real>::epsilon();
  int ne = 0;
  static const int n = 20;
  Real xs[n], ys[n];
  Options o;
  { // hfp int_-1^1 1/x^2 dx = -2
    const Real ds[] = {1};
    const Real a = -1, b = 1, c = 0;
    setup_problem(sizeof(ds)/sizeof(*ds), ds, a, b, c, n, xs, ys);
    const auto hfp = calc_hfp(o, xs, ys, n, a, b, c);
    if (std::abs(hfp + 2) > 1e4*eps) {
      printf("test_calc_hfp case 1: hfp %22.15e re %10.3e\n",
             hfp, std::abs(hfp + 2)/2);
      ++ne;
    }
  }
  { // cpv int_0^2 1/(x - 1) dx = 0
    const Real ds[] = {0, 1};
    const Real a = 0, b = 2, c = 1;
    setup_problem(sizeof(ds)/sizeof(*ds), ds, a, b, c, n, xs, ys);
    const auto hfp = calc_hfp(o, xs, ys, n, a, b, c);
    if (std::abs(hfp) > 1e4*eps) {
      printf("test_calc_hfp case 1: hfp %22.15e err %10.3e\n",
             hfp, std::abs(hfp));
      ++ne;
    }
  }
  { // cpv int_2^3 1/(x - 1) dx = log 2
    const Real ds[] = {0, 1};
    const Real a = 2, b = 3, c = 1;
    setup_problem(sizeof(ds)/sizeof(*ds), ds, a, b, c, n, xs, ys);
    const auto igal = calc_hfp(o, xs, ys, n, a, b, c);
    const auto igal_true = std::log(2);
    if (std::abs(igal - igal_true) > 0) {
      printf("test_calc_hfp case 1: igal %22.15e err %10.3e\n",
             igal, std::abs(igal - igal_true)/igal_true);
      ++ne;
    }
  }
  { // cpv int_0^3 1/(x - 1) dx = int_1^2 1/x dx = log(2)
    const Real ds[] = {0, 1};
    const Real a = 0, b = 3, c = 1;
    setup_problem(sizeof(ds)/sizeof(*ds), ds, a, b, c, n, xs, ys);
    const auto hfp = calc_hfp(o, xs, ys, n, a, b, c);
    const auto re = std::abs(hfp - std::log(2))/std::log(2);
    if (re > 1e4*eps) {
      printf("test_calc_hfp case 1: hfp %22.15e re %10.3e\n", hfp, re);
      ++ne;
    }
  }
  { // hfp int_0^1 1/x^2 dx = -1
    const Real ds[] = {1};
    const Real a = 0, b = 1, c = a;
    setup_problem(sizeof(ds)/sizeof(*ds), ds, a, b, c, n, xs, ys);
    const auto hfp = calc_hfp(o, xs, ys, n, a, b, c);
    if (std::abs(hfp + 1) > 1e4*eps) {
      printf("test_calc_hfp case 2: hfp %22.15e re %10.3e\n",
             hfp, std::abs(hfp + 1));
      ++ne;
    }
  }
  { // int_1^2 1/x^2 dx = 1/2
    const Real igal_true = 0.5;
    const Real ds[] = {1};
    const Real a = 1, b = 2, c = 0;
    setup_problem(sizeof(ds)/sizeof(*ds), ds, a, b, c, n, xs, ys);
    const auto igal = calc_hfp(o, xs, ys, n, a, b, c);
    if (std::abs(igal - igal_true) > 0) {
      printf("test_calc_hfp case 2: hfp %22.15e re %10.3e\n",
             igal, std::abs(igal - igal_true)/igal_true);
      ++ne;
    }
  }
  { // hfp int_0^1 1/x dx = 0
    const Real ds[] = {0, 1};
    const Real a = 0, b = 1, c = a;
    setup_problem(sizeof(ds)/sizeof(*ds), ds, a, b, c, n, xs, ys);
    const auto hfp = calc_hfp(o, xs, ys, n, a, b, c);
    if (std::abs(hfp) > 1e4*eps) {
      printf("test_calc_hfp case 2: hfp %22.15e err %10.3e\n",
             hfp, std::abs(hfp));
      ++ne;
    }
  }
  return ne;
}

int unittest () {
  int ne = 0;
  {
    const Real ds[] = {1, -1, 0.5, 2.1};
    ne += test_integrate_third_term(sizeof(ds)/sizeof(*ds), ds, -1, 1, 0);
    ne += test_integrate_third_term(sizeof(ds)/sizeof(*ds), ds, 0, 2, 1);
  }
  {
    const Real ds[] = {0.3, -0.4, 0.2, -1.3, 0.1, 0.05, -0.04};
    ne += test_integrate_third_term(sizeof(ds)/sizeof(*ds), ds, -1.1, 2.3, 0.5);
  }
  ne += test_calc_hfp();
  return ne;
}

} // namespace hfp
} // namespace acorn
} // namespace woodland
