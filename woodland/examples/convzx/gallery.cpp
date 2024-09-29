#include "woodland/examples/convzx/gallery.hpp"

#include "woodland/acorn/compose_quadrature.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/solver1d.hpp"
#include "woodland/acorn/dbg.hpp"

#include <cmath>

namespace woodland {
namespace examples {
namespace convzx {
namespace gallery {

typedef acorn::Matvec<2,Real> mv2;

ZxFn::Shape ZxFn::convert (const std::string& shape) {
  if (shape == "ramp") return Shape::ramp;
  if (shape == "quadratic") return Shape::quadratic;
  if (shape == "trig0") return Shape::trig0;
  if (shape == "trig1") return Shape::trig1;
  if (shape == "steep") return Shape::steep;
  if (shape == "zero") return Shape::zero;
  fprintf(stderr, "ZxFn::convert: invalid shape %s\n", shape.c_str());
  return Shape::zero;
}

std::string ZxFn::convert (const Shape shape) {
  switch (shape) {
  case Shape::ramp: return "ramp";
  case Shape::quadratic: return "quadratic";
  case Shape::trig0: return "trig0";
  case Shape::trig1: return "trig1";
  case Shape::steep: return "steep";
  case Shape::zero:
  default:
    return "zero";
  }  
}

void eval (const ZxFn::Shape shape, const Real x, Real& f, Real& g) {
  switch (shape) {
  case ZxFn::Shape::ramp: {
    const Real a = 0.1;
    f = a*x - 0.11;
    g = a;
  } break;
  case ZxFn::Shape::quadratic: {
    const Real a = 1.2;
    f = a*x*(1 - x) - 0.33;
    g = a*(1 - 2*x);
  } break;
  case ZxFn::Shape::trig0: {
    const Real a = 0.15, b = 2*M_PI;
    f = a*(1 + std::cos(b*(x - 0.5))) - 0.33;
    g = -a*b*std::sin(b*(x - 0.5));
  } break;
  case ZxFn::Shape::trig1:
  case ZxFn::Shape::steep: {
    const Real a = shape == ZxFn::Shape::steep ? 0.9 : 0.3, b = 2*M_PI;
    f = a*(x*std::sin(b*x) - 0.4);
    g = a*(std::sin(b*x) + x*b*std::cos(b*x));
  } break;
  case ZxFn::Shape::zero:
  default:
    f = -0.1; g = 0; break;
  }
}

static Real calc_arclength (const ZxFn::Shape shape, const Real x) {
  acorn::Quadrature q(20, acorn::Quadrature::gl);
  const Real* qx, * qw;
  q.get_xw(qx, qw);
  Real al = 0;
  for (int i = 0; i < q.nq; ++i) {
    Real f, g;
    eval(shape, x*(qx[i] + 1)/2, f, g);
    al += qw[i]*std::sqrt(1 + g*g);
  }
  return (al/2)*x;
}

Real ZxFn::calc_arclength (const Real x) const {
  return gallery::calc_arclength(shape, x);
}

typedef acorn::Solver1d<Real> s1d;

struct ArclengthFn : public s1d::Fn {
  ArclengthFn (const ZxFn::Shape shape_, const Real Ltgt_)
    : shape(shape_), Ltgt(Ltgt_)
  {}
  
  void eval (const Real& x, Real& resid, Real& resid_x) const override {
    resid = calc_arclength(shape, x) - Ltgt;
    Real f, g;
    gallery::eval(shape, x, f, g);
    resid_x = std::sqrt(1 + g*g);
  }

private:
  ZxFn::Shape shape;
  Real Ltgt;
};

static void print (const s1d::Info& i) {
  printf("x %12.5e f %9.2e fp %9.2e neval %2d nbisect %2d result %2d\n",
         i.x, i.f, i.fp, i.neval, i.nbisect, i.result);
}

static void
calc_points_by_arclength (const ZxFn::Shape shape, const int nx, const Real L,
                          std::vector<Real>& xs, std::vector<Real>& zs) {
  // Solve for each x in xs to make the arc length difference uniform.
  xs.resize(nx+1);
  xs[0] = 0;
  xs[nx] = 1;
  s1d::Tols tols;
  tols.xtol = tols.ftol = 1e-10;
  const Real Lfac = L/nx;
  for (int i = 1; i < nx; ++i) {
    const Real Ltgt = Real(i)*Lfac;
    s1d::Info info;
    ArclengthFn fn(shape, Ltgt);
    Real flo, fhi, unused;
    fn.eval(xs[i-1], flo, unused);
    fhi = L - Ltgt;
    s1d::fzero(tols, fn, xs[i-1], 1, flo, fhi, info);
    xs[i] = info.x;
    if (0) print(info);
    assert(info.success());
    assert(xs[i] > xs[i-1]);
    assert(xs[i] < 1);
  }

  // Evaluate the shape at each x in xs.
  zs.resize(nx+1);
  for (int i = 0; i <= nx; ++i) {
    Real unused;
    eval(shape, xs[i], zs[i], unused);
  }
}

void ZxFn::set_nx (int nx) {
  assert(nx > 1);

  // Set ny to make roughly uniform squares.
  const auto L = calc_arclength(1);
  ny = std::max(1, int(std::floor(nx/L)));

  calc_points_by_arclength(shape, nx, L, xs, zs);
}

void ZxFn::set_ny (int ny_) { ny = ny_; }

Disloc::Shape Disloc::convert (const std::string& shape) {
  if (shape == "pcosbell") return Shape::pcosbell;
  if (shape == "tapered") return Shape::tapered;
  if (shape == "stapered") return Shape::stapered;
  if (shape == "zero") return Shape::zero;
  fprintf(stderr, "ZxFn::convert: invalid shape %s\n", shape.c_str());
  return Shape::zero;  
}

std::string Disloc::convert (const Shape shape) {
  switch (shape) {
  case Shape::pcosbell: return "pcosbell";
  case Shape::tapered: return "tapered";
  case Shape::stapered: return "stapered";
  case Shape::zero:
  default:
    return "zero";
  }
}

void Disloc
::set (const int dim, const Shape shape, const Real amplitude,
       const Real ctr_x, const Real ctr_y,
       const Real axis_x, const Real axis_y,
       const Real length_x, const Real length_y) {
  assert(dim >= 0 && dim < 3);
  assert(length_x > 0 && length_y > 0);
  auto& d = ds[dim];
  d.shape = shape;
  d.amplitude = amplitude;
  d.ctr[0] = ctr_x; d.ctr[1] = ctr_y;
  d.axis[0][0] = axis_x; d.axis[0][1] = axis_y;
  mv2::normalize(d.axis[0]);
  mv2::perpccw2(d.axis[0], d.axis[1]);
  d.length[0] = length_x; d.length[1] = length_y;
}

static Real pcosbell (const Real x) {
  if (x > 0.5 || x < -0.5) return 0;
  return acorn::square(0.5*(1 + std::cos(2*M_PI*x)));
  //return std::pow(0.5*(1 + std::cos(2*M_PI*x)), 2);
}

static Real tapered (const Real x) {
  if (x > 0.5 || x < -0.5) return 0;
  return std::pow(1 - acorn::square(2*x), 1.5);
}

static Real stapered (const Real x) {
  if (x > 0.5 || x < -0.5) return 0;
  return std::pow(1 - acorn::square(2*x), 3);
}

void Disloc::eval (const Real xy[2], Real disloc[3]) const {
  for (int dim = 0; dim < 3; ++dim) {
    const auto& d = ds[dim];
    if (d.shape == Shape::zero) {
      disloc[dim] = 0;
      continue;
    }
    Real p[2];
    mv2::subtract(xy, d.ctr, p);
    const auto x = mv2::dot(d.axis[0], p), y = mv2::dot(d.axis[1], p);
    const auto r = std::sqrt(acorn::square(x/d.length[0]) +
                             acorn::square(y/d.length[1]));
    Real f;
    switch (d.shape) {
    case Shape::stapered:
      f = d.amplitude*stapered(r);
      break;
    case Shape::tapered:
      f = d.amplitude*tapered(r);
      break;
    case Shape::pcosbell:
    default:
      f = d.amplitude*pcosbell(r);
      break;
    }
    disloc[dim] = f;
  }
}

bool Disloc::is_boundary_zero (const Real threshold) const {
  const auto f = [&] (const Real x, const Real y) {
    const Real p[] = {x, y};
    Real disloc[3];
    eval(p, disloc);
    for (int i = 0; i < 3; ++i)
      if (std::abs(disloc[i]) > threshold) {
        pr(puf(i) pu(x) pu(y) pu(disloc[i]));
        return true;
      }
    return false;
  };

  const int n = 500;
  for (int i = 0; i <= n; ++i) if (f(Real(i)/n, 0)) return false;
  for (int i = 0; i <= n; ++i) if (f(Real(i)/n, 1)) return false;
  for (int i = 0; i <= n; ++i) if (f(0, Real(i)/n)) return false;
  for (int i = 0; i <= n; ++i) if (f(1, Real(i)/n)) return false;

  return true;
}

static int test_arclength () {
  int nerr = 0;
  ZxFn zxfn(ZxFn::Shape::trig1);
  const auto shape = zxfn.get_shape();
  const Real xend = 0.7;
  const auto L = zxfn.calc_arclength(xend);
  const int n = 100;
  Real Lsum = 0, xprev = 0, fprev = 0;
  for (int i = 0; i <= n; ++i) {
    const Real x = xend*Real(i)/n;
    Real f, g;
    eval(shape, x, f, g);
    if (i > 0)
      Lsum += std::sqrt(acorn::square(x - xprev) + acorn::square(f - fprev));
    xprev = x;
    fprev = f;
  }
  if (std::abs(L - Lsum) > 1e-4*Lsum) ++nerr;
  if (nerr) printf("convtest_zx: test_arclength failed\n");
  return nerr;
}

static int test_zxfn () {
  int nerr = 0;
  ZxFn zxfn(ZxFn::Shape::trig1);
  const int nx = 274;
  zxfn.set_nx(nx);
  const auto xs = zxfn.get_xbs();
  const auto zs = zxfn.get_zbs();
  Real segmin = 1, segmax = 0;
  for (int i = 0; i < nx; ++i) {
    const Real seg = std::sqrt(acorn::square(xs[i+1] - xs[i]) +
                               acorn::square(zs[i+1] - zs[i]));
    segmin = std::min(segmin, seg);
    segmax = std::max(segmax, seg);
  }
  if (segmax - segmin > 1e-3*segmax) ++nerr;
  if (nerr) printf("convtest_zx: test_zxfn failed\n");
  return nerr;  
}

int ZxFn::unittest () {
  const int nerr = test_arclength() + test_zxfn();
  if (nerr) printf("ZxFn::unittest failed\n");
  return nerr;
}

} // namespace gallery
} // namespace convzx
} // namespace examples
} // namespace woodland
