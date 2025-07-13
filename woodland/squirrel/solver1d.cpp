#include "woodland/squirrel/solver1d.hpp"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <limits>
#include <algorithm>

#include "woodland/acorn/acorn.hpp"
#include "woodland/acorn/dbg.hpp"

#include "woodland/squirrel/squirrel.hpp"

namespace woodland {
namespace squirrel {

template <typename Real>
Solver1d<Real>::Tols::Tols ()
  : xtol(std::sqrt(std::numeric_limits<Real>::epsilon())),
    ftol(std::sqrt(std::numeric_limits<Real>::epsilon()))
{}

template <typename Real>
Solver1d<Real>::Tols::Tols (const Real xtol_, const Real ftol_)
  : xtol(xtol_), ftol(ftol_)
{}

template <typename Real>
void Solver1d<Real>
::fzero (const Tols& tols, const Fn& fn, const Real& xlo_, const Real& xhi_,
         const Real& flo_, const Real& fhi_, Info& info) {
  static constexpr Real wall = 1e-4;

  Real xlo = xlo_, xhi = xhi_;
  Real flo = flo_, fhi = fhi_;

  assert(tols.xtol >= 0);
  assert(tols.ftol >= 0);

  if (fhi*flo >= 0) {
    info.result = Info::input_f_same_sign;
    return;
  }

#ifndef NDEBUG
  {
    Real f, fp;
    fn.eval(xlo, f, fp);
    assert(f == flo);
    fn.eval(xhi, f, fp);
    assert(f == fhi);
    info.neval += 2;
  }
#endif

  Real x = (xlo + xhi)/2;
  Real f, fp;
  for (;;) {
    fn.eval(x, f, fp);
    ++info.neval;
    if (std::abs(f) <= tols.ftol) {
      info.result = Info::success_on_ftol;
      break;
    }
    if (f*flo < 0) {
      xhi = x;
      fhi = f;
    } else {
      xlo = x;
      flo = f;
    }
    assert(flo*fhi < 0);
    const Real D = xhi - xlo;
    assert(D >= 0);
    if (D < tols.xtol) {
      info.result = Info::success_on_xtol;
      break;
    }
    const bool fp0 = fp == 0;
    if (not fp0)
      x += -f/fp;
    if (fp0 or std::min(x - xlo, xhi - x) < wall*D) {
      x = (xlo + xhi)/2;
      ++info.nbisect;
    }
  }
  info.x = x; info.f = f; info.fp = fp;
}

template <typename Real>
struct TestFn : Solver1d<Real>::Fn {
  TestFn (const int testno_)
    : testno(testno_)
  {}

  void eval (const Real& x, Real& f, Real& fp) const override {
    switch (testno) {
    case 0:
      f = -3*x + 11;
      fp = -3;
      break;
    case 1:
      f = std::cos(1.2*x);
      fp = -1.2*std::sin(1.2*x);
      break;
    default:
      f = fp = 0;
      assert(0);
    }
  }

private:
  int testno;
};

template <typename Real>
int Solver1d<Real>::unittest (const bool verbose) {
  int nerr = 0;
  Tols tols;
  for (int testno = 0; testno <= 1; ++testno) {
    TestFn<Real> fn(testno);
    const Real xlo = -2, xhi = 10;
    Real flo, fhi, unused;
    fn.eval(xlo, flo, unused);
    fn.eval(xhi, fhi, unused);
    Info info;
    fzero(tols, fn, xlo, xhi, flo, fhi, info);
    if (not info.success()) ++nerr;
    if (verbose)
      printf("fzero: x %10.3e f %10.3e fp %10.3e neval %2d nbisect %2d result %d\n",
             info.x, info.f, info.fp, info.neval, info.nbisect, int(info.result));
  }
  return nerr;
}

template struct Solver1d<Real>;

} // namespace squirrel
} // namespace woodland
