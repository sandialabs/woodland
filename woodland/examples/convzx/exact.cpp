#include "woodland/examples/convzx/exact.hpp"

#include "woodland/acorn/util.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/hs3d.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/elastostatics_integrals.hpp"
#include "woodland/acorn/dbg.hpp"

#include <cmath>

namespace woodland {
namespace examples {
namespace convzx {

typedef std::vector<Real> Ar;
typedef acorn::Matvec<3,Real> mv3;

static int
calc_segments (const Real xlo, const Real xhi, const Real x, const int nxs,
               Ar& xbs) {
  assert(xhi > xlo);
  assert(x >= xlo && x <= xhi);
  assert(nxs >= 0);
  const auto f = (x - xlo)/(xhi - xlo);
  const int nlo = std::round(f*(nxs - 1)), nhi = nxs - nlo - 1;
  assert(nlo >= 0 && nhi >= 0 && nlo+nhi+1 == nxs);
  const Real dx = (xhi - xlo)/nxs;
  Real xlim[2];
  if (nlo == 0) {
    xlim[0] = xlo;
    xlim[1] = xlo + dx;
  } else if (nhi == 0) {
    xlim[0] = xhi - dx;
    xlim[1] = xhi;
  } else {
    xlim[0] = x - dx/2;
    xlim[1] = x + dx/2;
  }
  xbs.resize(nxs+1);
  for (int i = 0; i <= nlo; ++i) {
    const Real a = nlo == 0 ? 0 : Real(i)/nlo;
    xbs[i] = (1-a)*xlo + a*xlim[0];
  }
  for (int i = 0; i <= nhi; ++i) {
    const Real a = nhi == 0 ? 1 : Real(i)/nhi;
    xbs[nlo+1+i] = (1-a)*xlim[1] + a*xhi;
  }
#ifndef NDEBUG
  for (int i = 0; i < nxs; ++i)
    assert(xbs[i] < xbs[i+1]);
  assert(xbs[  0] == xlo);
  assert(xbs[nxs] == xhi);
  assert(x >= xbs[nlo] && x <= xbs[nlo+1]);
#endif
  return nlo;
}

struct Exact::Integrands : public acorn::CallerIntegrands {
  enum { nint = 6 };
  
  Integrands (const Exact& e_, const Real x, const Real y,
              const Ar& xbs, const Ar& ybs, const int ix, const int iy)
    : e(e_)
  {
    rect[0] = xbs[ix  ]; rect[2] = xbs[ix+1];
    rect[4] = xbs[ix+1]; rect[6] = xbs[ix  ];
    rect[1] = ybs[iy  ]; rect[3] = ybs[iy  ];
    rect[5] = ybs[iy+1]; rect[7] = ybs[iy+1];
    rcv[0] = x;
    rcv[1] = y;
    e.d->get_surface(x, y, rcv[2]);
  }

  int nintegrands () const override { return nint; }

  void set_halfspace (const bool v) { halfspace_terms = v; }

  void eval (const int n, const Real* pts, RPtr integrands) const override {
    for (int k = 0; k < n; ++k) {
      const Real x = pts[2*k], y = pts[2*k+1];
      Real z, z_xy[2], lcs[9];
      e.d->get_surface(x, y, z, z_xy, lcs);
      const Real p[] = {x, y, z};
      const Real jacdet = std::sqrt(1 + acorn::square(z_xy[0]) +
                                    acorn::square(z_xy[1]));
      Real disloc_lcs[3], disloc_gcs[3];
      e.d->get_disloc(x, y, disloc_lcs);
      mv3::tmatvec(lcs, disloc_lcs, disloc_gcs);
      const auto integrand = &integrands[nint*k];
      if (halfspace_terms)
        acorn::hs3d::calc_sigma_point_halfspace_terms(
          e.lam, e.mu, p, &lcs[6], disloc_gcs, rcv, integrand);
      else
        acorn::fs3d::calc_sigma_point(e.lam, e.mu, p, &lcs[6], disloc_gcs, rcv,
                                      integrand);
      for (int i = 0; i < nint; ++i) integrand[i] *= jacdet;
    }
  }

  Real permitted_R_min (const Real R_max) const override {
    return e.o.permitted_R_min_factor*R_max;
  }

  acorn::plane::Polygon get_src_polygon () const {
    return acorn::plane::Polygon(rect, 4);
  }
  
private:
  const Exact& e;
  Real rect[8], rcv[3];
  bool halfspace_terms = false;
};

void Exact::init (const Description::CPtr& d_) {
  d = d_;
}

void Exact::set_lam_mu (const Real lam_, const Real mu_) {
  lam = lam_; mu = mu_;
}

void Exact::set_halfspace (const bool v) {
  halfspace = v;
}

void Exact::set_options (const Options& o_) {
  o = o_;
  assert(o.nxr > 0 && o.nyr > 0);
}

void Exact::calc_stress (const Real x, const Real y, Real sigma_lcs[6]) const {
  assert(d);
  Real xlo, xhi, ylo, yhi;
  d->get_rectangle_limits(xlo, xhi, ylo, yhi);
  Ar xbs, ybs;
  const int x_in = calc_segments(xlo, xhi, x, o.nxr, xbs);
  const int y_in = calc_segments(ylo, yhi, y, o.nyr, ybs);
  for (int i = 0; i < 6; ++i) sigma_lcs[i] = 0;
  acorn::integrals::Options io;
  io.np_radial = o.np_radial;
  io.np_angular = o.np_angular;
  const Real rcv[] = {x, y};
  Real sigma[6] = {0};
  for (int iy = 0; iy < o.nyr; ++iy)
    for (int ix = 0; ix < o.nxr; ++ix) {
      Integrands igs(*this, x, y, xbs, ybs, ix, iy);
      const bool
        hfp = ix == x_in && iy == y_in,
        tensor_quad = (not hfp &&
                       std::abs(ix - x_in) <= 1 && std::abs(iy - y_in) <= 1),
        tri_quad = not hfp && not tensor_quad;
      int triquad_order = o.triquad_order;
      int triquad_order_hs = o.triquad_order;
      if (o.triquad_order < 0 && tri_quad) {
        const Real
          L = std::max(xbs[ix+1] - xbs[ix], ybs[iy+1] - ybs[iy]),
          D = std::max(std::max(xbs[ix] - x, x - xbs[ix+1]),
                       std::max(ybs[iy] - y, y - ybs[iy+1]));
        const Real tol = 1e-16;
        triquad_order = acorn::get_triquad_order(L, D, tol);
        if (halfspace) {
          Real sz, rz;
          d->get_surface((xbs[ix] + xbs[ix+1])/2, (ybs[iy] + ybs[iy+1])/2, sz);
          d->get_surface(x, y, rz);
          triquad_order_hs = acorn::get_triquad_order(
            L, std::sqrt(acorn::square(D) + acorn::square(rz + sz)), tol);
        }
      }
      const auto poly = igs.get_src_polygon();
      if (hfp)
        acorn::integrals::calc_hfp(io, poly, rcv, igs, sigma);
      else if (tensor_quad)
        acorn::integrals::calc_integral_tensor_quadrature(
          io, poly, igs, rcv, true, sigma);
      else
        acorn::integrals::calc_integral(poly, igs, sigma, triquad_order);
      if (halfspace) {
        igs.set_halfspace(true);
        if (tri_quad) {
          assert(triquad_order_hs > 0);
          acorn::integrals::calc_integral(poly, igs, sigma, triquad_order_hs);
        } else {
          acorn::integrals::calc_integral_tensor_quadrature(
            io, poly, igs, rcv, not hfp, sigma);
        }
      }
    }
  {
    Real z, lcs[9];
    d->get_surface(x, y, z, nullptr, lcs);
    acorn::rotate_sym_tensor_3x3_RARt(lcs, sigma);
  }
  acorn::copy(6, sigma, sigma_lcs);
}

int Exact::unittest () {
  int nerr = 0;
  {
    const Real xlo = 0.2, xhi = 1.5;
    const int nxs = 7;
    Ar xbs;
    for (const Real x : {xlo + 0.01, (xlo + xhi)/3, xhi}) {
      const int x_in = calc_segments(xlo, xhi, x, nxs, xbs);
      if (not (x >= xbs[x_in] && x <= xbs[x_in+1])) ++nerr;
    }
  }
  return nerr;
}

} // namespace convzx
} // namespace examples
} // namespace woodland
