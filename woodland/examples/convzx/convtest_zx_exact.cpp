#include "woodland/examples/convzx/convtest_zx.hpp"

#include "woodland/acorn/matvec.hpp"

namespace woodland {
namespace examples {
namespace convzx {


typedef acorn::Matvec<3,Real> mv3;

struct ExactData : public Exact::Description {
  ExactData (const ZxFn::Shape zshape_, const Disloc::CPtr& disloc_)
    : zshape(zshape_), disloc(disloc_)
  {}

  void get_rectangle_limits (Real& xlo, Real& xhi, Real& ylo, Real& yhi)
    const override
  {
    xlo = ylo = 0;
    xhi = yhi = 1;
  }
  
  void get_surface (const Real x, const Real y,
                    Real& z, Real z_xy[2], Real lcs[9]) const override {
    Real g;
    eval(zshape, x, z, g);
    if (z_xy) {
      z_xy[0] = g;
      z_xy[1] = 0;
    }
    if (lcs) {
      for (int i = 0; i < 9; ++i) lcs[i] = 0;
      lcs[4] = 1;
      lcs[6] = -g;
      lcs[8] = 1;
      mv3::normalize(&lcs[6]);
      mv3::cross(&lcs[3], &lcs[6], &lcs[0]);
    }
  }

  void get_disloc (const Real x, const Real y, Real d[3]) const override {
    const Real p[] = {x, y};
    disloc->eval(p, d);
  }

private:
  ZxFn::Shape zshape;
  Disloc::CPtr disloc;
};

static void setup (ConvTest& ct, Exact::Ptr& exact, Exact::Options* co) {
  if (not ct.get_discretization()) ct.discretize();

  if (not exact) exact = std::make_shared<Exact>();
  const auto description = std::make_shared<ExactData>(
    ct.get_zxfn()->get_shape(), ct.get_disloc());
  exact->init(description);
  exact->set_lam_mu(ct.get_lam(), ct.get_mu());
  exact->set_halfspace(ct.get_use_halfspace());

  if (co)
    exact->set_options(*co);
  else {
    Exact::Options o;
    o.nxr = 11; o.nyr = 11;
    o.np_radial = 20; o.np_angular = 20;
    o.triquad_order = -1;
    exact->set_options(o);
  }
}

void ConvTest::eval_exact_at_tri_ctrs (RealArray& sigmas, Exact::Options* o) {
  setup(*this, exact, o);
  const auto& t = *d->get_triangulation();
  const auto& srf = *d->get_surface();
  const auto ncell = t.get_ntri();
  sigmas.resize(6*ncell);
  ompparfor for (Idx ir = 0; ir < ncell; ++ir) {
    if (verbosity > 0 && ir / (ncell/10) != (ir-1) / (ncell/10)) {
      printf(" %3.1f", 100*Real(ir)/ncell);
      fflush(stdout);
    }
    auto* const sigma = &sigmas[6*ir];
    Real xy[3];
    srf.tri_ctr_xyz(ir, xy);
    exact->calc_stress(xy[0], xy[1], sigma);
  }
  if (verbosity > 0) printf("\n");
}

void ConvTest
::eval_exact_at_rect_ctrs (const int nxr, const int nyr, RealArray& sigmas,
                           Exact::Options* o) {
  setup(*this, exact, o);
  ZxFn zxfno(zxfn->get_shape());
  zxfno.set_nx(nxr);
  const int ncell = nxr*nyr;
  sigmas.resize(6*ncell);
  ompparfor for (Idx ir = 0; ir < ncell; ++ir) {
    if (verbosity > 0 && ir / (ncell/10) != (ir-1) / (ncell/10)) {
      printf(" %3.1f", 100*Real(ir)/ncell);
      fflush(stdout);
    }
    auto* const sigma = &sigmas[6*ir];
    Real xy[2];
    calc_rect_ctr(zxfno, nyr, ir, xy);
    exact->calc_stress(xy[0], xy[1], sigma);
  }
  if (verbosity > 0) printf("\n");
}


} // namespace convzx
} // namespace examples
} // namespace woodland
