#include "woodland/acorn/vv.hpp"
#include "woodland/acorn/macros.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/caller_integrand.hpp"
#include "woodland/acorn/plane_geometry.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/elastostatics_integrals.hpp"
#include "woodland/acorn/compose_quadrature.hpp"
#include "woodland/acorn/compose_lagrange_polynomial.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace acorn {
namespace vv {
namespace convtest {
namespace flatstrip {

/* Some points:
     This test uses a 8x1 strip and refines only in the x direction. Thus, cells
   become increasingly skinny. In turn, calc_hfp is used for an increasing
   number of self-other calculations with refinement. Eventually this is a
   problem, but in practice one should never use such poor-quality cells, and
   doing so keeps this test fast.
     Self-other calc_hfp calls integrate over domains that are a superset of the
   src cell. All of the integral outside of the src cell cancels overall. But
   the disloc function must be smooth for the quadrature to work. Thus, the src
   cell disloc function is a high-order approximation to the true one that
   extends smoothly for some distance beyond the src cell. For high-quality
   cells, a polynomial fit is probably fine. For poor-quality cells like in this
   test, eventually the unstable extrapolation becomes a problem.
     Run example:
       OMP_NUM_THREADS=4 ./testmain conv_test_flat_strip
       hy vv.hy dev        # if python_outfn is defined
 */

using plane::Pt;
typedef Matvec<3,Real> mv3;

Config::Config () {
  width = 1;
  length = 8;
  lam = 1;
  mu = 1;
  n_cell_base = 8;
  n_refine = 8;
  const auto c = TestUvFn::Fn::constant;
  const TestUvFn::Fn ufns[] = {c, c, c};
  const Real uparms[] = {1, 0, 0, 0, 0, 0, 0, 0, 0};
#if 1
  // cosine_bell, tapered, hat
  const TestUvFn::Fn vfns[] = {TestUvFn::Fn::cosine_bell, c, c};
  const Real vparms[] = {1, M_PI, 1*M_PI, 0, 0, 0, 0, 0, 0};
#else
  const TestUvFn::Fn vfns[] = {TestUvFn::Fn::sin, c, c};
  const Real vparms[] = {1, 1, 0, 0, 0, 0, 0, 0, 0};
#endif
  fdisloc.init(ufns, uparms, vfns, vparms);
  // Include for output for `hy vv.hy dev`.
  python_outfn = "pyout_vv_flat_strip.py";
  run_fs = true;
}

static void
eval_sigma_at_ctrs (const Config& c, const int ncell,
                    const SigmaType::Enum stype, RPtr sigmas,
                    RPtr dislocs = nullptr,
                    const QuadratureParams qp = QuadratureParams()) {
  assert(SigmaType::is_valid(stype));
  Workspace w;
  ompparfor for (int ir = 0; ir < ncell; ++ir) {
    auto* const sigma = &sigmas[6*ir];
    for (int i = 0; i < 6; ++i) sigma[i] = 0;
    const Real rcv[] = {c.length*(ir + 0.5)/ncell, 0, 0};
    for (int is = 0; is < ncell; ++is) {
      Real s[6];
      const Real src[] = {c.length*(is + 0.5)/ncell, 0, 0};
      const Real nml[] = {0, 0, 1}, tan[] = {1, 0, 0};
      const Real suv[] = {0.5, 2*M_PI*(is + 0.5)/ncell};
      Real disloc[3];
      c.fdisloc.eval(suv, disloc);
      if (dislocs && is == ir) copy(3, disloc, &dislocs[3*is]);
      const Real xy_side_lens[] = {c.length/ncell, c.width};
      switch (stype) {
      case SigmaType::fs_okada:
#ifdef WOODLAND_ACORN_HAVE_DC3D
        fs3d::calc_sigma_const_disloc_rect_okada(
          c.lam, c.mu, src, nml, tan, xy_side_lens, disloc, rcv, s);
        break;
#endif
      case SigmaType::fs:
        fs3d::calc_sigma_const_disloc_rect(
          w, c.lam, c.mu, src, nml, tan, xy_side_lens, disloc, rcv, s,
          qp.np_radial, qp.np_angular, qp.triquad_order);
        break;
      case SigmaType::fs_exact_geometry_disloc:
      case SigmaType::invalid:
      default: assert(0);
      }
      for (int i = 0; i < 6; ++i) sigma[i] += s[i];
    }
  }
}

struct SigmaExactIntegrands : public CallerIntegrands {
  SigmaExactIntegrands (const Config& c, const int ncell_,
                        const int isrc_, const Pt rcv_)
    : ctfs(c), ncell(ncell_), isrc(isrc_), rcv{rcv_[0], rcv_[1], 0}
  {
    rect[0] = c.length* isrc   /ncell; rect[1] = -c.width/2;
    rect[2] = c.length*(isrc+1)/ncell; rect[3] = -c.width/2;
    rect[4] = c.length*(isrc+1)/ncell; rect[5] =  c.width/2;
    rect[6] = c.length* isrc   /ncell; rect[7] =  c.width/2;
    L = std::max(c.width, c.length/ncell);
    dist = std::max(rcv[0] - rect[2], rect[0] - rcv[0]);

    if (dist > 0) {
      Quadrature q(npt, Quadrature::gl);
      const Real* qx;
      q.get_x(qx);
      for (int i = 0; i < npt; ++i) {
        const auto a = (1 + qx[i])/2;
        x[i] = (1-a)*rect[0] + a*rect[2];
        const Pt uv = {0.5, 2*M_PI*x[i]/ctfs.length};
        Real disloc[3];
        ctfs.fdisloc.eval(uv, disloc);
        for (int j = 0; j < 3; ++j) f[j][i] = disloc[j];
      }
    }
  }

  int nintegrands () const override { return 6; }

  Real permitted_R_min (const Real R_max) const override {
    return 1e-3*R_max;
  }

  void eval (const int n, CRPtr ps, RPtr integrand) const override {
    for (int i = 0; i < n; ++i) {
      const auto p = &ps[2*i];
      const Real src[] = {p[0], p[1], 0};
      const Real nml[] = {0, 0, 1};
      Real disloc[3];
      if (dist > 0) {
        // Smoothly (but unstably) extrapolates; thus, suitable for self-other
        // call if the rcv isn't too far away form the src cell.
        for (int j = 0; j < 3; ++j)
          disloc[j] = eval_lagrange_poly(npt, x, f[j], p[0]);
      } else {
        // Doesn't work for calc_hfp self-other calls if disloc is not everywhere
        // smooth.
        const Real suv[] = {0.5*(1 + p[1]/ctfs.width), 2*M_PI*p[0]/ctfs.length};
        ctfs.fdisloc.eval(suv, disloc);
      }
      fs3d::calc_sigma_point(ctfs.lam, ctfs.mu, src, nml, disloc, rcv,
                             &integrand[6*i]);
      // jacdet = 1
    }
  }

  plane::Polygon get_src_polygon () const {
    return plane::Polygon(rect, 4);
  }

  const Config& ctfs;
  const int ncell;
  const int isrc;
  Real rcv[3], rect[8];
  Real dist, L;

private:
  // Not very stable extrapolation, but good enough for this test.
  enum : int { npt = 6 };
  Real x[npt], f[3][npt];
};

static void
eval_sigma_exact (Workspace& w, const SigmaExactIntegrands& s, RPtr sigma_out,
                  const QuadratureParams qp = QuadratureParams()) {
  Real sigma[6] = {0};
  const plane::Polygon p = s.get_src_polygon();
  if (s.dist < s.L/4) {
    integrals::Options io;
    io.np_radial = qp.np_radial;
    io.np_angular = qp.np_angular;
    integrals::calc_hfp(w, io, p, s.rcv, s, sigma);
  } else {
    const auto triquad_order = qp.triquad_order <= 0 ?
      get_triquad_order(s.L, s.dist) : qp.triquad_order;
    integrals::calc_integral(w, p, s, sigma, triquad_order);
  } 
  for (int i = 0; i < 6; ++i) sigma_out[i] = sigma[i];
}

static void eval_sigma_exact (
  Workspace& w, const Config& c, const int ncell, const int isrc,
  const int nrcv, CRPtr rcv_xy, RPtr sigmas,
  const QuadratureParams qp = QuadratureParams(),
  const bool thread = false)
{
  assert(isrc >= 0 && isrc < ncell);

  if (thread) {
    ompparfor for (int ir = 0; ir < nrcv; ++ir) {
      SigmaExactIntegrands sei(c, ncell, isrc, &rcv_xy[2*ir]);
      eval_sigma_exact(w, sei, &sigmas[6*ir], qp);
    }
  } else {
    for (int ir = 0; ir < nrcv; ++ir) {
      SigmaExactIntegrands sei(c, ncell, isrc, &rcv_xy[2*ir]);
      eval_sigma_exact(w, sei, &sigmas[6*ir], qp);
    }
  }
}

static void eval_sigma_exact_at_ctrs (
  const Config& c_src, const int ncell_src,
  const Config& c_rcv, const int ncell_rcv,
  RPtr sigmas, const QuadratureParams qp = QuadratureParams())
{
  Workspace w;
  ompparfor for (int ir = 0; ir < ncell_rcv; ++ir) {
    const Real rcv_xy[] = {c_rcv.length*(ir + 0.5)/ncell_rcv, 0};
    auto* const sigma = &sigmas[6*ir];
    for (int i = 0; i < 6; ++i) sigma[i] = 0;
    for (int is = 0; is < ncell_src; ++is) {
      Real s[6];
      eval_sigma_exact(w, c_src, ncell_src, is, 1, rcv_xy, s, qp);
      for (int i = 0; i < 6; ++i) sigma[i] += s[i];
    }
  }
}

Results run (const Config& c, const bool verbose) {
  using namespace flatstrip;

  const int max_ncell = c.n_cell_base*(1 << (c.n_refine-1));
  impl::Errors e(c.run_fs ? 3 : 2, c.n_refine, max_ncell);
  std::vector<Real> dislocs(3*max_ncell, 0);

  {
    int i = 0;
    if (c.run_fs)
      e.stypes[i++] = SigmaType::fs;
    e.stypes[i] = SigmaType::fs_exact_geometry_disloc;
  }

  QuadratureParams qpr, qpe;
  qpr.set_np_radial(20).set_np_angular(10);
  qpe.set_np_radial(8).set_np_angular(8);

  FILE* const fid = c.python_outfn.empty() ? nullptr :
    fopen(c.python_outfn.c_str(), "w");
  if (fid) impl::pyout::init(fid, e);
  
  for (int iref = 0; iref < c.n_refine; ++iref) {
    const int ncell = c.n_cell_base*(1 << iref);
    for (const auto stype : {SigmaType::fs_okada, SigmaType::fs}) {
      if (stype == SigmaType::fs && ! c.run_fs) continue;
      if (stype == SigmaType::fs && iref >= 5) {
        static bool first = true;
        if (first && verbose) printf("NOTE: skipping some fs runs\n");
        first = false;
        continue;
      }
      eval_sigma_at_ctrs(c, ncell, stype, e.sigmas[int(stype)].data(),
                         stype == SigmaType::fs_okada ? dislocs.data() : nullptr,
                         qpr);
    }
    {
      const int fac = 8;
      assert(ncell % fac == 0);
      eval_sigma_exact_at_ctrs(c, ncell/fac, c, ncell,
                               e.sigmas[e.stypes.size()].data(), qpe);
    }
    const Real dx = c.length/ncell;
    e.accum_errors(ncell, [&] (int) { return dx; }, iref);
    if (fid) impl::pyout::write(fid, e, dislocs, ncell);
  }

  if (fid) fclose(fid);

  if (verbose) printf("flatstrip\n");
  Results r;
  e.collect_errors(std::vector<int>({2}), &r, verbose);

  return r;
}

} // namespace flatstrip
} // namespace convtest
} // namespace vv
} // namespace acorn
} // namespace woodland
