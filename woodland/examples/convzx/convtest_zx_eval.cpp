#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/dbg.hpp"

#include "woodland/examples/convzx/convtest_zx.hpp"
#include "woodland/examples/convzx/convtest_zx_hmmvp.hpp"

#include <set>

#define USE_CALC_INTEGRAL_TENSOR_QUADRATURE

namespace woodland {
namespace examples {
namespace convzx {

typedef acorn::Matvec<3,Real> mv3;

static void
get_disloc_comps (const RealArray& dislocs, bool disloc_comps[3]) {
  for (int d = 0; d < 3; ++d)
    disloc_comps[d] = false;
  const auto ncell = dislocs.size() / 3;
  for (auto i = acorn::zero(ncell); i < ncell; ++i) {
    for (int d = 0; d < 3; ++d)
      if (dislocs[3*i+d] != 0)
        disloc_comps[d] = true;
    if (disloc_comps[0] && disloc_comps[1] && disloc_comps[2])
      break;
  }
}

// Determine the sign to apply to the Green's functions when isrc > ircv, for
// GFs obtained in the case isrc <= ircv. isrc and ircv are indices in the y
// direction. ic is the index in the polynomial
//    f(x,y) = c[0] + c[1] x + c[2] y + c[3] x y + c[4] x^2 + c[5] y^2.
inline int
calc_sign (const int isigma, const int idisloc, const int isrc, const int ircv,
           const int ic = 0) {
  return (isrc > ircv ?
          //  sigma:   xy             yz
          ((    isigma == 1 || isigma == 4 ? -1 : 1)*
           // disloc:   y
           (   idisloc == 1                ? -1 : 1)*
           // poly:     y             xy
           (        ic == 2 ||     ic == 3 ? -1 : 1)*
           // poly, cubic terms:
           //         x^2 y          y^3
           (        ic == 7 ||     ic == 9 ? -1 : 1)) :
          1);
}

static bool
are_nbrs (const TriangulationRelations& tr, const Idx t1, const Idx t2) {
  const auto& t2tsi = tr.get_t2etsi();
  const auto& t2ts = tr.get_t2ets();
  for (int i = t2tsi[t1]; i < t2tsi[t1+1]; ++i)
    if (t2ts[i] == t2)
      return true;
  return false;
}

void ConvTest::eval (RealArray& dislocs, RealArray& sigmas,
                     const EvalMethod eval_method) {
  if (not t) discretize();
  const auto ntri = t->get_ntri();
  dislocs.resize(3*ntri);
  sigmas.resize(6*ntri);
  const auto& srf = *d->get_surface();
  ompparfor for (Idx ti = 0; ti < ntri; ++ti) {
    Real p[3];
    srf.tri_ctr_xyz(ti, p);
    disloc->eval(p, &dislocs[3*ti]);
  }
  if (is_fast(eval_method))
    eval_fast(dislocs, sigmas);
  else
    eval_direct(dislocs, sigmas);
}

void ConvTest::eval_direct (const RealArray& dislocs, RealArray& sigmas) const {
  Stress::Options o;
  o.hfp_dist_fac = 0;
#ifdef USE_CALC_INTEGRAL_TENSOR_QUADRATURE
  Stress::Options o1 = o;
  o1.qp.np_radial = 40;
  o1.qp.np_angular = 40;
  o1.use_calc_integral_tensor_quadrature = true;
#endif

  const auto& t = *d->get_triangulation();
  const auto& tr = *d->get_triangulation_relations();
  const auto& srf = *d->get_surface();
  const auto ncell = t.get_ntri();

  Stress stress;
  stress.init(d);
  stress.init_dislocs(dislocs.data());

  ompparfor for (Idx ir = 0; ir < ncell; ++ir) {
    if (verbosity > 0 && ir / (ncell/10) != (ir-1) / (ncell/10)) {
      printf(" %3.1f", 100*Real(ir)/ncell);
      fflush(stdout);
    }
    auto* const sigma = &sigmas[6*ir];
    for (int i = 0; i < 6; ++i) sigma[i] = 0;
    for (Idx is = 0; is < ncell; ++is) {
      Real s[6];
      const bool nbrs = are_nbrs(tr, ir, is);
#ifdef USE_CALC_INTEGRAL_TENSOR_QUADRATURE
      stress.calc_s1_r1(lam, mu, is, ir, is == ir, s, nbrs ? o1 : o);
#else
      stress.calc_s1_r1(lam, mu, is, ir, is == ir || nbrs, s, o);
#endif
      for (int i = 0; i < 6; ++i) sigma[i] += s[i];
    }
    Real lcs[9];
    srf.tri_ctr_lcs(ir, lcs);
    acorn::rotate_sym_tensor_3x3_RARt(lcs, sigma);
  }
  if (verbosity > 0) printf("\n");
}

// For ntri/rect = 4.
inline int mirror_tri (const int ix) {
  const int tri = ix % 4;
  if (tri == 2) return ix+1;
  if (tri == 3) return ix-1;
  return ix;
}

void ConvTest::eval_fast (const RealArray& dislocs, RealArray& sigmas) const {
  const auto& t = *d->get_triangulation();
  const auto& tr = *d->get_triangulation_relations();
  const auto& srf = *d->get_surface();
  const auto ntri = t.get_ntri();

  const auto nxrect = zxfn->get_nx(), nyrect = zxfn->get_ny();
  const int ntri_per_rect = t.get_ntri() / (nxrect*nyrect);
  assert(ntri_per_rect == 2 || ntri_per_rect == 4);
  // There are y-strips of ntri_per_rect types, each type having a different
  // orientation of the triangles, and thus the following number of y-strips
  // total.
  const auto nystrip = nxrect*ntri_per_rect;

  bool disloc_comps[3];
  get_disloc_comps(dislocs, disloc_comps);

  Stress stress;
  stress.init(d);
  stress.init_dislocs(dislocs.data());
  const auto& dislocs_coefs = stress.get_dislocs_coefs();
  const int ncoef_per_dim = d->get_disloc_ncoef() + 1;

  const auto nthr = acorn::get_max_threads();
  const int ngfs = ncoef_per_dim*3*6*(2*nyrect-1);
  std::vector<Real> wrk(nthr*ngfs);

  Stress::Options options;
  // Use nbr relations rather than geometry to determine when to call calc_hfp
  // on other-interactions.
  options.hfp_dist_fac = 0;
#ifdef USE_CALC_INTEGRAL_TENSOR_QUADRATURE
  Stress::Options o1 = options;
  o1.qp.np_radial = 40;
  o1.qp.np_angular = 40;
  o1.use_calc_integral_tensor_quadrature = true;
#endif

  sigmas.clear();
  sigmas.resize(6*ntri, 0);

  // Provide as much outer-loop ||ism as possible. In the case of ntri/rect = 4,
  // we need to write to different y-strips in the same thread.
  const bool have_symmetry = ntri_per_rect == 4;
  const int stride = have_symmetry ? 2 : 1;
  const int outer_end = nystrip / stride;
  ompparfor for (int outer = 0; outer < outer_end; ++outer) {
    for (int inner = 0; inner < stride; ++inner) {
      // rx is a rcv y-strip.
      const int rx = stride*outer + inner;

      const auto tid = acorn::get_thread_num();
      if (verbosity > 0 && rx / (nystrip/10) != (rx-1) / (nystrip/10)) {
        printf(" %3.1f", 100*Real(rx)/nystrip);
        fflush(stdout);
      }

      // LCS is the same for every rcv in this y-strip.
      Real lcs[9];
      srf.tri_ctr_lcs(rx, lcs);

      for (int sx = 0; sx < nystrip; ++sx) {
        // sx is the src along the source x-strip. There are ntri_per_rect
        // sources per rectangle.

        // Compute Green's functions for this y-strip.
        Real* const gfs = &wrk[ngfs*tid];
        const int smr_end = have_symmetry ? 0 : (nyrect-1);
        for (int smr = -(nyrect-1); smr <= smr_end; ++smr) {
          const int siy = smr > 0 ? nyrect-1 : 0;
          const int is = siy*nystrip + sx;
          const int ir = (siy - smr)*nystrip + rx;
          const bool hfp = is == ir;
          const bool nbrs = are_nbrs(tr, ir, is);
          assert(is >= 0 && is < ntri);
          assert(ir >= 0 && ir < ntri);
          Real disloc[3] = {0}, coef[3*Discretization::reconstruct_ncoef] = {0};
          for (int id = 0; id < 3; ++id) {
            if (not disloc_comps[id]) continue;
            for (int ic = 0; ic < ncoef_per_dim; ++ic) {
              if (ic == 0) disloc[id] = 1;
              else coef[id*Discretization::reconstruct_ncoef + ic - 1] = 1;
              Real s[6];
#ifdef USE_CALC_INTEGRAL_TENSOR_QUADRATURE
              // Use tensor quadrature for adjacent-element other-
              // interactions. To switch to using calc_hfp, hfp -> hfp || nbrs.
              stress.calc_s1_r1(disloc, coef, lam, mu, is, ir, hfp, s,
                                nbrs ? o1 : options);
#else
              stress.calc_s1_r1(disloc, coef, lam, mu, is, ir, hfp || nbrs, s,
                                options);
#endif
              acorn::rotate_sym_tensor_3x3_RARt(lcs, s);
              acorn::copy(6, s,
                          &gfs[6*(ncoef_per_dim*(3*(nyrect-1 + smr) + id) + ic)]);
              if (ic == 0) disloc[id] = 0;
              else coef[id*Discretization::reconstruct_ncoef + ic - 1] = 0;
            }
          }
        }

        // Apply GFs to this y-strip.
        for (int riy = 0; riy < nyrect; ++riy) {
          for (int siy = 0; siy < nyrect; ++siy) {
            int ir, is;
            if (have_symmetry && siy > riy) {
              ir = nystrip*riy + mirror_tri(rx);
              is = nystrip*siy + mirror_tri(sx);
            } else {
              ir = nystrip*riy + rx;
              is = nystrip*siy + sx;
            }
            const int smr = have_symmetry ? -std::abs(siy - riy) : siy - riy;
            const auto sigma = &sigmas[6*ir];
            for (int id = 0; id < 3; ++id) {
              if (not disloc_comps[id]) continue;
              const auto disloc_idx = 3*is + id;
              const auto c0 = dislocs[disloc_idx];
              const auto coefs =
                &dislocs_coefs[disloc_idx*Discretization::reconstruct_ncoef];
              for (int ic = 0; ic < ncoef_per_dim; ++ic) {
                const auto disloc = ic == 0 ? c0 : coefs[ic-1];
                const auto gf =
                  &gfs[6*(ncoef_per_dim*(3*(nyrect-1 + smr) + id) + ic)];
                for (int isig = 0; isig < 6; ++isig) {
                  const int sign = (have_symmetry ?
                                    calc_sign(isig, id, siy, riy, ic) :
                                    1);
                  sigma[isig] += sign*disloc*gf[isig];
                }
              }
            }
          }
        }
      }
    }
  }
  
  if (verbosity > 0) printf("\n");
}

void ConvTest
::calc_rect_ctr (const ZxFn& zxfn, const int nyr, const int irect,
                 Real p[3], RPtr nml, RPtr lengths) {
  const auto xs = zxfn.get_xbs();
  const auto zs = zxfn.get_zbs();
  const auto nxr = zxfn.get_nx();
  const auto ix = irect % nxr, iy = irect / nxr;
  assert(ix < zxfn.get_nx());
  assert(iy < nyr);
  p[0] = (xs[ix] + xs[ix+1])/2;
  p[1] = (iy + 0.5)/nyr;
  p[2] = (zs[ix] + zs[ix+1])/2;
  if (nml) {
    nml[0] = zs[ix] - zs[ix+1];
    nml[1] = 0;
    nml[2] = xs[ix+1] - xs[ix];
    mv3::normalize(nml);
  }
  if (lengths) {
    lengths[0] = std::sqrt(acorn::square(xs[ix+1] - xs[ix]) +
                           acorn::square(zs[ix+1] - zs[ix]));
    lengths[1] = 1.0/nyr;
  }
}

void ConvTest
::fill_dislocs (const ZxFn& zxfn, const int nyr, const Disloc& d,
                RealArray& dislocs) {
  const int nxr = zxfn.get_nx();
  const int nrect = nxr*nyr;
  dislocs.resize(3*nrect);
  ompparfor for (int i = 0; i < nrect; ++i) {
    Real p[3];
    calc_rect_ctr(zxfn, nyr, i, p);
    d.eval(p, &dislocs[3*i]);
  }
}

// This fast method is based on the y-symmetry. First, in the y direction,
// z(x=const,y) = const. Second, the rectangles are uniform in size for
// x=const. Third, we're using an elastic full-space. Thus, we can evaluate
// Green's functions for a y-strip, then use these for the specific points in
// the y-strip, reducing the work by O(ny) over the direct method. If, in
// addition, supports are provided, then any y-strip with no support points in
// it is skipped.
void ConvTest
::eval_okada_fast (const int nxr, int& nyr, const SupportPoints& supports,
                   RealArray& dislocs, RealArray& sigmas) const {
  ZxFn zxfno(zxfn->get_shape());
  zxfno.set_nx(nxr);
  if (nyr < 1) nyr = zxfno.get_ny();

  const int nrect = nxr*nyr;
  const auto nthr = acorn::get_max_threads();

  fill_dislocs(zxfno, nyr, *disloc, dislocs);
  bool disloc_comps[3];
  get_disloc_comps(dislocs, disloc_comps);

  const int nsupport = int(supports.size());
  const bool use_sp = nsupport > 0;
  sigmas.clear();
  sigmas.resize(6*(use_sp ? nsupport : nrect), 0);

  std::vector<int> ix_ystrips;
  std::vector<std::vector<int>> irss(nthr), ir_slotss(nthr);
  if (use_sp) {
    for (int ix_ystrip = 0; ix_ystrip < nxr; ++ix_ystrip) {
      const auto it = std::lower_bound(supports.begin(), supports.end(),
                                       SupportPoint{ix_ystrip, -1});
      // If none, we can skip this y-strip.
      if (it == supports.end() || it->ix > ix_ystrip) continue;
      ix_ystrips.push_back(ix_ystrip);
    }
    for (int tid = 0; tid < nthr; ++tid) {
      irss[tid].resize(nyr);
      ir_slotss[tid].resize(nyr);
    }
  } else {
    ix_ystrips.resize(nxr);
    for (int ix_ystrip = 0; ix_ystrip < nxr; ++ix_ystrip)
      ix_ystrips[ix_ystrip] = ix_ystrip;
    for (int tid = 0; tid < nthr; ++tid) {
      irss[tid].resize(nyr);
      for (int i = 0; i < nyr; ++i)
        irss[tid][i] = i;
    }
  }
  const int nstrips = int(ix_ystrips.size());

  const Real yhat[] = {0, 1, 0};
  const int ngfs = 6*3*nyr;
  std::vector<Real> wrk(nthr*ngfs);
  ompparfor for (int idx = 0; idx < nstrips; ++idx) {
    // ix_ystrip is a set of receivers in a y-strip at the ix_ystrip'th x
    // position.
    const int ix_ystrip = ix_ystrips[idx];
    const auto tid = acorn::get_thread_num();
    if (verbosity > 0 && idx / (nstrips/10) != (idx-1) / (nstrips/10)) {
      printf(" %3.1f", 100*Real(idx)/nstrips);
      fflush(stdout);
    }

    // Get the y-strip's points.
    SupportPoints::const_iterator it;
    if (use_sp) {
      it = std::lower_bound(supports.begin(), supports.end(),
                            SupportPoint{ix_ystrip, -1});
      assert(it != supports.end());
    }

    // Get geometry data for each member of this rcv y-strip.
    Real unused[3], rnml[3], rxhat[3];
    calc_rect_ctr(zxfno, nyr, ix_ystrip, unused, rnml);
    mv3::cross(yhat, rnml, rxhat);

    // Gather the rcv y points that are needed.
    int nir;
    int* irs = irss[tid].data();
    int* ir_slots = ir_slotss[tid].data();
    if (use_sp) {
      nir = 0;
      for (int ir = 0; ir < nyr; ++ir) {
        const auto yit = std::lower_bound(it, supports.end(),
                                          SupportPoint{ix_ystrip, ir});
        if (yit == supports.end() || yit->ix != ix_ystrip || yit->iy != ir)
          continue;
        irs[nir] = ir;
        ir_slots[nir] = int(yit - supports.begin());
        assert(ir_slots[nir] >= 0 && ir_slots[nir] < nsupport);
        ++nir;
      }
    } else {
      nir = nyr;
    }

    Real* const gfs = &wrk[tid*ngfs];
    for (int ix = 0; ix < nxr; ++ix) {
      // ix is the src along the min-y x-strip.

      // Get geometry data for each src in this src y-strip.
      Real src[3], snml[3], slengths[2], sxhat[3];
      calc_rect_ctr(zxfno, nyr, ix, src, snml, slengths);
      mv3::cross(yhat, snml, sxhat);

      // Compute Green's functions.
      for (int iy = 0; iy < nyr; ++iy) {
        // iy is the rcv.
        Real rcv[3];
        calc_rect_ctr(zxfno, nyr, iy*nxr + ix_ystrip, rcv);

        Real disloc[3] = {0};
        for (int id = 0; id < 3; ++id) {
          if (not disloc_comps[id]) continue;
          disloc[id] = 1;
          Real dglbl[3], s[6];
          acorn::tmatvec(sxhat, yhat, snml, disloc, dglbl);
          if (use_woodland_rg0c0)
            acorn::fs3d::calc_sigma_const_disloc_rect(
              lam, mu, src, snml, sxhat, slengths, dglbl, rcv, s,
              20, 20);
          else
            acorn::fs3d::calc_sigma_const_disloc_rect_okada(
              lam, mu, src, snml, sxhat, slengths, dglbl, rcv, s);
          acorn::rotate_sym_tensor_3x3_RARt(rxhat, yhat, rnml, s);
          acorn::copy(6, s, &gfs[6*(3*iy + id)]);
          disloc[id] = 0;
        }
      }

      // Apply GFs to this y-strip.
      for (int ir_idx = 0; ir_idx < nir; ++ir_idx) {
        const int ir = irs[ir_idx];
        const int slot = use_sp ? ir_slots[ir_idx] : ir*nxr + ix_ystrip;        
        Real* const sigma = &sigmas[6*slot];
        for (int is = 0; is < nyr; ++is) {
          const int igf = std::abs(is - ir);
          for (int id = 0; id < 3; ++id) {
            if (not disloc_comps[id]) continue;
            const Real disloc = dislocs[3*(is*nxr + ix) + id];
            const Real* const gf = &gfs[6*(3*igf + id)];
            for (int isig = 0; isig < 6; ++isig)
              sigma[isig] += calc_sign(isig, id, is, ir)*disloc*gf[isig];
          }
        }
      }
    }
  }

  if (verbosity > 0) printf("\nnxr %d xeval %d\n", nxr, nstrips);
}

void ConvTest
::init_okada_hmmvp (const EvalMethod eval_method, const ZxFn& zxfn,
                    const int nyr) const {
  if (not hmmvp) hmmvp = std::make_shared<Hmmvp>();
  hmmvp->set_lam_mu(lam, mu);
  hmmvp->set_tol(hmmvp_tol);
  hmmvp->compress_okada_if_not(eval_method, zxfn, nyr);  
}

void ConvTest
::eval_okada_hmmvp (const int nxr, int& nyr, RealArray& dislocs,
                    RealArray& sigmas) const {
  ZxFn zxfno(zxfn->get_shape());
  zxfno.set_nx(nxr);
  if (nyr < 1) nyr = zxfno.get_ny();
  init_okada_hmmvp(EvalMethod::direct_hmmvp, zxfno, nyr);
  hmmvp->eval_okada(*disloc, dislocs, sigmas);
}

void ConvTest
::eval_okada_fast_hmmvp (const int nxr, int& nyr, const SupportPoints& supports,
                         RealArray& dislocs, RealArray& sigmas) const {
  throw_if(true, "eval_okada_fast_hmmvp is not implemented");
  ZxFn zxfno(zxfn->get_shape());
  zxfno.set_nx(nxr);
  if (nyr < 1) nyr = zxfno.get_ny();
  init_okada_hmmvp(EvalMethod::fast_hmmvp, zxfno, nyr);
  hmmvp->eval_okada_fast(supports, *disloc, dislocs, sigmas);
}

void ConvTest
::eval_okada (const int nxr, int& nyr, RealArray& dislocs, RealArray& sigmas,
              const EvalMethod eval_method) const {
  if (eval_method == EvalMethod::direct_hmmvp) {
    eval_okada_hmmvp(nxr, nyr, dislocs, sigmas);
    return;    
  }
  if (eval_method == EvalMethod::fast_hmmvp) {
    SupportPoints unused;
    eval_okada_fast_hmmvp(nxr, nyr, unused, dislocs, sigmas);
    return;
  }    
  if (eval_method == EvalMethod::fast) {
    SupportPoints unused;
    eval_okada_fast(nxr, nyr, unused, dislocs, sigmas);
    return;
  }

  ZxFn zxfno(zxfn->get_shape());
  zxfno.set_nx(nxr);
  if (nyr < 1) nyr = zxfno.get_ny();
  const int nrect = nxr*nyr;
  
  fill_dislocs(zxfno, nyr, *disloc, dislocs);

  sigmas.clear();
  sigmas.resize(6*nrect, 0);
  const Real yhat[] = {0, 1, 0};
  ompparfor for (int ir = 0; ir < nrect; ++ir) {
    if (verbosity > 0 && ir / (nrect/10) != (ir-1) / (nrect/10)) {
      printf(" %3.1f", 100*Real(ir)/nrect);
      fflush(stdout);
    }
    Real rcv[3], rnml[3];
    calc_rect_ctr(zxfno, nyr, ir, rcv, rnml);
    Real* const sigma = &sigmas[6*ir];
    for (int is = 0; is < nrect; ++is) {
      Real src[3], nml[3], lengths[2], s[6], xhat[3], dglbl[3];
      calc_rect_ctr(zxfno, nyr, is, src, nml, lengths);
      mv3::cross(yhat, nml, xhat);
      acorn::tmatvec(xhat, yhat, nml, &dislocs[3*is], dglbl);
      if (use_woodland_rg0c0)
        acorn::fs3d::calc_sigma_const_disloc_rect(
          lam, mu, src, nml, xhat, lengths, dglbl, rcv, s,
          20, 20);
      else
        acorn::fs3d::calc_sigma_const_disloc_rect_okada(
          lam, mu, src, nml, xhat, lengths, dglbl, rcv, s);
      for (int i = 0; i < 6; ++i) sigma[i] += s[i];
    }
    Real xhat[3];
    mv3::cross(yhat, rnml, xhat);
    acorn::rotate_sym_tensor_3x3_RARt(xhat, yhat, rnml, sigma);
  }

  if (verbosity > 0) printf("\n");
}

static void
calc_bilinear_interp_support (const int nxr, const int nyr, CRPtr xs,
                              const Real p[2], int support[2]) {
  assert(p[0] > xs[0]); assert(p[0] < xs[nxr]);
  assert(p[1] > 0); assert(p[1] < 1);
  // Find the rect p is in.
  int ix = std::lower_bound(xs, xs + nxr + 1, p[0]) - xs - 1;
  assert(ix >= 0); assert(ix < nxr);
  int iy = std::floor(p[1]*nyr);
  // Now determine the quadrant, which in turn determines the support since the
  // support points are rect centers.
  if (ix > 0 && p[0] - xs[ix] < 0.5*(xs[ix+1] - xs[ix])) --ix;
  if (ix == nxr-1) --ix; // handle right-side extrapolation
  const Real ylo = Real(iy)/nyr, yhi = Real(iy+1)/nyr;
  if (iy > 0 && p[1] - ylo < 0.5*(yhi - ylo)) --iy;
  if (iy == nyr-1) --iy;
  support[0] = ix; support[1] = iy;
}

static void
calc_bilinear_interp (const int nxr, const int nyr, CRPtr xs, const Real p[2],
                      const int support[2], const int sigmas_idxs[4],
                      const RealArray& sigmas, Real isigmas[6]) {
  const auto ix = support[0], iy = support[1];
  assert(ix >= 0 && ix < nxr-1);
  assert(iy >= 0 && iy < nyr-1);
  const Real
    xlo = (xs[ix  ] + xs[ix+1])/2,
    xhi = (xs[ix+1] + xs[ix+2])/2,
    a = (p[0] - xlo)/(xhi - xlo),
    ylo = Real(iy + 0.5)/nyr,
    yhi = Real(iy + 1.5)/nyr,
    b = (p[1] - ylo)/(yhi - ylo);
  assert((ix == 0 || ix == nxr-2) ||
         (p[0] >= xlo && p[0] <= xhi));
  assert((iy == 0 || iy == nyr-2) ||
         (p[1] >= ylo && p[1] <= yhi));
  for (int i = 0; i < 6; ++i) isigmas[i] = 0;
  const Real coefs[] = {(1-a)*(1-b), a*(1-b), (1-a)*b, a*b};
  for (int c = 0; c < 4; ++c) {
    const auto coef = coefs[c];
    CRPtr sigma = &sigmas[6*sigmas_idxs[c]];
    for (int i = 0; i < 6; ++i)
      isigmas[i] += coef*sigma[i];
  }
}

static void
calc_bilinear_interp (const int nxr, const int nyr, CRPtr xs, const Real p[2],
                      const int support[2], const RealArray& sigmas,
                      Real isigmas[6]) {
  const int ix = support[0], iy = support[1];
  const int sigmas_idxs[] = {nxr* iy    + ix, nxr* iy    + ix+1,                    
                             nxr*(iy+1) + ix, nxr*(iy+1) + ix+1};
  calc_bilinear_interp(nxr, nyr, xs, p, support, sigmas_idxs, sigmas, isigmas);
}

void ConvTest
::interp_okada (const int nxr, const int nyr, const RealArray& sigmas,
                RealArray& isigmas) const {
  assert(int(sigmas.size()/6) == nxr*nyr);
  ZxFn zxfno(zxfn->get_shape());
  zxfno.set_nx(nxr);
  const auto xs = zxfno.get_xbs();
  const auto ntri = t->get_ntri();
  const auto& srf = *d->get_surface();
  isigmas.resize(6*ntri);
  ompparfor for (int ti = 0; ti < ntri; ++ti) {
    Real p[3];
    srf.tri_ctr_xyz(ti, p);
    int support[2];
    calc_bilinear_interp_support(nxr, nyr, xs, p, support);
    calc_bilinear_interp(nxr, nyr, xs, p, support, sigmas, &isigmas[6*ti]);
  }
}

void ConvTest
::collect_rect_support_points (const int nxr, const int nyr,
                               SupportPoints& supports) {
  std::set<SupportPoint> ssp;
  ZxFn zxfno(zxfn->get_shape());
  zxfno.set_nx(nxr);
  const auto xs = zxfno.get_xbs();
  const auto ntri = t->get_ntri();
  const auto& srf = *d->get_surface();
  for (int ti = 0; ti < ntri; ++ti) {
    Real p[3];
    srf.tri_ctr_xyz(ti, p);
    int sp[2];
    calc_bilinear_interp_support(nxr, nyr, xs, p, sp);
    ssp.insert(SupportPoint{sp[0]  , sp[1]  });
    ssp.insert(SupportPoint{sp[0]+1, sp[1]  });
    ssp.insert(SupportPoint{sp[0]  , sp[1]+1});
    ssp.insert(SupportPoint{sp[0]+1, sp[1]+1});
  }
  supports.resize(ssp.size());
  int i = 0;
  for (const auto& sp : ssp) {
    assert(i == 0 ||
           sp.ix > supports[i].ix ||
           (sp.ix == supports[i].ix && sp.iy > supports[i].iy));
    supports[i++] = sp;
  }
}

void ConvTest
::eval_okada (const int nxr, const int nyr, const SupportPoints& supports,
              RealArray& sigmas, const bool use_hmmvp) {
  RealArray dislocs;
  int inyr = nyr;
  if (use_hmmvp)
    eval_okada_fast_hmmvp(nxr, inyr, supports, dislocs, sigmas);
  else
    eval_okada_fast(nxr, inyr, supports, dislocs, sigmas);
  assert(inyr == nyr);
}

void ConvTest
::interp_okada (const int nxr, const int nyr, const SupportPoints& supports,
                const RealArray& sigmas, RealArray& isigmas) {
  assert(sigmas.size() == 6*supports.size());
  ZxFn zxfno(zxfn->get_shape());
  zxfno.set_nx(nxr);
  const auto xs = zxfno.get_xbs();
  const auto ntri = t->get_ntri();
  const auto& srf = *d->get_surface();
  isigmas.resize(6*ntri);
  const auto b = supports.begin();
  const int nsupport = int(supports.size());
  ompparfor for (int ti = 0; ti < ntri; ++ti) {
    Real p[3];
    srf.tri_ctr_xyz(ti, p);
    int sp[2];
    calc_bilinear_interp_support(nxr, nyr, xs, p, sp);
    const auto it00 = std::lower_bound(supports.begin(), supports.end(),
                                       SupportPoint{sp[0], sp[1]});
    assert(it00->ix == sp[0] && it00->iy == sp[1]);
    auto it01 = it00;
    ++it01;
    assert(it01->ix == sp[0] && it01->iy == sp[1]+1);
    const auto it10 = std::lower_bound(it01, supports.end(),
                                       SupportPoint{sp[0]+1, sp[1]});
    assert(it10->ix == sp[0]+1 && it10->iy == sp[1]);
    auto it11 = it10;
    ++it11;
    assert(it11->ix == sp[0]+1 && it11->iy == sp[1]+1);
    const int sigmas_idxs[] = { int(it00 - b), int(it10 - b),
                                int(it01 - b), int(it11 - b)};
    for (int i = 0; i < 4; ++i)
      assert(sigmas_idxs[i] >= 0 && sigmas_idxs[i] < nsupport);
    calc_bilinear_interp(nxr, nyr, xs, p, sp, sigmas_idxs, sigmas,
                         &isigmas[6*ti]);
  }
}

void ConvTest::interp_okada_to_okada (
  const ZxFn& zxfn_src, const int nyr_src, const RealArray& sigmas_src,
  const ZxFn& zxfn_dst, const int nyr_dst, RealArray& sigmas_dst)
{
  const auto nxr_src = zxfn_src.get_nx(), nxr_dst = zxfn_dst.get_nx();
  assert(int(sigmas_src.size()) == 6*nxr_src*nyr_src);
  const auto nrect = nxr_dst*nyr_dst;
  const auto xs_src = zxfn_src.get_xbs();
  sigmas_dst.resize(6*nxr_dst*nyr_dst);
  ompparfor for (int ir = 0; ir < nrect; ++ir) {
    Real p[3];
    calc_rect_ctr(zxfn_dst, nyr_dst, ir, p);
    int support[2];
    calc_bilinear_interp_support(nxr_src, nyr_src, xs_src, p, support);
    calc_bilinear_interp(nxr_src, nyr_src, xs_src, p, support, sigmas_src,
                         &sigmas_dst[6*ir]);
  }
}

} // namespace convzx
} // namespace examples
} // namespace woodland
