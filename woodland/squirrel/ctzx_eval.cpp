#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/hs3d.hpp"
#include "woodland/acorn/dbg.hpp"

#include "woodland/squirrel/ctzx.hpp"

#include <set>

namespace woodland {
namespace squirrel {
namespace ctzx {

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

void ConvTest
::eval (RealArray& dislocs, RealArray& sigmas, const EvalMethod eval_method) {
  if (not d) discretize();
  const auto nc = d->get_mesh()->get_ncell();
  dislocs.resize(3*nc);
  sigmas.resize(6*nc);
  const auto& srf = *d->get_surface();
  ompparfor for (Idx ci = 0; ci < nc; ++ci) {
    Real p[3];
    srf.cell_ctr_xyz(ci, p);
    disloc->eval(p, &dislocs[3*ci]);
  }
  if (is_fast(eval_method)) {
    if (t) eval_tri_fast(dislocs, sigmas);
    else eval_nonunirect_fast(dislocs, sigmas);
  } else if (is_hmmvp(eval_method)) {
    eval_mesh_hmmvp(dislocs, sigmas);
  } else {
    eval_mesh_direct(dislocs, sigmas);
  }
}

void ConvTest
::eval_mesh_direct (const RealArray& dislocs, RealArray& sigmas) const {
  const auto& m = *d->get_mesh();
  const auto& re = *m.get_relations();
  const auto& srf = *d->get_surface();
  const auto ncell = m.get_ncell();

  Stress::Options o, onbr;
  set_standard_options(o, onbr);

  Stress stress(d);
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
      const bool nbrs = mesh::are_nbrs(re, ir, is);
      stress.calc_s1_r1(gfp, is, ir, s, nbrs ? onbr : o);
      for (int i = 0; i < 6; ++i) sigma[i] += s[i];
    }
    Real lcs[9];
    srf.cell_ctr_lcs(ir, lcs);
    acorn::rotate_sym_tensor_3x3_RARt(lcs, sigma);
  }
  if (verbosity > 0) printf("\n");
}

void ConvTest
::eval_mesh_hmmvp (const RealArray& dislocs, RealArray& sigmas) const {
#ifdef WOODLAND_HAVE_HMMVP
  bool disloc_comps[3];
  get_disloc_comps(dislocs, disloc_comps);
  const auto ncell = d->get_mesh()->get_ncell();
  RealArray x(ncell), y(ncell);
  auto& h = *hmats;
  h.set_greens_function_params(gfp);
  h.set_tol(1e-6);
  h.set_verbosity(1);
  for (int dc = 0; dc < h.ndisloc; ++dc) {
    if (not disloc_comps[dc]) continue;
    for (Idx i = 0; i < ncell; ++i)
      x[i] = dislocs[h.ndisloc*i + dc];
    for (int sc = 0; sc < h.nsigma; ++sc) {
      h.set_filename(dc, sc, "memory");
      h.compress(dc, sc);
      const auto hmat = h.get_hmatrix(dc, sc);
      hmat->get()->Mvp(x.data(), y.data(), 1);
      for (Idx i = 0; i < ncell; ++i)
        sigmas[h.nsigma*i + sc] += y[i];
    }
  }
#endif
}

// For ntri/rect = 4.
inline int mirror_tri (const int ix) {
  const int tri = ix % 4;
  if (tri == 2) return ix+1;
  if (tri == 3) return ix-1;
  return ix;
}

static void calc_gfs (
  const mesh::Relations& re, const GreensFnParams& gfp, const Stress& stress,
  const bool disloc_comps[3], const Stress::Options& options, const Stress::Options& onbr,
  const Real lcs[9], const Idx is, const Idx ir,
  const int ncoef_per_dim, const int gfs_idx, RPtr gfs)
{
  const bool nbrs = mesh::are_nbrs(re, ir, is);
  Real disloc[3] = {0}, coef[3*Discretization::reconstruct_ncoef] = {0};
  for (int id = 0; id < 3; ++id) {
    if (not disloc_comps[id]) continue;
    for (int ic = 0; ic < ncoef_per_dim; ++ic) {
      if (ic == 0) disloc[id] = 1;
      else coef[id*Discretization::reconstruct_ncoef + ic - 1] = 1;
      Real s[6];
      stress.calc_s1_r1(disloc, coef, gfp, is, ir, s,
                        nbrs ? onbr : options);
      acorn::rotate_sym_tensor_3x3_RARt(lcs, s);
      acorn::copy(6, s,
                  &gfs[6*(ncoef_per_dim*(3*gfs_idx + id) + ic)]);
      if (ic == 0) disloc[id] = 0;
      else coef[id*Discretization::reconstruct_ncoef + ic - 1] = 0;
    }
  }
}

static void apply_gfs (
  const bool disloc_comps[3], const RealArray& dislocs_coefs, const RealArray& dislocs,
  const Idx is, const Idx ir, const int siy, const int riy,
  const bool have_symmetry, const int ncoef_per_dim, const int gfs_idx, CRPtr gfs,
  RealArray& sigmas)
{
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
        &gfs[6*(ncoef_per_dim*(3*gfs_idx + id) + ic)];
      for (int isig = 0; isig < 6; ++isig) {
        const int sign = (have_symmetry ?
                          calc_sign(isig, id, siy, riy, ic) :
                          1);
        sigma[isig] += sign*disloc*gf[isig];
      }
    }
  }
}

void ConvTest
::eval_tri_fast (const RealArray& dislocs, RealArray& sigmas) const {
  const auto& m = *d->get_mesh();
  const auto& re = *m.get_relations();
  const auto& srf = *d->get_surface();
  const auto ntri = m.get_ncell();

  const auto nxrect = zxfn->get_nx(), nyrect = zxfn->get_ny();
  const int ntri_per_rect = ntri / (nxrect*nyrect);
  assert(ntri_per_rect == 2 || ntri_per_rect == 4);
  // There are y-strips of ntri_per_rect types, each type having a different
  // orientation of the triangles, and thus the following number of y-strips
  // total.
  const auto nystrip = nxrect*ntri_per_rect;

  bool disloc_comps[3];
  get_disloc_comps(dislocs, disloc_comps);

  Stress stress(d);
  stress.init_dislocs(dislocs.data());
  const auto& dislocs_coefs = stress.get_dislocs_coefs();
  const int ncoef_per_dim = d->get_disloc_ncoef() + 1;

  const auto nthr = acorn::get_max_threads();
  const int ngfs = ncoef_per_dim*3*6*(2*nyrect-1);
  std::vector<Real> wrk(nthr*ngfs);

  Stress::Options options, onbr;
  set_standard_options(options, onbr);

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
      if (verbosity > 0 && (10*rx)/nystrip != (10*(rx-1))/nystrip) {
        printf(" %3.1f", 100*Real(rx)/nystrip);
        fflush(stdout);
      }

      // LCS is the same for every rcv in this y-strip.
      Real lcs[9];
      srf.cell_ctr_lcs(rx, lcs);

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
          assert(is >= 0 && is < ntri);
          assert(ir >= 0 && ir < ntri);
          calc_gfs(re, gfp, stress, disloc_comps, options, onbr,
                   lcs, is, ir, ncoef_per_dim, nyrect-1 + smr, gfs);
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
            apply_gfs(disloc_comps, dislocs_coefs, dislocs, is, ir, siy, riy,
                      have_symmetry, ncoef_per_dim, nyrect-1 + smr, gfs, sigmas);
          }
        }
      }
    }
  }

  if (verbosity > 0) printf("\n");
}

static void
collect_nonunirect_y_strips (const mesh::Mesh& m, std::vector<Idx>& y_strip_ids) {
  const auto& v = *m.get_vtxs();
  const auto& to = *m.get_topo();
  const auto ncell = m.get_ncell();
  const auto get_y_min = [&] (const Idx ci) -> Real {
    const auto c = to.get_cell(ci);
    Real ymin = 1e3;
    for (int i = 0; i < c.n; ++i)
      ymin = std::min(ymin, v.get_vtx(c[i])[1]);
    return ymin;
  };
  y_strip_ids.push_back(0);
  Real ymin = get_y_min(0);
  for (Idx ci = 1; ci < ncell; ++ci) {
    const Real ci_ymin = get_y_min(ci);
    if (ci_ymin < ymin) y_strip_ids.push_back(ci);
    ymin = ci_ymin;
  }
  y_strip_ids.push_back(ncell);
}

void ConvTest
::eval_nonunirect_fast (const RealArray& dislocs, RealArray& sigmas) const {
  const auto& m = *d->get_mesh();
  const auto& re = *m.get_relations();
  const auto& srf = *d->get_surface();
  const auto ncell = m.get_ncell();
  const auto nyrect = zxfn->get_ny();

  std::vector<Idx> y_strip_ids;
  collect_nonunirect_y_strips(m, y_strip_ids);
  const int nystrip = static_cast<int>(y_strip_ids.size() - 1);

  bool disloc_comps[3];
  get_disloc_comps(dislocs, disloc_comps);

  Stress stress(d);
  stress.init_dislocs(dislocs.data());
  const auto& dislocs_coefs = stress.get_dislocs_coefs();
  const int ncoef_per_dim = d->get_disloc_ncoef() + 1;

  const auto nthr = acorn::get_max_threads();
  const int ngfs = ncoef_per_dim*3*6*(2*(2*nyrect)-1);
  std::vector<Real> wrk(nthr*ngfs);

  Stress::Options options, onbr;
  set_standard_options(options, onbr);

  sigmas.clear();
  sigmas.resize(6*ncell, 0);

  ompparfor for (int rx = 0; rx < nystrip; ++rx) {
    const auto tid = acorn::get_thread_num();
    if (verbosity > 0 and (10*rx)/nystrip != (10*(rx-1))/nystrip) {
      printf(" %3.1f", 100*Real(rx)/nystrip);
      fflush(stdout);
    }
    const Idx rcv_ids[] = {y_strip_ids[rx], y_strip_ids[rx+1]-1};
    const Idx nrcv = rcv_ids[1] - rcv_ids[0] + 1;

    // LCS is the same for every rcv in this y-strip.
    Real lcs[9];
    srf.cell_ctr_lcs(rcv_ids[0], lcs);

    for (int sx = 0; sx < nystrip; ++sx) {
      // sx is the src along the source min-y x-strip.
      const Idx src_ids[] = {y_strip_ids[sx], y_strip_ids[sx+1]-1};
      const Idx nsrc = src_ids[1] - src_ids[0] + 1;
      const bool have_symmetry = nsrc == nrcv;
      Real* const gfs = &wrk[ngfs*tid];

      // Compute Green's functions for this y-strip.
      const int nmax = std::max(nsrc, nrcv);
      const bool unequal = nsrc != nrcv;
      const bool nsrc_gt = nsrc > nrcv, nsrc_lt = nsrc < nrcv;
      const int smr_beg = -(nmax-1) + (nsrc_gt ? 1 : 0);
      const int smr_end = have_symmetry ? 0 : ((nmax-1) - (nsrc_lt ? 1 : 0));
      for (int smr = smr_beg; smr <= smr_end; ++smr) {
        // The following indexing is complicated because there are three cases:
        // nsrc = nrcv, and the two unequal cases. It might be better to rewrite
        // this as three separate blocks.
        //   src is either bottom or top cell in src y-strip, with adjustment
        // for unequal case.
        const Idx is =
          (smr <= 0 ?
           src_ids[0] + (nsrc_gt ? (smr % 2 == 0 ?  0 : 1) : 0) :
           src_ids[1] + (nsrc_gt ? (smr % 2 == 0 ? -1 : 0) : 0));
        //   rcv goes from top to bottom twice.
        const Idx ir = ((smr <= 0 ? rcv_ids[0] : rcv_ids[1])
                        - (nsrc_gt ? (smr <= 0 ? (smr-1)/2 : smr/2) : smr)
                        - (smr > 0 and nsrc_lt ? 1 : 0));
        assert(is >= src_ids[0] and is <= src_ids[1]);
        assert(ir >= rcv_ids[0] and ir <= rcv_ids[1]);
        calc_gfs(re, gfp, stress, disloc_comps, options, onbr,
                 lcs, is, ir, ncoef_per_dim, nmax-1 + smr, gfs);
      }

      // Apply GFs to this y-strip.
      for (int riy = 0; riy < nrcv; ++riy) {
        for (int siy = 0; siy < nsrc; ++siy) {
          const auto ir = rcv_ids[0] + riy, is = src_ids[0] + siy;
          int smr = (unequal ?
                     (nsrc_gt ? (siy - 2*riy) : (2*siy - riy)) :
                     (siy - riy));
          if (have_symmetry) smr = -std::abs(smr);
          apply_gfs(disloc_comps, dislocs_coefs, dislocs, is, ir, siy, riy,
                    have_symmetry, ncoef_per_dim, nmax-1 + smr, gfs, sigmas);
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
// x=const. Third, we're using an elastic full-space or a half-space with
// sufficient fault symmetries. Thus, we can evaluate Green's functions for a
// y-strip, then use these for the specific points in the y-strip, reducing the
// work by O(ny) over the direct method. If, in addition, supports are provided,
// then any y-strip with no support points in it is skipped.
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
    if (verbosity > 0 && (10*idx)/nstrips != (10*(idx-1))/nstrips) {
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
      std::swap(slengths[0], slengths[1]); // see [yhat] comment below
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
          acorn::hs3d::calc_sigma_const_disloc_rect(
            w, gfp.lam, gfp.mu, src, snml, yhat, slengths, dglbl, rcv, s,
            gfp.halfspace, not use_woodland_rg0c0, 20, 20);
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
::eval_okada (const int nxr, int& nyr, RealArray& dislocs, RealArray& sigmas,
              const EvalMethod eval_method) const {
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
      // [yhat] To handle the halfspace case with Okada's code, need the yhat
      // direction, which is horizontal, to be the "strike" direction. Swap
      // lengths to accommodate this.
      std::swap(lengths[0], lengths[1]);
      acorn::hs3d::calc_sigma_const_disloc_rect(
        w, gfp.lam, gfp.mu, src, nml, yhat, lengths, dglbl, rcv, s,
        gfp.halfspace, not use_woodland_rg0c0, 20, 20);
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
  const auto nc = d->get_mesh()->get_ncell();
  const auto& srf = *d->get_surface();
  isigmas.resize(6*nc);
  ompparfor for (int ti = 0; ti < nc; ++ti) {
    Real p[3];
    srf.cell_ctr_xyz(ti, p);
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
  const auto nc = d->get_mesh()->get_ncell();
  const auto& srf = *d->get_surface();
  for (int ti = 0; ti < nc; ++ti) {
    Real p[3];
    srf.cell_ctr_xyz(ti, p);
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
              RealArray& sigmas) {
  RealArray dislocs;
  int inyr = nyr;
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
  const auto nc = d->get_mesh()->get_ncell();
  const auto& srf = *d->get_surface();
  isigmas.resize(6*nc);
  const auto b = supports.begin();
  const int nsupport = int(supports.size());
  ompparfor for (int ti = 0; ti < nc; ++ti) {
    Real p[3];
    srf.cell_ctr_xyz(ti, p);
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

} // namespace ctzx
} // namespace squirrel
} // namespace woodland
