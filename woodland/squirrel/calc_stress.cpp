#include "woodland/squirrel/calc_stress.hpp"

#include "woodland/acorn/util.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/elastostatics_integrals.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/hs3d.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace squirrel {

static const int max_nvtx = 4;

typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;

static Real
calc_distance (const mesh::Mesh& m, const Idx isrc, const Real rcv[3],
               const bool mirror = false) {
  const auto& to = *m.get_topo();
  const auto& vtxs = *m.get_vtxs();
  const auto c = to.get_cell(isrc);
  assert(c.n <= max_nvtx);
  Real dist;
  for (int i = 0; i < c.n; ++i) {
    Real v[3];
    mv3::copy(vtxs.get_vtx(c[i]), v);
    if (mirror) v[2] *= -1;
    const auto idist = mv3::distance2(v, rcv);
    if (i == 0 or idist < dist) dist = idist;
  }
  dist = std::sqrt(dist);
  return dist;
}

struct Integrands : public acorn::CallerIntegrands {
  enum : int { ndisloc = 3, nsigma = 6 };
  
  Integrands (const Discretization::CPtr& d_, CRPtr disloc_ctr_lcs_,
              CRPtr coefs_, const GreensFnParams& gfp_, const Idx isrc_)
    : d(d_), srf(d->get_surface()), disloc_ctr_lcs(disloc_ctr_lcs_),
      coefs(coefs_), gfp(gfp_), isrc(isrc_)
  {
    const auto& m = *d->get_mesh();
    const auto& to = *m.get_topo();
    const auto& v = *m.get_vtxs();
    const auto c = to.get_cell(isrc);

    assert(c.n <= max_nvtx);
    nvtx = c.n;
    srf->cell_vtxs_uv(isrc, t_uv);

    //todo Add L to VtxLcs to avoid redundant calcs.
    assert(v.get_ndim() == 3);
    L = 0;
    for (int i = 0; i < c.n-1; ++i)
      for (int j = i+1; j < c.n; ++j)
        L = std::max(L, mv3::distance2(v.get_vtx(c[i]), v.get_vtx(c[j])));
    L = std::sqrt(L);
  }

  int nintegrands () const override { return nsigma; }

  acorn::plane::Polygon get_src_polygon () const {
    return acorn::plane::Polygon(t_uv, nvtx);
  }

  // rcv_uv makes sense and is only used if isrc == ircv.
  void init_rcv (const Real rcv_uv_[2], const Real rcv_gcs_[3]) {
    mv2::copy(rcv_uv_, rcv_uv);
    mv3::copy(rcv_gcs_, rcv_gcs);
    dist = calc_distance(*d->get_mesh(), isrc, rcv_gcs);
    dist_mirror = -1;
  }

  bool singular_pt (Real p[2]) const override {
    mv2::copy(rcv_uv, p);
    return true;
  }

  Real get_src_rcv_dist () const { return dist; }
  Real get_src_length () const { return L; }

  Real get_src_rcv_dist_mirror () const {
    if (dist_mirror < 0)
      dist_mirror = calc_distance(*d->get_mesh(), isrc, rcv_gcs, true);
    assert(dist_mirror >= dist);
    return dist_mirror;
  }

  void set_halfspace_terms (const bool v) { halfspace_terms = v; }

  Real permitted_r_min (const Real r_max) const override {
    assert(not srf->supports_J());
    return 1e-3*r_max;
  }

  bool supports_mult_by_R3 () const override {
    return not halfspace_terms and srf->supports_J();
  }

  void eval (const int n, CRPtr xys_lcs, RPtr integrands) const override {
    eval(n, xys_lcs, nullptr, integrands);
  }

  void eval_mult_by_R3 (const int n, CRPtr xys_lcs, CRPtr J_times_pdirs,
                        RPtr integrands) const override {
    eval(n, xys_lcs, J_times_pdirs, integrands);
  }

  void eval_shape_J(const Real xy_lcs[2], Real J[6]) const override {
    assert(srf->supports_J());
    Real p_gcs[3];
    srf->cell_position1(isrc, xy_lcs, p_gcs, nullptr, nullptr, J);
  }

  void eval (const int n, CRPtr xys_lcs, CRPtr J_times_pdirs,
             RPtr integrands) const {
    const int nper = 16;
    if (wrk.size() < size_t(n*nper)) wrk.resize(n*nper);
    RPtr p_gcs = wrk.data();
    RPtr lcs = p_gcs + 3*n;
    RPtr jacdet = lcs + 9*n;
    RPtr disloc_lcs = jacdet + n;

    const bool mult_by_R3 = J_times_pdirs;

    srf->cell_position(isrc, n, xys_lcs, p_gcs, lcs, jacdet);
    d->reconstruct_ci(isrc, ndisloc, coefs, n, disloc_ctr_lcs, xys_lcs,
                      disloc_lcs);

    for (int k = 0; k < n; ++k) {
      Real disloc_gcs[3];
      mv3::tmatvec(&lcs[9*k], &disloc_lcs[3*k], disloc_gcs);

      const auto integrand = &integrands[nsigma*k];
      if (halfspace_terms)
        acorn::hs3d::calc_sigma_point_halfspace_terms(
          gfp.lam, gfp.mu, &p_gcs[3*k], &lcs[9*k+6], disloc_gcs, rcv_gcs,
          integrand);
      else {
        const auto dir = mult_by_R3 ? &J_times_pdirs[3*k] : nullptr;
        acorn::fs3d::calc_sigma_point(gfp.lam, gfp.mu, &p_gcs[3*k], &lcs[9*k+6],
                                      disloc_gcs, rcv_gcs, integrand,
                                      mult_by_R3, dir);
      }

      for (int i = 0; i < nsigma; ++i) integrand[i] *= jacdet[k];
    }
  }

  Discretization::CPtr d;
  Surface::CPtr srf;
  // source data
  CRPtr disloc_ctr_lcs, coefs;
  GreensFnParams gfp;
  Idx isrc;
  int nvtx;
  Real t_uv[2*max_nvtx];
  // receiver data
  Real rcv_uv[2], rcv_gcs[3];
  Real dist, L;
  mutable Real dist_mirror;
  mutable std::vector<Real> wrk;
  bool halfspace_terms = false;
};

Stress::Stress (const Discretization::CPtr& d_) {
  d = d_;
}

void Stress::init_dislocs (CRPtr dislocs_lcs_) {
  assert(d);
  const auto nc = d->get_mesh()->get_ncell();
  dislocs_lcs = dislocs_lcs_;
  coefs.resize(nc*ncoef);
  ompparfor for (int ci = 0; ci < nc; ++ci)
    d->init_reconstruction(ci, ndisloc, dislocs_lcs, &coefs[ci*ncoef]);
}

void Stress
::calc_s1_r1 (const GreensFnParams& gfp, const int isrc, const int ircv,
              Real sigma_out[6], const Options o, Info* info) const {
  assert(dislocs_lcs);
  calc_s1_r1(&dislocs_lcs[isrc*ndisloc], &coefs[isrc*ncoef], gfp,
             isrc, ircv, sigma_out, o, info);
}

void Stress
::calc_s1_r1 (const Real isrc_disloc[3], CRPtr isrc_coefs,
              const GreensFnParams& gfp, const int isrc, const int ircv,
              Real sigma_out[6], const Options o, Info* info) const {
  assert(d);
  assert(isrc >= 0 && isrc < d->get_mesh()->get_ncell());
  assert(not (o.force_hfp and o.force_tensor_quadrature));
  assert(gfp.lam >= 0 and gfp.mu >= 0);

  Integrands igs(d, isrc_disloc, isrc_coefs, gfp, isrc);
  const auto src_polygon = igs.get_src_polygon();

  const auto srf = d->get_surface();
  Real rcv_uv[2], rcv_gcs[3];
  srf->cell_ctr_uv (ircv, rcv_uv );
  srf->cell_ctr_xyz(ircv, rcv_gcs);
  igs.init_rcv(rcv_uv, rcv_gcs);
  const auto dist = isrc == ircv ? 0 : igs.get_src_rcv_dist();
  const auto L = igs.get_src_length();

  Real sigma[6] = {0};
  if (info) {
    info->hfp = false;
    info->triquad_order = -1;
    info->triquad_order_halfspace = -1;
  }
  const bool near = dist < o.near_dist_fac*L;
  const bool use_hfp = (o.force_hfp or isrc == ircv or
                        (o.use_hfp_if_near and near));
  const bool use_tq = (o.force_tensor_quadrature or
                       (not o.use_hfp_if_near and near));

  igs.set_halfspace_terms(false);
  if (use_hfp) {
    if (isrc != ircv && not srf->can_interpolate_outside_cell()) {
      printf("Stress::calc_s1_r1: H.f.p. case but can't interp outside element:"
             "isrc %d ircv %d force_hfp %d dist %1.3e near_dist_fac %1.3e L "
             "%1.3e\n", isrc, ircv, int(o.force_hfp), dist, o.near_dist_fac, L);
      throw_if(true, "Stress::calc_s1_r1: hfp case");
    }

    acorn::integrals::Options io;
    io.np_radial = o.qp_hfp.np_radial;
    io.np_angular = o.qp_hfp.np_angular;
    acorn::integrals::calc_hfp(w, io, src_polygon, igs, sigma);
    if (info) info->hfp = true;
  } else {
    if (use_tq) {
      Real pt[2];
      const bool nearest_pt = srf->calc_cell_lcs(ircv, rcv_gcs, pt);
      if (not nearest_pt) {
        // Use the cell's center for the split.
        srf->cell_ctr_uv(isrc, pt);
      }
      acorn::integrals::Options io;
      io.np_radial = o.qp_tq.np_radial;
      io.np_angular = o.qp_tq.np_angular;
      acorn::integrals::calc_integral_tensor_quadrature(
        w, io, src_polygon, igs, pt, nearest_pt, sigma);
    } else {
      const int triquad_order =
        (o.triquad_order <= 0 ?
         acorn::get_triquad_order(L, dist, o.triquad_order_tol) :
         o.triquad_order);
      acorn::integrals::calc_integral(w, src_polygon, igs, sigma,
                                      triquad_order);
      if (info) info->triquad_order = triquad_order;
    }
  }

  if (gfp.halfspace) {
    igs.set_halfspace_terms(true);
    if (use_hfp or use_tq) {
      acorn::integrals::Options io;
      io.np_radial = o.qp_tq.np_radial;
      io.np_angular = o.qp_tq.np_angular;
      if (dist == 0) {
        Real pt[2];
        srf->cell_ctr_uv(isrc, pt);
        acorn::integrals::calc_integral_tensor_quadrature(
          w, io, src_polygon, igs, pt, false, sigma);
      } else {
        Real pt[2];
        const bool nearest_pt = srf->calc_cell_lcs(ircv, rcv_gcs, pt);
        if (not nearest_pt) {
          // Use the cell's center for the split.
          srf->cell_ctr_uv(isrc, pt);
        }
        acorn::integrals::calc_integral_tensor_quadrature(
          w, io, src_polygon, igs, pt, nearest_pt, sigma);        
      }
    } else {
      const auto dist_mirror = igs.get_src_rcv_dist_mirror();
      const int triquad_order =
        (o.triquad_order <= 0 ?
         acorn::get_triquad_order(L, dist_mirror, o.triquad_order_tol) :
         o.triquad_order);
      acorn::integrals::calc_integral(w, src_polygon, igs, sigma,
                                      triquad_order);
      if (info) info->triquad_order_halfspace = triquad_order;              
    }
  }
  
  acorn::copy(6, sigma, sigma_out);
}

void Stress
::calc_matrix_entries (const GreensFnParams& gfp, const Idx si, const Idx ri,
                       const int disloc_comp, Real sigma[6],
                       const Options o, const Options o_nbrs) const {
  const auto& m = *d->get_mesh();
  const auto& srf = *d->get_surface();
  const auto& re = *m.get_relations();
  const auto& c2csi = re.get_c2csi();
  const auto& c2cs = re.get_c2cs();
  const int nnbr = c2csi[si+1] - c2csi[si];
  for (int i = 0; i < 6; ++i) sigma[i] = 0;
  for (int ni = 0; ni <= nnbr; ++ni) {
    const auto sj = ni < nnbr ? c2cs[c2csi[si]+ni] : si;
    Real sj_disloc[3] = {0};
    if (sj == si) sj_disloc[disloc_comp] = 1;
    Real coefs[Stress::ncoef] = {0};
    const int ncoef_per_dim = Discretization::reconstruct_ncoef;
    d->init_reconstruction(sj, si, &coefs[disloc_comp*ncoef_per_dim]);
    const bool nbrs = mesh::are_nbrs(re, ri, sj);
    Real s[6];
    calc_s1_r1(sj_disloc, coefs, gfp, sj, ri, s, nbrs ? o_nbrs : o);
    Real lcs[9];
    srf.cell_ctr_lcs(ri, lcs);
    acorn::rotate_sym_tensor_3x3_RARt(lcs, s);
    for (int i = 0; i < 6; ++i)
      sigma[i] += s[i];
  }
}

void calc_dist_over_L (const Discretization::CPtr& d, const int isrc, RPtr doL) {
  Integrands igs(d, nullptr, nullptr, GreensFnParams(1, 1, false), isrc);
  const auto nc = d->get_mesh()->get_ncell();
  const auto srf = d->get_surface();
  for (int ir = 0; ir < nc; ++ir) {
    Real rcv_uv[2], rcv_gcs[3];
    srf->cell_ctr_uv(ir, rcv_uv);
    srf->cell_ctr_xyz(ir, rcv_gcs);
    igs.init_rcv(rcv_uv, rcv_gcs);
    const auto dist = igs.get_src_rcv_dist();
    const auto L = igs.get_src_length();
    doL[ir] = dist/L;
  }
}

} // namespace squirrel
} // namespace woodland
