#include "woodland/examples/convzx/calc_stress.hpp"

#include "woodland/acorn/util.hpp"
#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/elastostatics_integrals.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/hs3d.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace examples {
namespace convzx {

typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;
typedef acorn::Triangle2D t2d;

Stress::QuadratureParams& Stress::QuadratureParams::set_np_radial (int v)
{ np_radial = v; return *this; }
Stress::QuadratureParams& Stress::QuadratureParams::set_np_angular (int v)
{ np_angular = v; return *this; }
Stress::QuadratureParams& Stress::QuadratureParams::set_triquad_order (int v)
{ triquad_order = v; return *this; }

Stress::Options& Stress::Options::set_hfp_dist_fac (Real v) { hfp_dist_fac = v; return *this; }

static Real
calc_distance (const Triangulation& t, const Idx isrc, const Real rcv[3],
               const bool mirror = false) {
  Real ds[4];
  const auto tri = t.get_tri(isrc);
  for (int i = 0; i < 3; ++i) {
    Real v[3];
    mv3::copy(t.get_vtx(tri[i]), v);
    if (mirror) v[2] *= -1;
    ds[i] = mv3::distance2(v, rcv);
  }
  Real dist = ds[0];
  for (int i = 1; i < 3; ++i)
    dist = std::min(dist, ds[i]);
  dist = std::sqrt(dist);
  return dist;
}


struct Integrands : public acorn::CallerIntegrands {
  enum : int { ndisloc = 3, nsigma = 6 };
  
  Integrands (const Discretization::CPtr& d_, CRPtr disloc_ctr_lcs_,
              CRPtr coefs_, const Real lam_, const Real mu_, const Idx isrc_)
    : d(d_), srf(d->get_surface()), disloc_ctr_lcs(disloc_ctr_lcs_),
      coefs(coefs_), lam(lam_), mu(mu_), isrc(isrc_)
  {
    srf->tri_vtxs_uv(isrc, t_uv);

    //todo Add L to VtxLcs to avoid redundant calcs.
    const auto t = d->get_triangulation();
    const auto tri = t->get_tri(isrc);
    L = mv3::distance2(t->get_vtx(tri[0]), t->get_vtx(tri[1]));
    L = std::max(L, mv3::distance2(t->get_vtx(tri[0]), t->get_vtx(tri[2])));
    L = std::max(L, mv3::distance2(t->get_vtx(tri[1]), t->get_vtx(tri[2])));
    L = std::sqrt(L);
  }

  int nintegrands () const override { return nsigma; }

  Real permitted_R_min (const Real R_max) const override {
    return 1e-3*R_max;
  }

  acorn::plane::Polygon get_src_polygon () const {
    return acorn::plane::Polygon(t_uv, 3);
  }

  void init_rcv_gcs (const Real rcv_gcs_[3]) {
    mv3::copy(rcv_gcs_, rcv_gcs);
    dist = calc_distance(*d->get_triangulation(), isrc, rcv_gcs);
    dist_mirror = -1;
  }

  Real get_src_rcv_dist () const { return dist; }
  Real get_src_length () const { return L; }

  Real get_src_rcv_dist_mirror () const {
    if (dist_mirror < 0)
      dist_mirror = calc_distance(*d->get_triangulation(), isrc, rcv_gcs, true);
    assert(dist_mirror >= dist);
    return dist_mirror;
  }

  void set_halfspace (const bool v) { halfspace_terms = v; }

  void eval (const int n, CRPtr xys_lcs, RPtr integrands) const override {
    const int nper = 16;
    if (wrk.size() < size_t(n*nper)) wrk.resize(n*nper);
    RPtr p_gcs = wrk.data();
    RPtr lcs = p_gcs + 3*n;
    RPtr jacdet = lcs + 9*n;
    RPtr disloc_lcs = jacdet + n;

    srf->tri_position(isrc, n, xys_lcs, p_gcs, lcs, jacdet);
    d->tri_reconstruct_ti(isrc, ndisloc, coefs, n, disloc_ctr_lcs, xys_lcs,
                          disloc_lcs);

    for (int k = 0; k < n; ++k) {
      Real disloc_gcs[3];
      mv3::tmatvec(&lcs[9*k], &disloc_lcs[3*k], disloc_gcs);

      const auto integrand = &integrands[nsigma*k];
      if (halfspace_terms)
        acorn::hs3d::calc_sigma_point_halfspace_terms(
          lam, mu, &p_gcs[3*k], &lcs[9*k+6], disloc_gcs, rcv_gcs, integrand);
      else
        acorn::fs3d::calc_sigma_point(lam, mu, &p_gcs[3*k], &lcs[9*k+6],
                                      disloc_gcs, rcv_gcs, integrand);

      for (int i = 0; i < nsigma; ++i) integrand[i] *= jacdet[k];
    }
  }

  Discretization::CPtr d;
  Surface::CPtr srf;
  // source data
  CRPtr disloc_ctr_lcs, coefs;
  Real lam, mu;
  Idx isrc;
  // receiver data
  Real rcv_gcs[3];
  Real t_uv[6];
  Real dist, L;
  mutable Real dist_mirror;
  mutable std::vector<Real> wrk;
  bool halfspace_terms = false;
};

void Stress::init (const Discretization::CPtr& d_) {
  d = d_;
  dislocs_lcs = nullptr;
}

void Stress::init_dislocs (CRPtr dislocs_lcs_) {
  assert(d);
  const auto ntri = d->get_triangulation()->get_ntri();
  dislocs_lcs = dislocs_lcs_;
  coefs.resize(ntri*ncoef);
  ompparfor for (int ti = 0; ti < ntri; ++ti)
    d->tri_reconstruct_fit(ti, ndisloc, dislocs_lcs, &coefs[ti*ncoef]);
}

void Stress
::calc_s1_r1 (const Real lam, const Real mu, const int isrc, const int ircv,
              const bool self_interaction, Real sigma_out[6],
              const Options o, Info* info) const {
  assert(dislocs_lcs);
  calc_s1_r1(&dislocs_lcs[isrc*ndisloc], &coefs[isrc*ncoef],
             lam, mu, isrc, ircv, self_interaction,
             sigma_out, o, info);
}

void Stress
::calc_s1_r1 (const Real isrc_disloc[3], CRPtr isrc_coefs,
              const Real lam, const Real mu,
              const int isrc, const int ircv,
              const bool hfp, Real sigma_out[6],
              const Options o, Info* info) const {
  assert(d);
  assert(isrc >= 0 && isrc < d->get_triangulation()->get_ntri());

  Integrands igs(d, isrc_disloc, isrc_coefs, lam, mu, isrc);
  const auto src_polygon = igs.get_src_polygon();

  const auto srf = d->get_surface();
  Real rcv_gcs[3];
  srf->tri_ctr_xyz(ircv, rcv_gcs);
  igs.init_rcv_gcs(rcv_gcs);
  const auto dist = isrc == ircv ? 0 : igs.get_src_rcv_dist();
  const auto L = igs.get_src_length();

  Real sigma[6] = {0};
  if (info) {
    info->hfp = false;
    info->triquad_order = -1;
    info->triquad_order_halfspace = -1;
  }
  const bool calc_hfp = dist == 0 || hfp || dist < o.hfp_dist_fac*L;

  acorn::integrals::Options io;
  io.np_radial = o.qp.np_radial;
  io.np_angular = o.qp.np_angular;

  igs.set_halfspace(false);
  if (calc_hfp) {
    if (isrc != ircv && not srf->can_interpolate_outside_element()) {
      printf("Stress::calc_s1_r1: H.f.p. case but can't interp outside element:"
             "isrc %d ircv %d hfp %d dist %1.3e hfp_dist_fac %1.3e L %1.3e\n",
             isrc, ircv, int(hfp), dist, o.hfp_dist_fac, L);
      throw_if(true, "Stress::calc_s1_r1: hfp case");
    }

    Real rcv_lcs[2];
    srf->tri_ctr_uv(ircv, rcv_lcs);
    acorn::integrals::calc_hfp(w, io, src_polygon, rcv_lcs, igs, sigma);
    if (info) info->hfp = true;
  } else {
    if (o.use_calc_integral_tensor_quadrature) {
      Real pt[2];
      const bool nearest_pt = srf->calc_tri_lcs(ircv, rcv_gcs, pt);
      if (not nearest_pt) {
        // Use the cell's center for the split.
        srf->tri_ctr_uv(isrc, pt);
      }
      acorn::integrals::calc_integral_tensor_quadrature(
        w, io, src_polygon, igs, pt, nearest_pt, sigma);
    } else {
      const int triquad_order = (o.qp.triquad_order <= 0 ?
                                 acorn::get_triquad_order(L, dist) :
                                 o.qp.triquad_order);
      acorn::integrals::calc_integral(w, src_polygon, igs, sigma,
                                      triquad_order);
      if (info) info->triquad_order = triquad_order;
    }
  }

  if (o.halfspace) {
    igs.set_halfspace(true);
    if (calc_hfp || o.use_calc_integral_tensor_quadrature) {
      if (dist == 0) {
        Real pt[2];
        srf->tri_ctr_uv(isrc, pt);
        acorn::integrals::calc_integral_tensor_quadrature(
          w, io, src_polygon, igs, pt, false, sigma);
      } else {
        Real pt[2];
        const bool nearest_pt = srf->calc_tri_lcs(ircv, rcv_gcs, pt);
        if (not nearest_pt) {
          // Use the cell's center for the split.
          srf->tri_ctr_uv(isrc, pt);
        }
        acorn::integrals::calc_integral_tensor_quadrature(
          w, io, src_polygon, igs, pt, nearest_pt, sigma);        
      }
    } else {
      const auto dist_mirror = igs.get_src_rcv_dist_mirror();
      const int triquad_order = (o.qp.triquad_order <= 0 ?
                                 acorn::get_triquad_order(L, dist_mirror) :
                                 o.qp.triquad_order);
      acorn::integrals::calc_integral(w, src_polygon, igs, sigma,
                                      triquad_order);
      if (info) info->triquad_order_halfspace = triquad_order;              
    }
  }
  
  acorn::copy(6, sigma, sigma_out);
}

void calc_dist_over_L (const Discretization::CPtr& d, const int isrc, RPtr doL) {
  Integrands igs(d, nullptr, nullptr, 1, 1, isrc);
  const auto ntri = d->get_triangulation()->get_ntri();
  const auto srf = d->get_surface();
  for (int ir = 0; ir < ntri; ++ir) {
    Real rcv_gcs[3];
    srf->tri_ctr_xyz(ir, rcv_gcs);
    igs.init_rcv_gcs(rcv_gcs);
    const auto dist = igs.get_src_rcv_dist();
    const auto L = igs.get_src_length();
    doL[ir] = dist/L;
  }
}

} // namespace convzx
} // namespace examples
} // namespace woodland
