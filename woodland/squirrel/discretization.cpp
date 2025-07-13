#include "woodland/squirrel/discretization.hpp"

#include "woodland/acorn/compose_triquad.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/linalg.hpp"
#include "woodland/acorn/dbg.hpp"
#include "woodland/squirrel/ctzx_triangulation.hpp"

#include <set>
#include <functional>

namespace woodland {
namespace squirrel {

typedef Discretization::RealArray RealArray;
typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;
typedef acorn::Triangle2D t2d;
using acorn::prarr;
using mesh::Mesh;

static bool
on_boundary (const mesh::BoolArray& v2bdy, const mesh::Cell& c) {
  for (int i = 0; i < c.n; ++i)
    if (v2bdy[c[i]])
      return true;
  return false;
}

struct Polynomial2D {
  enum : int { ncoef = Discretization::reconstruct_ncoef };

  static int order2nc (const int order) {
    if (order == 0) return 0;
    if (order == 1) return 2;
    if (order == 2) return 5;
    if (order == 3) return 9;
    return -1;
  }

  // Fit
  //     z(x,y) = c[0] x + c[1] y + c[2] x y + c[3] x^2 + c[4] y^2,
  // dropping terms as needed for small nv. N.B. that we assume z(x,y) = 0.
  //   rwrk has size >= nv (Polynomial2D::ncoef + 1).
  //   For nv large enough, add cubic terms, either all or none.
  static void fit (const int nv, CRPtr xyzs, RPtr wrk, Real coef[ncoef],
                   const int nc_max = 9) {
    static_assert(ncoef == 9, "Check ncoef consistency.");
    assert(nv >= 2);
    const int nv_cubic = 11;
    Real R[(ncoef*(ncoef+1))/2], rwrk[ncoef];
    Real* const A = wrk;
    Real* const bx = A + nv*ncoef;
    for (int i = 0; i < nv; ++i) {
      const Real* p = &xyzs[3*i];
      A[0*nv+i] = p[0];
      A[1*nv+i] = p[1];
      A[2*nv+i] = p[0]*p[1];
      A[3*nv+i] = p[0]*p[0];
      A[4*nv+i] = p[1]*p[1];
      if (nv < nv_cubic) continue;
      A[5*nv+i] = p[0]*p[0]*p[0];
      A[6*nv+i] = p[0]*p[0]*p[1];
      A[7*nv+i] = p[0]*p[1]*p[1];
      A[8*nv+i] = p[1]*p[1]*p[1];
    }
    const int nc = std::min(nc_max,
                            nv < 3 ? 2 : nv < 5 ? 3 : nv < nv_cubic ? 5 : 9);
    int iwrk[ncoef];
    acorn::linalg::qr_fac(nv, nc, A, R, iwrk);
    for (int i = 0; i < nv; ++i)
      bx[i] = xyzs[3*i+2];
    acorn::linalg::ls_slv(nv, nc, A, R, iwrk, 1, bx, rwrk);
    for (int i = 0; i < ncoef; ++i)
      coef[i] = 0;
    acorn::copy(nc, bx, coef);
  }

  static Real eval (const Real c[ncoef], const Real x, const Real y,
                    const int nc_max = 9) {
    assert(nc_max == 0 || nc_max == 2 || nc_max == 3 || nc_max == 5 ||
           nc_max == 9);
    Real z = 0;
    if (nc_max > 0)
      z += c[0]*x + c[1]*y;
    if (nc_max > 2)
      z += c[2]*x*y;
    if (nc_max > 3)
      z += c[3]*x*x + c[4]*y*y;
    if (nc_max > 5)
      z += c[5]*x*x*x + c[6]*x*x*y + c[7]*x*y*y + c[8]*y*y*y;
    return z;
  }

  static void eval_grad_at_0 (const Real coef[ncoef], Real grad[2]) {
    grad[0] = coef[0];
    grad[1] = coef[1];
  }

  static int unittest () {
    int nerr = 0;
    const auto z_fn =
      [&] (const Real x, const Real y) {
        return 0.5*x - 1.1*y + 0.2*x*x - 0.3*y*y + 0.3*x*y + 0.1*x*x*y;
      };
    const int nv = 11;
    Real xyzs[3*nv];
    for (int vi = 0; vi < nv; ++vi) {
      const auto x = (acorn::urand() - 0.5);
      const auto y = (acorn::urand() - 0.5);
      xyzs[3*vi+0] = x;
      xyzs[3*vi+1] = y;
      xyzs[3*vi+2] = z_fn(x, y);
    }
    Real c[ncoef], rwrk[nv*(ncoef+1)];
    Polynomial2D::fit(nv, xyzs, rwrk, c);
    for (int vi = 0; vi < nv; ++vi) {
      const auto x = xyzs[3*vi+0];
      const auto y = xyzs[3*vi+1];
      const auto z = xyzs[3*vi+2];
      const auto ze = Polynomial2D::eval(c, x, y);
      if (std::abs(ze - z) > 100*mv2::eps) ++nerr;
    }
    return nerr;
  }
};

static void
calc_reconstruction_coefs (const Mesh& m, const Surface& srf, const int order,
                           RealArray& rcs) {
  const auto to = *m.get_topo();
  const auto re = *m.get_relations();
  const auto& c2csi = re.get_c2csi();
  const auto& c2cs = re.get_c2cs();
  const auto& v2bdy = re.get_v2bdy();
  rcs.resize(Polynomial2D::ncoef*c2cs.size());
  const auto nc = m.get_ncell();
  std::vector<std::vector<Real>> wrk(acorn::get_max_threads());
  const auto nc_linear = Polynomial2D::order2nc(1);
  const auto nc_max = std::max(nc_linear, Polynomial2D::order2nc(order));
  ompparfor for (int ci = 0; ci < nc; ++ci) {
    const auto tid = acorn::get_thread_num();
    const Idx os = c2csi[ci];
    Real* const rc = &rcs[Polynomial2D::ncoef*os];
    const Idx nnbr = c2csi[ci+1] - os;
    const size_t nwrk = nnbr*(Polynomial2D::ncoef + 4);
    if (wrk[tid].size() < nwrk) wrk[tid].resize(nwrk);
    Real* const ctrs = wrk[tid].data();
    Real* const rwrk = ctrs + 3*nnbr;
    Real ci_ctr[3], lcs[9];
    srf.cell_ctr_xyz(ci, ci_ctr);
    srf.cell_ctr_lcs(ci, lcs);
    for (int j = 0; j < nnbr; ++j) {
      const Idx cj = c2cs[os+j];
      Real cj_ctr[3], tmp[3];
      srf.cell_ctr_xyz(cj, cj_ctr);
      mv3::subtract(cj_ctr, ci_ctr, tmp);
      mv3::matvec(lcs, tmp, cj_ctr);
      ctrs[3*j+0] = cj_ctr[0];
      ctrs[3*j+1] = cj_ctr[1];
      // Ignore the LCS z component. This will instead be disloc(x,y).
      ctrs[3*j+2] = 0;
    }
    bool on_bdy = on_boundary(v2bdy, to.get_cell(ci));
    for (int j = 0; j < nnbr; ++j) {
      ctrs[3*j+2] = 1;
      Real coef[Polynomial2D::ncoef];
      Polynomial2D::fit(nnbr, ctrs, rwrk, coef, on_bdy ? nc_linear : nc_max);
      acorn::copy(Polynomial2D::ncoef, coef, &rc[Polynomial2D::ncoef*j]);
      ctrs[3*j+2] = 0;
    }
  }
}

struct Area : public acorn::CallerIntegrands {
  Area (const Surface& srf_, const Idx ci_) : srf(srf_), ci(ci_) {}
  int nintegrands () const override { return 1; }
  void eval (const int n, CRPtr xys_lcs, RPtr integrand) const override {
    for (int i = 0; i < n; ++i) {
      Real p_gcs[3], jacdet;
      const auto xy_lcs = &xys_lcs[2*i];
      srf.cell_position1(ci, xy_lcs, p_gcs, nullptr, &jacdet);
      integrand[i] = jacdet;
    }
  }
private:
  const Surface& srf;
  const Idx ci;
};

Real cell_surface_area (Workspace& w, const Surface& srf, const Idx ci) {
  const int nv_max = 8;
  const auto nv = srf.cell_nvtx(ci);
  assert(nv <= nv_max);
  Real vtxs_lcs[2*nv_max], area = 0;
  srf.cell_vtxs_uv(ci, vtxs_lcs);
  const acorn::integrals::Polygon polygon(vtxs_lcs, nv);
  acorn::integrals::calc_integral(w, polygon, Area(srf, ci), &area);
  return area;
}

// Check that m->get_vtxs() 3D geometry matches srf 3D geometry.
static void check_vtxs (const Mesh& m, const Surface& srf) {
  const auto& v = *m.get_vtxs();
  const auto nv = v.get_nvtx();
  const auto nc = m.get_ncell();
  const auto& to = *m.get_topo();
  std::vector<std::pair<Idx,int>> vi2ci(nv);
  for (Idx ci = 0; ci < nc; ++ci) {
    const auto c = to.get_cell(ci);
    for (int j = 0; j < c.n; ++j)
      vi2ci[c[j]] = std::make_pair(ci, j);
  }
  const auto f = [&] (const Idx vi, Real p[3], const Real*& q) {
    const auto& e = vi2ci[vi];
    const auto ci = e.first;
    const auto j = e.second;
    const auto nv = to.get_cell_nvtx(ci);
    assert(nv <= 8);
    Real uvs[2*8];
    srf.cell_vtxs_uv(ci, uvs);
    srf.cell_position1(ci, &uvs[2*j], p);
    q = v.get_vtx(vi);
  };
  std::vector<Real> dists(nv), norms(nv);
  ompparfor for (Idx vi = 0; vi < nv; ++vi) {
    Real p[3];
    const Real* q;
    f(vi, p, q);
    dists[vi] = mv3::distance2(p, q);
    norms[vi] = mv3::norm22(q);
  }
  Real den = 0;
  for (Idx vi = 0; vi < nv; ++vi)
    den = std::max(den, norms[vi]);
  Real num = 0;
  int vi_max = 0;
  for (Idx vi = 0; vi < nv; ++vi)
    if (dists[vi] > num) {
      num = dists[vi];
      vi_max = vi;
    }
  if (num > acorn::square(1e4*mv3::eps)*den) {
    const Real err = std::sqrt(num/den);
    Real p[3];
    const Real* q;
    f(vi_max, p, q);    
    printf("WARNING check_vtxs:\n"
           "  vi %d vtx %1.5e %1.5e %1.5e srf %1.5e %1.5e %1.5e err %1.3e\n",
           int(vi_max), q[0], q[1], q[2], p[0], p[1], p[2], err);
  }
}

Discretization
::Discretization (const Mesh::CPtr& m_, const Surface::CPtr& srf_,
                  const bool check_consistency) {
  m = m_;
  srf = srf_;
  disloc_order = -1;
  set_disloc_order(2);
  if (check_consistency) {
    throw_if_nomsg(disloc_order < 0);
    throw_if_nomsg(m->get_vtxs()->get_ndim() != 3);
    check_vtxs(*m, *srf);
  }
}

void Discretization::set_disloc_order (const int order) {
  assert(order >= 0 && order <= 3);
  if (disloc_order != order)
    calc_reconstruction_coefs(*m, *srf, order, reconstruct_coefs);
  disloc_order = order;
}

int Discretization::get_disloc_ncoef () const {
  return Polynomial2D::order2nc(disloc_order);
}

void Discretization
::init_reconstruction (const Idx ci, const int nfn, CRPtr fns_ctr,
                       RPtr coefs) const {
  const auto& re = *m->get_relations();
  const auto& c2csi = re.get_c2csi();
  const auto& c2cs = re.get_c2cs();
  const Idx os = c2csi[ci];
  const Real* const rc = &reconstruct_coefs[Polynomial2D::ncoef*os];
  const Idx nc = c2csi[ci+1] - os;
  const int nc_max = Polynomial2D::order2nc(disloc_order);
  for (int i = 0; i < nfn*nc_max; ++i)
    coefs[i] = 0;
  for (int fi = 0; fi < nfn; ++fi) {
    const auto f_ci = fns_ctr[nfn*ci + fi];
    auto coef = &coefs[Polynomial2D::ncoef*fi];
    for (Idx j = 0; j < nc; ++j) {
      const auto cj = c2cs[os+j];
      const auto f_cj = fns_ctr[nfn*cj + fi];
      const auto df = f_cj - f_ci;
      for (int i = 0; i < nc_max; ++i)
        coef[i] += df*rc[j*Polynomial2D::ncoef+i];
    }
  }
}

void Discretization
::init_reconstruction (const Idx ci, const Idx fn_one_idx, RPtr coefs) const {
  const int nfn = 1;
  const auto& re = *m->get_relations();
  const auto& c2csi = re.get_c2csi();
  const auto& c2cs = re.get_c2cs();
  const Idx os = c2csi[ci];
  const Real* const rc = &reconstruct_coefs[Polynomial2D::ncoef*os];
  const Idx nc = c2csi[ci+1] - os;
  const int nc_max = Polynomial2D::order2nc(disloc_order);
  for (int i = 0; i < nfn*nc_max; ++i)
    coefs[i] = 0;
  for (int fi = 0; fi < nfn; ++fi) {
    const auto f_ci = fn_one_idx == ci ? 1 : 0;
    auto coef = &coefs[Polynomial2D::ncoef*fi];
    for (Idx j = 0; j < nc; ++j) {
      const auto cj = c2cs[os+j];
      const auto f_cj = fn_one_idx == cj ? 1 : 0;
      const auto df = f_cj - f_ci;
      for (int i = 0; i < nc_max; ++i)
        coef[i] += df*rc[j*Polynomial2D::ncoef+i];
    }
  }  
}

void Discretization
::reconstruct (const Idx ci, const int nfn, CRPtr coefs, const int n,
               CRPtr fns_ctr, CRPtr uv, RPtr fns) const {
  reconstruct_ci(ci, nfn, coefs, n, &fns_ctr[nfn*ci], uv, fns);
}

void Discretization
::reconstruct_ci (const Idx ci, const int nfn, CRPtr coefs, const int n,
                  CRPtr fns_ctr, CRPtr uv, RPtr fns) const {
  Real ctr[3], lcs[9];
  srf->cell_ctr_xyz(ci, ctr);
  srf->cell_ctr_lcs(ci, lcs);
  const int nc_max = Polynomial2D::order2nc(disloc_order);
  for (int k = 0; k < n; ++k) {
    Real p[3];
    srf->cell_position1(ci, &uv[2*k], p);
    Real pmc[3];
    mv3::subtract(p, ctr, pmc);
    Real p_lcs[3];
    mv3::matvec(lcs, pmc, p_lcs);
    for (int i = 0; i < nfn; ++i)
      fns[nfn*k + i] = (fns_ctr[i] +
                        Polynomial2D::eval(&coefs[Polynomial2D::ncoef*i],
                                           p_lcs[0], p_lcs[1],
                                           nc_max));
  }
}

int Discretization::unittest () {
  int nerr = 0, ne;
  rununittest(Polynomial2D::unittest);
  return nerr;
}

} // namespace squirrel
} // namespace woodland
