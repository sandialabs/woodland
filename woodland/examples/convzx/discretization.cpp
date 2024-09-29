#include "woodland/examples/convzx/discretization.hpp"

#include "woodland/acorn/compose_triquad.hpp"
#include "woodland/acorn/compose_lagrange_polynomial.hpp"
#include "woodland/acorn/bezier_cubic.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/linalg.hpp"
#include "woodland/acorn/dbg.hpp"

#include "woodland/examples/convzx/convtest_zx.hpp"

#include <set>
#include <functional>

namespace woodland {
namespace examples {
namespace convzx {

typedef TriangulationRelations::IdxArray IdxArray;
typedef Discretization::RealArray RealArray;
typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;
typedef acorn::Triangle2D t2d;
typedef acorn::BezierCubic<Real> bc;
using acorn::prarr;

static void asserti (const Triangulation& t, const Idx ti) {
  assert(ti >= 0);
  assert(ti < t.get_ntri());
}

static void init_v2ts (const Triangulation& t, IdxArray& v2tsi, IdxArray& v2ts) {
  const auto nvtx = t.get_nvtx();
  const auto ntri = t.get_ntri();

  std::vector<std::set<Idx>> d(nvtx);
  for (Idx ti = 0; ti < ntri; ++ti) {
    const auto tri = t.get_tri(ti);
    for (int i = 0; i < 3; ++i)
      d[tri[i]].insert(ti);
  }

  v2tsi.resize(nvtx+1);
  v2tsi[0] = 0;
  for (Idx vi = 0; vi < nvtx; ++vi) {
    v2tsi[vi+1] = v2tsi[vi] + d[vi].size();
    for (const auto& e : d[vi])
      v2ts.push_back(e);
  }
}

static void
init_t2ts (const Triangulation& t, const IdxArray& v2tsi, const IdxArray& v2ts,
           IdxArray& t2tsi, IdxArray& t2ts, const int halo = 1) {
  assert(halo == 1 || halo == 2);
  const auto ntri = t.get_ntri();
  t2tsi.resize(ntri+1);
  t2tsi[0] = 0;
  std::set<Idx> d;
  for (Idx ti = 0; ti < ntri; ++ti) {
    d.clear();
    const auto tri = t.get_tri(ti);
    for (int i = 0; i < 3; ++i) {
      const auto vi = tri[i];
      for (Idx iidx = v2tsi[vi]; iidx < v2tsi[vi+1]; ++iidx) {
        if (v2ts[iidx] == ti) continue;
        d.insert(v2ts[iidx]);
        if (halo > 1) {
          const auto trj = t.get_tri(v2ts[iidx]);
          for (int j = 0; j < 3; ++j) {
            const auto vj = trj[j];
            for (Idx jidx = v2tsi[vj]; jidx < v2tsi[vj+1]; ++jidx) {
              if ((v2ts[jidx] == ti)) continue;
              d.insert(v2ts[jidx]);
            }
          }
        }
      }
    }
    t2tsi[ti+1] = t2tsi[ti] + d.size();
    for (const auto& e : d)
      t2ts.push_back(e);
  }
}

static void
init_t2ets (const Triangulation& t, const IdxArray& t2tsi, const IdxArray& t2ts,
            IdxArray& t2etsi, IdxArray& t2ets_) {
  const auto ntri = t.get_ntri();
  t2etsi.resize(ntri+1);
  t2etsi[0] = 0;
  std::vector<Idx> t2ets;
  for (Idx ti = 0; ti < ntri; ++ti) {
    t2etsi[ti+1] = t2etsi[ti];
    const auto tri = t.get_tri(ti);
    for (Idx j = t2tsi[ti]; j < t2tsi[ti+1]; ++j) {
      const auto tj = t2ts[j];
      const auto trj = t.get_tri(tj);
      int n = 0;
      for (int vi = 0; vi < 3; ++vi)
        for (int vj = 0; vj < 3; ++vj)
          if (tri[vi] == trj[vj]) ++n;
      assert(n <= 2);
      if (n == 2) {
        t2ets.push_back(tj);
        ++t2etsi[ti+1];
      }
    }
    assert(t2etsi[ti+1] - t2etsi[ti] >= 0 &&
           t2etsi[ti+1] - t2etsi[ti] <= 3);
  }
  t2ets_.resize(t2ets.size());
  std::copy(t2ets.begin(), t2ets.end(), t2ets_.begin());
}

static void calc_cross (const Triangulation& t, const Idx ti, RPtr nml) {
  const auto tri = t.get_tri(ti);
  const auto a = t.get_vtx(tri[0]);
  Real v1[3], v2[3];
  mv3::subtract(t.get_vtx(tri[1]), a, v1);
  mv3::subtract(t.get_vtx(tri[2]), a, v2);
  mv3::cross(v1, v2, nml);
}

static void
init_xhat_from_primary (const Real zhat[3], const Real primary[3],
                        Real xhat[3]) {
  const auto alpha = mv3::dot(zhat, primary);
  mv3::axpbyz(1, primary, -alpha, zhat, xhat);
  mv3::normalize(xhat);
}

// Collect all vtxs near, but not including, vi.
static Size collect_vtxs (
  const Triangulation& t, const IdxArray& v2tsi, const IdxArray& v2ts,
  const Idx vi, std::vector<Real>& w, const int halo = 1)
{
  assert(halo > 0 && halo < 3);
  std::set<Idx> vs;
  for (Idx j = v2tsi[vi]; j < v2tsi[vi+1]; ++j) {
    const auto tri = t.get_tri(v2ts[j]);
    for (int i = 0; i < 3; ++i) {
      const int vo = tri[i];
      if (vo != vi) {
        vs.insert(vo);
        if (halo > 1) {
          for (Idx k = v2tsi[vo]; k < v2tsi[vo+1]; ++k) {
            const auto trio = t.get_tri(v2ts[k]);
            for (int l = 0; l < 3; ++l)
              if (trio[l] != vi) vs.insert(trio[l]);
          }
        }
      }
    }
  }
  const auto nv = vs.size();
  const auto n = 3*nv;
  if (w.size() < n) w.resize(n);
  Idx i = 0;
  for (const auto e : vs) {
    mv3::copy(t.get_vtx(e), &w[3*i]);
    ++i;
  }
  return nv;
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

void estimate_vertex_normals (const Triangulation& t, const IdxArray& v2tsi,
                              const IdxArray& v2ts, RPtr vtx_nmls,
                              const bool simple) {
  const auto nvtx = t.get_nvtx();
  const auto nthr = acorn::get_max_threads();
  std::vector<std::vector<Real>> w(nthr), rwrk(nthr);
  ompparfor for (Idx vi = 0; vi < nvtx; ++vi) {
    const auto tid = acorn::get_thread_num();
    Real vnml[3] = {0};
    // Initial vtx nml based on weighted average of adjacent triangle face nmls.
    for (Idx j = v2tsi[vi]; j < v2tsi[vi+1]; ++j) {
      Real tnml[3];
      calc_cross(t, v2ts[j], tnml);
      mv3::axpy(1, tnml, vnml);
    }
    mv3::normalize(vnml);
    if (simple) {
      mv3::copy(vnml, &vtx_nmls[3*vi]);
      continue;
    }
    // Use this initial nml to project adjacent vertices onto a plane.
    const auto nv = collect_vtxs(t, v2tsi, v2ts, vi, w[tid], 2);
    const auto vi_vtx = t.get_vtx(vi);
    // Transform vtxs to LCS.
    Real xhat[3] = {1,0,0}, yhat[3];
    if (mv3::dot(vnml, xhat) < 0.1) { xhat[0] = 0; xhat[1] = 1; }
    mv3::cross(vnml, xhat, yhat);
    mv3::normalize(yhat);
    mv3::cross(yhat, vnml, xhat);
    mv3::normalize(xhat);
    for (Idx i = 0; i < nv; ++i) {
      Real tmp[3];
      Real* const vtx = &w[tid][3*i];
      mv3::subtract(vtx, vi_vtx, tmp);
      acorn::matvec(xhat, yhat, vnml, tmp, vtx);
    }
    // LS fit a quadratic to z(x,y) in the temporary LCS.
    Real coef[Polynomial2D::ncoef];
    const size_t rwrk_sz = nv*(Polynomial2D::ncoef+1);
    if (rwrk[tid].size() < rwrk_sz) rwrk[tid].resize(rwrk_sz);
    Polynomial2D::fit(nv, w[tid].data(), rwrk[tid].data(), coef);
    // Calculate the vtx nml in the temporary LCS from the quadratic fit.
    Real g[2];
    Polynomial2D::eval_grad_at_0(coef, g);
    Real nml_lcs[] = {-g[0], -g[1], 1};
    mv3::normalize(nml_lcs);
    // Transform nml to GCS.
    acorn::tmatvec(xhat, yhat, vnml, nml_lcs, vnml);
    mv3::normalize(vnml);
    mv3::copy(vnml, &vtx_nmls[3*vi]);
  }
}

static void init_reconstruction (
  const Triangulation& t, const Surface& srf, const IdxArray& t2tsi,
  const IdxArray& t2ts, const IdxArray& t2etsi, const IdxArray& t2ets,
  const int order, RealArray& rcs)
{
  rcs.resize(Polynomial2D::ncoef*t2ts.size());
  const auto ntri = t.get_ntri();
  std::vector<std::vector<Real>> wrk(acorn::get_max_threads());
  const auto nc_linear = Polynomial2D::order2nc(1);
  const auto nc_max = std::max(nc_linear, Polynomial2D::order2nc(order));
  ompparfor for (int ti = 0; ti < ntri; ++ti) {
    const auto tid = acorn::get_thread_num();
    const Idx os = t2tsi[ti];
    Real* const rc = &rcs[Polynomial2D::ncoef*os];
    const Size nt = t2tsi[ti+1] - os;
    const size_t nwrk = nt*(Polynomial2D::ncoef + 4);
    if (wrk[tid].size() < nwrk) wrk[tid].resize(nwrk);
    Real* const ctrs = wrk[tid].data();
    Real* const rwrk = ctrs + 3*nt;
    Real ti_ctr[3], lcs[9];
    srf.tri_ctr_xyz(ti, ti_ctr);
    srf.tri_ctr_lcs(ti, lcs);
    for (int j = 0; j < nt; ++j) {
      const Idx tj = t2ts[os+j];
      Real tj_ctr[3], tmp[3];
      srf.tri_ctr_xyz(tj, tj_ctr);
      mv3::subtract(tj_ctr, ti_ctr, tmp);
      mv3::matvec(lcs, tmp, tj_ctr);
      ctrs[3*j+0] = tj_ctr[0];
      ctrs[3*j+1] = tj_ctr[1];
      // Ignore the LCS z component. This will instead be disloc(x,y).
      ctrs[3*j+2] = 0;
    }
    bool on_bdy = t2etsi[ti+1] - t2etsi[ti] < 3;
    if (not on_bdy) {
      for (int j = t2etsi[ti]; j < t2etsi[ti+1]; ++j) {
        const auto tj = t2ets[j];
        if (t2etsi[tj+1] - t2etsi[tj] < 3) {
          on_bdy = true;
          break;
        }
      }
    }
    for (int j = 0; j < nt; ++j) {
      ctrs[3*j+2] = 1;
      Real coef[Polynomial2D::ncoef];
      Polynomial2D::fit(nt, ctrs, rwrk, coef, on_bdy ? nc_linear : nc_max);
      acorn::copy(Polynomial2D::ncoef, coef, &rc[Polynomial2D::ncoef*j]);
      ctrs[3*j+2] = 0;
    }
  }
}

TriangulationRelations::TriangulationRelations (const Triangulation::CPtr& t_) {
  t = t_;
  init_v2ts(*t, v2tsi, v2ts);
  init_t2ts(*t, v2tsi, v2ts, t2tsi, t2ts);
  init_t2ets(*t, t2tsi, t2ts, t2etsi, t2ets);
}

void TriangulationRelations::make_t2ts2 () {
  if (t2ts2.empty())
    init_t2ts(*t, v2tsi, v2ts, t2ts2i, t2ts2, 2);
}

struct Area : public acorn::CallerIntegrands {
  Area (const Surface& srf_, const Idx ti_) : srf(srf_), ti(ti_) {}
  int nintegrands () const override { return 1; }
  void eval (const int n, CRPtr xys_lcs, RPtr integrand) const override {
    for (int i = 0; i < n; ++i) {
      Real p_gcs[3], jacdet;
      const auto xy_lcs = &xys_lcs[2*i];
      srf.tri_position1(ti, xy_lcs, p_gcs, nullptr, &jacdet);
      integrand[i] = jacdet;
    }
  }
private:
  const Surface& srf;
  const Idx ti;
};

Real tri_surface_area (Workspace& w, const Surface& srf, const Idx ti) {
  Real vtxs_lcs[6], area = 0;
  srf.tri_vtxs_uv(ti, vtxs_lcs);
  const acorn::integrals::Polygon polygon(vtxs_lcs, 3);
  acorn::integrals::calc_integral(w, polygon, Area(srf, ti), &area);
  return area;
}

static Real eval_f (const ExactParam2DSurface::Name name,
                    const ExactParam2DSurface::UserFn::CPtr& ufn,
                    const Real x, const Real y) {
  using N = ExactParam2DSurface::Name;
  switch (name) {
  case N::zero: return 0;
  case N::trig1:
    return std::sin(x)*std::sin(y);
  case N::trig2:
    return 0.3*std::cos(M_PI*(x - 0.5))*std::cos(M_PI*(y - 0.5));
  case N::user: {
    Real f;
    ufn->eval(x, y, f, nullptr);
    return f;
  } break;
  default:
    assert(0);
    return 0;
  }
}

static void eval_g (const ExactParam2DSurface::Name name,
                    const ExactParam2DSurface::UserFn::CPtr& ufn,
                    const Real x, const Real y, Real g[2]) {
  using N = ExactParam2DSurface::Name;
  switch (name) {
  case N::zero: g[0] = g[1] = 0; break;
  case N::trig1:
    g[0] = std::cos(x)*std::sin(y);
    g[1] = std::sin(x)*std::cos(y);
    break;
  case N::trig2:
    g[0] = -0.3*M_PI*std::sin(M_PI*(x - 0.5))*std::cos(M_PI*(y - 0.5));
    g[1] = -0.3*M_PI*std::cos(M_PI*(x - 0.5))*std::sin(M_PI*(y - 0.5));
    break;
  case N::user: {
    Real f;
    ufn->eval(x, y, f, g);
  } break;
  default:
    assert(0);
    g[0] = g[1] = 0;
  }
}

ExactParam2DSurface
::ExactParam2DSurface (const Triangulation::Ptr& t, const Name name) {
  assert(name != Name::user);
  init(t, name, nullptr);
}

ExactParam2DSurface
::ExactParam2DSurface (const Triangulation::Ptr& t, const UserFn::CPtr& ufn) {
  init(t, Name::user, ufn);
}

void ExactParam2DSurface::init (const Triangulation::Ptr& t_, const Name name_,
                                const UserFn::CPtr ufn_) {
  t = t_;
  name = name_;
  ufn = ufn_;
  const auto ntri = t->get_ntri();
  xy_tris.resize(ntri);
  tri_ctrs_gcs.resize(3*ntri);
  tri_lcs_at_ctrs.resize(9*ntri);
  ompparfor for (Idx ti = 0; ti < ntri; ++ti) {
    const auto tri = t->get_tri(ti);
    auto& xy_tri = xy_tris[ti];
    auto* vtxs = xy_tri.vtxs();
    for (int i = 0; i < 3; ++i)
      mv2::copy(t->get_vtx(tri[i]), &vtxs[2*i]);
    const Real lam[] = {1.0/3.0, 1.0/3.0, 1.0/3.0};
    t2d::barycentric_to_xy(vtxs, lam, xy_tri.ctr);
    auto lcs = &tri_lcs_at_ctrs[9*ti];
    lcs[0] = 1; lcs[1] = 0; lcs[2] = 0;
    tri_position1(ti, xy_tri.ctr, &tri_ctrs_gcs[3*ti], lcs);
  }
}

void ExactParam2DSurface::reset_triangulation_z () {
  const auto nvtx = t->get_nvtx();
  ompparfor for (Idx vi = 0; vi < nvtx; ++vi) {
    auto vtx = t->get_vtx(vi);
    vtx[2] = eval_f(name, ufn, vtx[0], vtx[1]);
  }
}

Size ExactParam2DSurface::get_ntri () const { return t->get_ntri(); }

void ExactParam2DSurface::tri_vtxs_uv (const Idx ti, Real t_lcs[6]) const {
  acorn::copy(6, xy_tris[ti].vtxs(), t_lcs);
}

void ExactParam2DSurface::tri_ctr_uv (const Idx ti, Real p[2]) const {
  asserti(*t, ti);
  mv2::copy(xy_tris[ti].ctr, p);
}

void ExactParam2DSurface::tri_ctr_xyz (const Idx ti, Real p[3]) const {
  asserti(*t, ti);
  mv3::copy(&tri_ctrs_gcs[3*ti], p);
}

void ExactParam2DSurface::tri_ctr_lcs (const Idx ti, Real lcs[9]) const {
  asserti(*t, ti);
  acorn::copy(9, &tri_lcs_at_ctrs[9*ti], lcs);
}

bool ExactParam2DSurface
::calc_tri_lcs (const Idx ti, const Real xyz[3], Real uv[2]) const {
  mv2::copy(xyz, uv);
  return true;
}

void ExactParam2DSurface
::tri_position (const Idx ti, const int n, CRPtr xy, RPtr p_gcs_, RPtr lcs_,
                RPtr jacdet) const {
  for (int i = 0; i < n; ++i) {
    Real p_gcs[3];
    mv2::copy(&xy[2*i], p_gcs);
    p_gcs[2] = eval_f(name, ufn, p_gcs[0], p_gcs[1]);
    mv3::copy(p_gcs, &p_gcs_[3*i]);
    Real grad[2];
    if (jacdet || lcs_) eval_g(name, ufn, p_gcs[0], p_gcs[1], grad);
    if (jacdet) jacdet[i] = std::sqrt(1 + mv2::norm22(grad));
    if (lcs_) {
      Real zhat[] = {-grad[0], -grad[1], 1};
      mv3::normalize(zhat);
      Real xhat[3], yhat[3];
      init_xhat_from_primary(zhat, &tri_lcs_at_ctrs[9*ti], xhat);
      mv3::cross(zhat, xhat, yhat);
      mv3::normalize(yhat);
      mv3::copy(xhat, &lcs_[9*i  ]);
      mv3::copy(yhat, &lcs_[9*i+3]);
      mv3::copy(zhat, &lcs_[9*i+6]);
    }
  }
}

Discretization
::Discretization (const Triangulation::CPtr& t_,
                  const TriangulationRelations::Ptr& tr_,
                  const Surface::CPtr& srf_) {
  t = t_;
  tr = tr_;
  srf = srf_;
  disloc_order = -1;
  set_disloc_order(2);
  assert(disloc_order >= 0);
}

void Discretization::set_disloc_order (const int order) {
  assert(order >= 0 && order <= 3);
  const bool two_halo = order == 3;
  if (two_halo) tr->make_t2ts2();
  const auto& t2tsi = two_halo ? tr->get_t2ts2i() : tr->get_t2tsi();
  const auto& t2ts = two_halo ? tr->get_t2ts2() : tr->get_t2ts();
  if (disloc_order != order)
    init_reconstruction(*t, *srf, t2tsi, t2ts,
                        tr->get_t2etsi(), tr->get_t2ets(), order,
                        reconstruct_coefs);
  disloc_order = order;
}

int Discretization::get_disloc_ncoef () const {
  return Polynomial2D::order2nc(disloc_order);
}

void Discretization
::tri_reconstruct_fit (const Idx ti, const int nfn, CRPtr fns_ctr,
                       RPtr coefs) const {
  const bool two_halo = disloc_order == 3;
  if (two_halo) tr->make_t2ts2();
  const auto& t2tsi = two_halo ? tr->get_t2ts2i() : tr->get_t2tsi();
  const auto& t2ts = two_halo ? tr->get_t2ts2() : tr->get_t2ts();
  const Idx os = t2tsi[ti];
  const Real* const rc = &reconstruct_coefs[Polynomial2D::ncoef*os];
  const Size nt = t2tsi[ti+1] - os;
  for (int i = 0; i < nfn*Polynomial2D::ncoef; ++i)
    coefs[i] = 0;
  const int nc_max = Polynomial2D::order2nc(disloc_order);
  for (int fi = 0; fi < nfn; ++fi) {
    const auto f_ti = fns_ctr[nfn*ti + fi];
    auto coef = &coefs[Polynomial2D::ncoef*fi];
    for (Idx j = 0; j < nt; ++j) {
      const auto tj = t2ts[os+j];
      for (int i = 0; i < nc_max; ++i)
        coef[i] += (fns_ctr[nfn*tj + fi] - f_ti)*rc[j*Polynomial2D::ncoef+i];
    }
  }
}

void Discretization
::tri_reconstruct (const Idx ti, const int nfn, CRPtr coefs, const int n,
                   CRPtr fns_ctr, CRPtr uv, RPtr fns) const {
  tri_reconstruct_ti(ti, nfn, coefs, n, &fns_ctr[nfn*ti], uv, fns);
}

void Discretization
::tri_reconstruct_ti (const Idx ti, const int nfn, CRPtr coefs, const int n,
                      CRPtr fns_ctr, CRPtr uv, RPtr fns) const {
  Real ctr[3], lcs[9];
  srf->tri_ctr_xyz(ti, ctr);
  srf->tri_ctr_lcs(ti, lcs);
  const int nc_max = Polynomial2D::order2nc(disloc_order);
  for (int k = 0; k < n; ++k) {
    Real p[3];
    srf->tri_position1(ti, &uv[2*k], p);
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
  int nerr = 0;
  nerr += Polynomial2D::unittest();
  nerr += ExtrudedCubicSplineSurface::unittest();
  nerr += FlatElementSurface::unittest();
  return nerr;
}

} // namespace convzx
} // namespace examples
} // namespace woodland
