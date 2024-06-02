#include "woodland/examples/convzx/convtest_zx_surface.hpp"
#include "woodland/examples/convzx/convtest_zx.hpp"

#include "woodland/acorn/bezier_cubic.hpp"
#include "woodland/acorn/compose_lagrange_polynomial.hpp"
#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace examples {
namespace convzx {


typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;
typedef acorn::Triangle2D t2d;
typedef acorn::BezierCubic<Real> bc;
typedef std::vector<Real> RealArray;

static void asserti (const Triangulation& t, const Idx ti) {
  assert(ti >= 0);
  assert(ti < t.get_ntri());
}

static void
init_xhat_from_primary (const Real zhat[3], const Real primary[3],
                        Real xhat[3]) {
  const auto alpha = mv3::dot(zhat, primary);
  mv3::axpbyz(1, primary, -alpha, zhat, xhat);
  mv3::normalize(xhat);
}

ExtrudedCubicSplineSurface::ExtrudedCubicSplineSurface(
  const gallery::ZxFn::Shape zshape, const Triangulation::CPtr& t,
  const bool use_exact_normals, const int nml_recon_order)
{
  init(zshape, t, use_exact_normals, nml_recon_order);
}

static void
init_splines (const gallery::ZxFn::Shape zshape, const bool use_exact_normals,
              const int nml_recon_order, const int nseg, const RealArray& xs,
              const RealArray& zs, const RealArray& ps, RealArray& cs) {
  // Normals.
  RealArray nmls(2*(nseg+1));
  if (use_exact_normals) {
    for (int i = 0; i <= nseg; ++i) {
      Real z, zg;
      gallery::eval(zshape, ps[2*i], z, zg);
      Real nml[] = {-zg, 1};
      mv2::normalize(nml);
      mv2::copy(nml, &nmls[2*i]);
    }
  } else {
    assert(nml_recon_order == 2 || nml_recon_order == 4);
    const int npt = nml_recon_order/2;
    for (int ix = 0; ix <= nseg; ++ix) {
      int s0, s1;
      if (ix - npt < 0) {
        s0 = 0;
        s1 = std::min(nseg, s0 + 2*npt);
      } else if (ix + npt > nseg) {
        s1 = nseg;
        s0 = std::max(0, s1 - 2*npt);
      } else {
        s0 = ix - npt;
        s1 = ix + npt;
      }
      const Real z_x =
        acorn::eval_lagrange_poly_derivative(s1-s0+1, &xs[s0], &zs[s0], xs[ix]);
      Real nml[] = {-z_x, 1};
      mv2::normalize(nml);
      mv2::copy(nml, &nmls[2*ix]);
    }
  }

  // Init the splines.
  cs.resize(8*nseg);
  bc::init_from_nml(nseg, &ps[0], &nmls[0], &cs[0]);
}

void ExtrudedCubicSplineSurface
::init (const gallery::ZxFn::Shape zshape_, const Triangulation::CPtr& t_,
        const bool use_exact_normals, const int nml_recon_order) {
  zshape = zshape_;
  t = t_;

  // Collect (x,z).
  RealArray ps, xs, zs;
  const auto nvtx = t->get_nvtx();
  for (int i = 0; i < nvtx; ++i) {
    const auto vtx = t->get_vtx(i);
    Real z, zg;
    gallery::eval(zshape, vtx[0], z, zg);
    // We expect the vertices to have exact z(x) values for this zshape.
    assert(vtx[2] == z);
    if (vtx[1] == 0) {
      xs.push_back(vtx[0]);
      zs.push_back(vtx[2]);
      ps.push_back(vtx[0]);
      ps.push_back(vtx[2]);
      // We expect the first entries in the vtxs list to be the y=0 strip ...
      assert(ps.size() == 2*size_t(i+1));
      // ... with x increasing.
      assert(i == 0 || xs[i] > xs[i-1]);
    }
  }
  const int nseg = int(ps.size()/2)-1;
  assert(t->get_nvtx() % (nseg+1) == 0); // only ntri_per_rect=2 is allowed

  init_splines(zshape, use_exact_normals, nml_recon_order, nseg, xs, zs, ps, cs);

  // Fill triangle data.
  const auto ntri = t->get_ntri();
  lcl_tris.resize(ntri);
  tri_ctrs_gcs.resize(3*ntri);
  tri_lcs_at_ctrs.resize(9*ntri);
  ompparfor for (Idx ti = 0; ti < ntri; ++ti) {
    const auto tri = t->get_tri(ti);
    auto& lt = lcl_tris[ti];

    Real xy_ctr[2] = {0};
    for (int i = 0; i < 3; ++i)
      mv2::axpy(1.0/3.0, t->get_vtx(tri[i]), xy_ctr);

    lt.seg = std::lower_bound(xs.begin(), xs.end(), xy_ctr[0]) - xs.begin() - 1;
    assert(lt.seg >= 0 && lt.seg <= nseg);

    const auto c = &cs[8*lt.seg];
    assert(xy_ctr[0] >= c[0] && xy_ctr[0] <= c[6]);
    lt.ctr[0] = bc::calc_t_for_x(c, xy_ctr[0]);
    assert(lt.ctr[0] > 0 && lt.ctr[0] < 1);
    lt.ctr[1] = xy_ctr[1];

    Real xmin = 1e20, xmax = -1e20;
    for (int i = 0; i < 3; ++i) {
      const auto xy = t->get_vtx(tri[i]);
      auto v = &lt.vtxs()[2*i];
      xmin = std::min(xmin, xy[0]);
      xmax = std::max(xmax, xy[0]);
      assert(xy[0] >= c[0] && xy[0] <= c[6]);
      v[0] = bc::calc_t_for_x(c, xy[0]);
      assert(v[0] >= 0 && v[0] <= 1);
      v[1] = xy[1];
    }

    // We want to keep the polynomial scaled about the same in (x,y) so that
    // dist, L, acorn::integrals R0, and similar calculations are fine.
    //   Thus, map t in [0,1] to another local coordinate [xmin,xmax].
    //   If this is not done, in refinement tests, the x dim of a cell remains
    // O(1) while the y dim shrinks under refinement.
    lt.xmin = xmin;
    lt.xmax = xmax;
    for (int i = 0; i < 3; ++i) {
      auto v = &lt.vtxs()[2*i];
      v[0] = xmin + (xmax - xmin)*v[0];
    }
    lt.ctr[0] = xmin + (xmax - xmin)*lt.ctr[0];

    auto lcs = &tri_lcs_at_ctrs[9*ti];
    lcs[0] = 1; lcs[1] = 0; lcs[2] = 0;
    tri_position1(ti, lt.ctr, &tri_ctrs_gcs[3*ti], lcs);
  }
}

Size ExtrudedCubicSplineSurface::get_ntri () const { return t->get_ntri(); }

void ExtrudedCubicSplineSurface::tri_vtxs_uv (const Idx ti, Real t_lcs[6]) const {
  acorn::copy(6, lcl_tris[ti].vtxs(), t_lcs);
}

void ExtrudedCubicSplineSurface::tri_ctr_uv (const Idx ti, Real p[2]) const {
  asserti(*t, ti);
  mv2::copy(lcl_tris[ti].ctr, p);
}

void ExtrudedCubicSplineSurface::tri_ctr_xyz (const Idx ti, Real p[3]) const {
  asserti(*t, ti);
  mv3::copy(&tri_ctrs_gcs[3*ti], p);
}

void ExtrudedCubicSplineSurface::tri_ctr_lcs (const Idx ti, Real lcs[9]) const {
  asserti(*t, ti);
  acorn::copy(9, &tri_lcs_at_ctrs[9*ti], lcs);
}

bool ExtrudedCubicSplineSurface
::calc_tri_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const {
  const auto& lt = lcl_tris[ti];
  if (xyz[0] >= lt.xmin && xyz[0] <= lt.xmax)
    uv[0] = bc::calc_t_for_x(&cs[8*lt.seg], xyz[0]);
  else if (xyz[0] < lt.xmin) uv[0] = 0;
  else uv[0] = 1;
  uv[1] = xyz[1];
  return true;
}

void ExtrudedCubicSplineSurface
::tri_position (const Idx ti, const int n, CRPtr ty, RPtr p_gcs_, RPtr lcs_,
                Real* jacdet) const {
  const auto& lt = lcl_tris[ti];
  for (int i = 0; i < n; ++i) {
    const Real t = (ty[2*i] - lt.xmin)/(lt.xmax - lt.xmin);
    Real p[2], pt[2];
    bc::eval_p (&cs[8*lt.seg], t, p );
    bc::eval_pt(&cs[8*lt.seg], t, pt);
    Real p_gcs[3];
    p_gcs[0] = p[0];
    p_gcs[1] = ty[2*i+1];
    p_gcs[2] = p[1];
    mv3::copy(p_gcs, &p_gcs_[3*i]);
    // norm(cross((pt[0], 0, pt[1]), (0, 1, 0)), 2)
    if (jacdet) jacdet[i] = mv2::norm2(pt)/(lt.xmax - lt.xmin);
    if (lcs_) {
      // dz/dx = dz/dt / dx/dt
      Real zhat[] = {-pt[1]/pt[0], 0, 1};
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

int ExtrudedCubicSplineSurface::unittest () {
  int nerr = 0;
  const auto zshape = gallery::ZxFn::Shape::trig1;
  const int nseg = 11;
  RealArray xs(nseg+1), zs(nseg+1), ps(2*(nseg+1)), cs;

  for (int i = 0; i <= nseg; ++i) {
    xs[i] = ps[2*i] = Real(i)/nseg;
    Real unused;
    gallery::eval(zshape, xs[i], zs[i], unused);
    ps[2*i+1] = zs[i];
  }

  init_splines(zshape, true, -1, nseg, xs, zs, ps, cs);

  //todo

  if (nerr) printf("ExtrudedCubicSplineSurface::unittest failed\n");
  return nerr;
}

FlatElementSurface
::FlatElementSurface (const Triangulation::CPtr& t, const Real primary[3]) {
  init(t, primary);
}

static void calc_cross (const Triangulation& t, const Idx ti, RPtr nml) {
  const auto tri = t.get_tri(ti);
  const auto a = t.get_vtx(tri[0]);
  Real v1[3], v2[3];
  mv3::subtract(t.get_vtx(tri[1]), a, v1);
  mv3::subtract(t.get_vtx(tri[2]), a, v2);
  mv3::cross(v1, v2, nml);
}

void FlatElementSurface
::init (const Triangulation::CPtr& t_, const Real primary[3]) {
  t = t_;
  const auto& rt = *t;
  const auto ntri = rt.get_ntri();
  lcss.resize(ntri);
  vtx_lcss.resize(ntri);
  tri_ctrs_gcs.resize(3*ntri);
  ompparfor for (Idx ti = 0; ti < ntri; ++ti) {
    const auto tri = rt.get_tri(ti);
    auto& lcs = lcss[ti];
    auto xhat = lcs.xhat;
    auto yhat = lcs.yhat;
    auto zhat = lcs.zhat;
    auto& vtx_lcs = vtx_lcss[ti];
    auto* const ltri = vtx_lcs.vtxs();

    // LCS
    calc_cross(rt, ti, zhat);
    mv3::normalize(zhat);
    init_xhat_from_primary(zhat, primary, xhat);
    assert(mv3::norm22(xhat) != 0);
    mv3::cross(zhat, xhat, yhat);
    mv3::normalize(yhat); // for numerics

    // LCS vtxs
    for (int i = 0; i < 3; ++i) {
      const auto vi = tri[i];
      const auto vtx = rt.get_vtx(vi);
      ltri[2*i  ] = mv3::dot(xhat, vtx);
      ltri[2*i+1] = mv3::dot(yhat, vtx);
    }
    // other LCS triangle data
    vtx_lcs.z = mv3::dot(zhat, rt.get_vtx(tri[0]));
    t2d::calc_barycentric_matrix(ltri, vtx_lcs.b);
    vtx_lcs.ctr[0] = vtx_lcs.ctr[1] = 0;
    for (int i = 0; i < 3; ++i)
      mv2::axpy(1.0/3.0, &ltri[2*i], vtx_lcs.ctr);

    // GCS data
    tri_position1(ti, vtx_lcs.ctr, &tri_ctrs_gcs[3*ti]);
  }
}

Size FlatElementSurface::get_ntri () const { return t->get_ntri(); }

void FlatElementSurface::tri_vtxs_uv (const Idx ti, Real t_lcs[6]) const {
  asserti(*t, ti);
  acorn::copy(6, vtx_lcss[ti].vtxs(), t_lcs);
}

void FlatElementSurface::tri_ctr_uv (const Idx ti, Real p[2]) const {
  asserti(*t, ti);
  mv2::copy(vtx_lcss[ti].ctr, p);
}

void FlatElementSurface::tri_ctr_xyz (const Idx ti, Real p[3]) const {
  asserti(*t, ti);
  mv3::copy(&tri_ctrs_gcs[3*ti], p);
}

void FlatElementSurface::tri_ctr_lcs (const Idx ti, Real lcs[9]) const {
  asserti(*t, ti);
  acorn::copy(9, lcss[ti].data(), lcs);
}

bool FlatElementSurface
::calc_tri_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const {
  const auto& lcs = lcss[ti];
  uv[0] = mv3::dot(lcs.xhat, xyz);
  uv[1] = mv3::dot(lcs.yhat, xyz);
  return true;
}

void FlatElementSurface
::tri_position (const Idx ti, const int n, CRPtr xy_lcs, RPtr p_gcs, RPtr lcs,
                RPtr jacdet) const {
  const auto& vtx_lcs = vtx_lcss[ti];
  const auto& rt = *t;
  const auto tri = rt.get_tri(ti);
  for (int i = 0; i < n; ++i) {
    Real lam[3];
    t2d::xy_to_barycentric(vtx_lcs.v2, vtx_lcs.b, &xy_lcs[2*i], lam);
    mv3::sum3(lam[0], rt.get_vtx(tri[0]), lam[1], rt.get_vtx(tri[1]),
              lam[2], rt.get_vtx(tri[2]), &p_gcs[3*i]);
    if (jacdet) jacdet[i] = 1;
    if (lcs) tri_ctr_lcs(ti, &lcs[9*i]);
  }
}

int FlatElementSurface::unittest () {
  int nerr = 0;

  for (const int ntriper : {2, 4}) {
    gallery::ZxFn zxfn(gallery::ZxFn::Shape::trig1);
    zxfn.set_nx(5);
    const auto t = ConvTest::triangulate(zxfn, ntriper);
    const Real primary[] = {1,1,1};
    FlatElementSurface s(t, primary);
    const int ntri = t->get_ntri();
    for (int ti = 0; ti < ntri; ++ti) {
      const auto tri = t->get_tri(ti);
      for (int trial = 0; trial < 5; ++trial) {
        Real lam[] = {acorn::urand(), acorn::urand(), acorn::urand()};
        Real den = 0;
        for (int i = 0; i < 3; ++i) den += lam[i];
        mv3::divide(den, lam);
        Real p_gcs[3] = {0};
        for (int i = 0; i < 3; ++i)
          mv3::axpy(lam[i], t->get_vtx(tri[i]), p_gcs);
        Real p_lcs[2];
        s.calc_tri_lcs(ti, p_gcs, p_lcs);
        Real p_gcs1[3], lcs[9], jacdet;
        s.tri_position1(ti, p_lcs, p_gcs1, lcs, &jacdet);
        if (mv3::distance(p_gcs, p_gcs1) > 1e2*mv3::eps*mv3::norm2(p_gcs)) {
          acorn::prarr("lam",lam,3);
          pr(puf(mv3::distance(p_gcs, p_gcs1)) pu(mv3::norm2(p_gcs)));
          ++nerr;
        }
        if (jacdet != 1) ++nerr;
      }
    }
  }

  if (nerr) printf("FlatElementSurface::unittest failed\n");
  return nerr;
}


} // namespace convzx
} // namespace examples
} // namespace woodland
