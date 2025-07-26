#include "woodland/squirrel/ctzx_surface.hpp"
#include "woodland/squirrel/ctzx.hpp"
#include "woodland/squirrel/bezier_cubic.hpp"

#include "woodland/acorn/compose_lagrange_polynomial.hpp"
#include "woodland/acorn/plane_geometry.hpp"
#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/linalg.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace squirrel {
namespace ctzx {

typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;
typedef acorn::Triangle2D t2d;
typedef squirrel::BezierCubic<Real> bc;
typedef std::vector<Real> RealArray;
using mesh::Mesh;

static void asserti (const Triangulation& t, const Idx ti) {
  assert(ti >= 0);
  assert(ti < t.get_ncell());
}

namespace {

struct C2CubicSplineCoefs {
  int nseg, n;
  CRPtr xs, ys;
  Real y_x_beg, y_x_end;
  RPtr cs;
  RealArray A, cpolys;
  bool out = false;
  
  C2CubicSplineCoefs (const int nseg_, CRPtr xs_, CRPtr ys_,
                      const Real y_x_beg_, const Real y_x_end_,
                      RPtr cs_)
    : nseg(nseg_), xs(xs_), ys(ys_), y_x_beg(y_x_beg_), y_x_end(y_x_end_),
      cs(cs_)
  {
    run();
  }

private:
  int idx (const int r, const int c) const { return c*n + r; }

  Real pow (const Real x, const int e) const {
    assert(e >= 0 && e <= 3);
    return e == 0 ? 1 : e == 1 ? x : e == 2 ? x*x : x*x*x;
  }

  void fill_y (const int r, const int c0, const Real x) {
    for (int i = 0; i < 4; ++i) A[idx(r,c0+i)] = pow(x, i);
  }

  void fill_y_x (const int sign, const int r, const int c0, const Real x) {
    A[idx(r,c0)] = 0;
    for (int i = 1; i < 4; ++i) A[idx(r,c0+i)] = sign*i*pow(x, i-1);
  }

  void fill_y_xx (const int sign, const int r, const int c0, const Real x) {
    A[idx(r,c0)] = 0;
    A[idx(r,c0+1)] = 0;
    for (int i = 2; i < 4; ++i) A[idx(r,c0+i)] = sign*i*(i-1)*pow(x, i-2);
  }

  void fill () {
    n = 4*nseg;
    A.resize(n*n, 0);
    cpolys.resize(n);
    RPtr b = &cpolys[0];

    int r = 0;
    for (int k = 0; k < nseg; ++k) {
      const int c = 4*k;
      if (k == 0) {
        fill_y_x ( 1, r  , c, xs[k]);
        b[r] = y_x_beg;
        ++r;
      } else {
        fill_y_x (-1, r-2, c, xs[k]);
        fill_y_xx(-1, r-1, c, xs[k]);
      }
      fill_y(r  , c, xs[k  ]);
      fill_y(r+1, c, xs[k+1]);
      b[r  ] = ys[k  ];
      b[r+1] = ys[k+1];
      if (k == nseg-1) {
        assert(r+2 == n-1);
        fill_y_x (1, r+2, c, xs[k+1]);
        b[r+2] = y_x_end;
      } else {
        fill_y_x (1, r+2, c, xs[k+1]);
        fill_y_xx(1, r+3, c, xs[k+1]);
        r += 4;
        assert(r < n);
      }
    }

    if (out)
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          if (A[idx(i,j)] == 0) printf("          ");
          else printf(" %9.2e", A[idx(i,j)]);
        }
        if (b[i] == 0) printf(" |\n");
        else printf(" | %9.2e\n", b[i]);
      }
  }

  // At the end of this routine, a coefficient vector c gives
  //     f(x) = sum_i=0^3 c[i] x^i.
  void solve () {
    RealArray R(n*n), rwrk(n);
    std::vector<int> iwrk(n);
    acorn::linalg::qr_fac(n, n, &A[0], &R[0], &iwrk[0]);
    acorn::linalg::ls_slv(n, n, &A[0], &R[0], &iwrk[0], 1, &cpolys[0], &rwrk[0]);
    if (out) {
      for (int i = 0; i < n; ++i) printf(" %9.2e", cpolys[i]);
      printf("\n");
    }
  }

  // Evaluate in cpoly form.
  static Real eval (CRPtr c, const Real x) {
    return ((c[3]*x + c[2])*x + c[1])*x + c[0];
  }
  static Real eval_deriv (CRPtr c, const Real x) {
    return (3*c[3]*x + 2*c[2])*x + c[1];
  }

  /* At the end of this routine, a coefficient vector c is in the form used in
     BezierCubic:
         t in [0,1], r = 1-t
         t(x) = (x - xbeg)/(xend - xbeg)
         f(x) = p(t(x)) = c[0] r^3 + c[1] r^2 t + c[2] r t^2 + c[3] t^3.
     Now,
         f'(x) = p'(t) t_x = p'(t)/(xend - xbeg)
         p(0) = c[0]
         p(1) = c[3]
         p'(t) = -3 c[0] r^2 + r (r - 2 t) c[1] + t (2 r - t) c[2] + 3 c[3] t^2
         p'(0) = -3 c[0] + c[1] => c[1] = f'(xbeg) (xend - xbeg) + 3 c[0]
         p'(1) =  3 c[3] - c[2] => c[2] = 3 c[3] - f'(xend) (xend - xbeg).
     BezierCubic is for a 2D parameterized curve, so we actually have c[0:7]
     with 2-vectors for each of the four positions giving the coefficients for
     (x,y), as in BezierCubic.
  */ 
  void convert () {
    for (int is = 0; is < nseg; ++is) {
      const Real d = xs[is+1] - xs[is];
      CRPtr cpoly = &cpolys[4*is];
      RPtr c = &cs[8*is];
      // x. x_x(x) = 1.
      c[0] = xs[is];
      c[2] = d + 3*c[0];
      c[6] = xs[is+1];
      c[4] = 3*c[6] - d;
      // y(x).
      c[1] = eval(cpoly, xs[is]);
      c[3] = d*eval_deriv(cpoly, xs[is]) + 3*c[1];
      c[7] = eval(cpoly, xs[is+1]);
      c[5] = 3*c[7] - d*eval_deriv(cpoly, xs[is+1]);
    }
    if (out) {
      for (int i = 0; i < n; ++i) printf(" %9.2e", cs[2*i+1]);
      printf("\n");
    }
  }

  void run () {
    fill();
    solve();
    convert();
  }

public:
  static int unittest () {
    int nerr = 0;
    const auto f = [=] (const Real x) { return 0.3*x*x*x + 0.1*x*x - 0.5*x + 0.6; };
    const auto g = [=] (const Real x) { return 3*0.3*x*x + 2*0.1*x - 0.5; };
    const Real xs[] = {-0.75, -0.6, -0.3, 0.1, 0.9};
    const int nseg = sizeof(xs)/sizeof(*xs) - 1;
    Real ys[nseg+1];
    for (int i = 0; i < nseg+1; ++i) ys[i] = f(xs[i]);
    Real cs[8*nseg];
    C2CubicSplineCoefs(nseg, xs, ys, g(xs[0]), g(xs[nseg]), cs);
    for (int is = 0; is < 2; ++is) {
      for (int i = 0; i < 11; ++i) {
        Real p[2];
        const Real
          a = Real(i)/10,
          x = (1-a)*xs[is] + a*xs[is+1],
          y = f(x),
          yp = g(x),
          t = (x - xs[is])/(xs[is+1] - xs[is]);
        bc::eval_p(&cs[8*is], t, p);
        if (std::abs(p[0] - x) > 10*mv2::eps) {
          pr(puf(is) pu(t) pu(x) pu(p[0] - x));
          ++nerr;; 
        }
        if (std::abs(p[1] - y) > 1e3*mv2::eps) {
          pr(puf(is) pu(x) pu(t) pu(y) pu(p[1] - y));
          ++nerr;
        }
        bc::eval_pt(&cs[8*is], t, p);
        const Real d = xs[is+1] - xs[is];
        if (std::abs(p[0] - d) > 10*mv2::eps) {
          pr("x_t" pu(is) pu(t) pu(d) pu(p[0] - d));
          ++nerr;
        }
        if (std::abs(p[1] - yp*d) > 1e3*mv2::eps) {
          pr("y_t" pu(is) pu(x) pu(t) pu(yp) pu(p[1] - yp*d));
          ++nerr;
        }
      }
    }
    return 0;
  }
};

static void calc_c2_cubic_spline_coefs (
  // xs[0:nseg-1] is sorted ascending.
  const int nseg, CRPtr xs, CRPtr ys,
  // y'(x) at xs[0] and xs[nseg-1].
  const Real y_x_beg, const Real y_x_end,
  // cs has size 8*nseg. On output, its data are in the format BezierCubic uses.
  RPtr cs)
{
  C2CubicSplineCoefs(nseg, xs, ys, y_x_beg, y_x_end, cs);
}

} // namespace

ExtrudedCubicSplineSurface::ExtrudedCubicSplineSurface(
  const gallery::ZxFn::Shape zshape, const Triangulation::CPtr& t,
  const Recon recon)
{
  init(zshape, t, recon);
}

static void
init_splines (const gallery::ZxFn::Shape zshape,
              const ExtrudedCubicSplineSurface::Recon recon,
              const int nseg, const RealArray& xs, const RealArray& zs,
              const RealArray& ps, RealArray& cs) {
  using Recon = ExtrudedCubicSplineSurface::Recon;

  // Tangents.
  RealArray tans(2*(nseg+1));
  if (recon == Recon::c1_tan_exact) {
    for (int i = 0; i <= nseg; ++i) {
      Real z, zg;
      gallery::eval(zshape, ps[2*i], z, zg);
      Real tan[] = {1, zg};
      mv2::copy(tan, &tans[2*i]);
    }
  } else if (recon == Recon::c1_tan_2 || recon == Recon::c1_tan_4) {
    const int tan_recon_order = recon == Recon::c1_tan_2 ? 2 : 4;
    const int npt = tan_recon_order/2;
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
      Real tan[] = {1, z_x};
      mv2::copy(tan, &tans[2*ix]);
    }
  }

  // Lengths.
  RealArray ds(nseg);
  for (int i = 0; i < nseg; ++i)
    ds[i] = ps[2*(i+1)] - ps[2*i];

  // Init the splines.
  cs.resize(8*nseg);
  if (recon == Recon::c2) {
    const int n = std::min(nseg+1, 4), e0 = nseg-n+1;
    const Real
      z_x_beg = acorn::eval_lagrange_poly_derivative(n, &xs[0], &zs[0], xs[0]),
      z_x_end = acorn::eval_lagrange_poly_derivative(n, &xs[e0], &zs[e0], xs[nseg]);
    cs.resize(8*nseg);
    calc_c2_cubic_spline_coefs(nseg, &xs[0], &zs[0], z_x_beg, z_x_end, &cs[0]);
  } else {
    bc::init_from_tan(nseg, &ps[0], &tans[0], &ds[0], &cs[0]);
  }
}

void ExtrudedCubicSplineSurface
::init (const gallery::ZxFn::Shape zshape_, const Triangulation::CPtr& t_,
        const Recon recon) {
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

  init_splines(zshape, recon, nseg, xs, zs, ps, cs);

  // Fill triangle data.
  const auto ntri = t->get_ncell();
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
    cell_position1(ti, lt.ctr, &tri_ctrs_gcs[3*ti], lcs);
  }
}

Idx ExtrudedCubicSplineSurface::get_ncell () const { return t->get_ncell(); }

void ExtrudedCubicSplineSurface::cell_vtxs_uv (const Idx ti, Real t_lcs[6]) const {
  acorn::copy(6, lcl_tris[ti].vtxs(), t_lcs);
}

void ExtrudedCubicSplineSurface::cell_ctr_uv (const Idx ti, Real p[2]) const {
  asserti(*t, ti);
  mv2::copy(lcl_tris[ti].ctr, p);
}

void ExtrudedCubicSplineSurface::cell_ctr_xyz (const Idx ti, Real p[3]) const {
  asserti(*t, ti);
  mv3::copy(&tri_ctrs_gcs[3*ti], p);
}

void ExtrudedCubicSplineSurface::cell_ctr_lcs (const Idx ti, Real lcs[9]) const {
  asserti(*t, ti);
  acorn::copy(9, &tri_lcs_at_ctrs[9*ti], lcs);
}

bool ExtrudedCubicSplineSurface
::calc_cell_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const {
  const auto& lt = lcl_tris[ti];
  if (xyz[0] >= lt.xmin && xyz[0] <= lt.xmax)
    uv[0] = bc::calc_t_for_x(&cs[8*lt.seg], xyz[0]);
  else if (xyz[0] < lt.xmin) uv[0] = 0;
  else uv[0] = 1;
  uv[1] = xyz[1];
  return true;
}

bool ExtrudedCubicSplineSurface::supports_J () const { return true; }

void ExtrudedCubicSplineSurface
::cell_position (const Idx ti, const int n, CRPtr ty, RPtr p_gcs_, RPtr lcs_,
                 RPtr jacdet, RPtr jac) const {
  const auto& lt = lcl_tris[ti];
  for (int i = 0; i < n; ++i) {
    const Real D = lt.xmax - lt.xmin;
    const Real t = (ty[2*i] - lt.xmin)/D;
    Real p[2], pt[2];
    bc::eval_p (&cs[8*lt.seg], t, p );
    bc::eval_pt(&cs[8*lt.seg], t, pt);
    Real p_gcs[3];
    p_gcs[0] = p[0];
    p_gcs[1] = ty[2*i+1];
    p_gcs[2] = p[1];
    mv3::copy(p_gcs, &p_gcs_[3*i]);
    // J = [p'(t)[0]/D 0; 0 1; p'(t)[1]/D 0]
    // sqrt(det(J'J)) = norm(pt, 2)/D
    if (jacdet) jacdet[i] = mv2::norm2(pt)/D;
    if (jac) {
      auto J = &jac[6*i];
      for (int k = 0; k < 6; ++k) J[k] = 0;
      J[0] = pt[0]/D; J[3] = 1; J[4] = pt[1]/D;
    }
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

  nerr += C2CubicSplineCoefs::unittest();

  const auto zshape = gallery::ZxFn::Shape::trig1;
  const int nseg = 11;
  RealArray xs(nseg+1), zs(nseg+1), ps(2*(nseg+1)), cs;

  for (int i = 0; i <= nseg; ++i) {
    xs[i] = ps[2*i] = Real(i)/nseg;
    Real unused;
    gallery::eval(zshape, xs[i], zs[i], unused);
    ps[2*i+1] = zs[i];
  }

  init_splines(zshape, Recon::c1_tan_exact, nseg, xs, zs, ps, cs);

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
  const auto ntri = rt.get_ncell();
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
    cell_position1(ti, vtx_lcs.ctr, &tri_ctrs_gcs[3*ti]);
  }
}

Idx FlatElementSurface::get_ncell () const { return t->get_ncell(); }

void FlatElementSurface::cell_vtxs_uv (const Idx ti, Real t_lcs[6]) const {
  asserti(*t, ti);
  acorn::copy(6, vtx_lcss[ti].vtxs(), t_lcs);
}

void FlatElementSurface::cell_ctr_uv (const Idx ti, Real p[2]) const {
  asserti(*t, ti);
  mv2::copy(vtx_lcss[ti].ctr, p);
}

void FlatElementSurface::cell_ctr_xyz (const Idx ti, Real p[3]) const {
  asserti(*t, ti);
  mv3::copy(&tri_ctrs_gcs[3*ti], p);
}

void FlatElementSurface::cell_ctr_lcs (const Idx ti, Real lcs[9]) const {
  asserti(*t, ti);
  acorn::copy(9, lcss[ti].data(), lcs);
}

bool FlatElementSurface
::calc_cell_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const {
  const auto& lcs = lcss[ti];
  uv[0] = mv3::dot(lcs.xhat, xyz);
  uv[1] = mv3::dot(lcs.yhat, xyz);
  return true;
}

bool FlatElementSurface::supports_J () const { return true; }

static void calc_flat_elem_jacobian (
  // Triangle GCS vertices.
  const Real v0[3], const Real v1[3], const Real v2[3],
  // Triangle2D::calc_barycentric_matrix data.
  const Real b[4],
  // Jacobian of p_gcs w.r.t. p_lcs.
  Real J[6])
{
  // lam = [b (p_lcs - v2_lcs);
  //        1 - e'B (p_lcs - v2_lcs)]
  // p_gcs = [v0,v1,v2] lam
  // => J = [v0,v1,v2] [B; -e'B]
  for (int k = 0; k < 2; ++k) {
    Real c[3];
    const int os = 2*k;
    mv3::sum3(b[os], v0, b[os+1], v1, -(b[os+0] + b[os+1]), v2, c);
    for (int d = 0; d < 3; ++d) J[2*d+k] = c[d];
  }
}

void FlatElementSurface
::cell_position (const Idx ti, const int n, CRPtr xy_lcs, RPtr p_gcs, RPtr lcs,
                 RPtr jacdet, RPtr J) const {
  const auto& vtx_lcs = vtx_lcss[ti];
  const auto& rt = *t;
  const auto tri = rt.get_tri(ti);
  for (int i = 0; i < n; ++i) {
    Real lam[3];
    t2d::xy_to_barycentric(vtx_lcs.v2, vtx_lcs.b, &xy_lcs[2*i], lam);
    mv3::sum3(lam[0], rt.get_vtx(tri[0]), lam[1], rt.get_vtx(tri[1]),
              lam[2], rt.get_vtx(tri[2]), &p_gcs[3*i]);
    if (jacdet) jacdet[i] = 1;
    if (J)
      calc_flat_elem_jacobian(
        rt.get_vtx(tri[0]), rt.get_vtx(tri[1]), rt.get_vtx(tri[2]),
        vtx_lcs.b, J);
    if (lcs) cell_ctr_lcs(ti, &lcs[9*i]);
  }
}

static int test_flat_elem_jacobian () {
  int ne = 0;
  const Real v_lcs[3][2] = {{0,0}, {1,0}, {0,1}};
  const Real v_gcs[3][3] = {{0.1, -0.1, 0.8}, {1.5, -0.2, -0.2},
                            {-0.1, 0.8, -1.6}};
  const Real p_lcs[2] = {0.5, 0.3};
  Real b[4];
  t2d::calc_barycentric_matrix(v_lcs[0], v_lcs[1], v_lcs[2], b);
  Real J[6];
  calc_flat_elem_jacobian(v_gcs[0], v_gcs[1], v_gcs[2], b, J);
  Real J_fd[6];
  const Real delta = 1e-6;
  for (int i = 0; i < 2; ++i) {
    Real gs[2][3];
    for (int j = 0; j < 2; ++j) {
      Real p[2];
      mv2::copy(p_lcs, p);
      const int sign = j == 0 ? -1 : 1;
      p[i] += sign*delta;
      Real lam[3];
      t2d::xy_to_barycentric(v_lcs[2], b, p, lam);
      mv3::sum3(lam[0], v_gcs[0], lam[1], v_gcs[1], lam[2], v_gcs[2], gs[j]);
    }
    Real d[3];
    mv3::subtract(gs[1], gs[0], d);
    for (int j = 0; j < 3; ++j) J_fd[2*j+i] = (gs[1][j] - gs[0][j])/(2*delta);
  }
  Real sum = 0;
  for (int i = 0; i < 6; ++i) sum += std::abs(J[i] - J_fd[i]);
  if (sum > 1e-9) {
    ++ne;
    prc(sum);
    acorn::prarr("J",J,6);
    acorn::prarr("J_fd",J_fd,6);
  }
  return ne;
}

int FlatElementSurface::unittest () {
  int nerr = 0, ne;

  for (const int ntriper : {2, 4}) {
    gallery::ZxFn zxfn(gallery::ZxFn::Shape::trig1);
    zxfn.set_nx(5);
    const auto t = ConvTest::triangulate(zxfn, ntriper);
    const Real primary[] = {1,1,1};
    FlatElementSurface s(t, primary);
    const int ntri = t->get_ncell();
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
        s.calc_cell_lcs(ti, p_gcs, p_lcs);
        Real p_gcs1[3], lcs[9], jacdet;
        s.cell_position1(ti, p_lcs, p_gcs1, lcs, &jacdet);
        if (mv3::distance(p_gcs, p_gcs1) > 1e2*mv3::eps*mv3::norm2(p_gcs)) {
          acorn::prarr("lam",lam,3);
          pr(puf(mv3::distance(p_gcs, p_gcs1)) pu(mv3::norm2(p_gcs)));
          ++nerr;
        }
        if (jacdet != 1) ++nerr;
      }
    }
  }

  rununittest(test_flat_elem_jacobian);
  
  return nerr;
}

} // namespace ctzx
} // namespace squirrel
} // namespace woodland
