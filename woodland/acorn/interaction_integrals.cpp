#include "woodland/acorn/interaction_integrals.hpp"

#include "woodland/acorn/compose_quadrature.hpp"
#include "woodland/acorn/compose_triquad.hpp"

#include "woodland/acorn/vectorization.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/hfp.hpp"
#include "woodland/acorn/dbg.hpp"

#include <cstdio>

#ifdef WOODLAND_ACORN_FIGURE
# include <memory>
#endif

namespace woodland {
namespace acorn {
namespace integrals {

#ifdef WOODLAND_ACORN_FIGURE
namespace {

struct Figure {
  typedef std::shared_ptr<Figure> Ptr;

  Figure (const std::string fname = "fig_acorn_integrals.py") {
    fid = fopen(fname.c_str(), "w");
    assert(fid);
    fprintf(fid, "import numpy as npy\n");
  }

  ~Figure () { if (fid) fclose(fid); }

  void calc_hfp_1 (const Polygon& p, const Pt cc, const CallerIntegrands& f,
                   const int ie_) {
    ie = ie_;
    if (ie != 0) return;
    {
      Sprinter s;
      s.add("cc = npy.array([%1.6e,%1.6e])\n", cc[0], cc[1]);
      s.add("poly = npy.array([");
      for (int i = 0; i < p.n; ++i)
        s.add("[%1.6e,%1.6e],", p.xys[2*i], p.xys[2*i+1]);
      s.add("])\n");
      s.out(fid, false);
    }
    {
      Real xmin = 1e20, xmax = -xmin, ymin = xmin, ymax = xmax;
      for (int i = 0; i < p.n; ++i) {
        xmin = std::min(xmin, p.xys[2*i]);
        xmax = std::max(xmax, p.xys[2*i]);
        ymin = std::min(ymin, p.xys[2*i+1]);
        ymax = std::max(ymax, p.xys[2*i+1]);
      }
      const Real e = 0.02, dx = xmax - xmin, dy = ymax - ymin;
      xmin -= e*dx; xmax += e*dx;
      ymin -= e*dy; ymax += e*dy;
      fprintf(fid, "extent = [%1.6e,%1.6e,%1.6e,%1.6e]\n",
              xmin, xmax, ymin, ymax);
      const int n = 128;
      fprintf(fid, "f = npy.array([");
      for (int i = 0; i < n; ++i) {
        fprintf(fid, "[");
        const Real a = (i + 0.5)/n, y = (1-a)*ymin + a*ymax;
        for (int j = 0; j < n; ++j) {
          const Real b = (j + 0.5)/n, x = (1-b)*xmin + b*xmax;
          const Pt pt = {x, y};
          if (not plane::is_point_inside_convex_polygon(p, pt)) {
            fprintf(fid, "npy.nan,");
            continue;
          }
          Real integrand[6];
          f.eval(1, pt, integrand);
          fprintf(fid, "%1.6e,", integrand[iidx]);
        }
        fprintf(fid, "],\n");
      }
      fprintf(fid, "])\n");
    }
  }

  void circular_sector_term_1 (
    const Options& o, const Real th0, const Real th1, const Pt cc, const Real R,
    CRPtr radii, CRPtr radial_integrand, const int integrand_idx)
  {
    if (integrand_idx != iidx) return;
    Sprinter s;
    if (ie == 0) s.add("th0 = []; th1 = []; R = []; radii = []\n");
    s.add("th0.append(%1.6e); th1.append(%1.6e); R.append(%1.6e)\n", th0, th1, R);
    if (ie == iie) {
      s.add("ie_circsec = %d\n", ie);
      s.add("radii = npy.array([");
      for (int i = 0; i < o.np_radial; ++i)
        s.add("%1.6e,", radii[i]);
      s.add("])\n");
      s.add("radial_integrand = npy.array([");
      for (int i = 0; i < o.np_radial; ++i)
        s.add("%1.6e,", radial_integrand[i]);
      s.add("])\n");
    }
    s.out(fid, false);
  }

private:
  FILE* fid = nullptr;
  const int iie = 1, iidx = 2;
  int ie;
};

Figure::Ptr g_fig;

} // namespace

void fig_init () { g_fig = std::make_shared<Figure>(); }
void fig_fin  () { g_fig = nullptr; }

#else

void fig_init () {}
void fig_fin  () {}

#endif // WOODLAND_ACORN_FIGURE

typedef Matvec2d<Real> mv2;

struct EvalAccumulator {
  enum : int { capacity = RealPack::n };

  EvalAccumulator (const CallerIntegrands& f_)
    : f(f_), nint(f.nintegrands()), n(0)
  {
    for (int i = 0; i < nint; ++i) integral[i] = 0;
  }
  
  void accum (const Pt pt, const Real wt) {
    pts[2*n  ] = pt[0];
    pts[2*n+1] = pt[1];
    wts[n] = wt;
    ++n;
    if (n == capacity) flush();
  }

  const Real* get_integral () {
    flush();
    return integral;
  }

private:
  const CallerIntegrands& f;
  const int nint;
  int n;
  Real pts[2*capacity], wts[capacity];
  Real integrand[CallerIntegrands::max_n_integrand*capacity];
  Real integral[CallerIntegrands::max_n_integrand];

  void flush () {
    if ( ! n) return;
#if 1
    f.eval(n, pts, integrand);
#else
    for (int k = 0; k < n; ++k)
      f.eval(1, &pts[2*k], &integrand[nint*k]);
#endif
    for (int k = 0; k < n; ++k)
      for (int i = 0; i < nint; ++i)
        integral[i] += wts[k]*integrand[nint*k+i];
    n = 0;
  }
};

// Angles associated with v0, v1, w.r.t. a circle center at (0,0).
static void calc_thetas(const Pt v0, const Pt v1, Real& th0, Real& th1) {
  th0 = std::atan2(v0[1], v0[0]);
  th1 = std::atan2(v1[1], v1[0]);
  if (std::abs(th1 - th0) > M_PI) {
    if (th1 > th0) th0 += 2*M_PI;
    else th1 += 2*M_PI;
  }
  assert(std::abs(th1 - th0) <= M_PI);
}

// Input values are overwritten.
static void calc_hfp_circle_arc_integral_times_rsquared (
  const Options& o, const Pt cc, const Real r, const Real th0, const Real th1,
  const CallerIntegrands& f, RPtr integral)
{
  assert(is_gll_supported(o.np_angular));
  const auto* const qth = get_x_gll(o.np_angular);
  const auto* const wth = get_w_gll(o.np_angular);
  EvalAccumulator ea(f);
  // Factor is as follows:
  //     (1/2 for [-1,1] quadrature)
  //   * (|th1 - th0| r for circle circumference)
  //   * (r^2 for singularity)
  const auto fac = 0.5*std::abs(th1 - th0)*cube(r);
  for (int igll = 0; igll < o.np_angular; ++igll) {
    const auto a = (1 + qth[igll])/2;
    const auto theta = (1-a)*th0 + a*th1;
    Pt pt{r*std::cos(theta), r*std::sin(theta)};
    mv2::axpy(1, cc, pt);
    ea.accum(pt, wth[igll]*fac);
  }
  const auto* eai = ea.get_integral();
  for (int i = 0, n = f.nintegrands(); i < n; ++i)
    integral[i] = eai[i];
}

struct RadialP : public hfp::CallerPolyP {
  RadialP(const int n, CRPtr xs, CRPtr ys);
  bool in_support(const Real x) const override;
};

RadialP::RadialP (const int n, CRPtr xs, CRPtr ys)
  : hfp::CallerPolyP(n, xs, ys)
{
  assert(n > 0);
  assert(xs[0] >= 0);
}

bool RadialP::in_support (const Real x) const {
  return x >= 0 && x <= xs[n-1];
}

static void
calc_hfp_circular_sector_term (const Options& o, const Pt cc, const Real R,
                               const Real th0, const Real th1,
                               const CallerIntegrands& f, RPtr hfps) {
  const auto R_min = f.permitted_R_min(R);
  assert(R_min < R);
  const int nintegrands = f.nintegrands();

  Real radii[Quadrature::max_nq];
  Real radial_integrand[CallerIntegrands::max_n_integrand][Quadrature::max_nq];

  for (int i = 0; i < nintegrands; ++i)
    for (int j = 0; j < o.np_radial; ++j)
      radial_integrand[i][j] = 0;
  
  { // Form the integrand in r from R_min to R.
    assert(is_gll_supported(o.np_radial));
    const auto* const qr = get_x_gll(o.np_radial);
    for (int j = 0; j < o.np_radial; ++j) {
      const auto a_gl = j == 0 ? 0 : j == o.np_radial-1 ? 1 : (1 + qr[j])/2;
      const auto r = a_gl*R + (1 - a_gl)*R_min;
      radii[j] = r;
      Real integral[CallerIntegrands::max_n_integrand];
      calc_hfp_circle_arc_integral_times_rsquared(o, cc, r, th0, th1, f,
                                                  integral);
      for (int i = 0; i < nintegrands; ++i)
        radial_integrand[i][j] = integral[i];
    }
  }

  { // Calculate hfp int_0^R ...
    hfp::Options o1d;
    o1d.gl_np = o.np_radial;
    for (int i = 0; i < nintegrands; ++i) {
      RadialP p(o.np_radial, radii, radial_integrand[i]);
      hfps[i] += hfp::calc_hfp(o1d, p, 0, R, 0);
#ifdef WOODLAND_ACORN_FIGURE
      if (g_fig)
        g_fig->circular_sector_term_1(o, th0, th1, cc, R, radii,
                                      radial_integrand[i], i);
#endif
    }
  }
}

// A radtricirc is a triangle having the following three edge types:
//   a circle-arc edge from p_base to p_top, for a circle having radius Rc;
//   a radial line segment from p_base to v = (Rv/Rc) p_base;
//   a line segment (v, p_top).
//
// This routine approximates an integral over a radtricirc using a
// tensor-product quadrature. One direction is orthogonal to the radial line
// segment. The other is parallel to the radial line segment. The motivation of
// this approach is to minimize the Jacobian determinant's deviation from
// polynomial.
//   The integral is as follows, using slightly different symbols than in the
// routine for convenience. In the discussion, picture the shape as oriented
// with the arc on the left, the radial base on the bottom and horizontal, the
// other edge on the right, and the circle centered at 0.
//   Let theta = t1 - t0.
//   The vertices in cartesian coords (x,y) are these:
//     bottom left  (Rc,0)                                     p_base
//     bottom right (Rv,0)                                     v
//     top          (Rc cos theta, Rc sin theta) := (xt, yt)   p_top
//   Given y, the point on the right edge is
//     alpha(y) := y / (Rc sin theta)
//     x_e(y) := (1 - alpha) Rv + alpha Rc cos theta
//             = Rv + (Rc cos theta - Rv) alpha
//             = Rv + y (xt - Rv)/yt,
// and the point on the arc is
//     x_c(y) := sqrt(Rc^2 - y^2).
//   The integral is then
//     int_0^yt int_[sqrt(Rc^2 - y^2)]^[Rv + y (xt - Rv)/yt] f(x,y) dx dy
// and the tensor quadrature is over x and y.
//
// In this routine, v and p_top are w.r.t. a center at 0. cc is used to
// translate to the caller's coordinates.
static void
calc_integral_radtricirc (const Options& o, const Pt cc, const Real R,
                          const Pt v, const Pt p_top, const CallerIntegrands& f,
                          RPtr integrals) {
  const int nint = f.nintegrands();

  Quadrature iq1d(o.np_radial, Quadrature::gl);
  const Real* iqx, * iqw;
  iq1d.get_xw(iqx, iqw);
  const int inq = iq1d.nq;
  Quadrature jq1d(o.np_angular, Quadrature::gl);
  const Real* jqx, * jqw;
  jq1d.get_xw(jqx, jqw);
  const int jnq = jq1d.nq;

  // Intersection of (zero, v) with circle.
  Pt p_base;
  const auto Rv = mv2::norm2(v);
  //todo If Rv is very close to R, then this integral doesn't contribute very
  // much. In the future, I'd like to filter out such cases. For now, assert
  // here so I know when this case is occurring at the machine precision level.
  assert(Rv > R);
  mv2::copy(v, p_base);
  mv2::scale(R/Rv, p_base);

  EvalAccumulator ea(f);
  Real t0, t1;
  calc_thetas(p_base, p_top, t0, t1);
  const auto ct0 = std::cos(t0), st0 = std::sin(t0);
  const auto theta = t1 - t0;
  assert(std::abs(theta) < M_PI/2);
  const auto xt = R*std::cos(theta);
  const auto yt = R*std::sin(theta);
  for (int iq = 0; iq < inq; ++iq) {
    const auto u = (1 + iqx[iq])/2;
    const auto y = u*yt;
    const auto xb = std::sqrt(square(R) - square(y));
    const auto xe = (1 - u)*Rv + u*xt;
    const auto i_factor = 0.25*iqw[iq]*std::abs(yt);
    for (int jq = 0; jq < jnq; ++jq) {
      const auto v = (1 + jqx[jq])/2;
      const auto x = (1 - v)*xb + v*xe;
      Pt xy_g{x, y};
      rotate_vector_2(ct0, st0, xy_g);
      mv2::axpy(1, cc, xy_g);
      ea.accum(xy_g, i_factor*jqw[jq]*(xe - xb));
    }
  }

  const auto* eai = ea.get_integral();
  for (int i = 0; i < nint; ++i)
    integrals[i] += eai[i];
}

// pt is the point at which the maximal circle touches the line segment (v1,v2).
// pt_at is -1 if pt is at v1, 0 if between v1 and v2, 1 if at v2.
static void decompose_tri (const Pt v1, const Pt v2,
                           Real& R02, Pt pt, int& pt_at) {
  Real a;
  const Pt zero = {0};
  plane::project_p_onto_line(zero, v1, v2, a, pt);
  pt_at = a <= -0.5 ? -1 : (a >= 0.5 ? 1 : 0);
  if (pt_at != 0) mv2::copy(pt_at == -1 ? v1 : v2, pt);
  R02 = mv2::norm22(pt);
}

bool calc_hfp (const Options& o, const Polygon& p, const Pt cc,
               const CallerIntegrands& f, RPtr hfps) {
  const int nint = f.nintegrands();
  Real ghfps[CallerIntegrands::max_n_integrand] = {0};
  Real area = 0; // signed area of p
  for (int ie = 0; ie < p.n; ++ie) {
    // Compute decomposition of triangle (cc, poly vtx ie, poly vtx (ie+1)%n)
    // into    
    //   1 circular sector with radius R0 and angles (th0,th1) and
    //   1 or 2 radtricircs with top point pt_top.
    Pt v0, v1; // centered at (0,0)
    mv2::subtract(&p.xys[2*ie], cc, v0);
    mv2::subtract(&p.xys[2*((ie+1) % p.n)], cc, v1);
    Real R02;
    Pt pt_top; // top of each radtricirc
    int pt_top_at;
    decompose_tri(v0, v1, R02, pt_top, pt_top_at);
    const auto R0 = std::sqrt(R02);
    Real th0, th1;
    calc_thetas(v0, v1, th0, th1);
#ifdef WOODLAND_ACORN_FIGURE
    if (g_fig) g_fig->calc_hfp_1(p, cc, f, ie);
#endif
    Real lhfps[CallerIntegrands::max_n_integrand] = {0};
    // H.f.p. over the circular sector.
    calc_hfp_circular_sector_term(o, cc, R0, th0, th1, f, lhfps);
    // Integrals over the radtricirc(s).
    if (pt_top_at >= 0)
      calc_integral_radtricirc(o, cc, R0, v0, pt_top, f, lhfps);
    if (pt_top_at <= 0)
      calc_integral_radtricirc(o, cc, R0, v1, pt_top, f, lhfps);
    // Multiply by -1 if the triangle is CW, 1 if CCW. At the end, we compensate
    // for the signed area of p so that p can be CW or CCW.
    //   The result of accounting for sign is that cc can be outside of p. If cc
    // is near p, calc_hfp is more accurate than calc_integral even though the
    // integral is proper.
    Pt zero = {0};
    const auto signed_area = Triangle2D::calc_signed_area(zero, v0, v1);
    area += signed_area;
    const auto sign = signed_area >= 0 ? 1 : -1;
    for (int i = 0; i < nint; ++i) ghfps[i] += sign*lhfps[i];
  }
  const auto sign = area >= 0 ? 1 : -1;
  for (int i = 0; i < nint; ++i) hfps[i] += sign*ghfps[i];
  return true;
}

static void
calc_integral (const int tq_order, const Pt a, const Pt b, const Pt c,
               const CallerIntegrands& f, RPtr integral) {
  const int nint = f.nintegrands();

  const Real* qx, * qw;
  int nq;
  TriangleQuadrature::get_coef(tq_order, qx, qw, nq);

  const Real tri_area = Triangle2D::calc_signed_area(a, b, c);

  EvalAccumulator ea(f);
  for (int iq = 0; iq < nq; ++iq) {
    Pt x;
    Triangle2D::barycentric_to_xy(a, b, c, &qx[3*iq], x);
    ea.accum(x, std::abs(tri_area) * qw[iq]);
  }
  const auto* eai = ea.get_integral();
  for (int i = 0; i < nint; ++i)
    integral[i] += eai[i];
}

bool calc_integral (const Polygon& p, const CallerIntegrands& f, RPtr integral,
                    const int tq_order) {
  assert(TriangleQuadrature::is_order_supported(tq_order));
  Pt centroid = {0};
  for (int i = 0; i < p.n; ++i)
    mv2::axpy(1, &p.xys[2*i], centroid);
  for (int d = 0; d < 2; ++d) centroid[d] /= p.n;
  for (int i = 0; i < p.n; ++i)
    calc_integral(tq_order, centroid, &p.xys[2*i], &p.xys[2*((i+1)%p.n)], f,
                  integral);
  return true;
}

bool calc_integral_tensor_quadrature (
  const Options& o, const Polygon& p, const CallerIntegrands& f,
  const Pt anchor, const bool nearest_bdy_pt_to_anchor,
  RPtr integral)
{
  const int nint = f.nintegrands();

  Pt p_common;
  if (nearest_bdy_pt_to_anchor) {
    int on_vertex, on_edge;
    plane::calc_nearest_point_on_polygon_boundary_to_point(
      p, anchor, p_common, on_vertex, on_edge);
  } else {
    mv2::copy(anchor, p_common);
  }
  
  Quadrature iq1d(o.np_radial, Quadrature::gl);
  const Real* iqx, * iqw;
  iq1d.get_xw(iqx, iqw);
  const int inq = iq1d.nq;
  Quadrature jq1d(o.np_angular, Quadrature::gl);
  const Real* jqx, * jqw;
  jq1d.get_xw(jqx, jqw);
  const int jnq = jq1d.nq;

  Real p_area = 0;
  for (int i = 1; i < p.n-1; ++i)
    p_area += std::abs(Triangle2D::calc_signed_area(
                         &p.xys[0], &p.xys[2*i], &p.xys[2*(i+1)]));

  EvalAccumulator ea(f);
  for (int tri = 0; tri < p.n; ++tri) {
    const auto v1 = &p.xys[2*tri], v2 = &p.xys[2*((tri+1)%p.n)];
    const Real area = std::abs(Triangle2D::calc_signed_area(p_common, v1, v2));
    if (area < o.relative_area_tol*p_area) continue;
    for (int i = 0; i < inq; ++i) {
      const Real b = (1 + iqx[i])/2;
      Pt j1;
      mv2::sum2(b, v1, -b, v2, j1);
      for (int j = 0; j < jnq; ++j) {
        const Real a = (1 + jqx[j])/2;
        Pt x, j2;
        // Degenerate side of quad for v = p_common:
        //   (1-a)*(1-b) v + a*(1-b) v = (1-b) v.
        mv2::sum3(1-b, p_common, a*b, v1, (1-a)*b, v2, x);
        mv2::sum3(-1, p_common, a, v1, 1-a, v2, j2);
        const Real jacdet = std::abs(j1[0]*j2[1] - j1[1]*j2[0]);
        ea.accum(x, jacdet*iqw[i]*jqw[j]/4);
      }
    }
  }
  const auto* eai = ea.get_integral();
  for (int i = 0; i < nint; ++i)
    integral[i] += eai[i];
  
  return true;
}

namespace {

struct One : public CallerIntegrands {
  int nintegrands () const override { return 1; }
  void eval (const int n, CRPtr p, RPtr integrand) const override {
    for (int i = 0; i < n; ++i) integrand[i] = 1;
  }
};

int test_area (const Polygon& p, const Pt cc) {
  int ne = 0;
  Real area_true = 0;
  {
    Real pin[2] = {0};
    for (int i = 0; i < p.n; ++i)
      mv2::axpy(1, &p.xys[2*i], pin);
    mv2::divide(p.n, pin);
    for (int i = 0; i < p.n; ++i)
      area_true += Triangle2D::calc_signed_area(pin, &p.xys[2*i],
                                                &p.xys[2*((i+1) % p.n)]);
  }
  Options o;
  o.np_angular = 9;
  Real area = 0;
  if ( ! calc_hfp(o, p, cc, One(), &area)) ++ne;
  if (std::abs(area - area_true) >= 1e3*area_true*mv2::eps) {
    printf("integrals::test_area (hfp) %1.3f %1.2e\n",
           area, std::abs(area - area_true)/area_true);
    ++ne;
  }
  area = 0;
  if ( ! calc_integral(p, One(), &area)) ++ne;
  if (std::abs(area - area_true) >= 1e1*area_true*mv2::eps) {
    printf("otherint::test_area (integral) %1.3f %1.2e\n",
           area, std::abs(area - area_true)/area_true);
    ++ne;
  }
  if (ne) printf("integrals::test_area failed %d\n", ne);
  return ne;
}

int test_area () {
  int ne = 0;
  {
    const Real xys[] = {-1,-1,1,-1,1,1,-1,1};
    const Polygon p(xys, sizeof(xys)/sizeof(*xys)/2);
    {
      const Real cc[2] = {0};
      ne += test_area(p, cc);
    }
    {
      const Real cc[] = {0.3, 0.05};
      ne += test_area(p, cc);
    }
  }
  {
    const Real xys[] = {-0.5,-0.6, 0.2,-0.5, 0.8,-0.2, 1.1,0.1, 0.5,0.9, -0.5,-0.5};
    const Real cc[] = {0.3, 0.05};
    ne += test_area(Polygon(xys, sizeof(xys)/sizeof(*xys)/2), cc);
  }
  return ne;
}

/* Tests:

   Entry 1: 1/r^3
     hfp int_0^R int_0^{2 pi} 1/r^3 r dt dr
       = hfp int_0^R int_0^{2 pi} 1/r^2 dt dr
       = hfp int_0^R 2 pi/r^2 dr
       = -2 pi/R.

     hfp int_-2^1 int_-1^2 1/sqrt(x^2 + y^2)^3 dx dy
       = sqrt(2)/2

   Entry 2: 1/r^2
     hfp int_0^R int_0^{2 pi} 1/r^2 r dt dr
       = hfp int_0^R int_0^{2 pi} 1/r dt dr
       = hfp int_0^R 2 pi/r dr
       = 2 pi log R.

   Entry 3. 1

   Entry 4. x*y
 */

struct SingularTestIntegrands : public CallerIntegrands {
  enum : int { nint = 4 };
  int nintegrands () const override { return nint; }
  Real permitted_R_min (const Real /*R_max*/) const override {
    return std::cbrt(std::numeric_limits<Real>::min());
  }
  void eval (const int n, CRPtr p, RPtr integrand) const override {
    for (int i = 0; i < n; ++i) {
      const auto d2 = mv2::norm22(&p[2*i]);
      const auto os = 4*i;
      integrand[os+0] = 1/(d2*std::sqrt(d2));
      integrand[os+1] = 1/d2;
      integrand[os+2] = 1;
      integrand[os+3] = p[2*i]*p[2*i+1];
    }
  }
  Real true_hfp_circle (const int idx, const Real R) const {
    assert(idx >= 0 && idx < nint);
    assert(R > 0);
    switch (idx) {
    case 0: return -2*M_PI/R;
    case 1: return 2*M_PI*std::log(R);
    case 2: return M_PI*square(R);
    case 3: return 0;
    default:
      assert(false);
      return 0;
    }
  }
};

int test_calc_hfp_circle () {
  int ne = 0;
  { // circle term
    const Pt cc = {0};
    const Real R = 0.5;
    Real hfps[SingularTestIntegrands::nint] = {0};
    SingularTestIntegrands f;
    calc_hfp_circular_sector_term(Options(), cc, R, 0, 2*M_PI, f, hfps);
    const int n = f.nintegrands();
    for (int i = 0; i < n; ++i) {
      const auto t = f.true_hfp_circle(i, R);
      if (std::abs(t - hfps[i]) > 5e4*std::max(1.0, std::abs(t))*mv2::eps) {
        ++ne;
        printf("integrals::test_calc_hfp_circle: %d %f %e %e %e\n",
               i, R, t, hfps[i],
               std::abs(t - hfps[i])/std::max(1.0, std::abs(t)));
      }
    }
  }
  return ne;
}

struct NonsingularTestIntegrands : public CallerIntegrands {
  enum : int { nint = 1 };
  Pt cc;
  NonsingularTestIntegrands () { cc[0] = 0.2; cc[1] = -0.3; }
  int nintegrands () const override { return nint; }
  void eval (const int n, CRPtr p, RPtr integrand) const override {
    for (int i = 0; i < n; ++i) {
      const auto x = p[2*i], y = p[2*i+1];
      integrand[i] = square(x - cc[0])*cube(y - cc[1]);
    }
  };
};

int test_integral () {
  int ne = 0;
  // Over [-2,1]x[-1,2].
  NonsingularTestIntegrands f;
  Real xys[] = {-2,-1, 1,-1, 1,2, -2,2};
  for (int i = 0; i < 4; ++i) mv2::axpy(1, f.cc, &xys[2*i]);
  const Polygon p(xys, 4);
  Real hfps[NonsingularTestIntegrands::nint] = {0};
  if ( ! calc_hfp(Options(), p, f.cc, f, hfps)) ++ne;
  if (std::abs(hfps[0] - 45.0/4) > 1e4*mv2::eps) {
    printf("integrals::test_integral: %f %e\n",
           hfps[0], std::abs(hfps[0] - 45.0/4));
    ++ne;
  }
  const Real value = 45.0/4;
  for (int i = 0; i < 2; ++i) {
    Real integrals[NonsingularTestIntegrands::nint] = {0};
    if (i == 0) {
      if ( ! calc_integral(p, f, integrals)) ++ne;
    } else {
      Pt pt{5,5};
      Options o;
      o.np_radial = o.np_angular = 12;
      if ( ! calc_integral_tensor_quadrature(o, p, f, pt, true, integrals))
        ++ne;
    }
    if (std::abs(integrals[0] - value) > 1e2*mv2::eps) {
      printf("otherint::test_integral: %f %e\n",
             integrals[0], std::abs(integrals[0] - value));
      ++ne;
    }
  }
  return ne;
}

int analyze_area () {
  int ne = 0;
  struct AnalyzeP : hfp::CallerP {
    Real eval (const Real x) const override { return 2*M_PI*cube(x); }
    void eval_p_c (const Real c, Real& pc, Real& pdc) const override
    { pc = 0; pdc = 0; }
  };
  Real xys[] = {-1,-1,1,-1,1,1,-1,1}, cc[2] = {0.1, -0.2};
  const int n = sizeof(xys)/sizeof(*xys)/2;
  for (int i = 0; i < n; ++i) mv2::axpy(1, cc, &xys[2*i]);
  const Polygon p(xys, n);
  Real area = 0;
  calc_hfp(Options(), p, cc, One(), &area);
  if (std::abs(area -4) > 1e4*mv2::eps) ++ne;
  area = hfp::calc_hfp(hfp::Options(), AnalyzeP(), 0, 1, 0);
  if (std::abs(area - M_PI) > 5e1*mv2::eps) ++ne;
  if (ne > 0) printf("integrals::analyze_area failed %d\n", ne);
  return ne;
}

} // namespace

int unittest () {
  int ne = 0;
  ne += analyze_area();
  ne += test_area();
  ne += test_calc_hfp_circle();
  ne += test_integral();
  return ne;
}
  
} // namespace integrals
} // namespace acorn
} // namespace woodland
