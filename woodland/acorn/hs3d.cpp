#include "woodland/acorn/hs3d.hpp"

#include "woodland/acorn/macros.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/elastostatics.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/elastostatics_integrals.hpp"
#include "woodland/acorn/vectorization.hpp"
#include "woodland/acorn/dbg.hpp"

#include <cmath>
#include <limits>

namespace woodland {
namespace acorn {
namespace hs3d {

typedef Matvec<2,Real> mv2;
typedef Matvec<3,Real> mv3;

// Canonical orientation: src y = (0,0,y3), rcv x = (x1,0,x3).
template <typename RealT>
void calc_halfspace_term_A(
  const Real lam, const Real mu, const RealT v[3], const RealT d[3],
  const RealT& y3, const RealT& x1, const RealT& x3, RealT sigma[6]);
template <typename RealT>
void calc_halfspace_term_B(
  const Real lam, const Real mu, const RealT v[3], const RealT d[3],
  const RealT& y3, const RealT& x1, const RealT& x3, RealT sigma[6]);
template <typename RealT>
void calc_halfspace_term_C(
  const Real lam, const Real mu, const RealT v[3], const RealT d[3],
  const RealT& y3, const RealT& x1, const RealT& x3, RealT sigma[6]);

// R = [c s 0; s -c 0; 0 0 1] such that (R x)[1] = 0.
static void form_rotation_to_canonical (const Real x[3], Real& c, Real& s) {
  if (x[1] == 0) {
    c = 1;
    s = 0;
    return;
  }
  c = x[0];
  s = x[1];
  const Real den = std::sqrt(c*c + s*s);
  c /= den;
  s /= den;
}

// General orientation.
template <typename RealT> void
calc_sigma_point_halfspace_terms (
  const Real lam, const Real mu, const RealT src[3], const RealT nml[3],
  const RealT disloc[3], const RealT rcv[3], RealT sigma[6])
{
  assert(src[2] <= 0);
  assert(rcv[2] <= 0);
  Real rcv_lcl[] = {rcv[0] - src[0], rcv[1] - src[1], rcv[2]};
  Real c, s;
  form_rotation_to_canonical(rcv_lcl, c, s);
  const Real
    x1 = c*rcv_lcl[0] + s*rcv_lcl[1],
    x3 = rcv_lcl[2];
  Real nml_lcl[3], disloc_lcl[3];
  mv3::copy(nml, nml_lcl);
  rotate_vector_2(c, -s, nml_lcl);
  mv3::copy(disloc, disloc_lcl);
  rotate_vector_2(c, -s, disloc_lcl);
  Real sigma_A[6], sigma_B[6], sigma_C[6];
  calc_halfspace_term_A(lam, mu, nml_lcl, disloc_lcl, src[2], x1, x3, sigma_A);
  calc_halfspace_term_B(lam, mu, nml_lcl, disloc_lcl, src[2], x1, x3, sigma_B);
  calc_halfspace_term_C(lam, mu, nml_lcl, disloc_lcl, src[2], x1, x3, sigma_C);
  for (int i = 0; i < 6; ++i) sigma[i] = sigma_A[i] + sigma_B[i] + sigma_C[i];
  const Real R[] = {c, s, 0, -s, c, 0, 0, 0, 1};
  rotate_sym_tensor_3x3_RtAR(R, sigma);
}

template void calc_sigma_point_halfspace_terms<Real>(
  const Real lam, const Real mu, const Real src[3], const Real nml[3],
  const Real disloc[3], const Real rcv[3], Real sigma[6]);

struct HalfspaceTermsIntegrands : public CallerIntegrands {
  HalfspaceTermsIntegrands (
    const Real lam_, const Real mu_, CRPtr src_, CRPtr tangent, CRPtr nml_,
    CRPtr disloc_, CRPtr rcv_)
    : lam(lam_), mu(mu_), src(src_), xhat(tangent), nml(nml_), disloc(disloc_),
      rcv(rcv_)
  {
    mv3::cross(nml, tangent, yhat);
  }
  
  int nintegrands () const override { return 6; }

  void eval (const int n, CRPtr p, RPtr integrands) const override {
    for (int k = 0; k < n; ++k) {
      Real pt_src[3];
      mv3::sum3(1, src, p[2*k], xhat, p[2*k+1], yhat, pt_src);
      calc_sigma_point_halfspace_terms(lam, mu, pt_src, nml, disloc, rcv,
                                       &integrands[6*k]);
    }
  }

  const Real lam, mu;
  CRPtr src, xhat, nml, disloc, rcv;
  Real yhat[3];
};

void calc_sigma_const_disloc_rect(
  Workspace& w, const Real lam, const Real mu, const Real src[3], const Real nml[3],
  const Real tangent[3], const Real xy_side_lens[2], const Real disloc[3],
  const Real rcv[3], Real sigma[6], const bool halfspace, const bool use_okada,
  const int hfp_np_radial, const int hfp_np_angular, int triquad_order,
  const Real triquad_tol)
{
  throw_if(halfspace && tangent[2] != 0,
           "calc_sigma_const_disloc_rect: if halfspace, tangent[2] must be 0");
  if (use_okada) {
    if (halfspace) {
      const Real
        z[3] = {0,0,1},
        s[] = {0, 0, src[2]};        
      Real R[9], d[3], n[3], tmp[3],
        r[] = {rcv[0] - src[0], rcv[1] - src[1], rcv[2]};
      mv3::copy(tangent, R);
      mv3::cross(z, tangent, R+3);
      mv3::copy(z, R+6);
      mv3::copy(r, tmp);
      mv3::matvec(R, tmp, r);
      mv3::matvec(R, disloc, d);
      mv3::matvec(R, nml, n);
      call_okada(false, false, lam, mu, 1, s, r, d, n, xy_side_lens, sigma);
      rotate_sym_tensor_3x3_RtAR(R, sigma);
    }
    else {
      fs3d::calc_sigma_const_disloc_rect_okada(
        lam, mu, src, nml, tangent, xy_side_lens, disloc, rcv, sigma);
    }
  } else {
    acorn::fs3d::RectInfo info;
    acorn::fs3d::calc_sigma_const_disloc_rect(
      w, lam, mu, src, nml, tangent, xy_side_lens, disloc, rcv, sigma,
      hfp_np_radial, hfp_np_angular, triquad_order, triquad_tol, &info);
    if (halfspace) {
      Real sigma_hs[6] = {0};
      HalfspaceTermsIntegrands hti(lam, mu, src, tangent, nml, disloc, rcv);
      const auto wh = xy_side_lens[0]/2, hh = xy_side_lens[1]/2;
      const Real xys[] = {-wh, -hh, wh, -hh, wh, hh, -wh, hh};
      plane::Polygon rect(xys, 4);
      triquad_order = info.triquad_order;
      integrals::calc_integral(w, rect, hti, sigma_hs, triquad_order);
      for (int i = 0; i < 6; ++i) sigma[i] += sigma_hs[i];
    }
  }
}

static bool cmp (const bool fullspace, const Real src[3], const Real rcv[3],
                 const Real disloc[3], const Real nml[3],
                 const Real sigma_o[6], const Real sigma_me[6], const Real tol) {
  Real diff[6], maxval = 0;
  for (int i = 0; i < 6; ++i) {
    diff[i] = sigma_me[i] - sigma_o[i];
    maxval = std::max(maxval, sigma_o[i]);
  }
  bool ok = true;
  for (int i = 0; i < 6; ++i)
    if (std::abs(diff[i]) > tol*maxval)
      ok = false;
  if (not ok) {
#ifdef WOODLAND_ACORN_HAVE_DC3D
    prc(fullspace);
    prarr("src",src,3);
    prarr("rcv",rcv,3);
    prarr("disloc",disloc,3);
    prarr("nml",nml,3);
    prarr("sigma_o",sigma_o,6);
    prarr("sigma_m",sigma_me,6);
    prarr("diff",diff,6);
#endif
  }
  return ok;
}

// Test point GF.
static int testpt (Real s1, Real s2, Real s3, Real n1, Real n2, Real n3,
                   Real d1, Real d2, Real d3, Real r1, Real r2, Real r3) {
  const Real lam = 0.95, mu = 1.1;
  const Real p1[] = {s1,s2,s3}, p2[] = {r1,r2,r3};
  const Real disloc[] = {d1,d2,d3};
  Real nml[] = {n1,n2,n3};
  mv3::normalize(nml);
  const bool point = true;
  bool all_ok = true;
  for (const bool fullspace : {true, false}) {
    for (int trial = 0; trial < 2; ++trial) {
      const Real* src = trial == 0 ? p1 : p2;
      const Real* rcv = trial == 0 ? p2 : p1;
      Real sigma_o[6], sigma_me[6];
      call_okada(fullspace, point, lam, mu, 1, src, rcv, disloc, nml, nullptr,
                 sigma_o);
      fs3d::calc_sigma_point(lam, mu, src, nml, disloc, rcv, sigma_me);
      if (not fullspace) {
        Real sigma_me_hs[6];
        calc_sigma_point_halfspace_terms(lam, mu, src, nml, disloc, rcv,
                                         sigma_me_hs);
        for (int i = 0; i < 6; ++i) sigma_me[i] += sigma_me_hs[i];
      }
      const bool ok = cmp(fullspace, src, rcv, disloc, nml,
                          sigma_o, sigma_me, 100*mv2::eps);
      if (not ok) all_ok = false;
    }
  }
  return all_ok ? 0 : 1;
}

// Test rectangular GF with Okada's constraints on geometry. The source is a
// rectangle having width w, height h, and face normal n with its w-axis
// parallel to the free surface.
static int testrect (
  Real s1, Real s2, Real s3, Real n1, Real n2, Real n3, Real w, Real h,
  Real d1, Real d2, Real d3, Real r1, Real r2, Real r3, Real tolmult = 1e4)
{
  Workspace ws;
  const Real lam = 0.95, mu = 1.1;
  const Real src[] = {s1,s2,s3}, rcv[] = {r1,r2,r3}, disloc[] = {d1,d2,d3},
    xy_side_lens[] = {w,h};
  Real nml[] = {n1,n2,n3};
  mv3::normalize(nml);

  Real tangent[3]; {
    const Real zhat[] = {0,0,1};
    mv3::cross(nml, zhat, tangent);
    const Real den = mv3::norm2(tangent);
    if (den == 0) {
      tangent[0] = 1;
      tangent[1] = tangent[2] = 0;
    } else {
      mv3::divide(den, tangent);
    }
    tangent[2] = 0; // numerics
  }

  bool all_ok = true;
  for (const bool fullspace : {true, false}) {
    Real sigma_o[6], sigma_me[6];
    calc_sigma_const_disloc_rect(ws, lam, mu, src, nml, tangent, xy_side_lens,
                                 disloc, rcv, sigma_o , not fullspace, true);
    calc_sigma_const_disloc_rect(ws, lam, mu, src, nml, tangent, xy_side_lens,
                                 disloc, rcv, sigma_me, not fullspace, false);
    const bool ok = cmp(fullspace, src, rcv, disloc, nml,
                        sigma_o, sigma_me, tolmult*mv2::eps);
    if (not ok) all_ok = false;
  }
  
  return all_ok ? 0 : 1;
}

int unittest () {
#ifndef WOODLAND_ACORN_HAVE_DC3D
  printf("hs3d::unittest: Because extern/dc3d.f is not available, "
         "this unit test will purposely fail.\n");
#endif
  int n = 0;
  n += testpt(0.1, -0.3, -1,  -1, 0.5, 0.4,  1, -0.7, 0.3,  -0.2, 0.1, -0.1);
  n += testpt(0.1, -0.3, -1,  -1, 0.5, 0.4,  1, -0.7, 0.3,  0, 0, 0);
  n += testpt(0, 0, -1,       -1, 0.5, 0.4,  1, -0.7, 0.3,  0, 0, 0);
  n += testpt(0, 0, -1,       0, 0, 1,       1, -0.7, 0.3,  0, 0, -0.1);
  n += testrect(0.1, -0.3, -1,    -1, 0.5, 0.4,  0.4, 0.2,  1, -0.7, 0.3,  -0.2, 0.1, -0.1);
  n += testrect(0.1, -0.3, -1,    -1, 0.5, 0.4,  0.2, 0.4,  1, -0.7, 0.3,  -0.2, 0.1, -0.1);
  n += testrect(-0.2, 0.1, -0.2,  -1, 0.5, 0.4,  0.4, 0.1,  1, -0.7, 0.3,  0.1, -0.3, -1);
  n += testrect(0, 0, -1,         1e-8, 0, 1,    0.4, 0.2,  1, -0.7, 0.3,  0, 0, -0.1);
  n += testrect(0.1, -0.3, -1,    -1, 0.5, 0.4,  0.4, 0.2,  1, -0.7, 0.3,  3.2, 0.1, -0.1, 1e5);
  n += testrect(0.1, -0.3, -1,    -1, 0.5, 0.4,  0.4, 0.2,  1, -0.7, 0.3,  3.2, -7.1, -1, 1e5);
  return n;
}

} // namespace hs3d
} // namespace acorn
} // namespace woodland
