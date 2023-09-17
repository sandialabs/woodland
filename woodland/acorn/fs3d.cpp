#include "woodland/acorn/fs3d.hpp"

#include "woodland/acorn/macros.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/elastostatics.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/elastostatics_integrals.hpp"
#include "woodland/acorn/vectorization.hpp"
#include "woodland/acorn/dbg.hpp"

#include <cmath>
#include <limits>

namespace woodland {
namespace acorn {
namespace fs3d {

typedef Matvec<2,Real> mv2;
typedef Matvec<3,Real> mv3;

// Solve for p such that
//     r ls= s + [xhat, yhat] p,
// where ls= means least-squares fit.
static void
calc_singular_point (const Real xhat[3], const Real yhat[3], const Real s[3],
                     const Real r[3], Real p[2]) {
  Real tmp[3];
  mv3::subtract(r, s, tmp);
  const Real rhs[] = {mv3::dot(xhat, tmp), mv3::dot(yhat, tmp)};
  const Real a = mv3::dot(xhat, xhat), b = mv3::dot(xhat, yhat),
    d = mv3::dot(yhat, yhat), det = a*d - b*b;
  p[0] = ( d*rhs[0] - b*rhs[1])/det;
  p[1] = (-b*rhs[0] + a*rhs[1])/det;
}

struct Integrands : public acorn::CallerIntegrands {
  Integrands (const Real lam_, const Real mu_, CRPtr src_, CRPtr zhat_,
              CRPtr xhat_, CRPtr disloc_, CRPtr rcv_)
    : lam(lam_), mu(mu_), src(src_), zhat(zhat_), xhat(xhat_), disloc(disloc_),
      rcv(rcv_)
  {
    mv3::cross(zhat, xhat, yhat);
    calc_singular_point(xhat, yhat, src, rcv, cc);
  }

  const Real* get_singular_point () const { return cc; }

  virtual int nintegrands () const override { return 6; }

  virtual Real permitted_R_min (const Real R_max) const override {
    return 1e-3*R_max;
  }

  virtual void eval (const int n, CRPtr p, RPtr integrand) const override {
    for (int i = 0; i < n; ++i) {
      Real pt_src[3];
      mv3::sum3(1, src, p[2*i], xhat, p[2*i+1], yhat, pt_src);
      calc_sigma_point(lam, mu, pt_src, zhat, disloc, rcv, &integrand[6*i]);
    }
  }

  const Real lam, mu;
  CRPtr src, zhat, xhat, disloc, rcv;
  Real yhat[3], cc[2];
};

static void rcv_distance (const Real xy_side_lens[2], const Integrands& f,
                          // (distance outside of rect, distance above/below rect)
                          Real dist[2]) {
  Real pg[3];
  mv3::subtract(f.rcv, f.src, pg);
  const Real pl[] = {mv3::dot(f.xhat, pg), mv3::dot(f.yhat, pg),
                     mv3::dot(f.zhat, pg)};
  const Real hx = xy_side_lens[0]/2, hy = xy_side_lens[1]/2;
  dist[1] = std::abs(pl[2]);
  const Real d4[] = {pl[0] - hx, -hx - pl[0],
                     pl[1] - hy, -hy - pl[1]};
  dist[0] = 0;
  for (int i = 0; i < 4; ++i) dist[0] = std::max(dist[0], d4[i]);
}

void calc_sigma_const_disloc_rect (
  const Real lam, const Real mu, const Real src[3], const Real nml[3],
  const Real tangent[3], const Real xy_side_lens[2], const Real disloc[3],
  const Real rcv[3], Real sigma[6],
  const int np_radial, const int np_angular, int triquad_order)
{
  using namespace acorn;
  integrals::Options io;
  io.np_radial  = np_radial;
  io.np_angular = np_angular;
  const Real hx = xy_side_lens[0]/2, hy = xy_side_lens[1]/2;
  const Real xys[] = {-hx,-hy, hx,-hy, hx,hy, -hx,hy};
  const plane::Polygon p(xys, 4);
  Integrands f(lam, mu, src, nml, tangent, disloc, rcv);
  for (int i = 0; i < 6; ++i) sigma[i] = 0;
  Real dist[2];
  rcv_distance(xy_side_lens, f, dist);
  const Real L = std::max(hx, hy);
  if (dist[1] < 1e-4*L && dist[0] < L)
    integrals::calc_hfp(io, p, f.get_singular_point(), f, sigma);
  else {
    if (triquad_order <= 0) triquad_order = get_triquad_order(L, dist[0]);
    integrals::calc_integral(p, f, sigma, triquad_order);
  }
}

void calc_sigma_const_disloc_rect_okada (
  const Real lam, const Real mu, const Real src[3], const Real nml[3],
  const Real tangent[3], const Real xy_side_lens[2], const Real disloc[3],
  const Real rcv[3], Real sigma[6])
{
  const Real n[] = {0, 0, 1}, s[3] = {0};
  Real R[9], r[3], d[3], tmp[3];
  mv3::copy(tangent, R);
  mv3::copy(nml, R+6);
  mv3::cross(nml, tangent, R+3);
  mv3::subtract(rcv, src, tmp);
  mv3::matvec(R, tmp, r);
  mv3::matvec(R, disloc, d);
  acorn::call_okada(true, false, lam, mu, 1, s, r, d, n, xy_side_lens, sigma);
  rotate_sym_tensor_3x3_RtAR(R, sigma);
}

bool time_calc_sigma_point (const int n, const bool verbose) {
  const int N = RealPack::n*n;
  const Real lam = 1, mu = 1, src[] = {0,0,0}, nml[] = {0,0,1}, rcv[] = {1,0,0};
  Real a1 = 0, a2 = 0;
  {
    const auto t1 = dbg::gettime();
    Real disloc[3] = {0}, sigma[6];
    for (int i = 0; i < N; ++i) {
      disloc[0] = Real(i)/N;
      fs3d::calc_sigma_point(lam, mu, src, nml, disloc, rcv, sigma);
      a1 += sigma[2];
    }
    const auto t2 = dbg::gettime();
    if (verbose) printf("scalar %1.3e\n", t2 - t1);
  }
#ifdef WOODLAND_ACORN_VECTORIZE
  {
    RealPack psrc[3], pnml[3], prcv[3], pdisloc[3];
    for (int i = 0; i < RealPack::n; ++i) {
      for (int d = 0; d < 3; ++d) psrc[d][i] = src[d];
      for (int d = 0; d < 3; ++d) pnml[d][i] = nml[d];
      for (int d = 0; d < 3; ++d) prcv[d][i] = rcv[d];
      for (int d = 0; d < 3; ++d) pdisloc[d][i] = 0;
    }
    const auto t1 = dbg::gettime();
    RealPack psigma[6];
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < RealPack::n; ++j)
        pdisloc[0][j] = Real(RealPack::n*i+j)/N;
      fs3d::calc_sigma_point(lam, mu, psrc, pnml, pdisloc, prcv, psigma);
      for (int j = 0; j < RealPack::n; ++j)
        a2 += psigma[2][j];
    }
    const auto t2 = dbg::gettime();
    if (verbose) printf("vector %1.3e\n", t2 - t1);
  }
#endif
  if (verbose) printf("%23.15e %23.15e\n", a1, a2-a1);
  return a1 == a2;
}

// Just run the simplest problem. test_*space_against_okada in unittest.cpp
// tests the routine more extensively.
int unittest () {
  int ne = 0;
  const Real lam = 1, mu = 1, src[] = {0,0,0}, nml[] = {0,0,1},
    tangent[] = {1,0,0}, xy_side_lens[] = {1,1}, disloc[] = {1,1,1},
    rcv[] = {0,0,0};
  Real sigma[6];
  calc_sigma_const_disloc_rect(lam, mu, src, nml, tangent, xy_side_lens, disloc,
                               rcv, sigma);
  return ne;
}

} // namespace fs3d
} // namespace acorn
} // namespace woodland
