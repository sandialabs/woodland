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

  bool singular_pt (Real p[2]) const override {
    mv2::copy(cc, p);
    return true;
  }

  int nintegrands () const override { return 6; }

  void eval (const int n, CRPtr p, RPtr integrand) const override {
    return eval(n, p, nullptr, false, integrand);
  }

  //Real permitted_r_min (const Real r_max) const override { return 1e-3*r_max; }

  bool supports_mult_by_R3 () const override { return true; }

  void eval_mult_by_R3 (const int n, CRPtr p, CRPtr J_times_pdir,
                        RPtr integrand) const override {
    eval(n, p, J_times_pdir, true, integrand);
  }

  void eval_shape_J(const Real p[2], Real J[6]) const override {
    for (int d = 0; d < 3; ++d) J[2*d  ] = xhat[d];
    for (int d = 0; d < 3; ++d) J[2*d+1] = yhat[d];
  }

  const Real lam, mu;
  CRPtr src, zhat, xhat, disloc, rcv;
  Real yhat[3], cc[2];

private:
  void eval (const int n, CRPtr p, CRPtr J_times_pdir, const bool mult_by_R3,
             RPtr integrand) const {
    for (int i = 0; i < n; ++i) {
      Real pt_src[3];
      mv3::sum3(1, src, p[2*i], xhat, p[2*i+1], yhat, pt_src);
      calc_sigma_point(lam, mu, pt_src, zhat, disloc, rcv, &integrand[6*i],
                       mult_by_R3, mult_by_R3 ? &J_times_pdir[3*i] : nullptr);
    }
  }
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
  Workspace& w, const Real lam, const Real mu, const Real src[3],
  const Real nml[3], const Real tangent[3], const Real xy_side_lens[2],
  const Real disloc[3], const Real rcv[3], Real sigma[6],
  const int np_radial, const int np_angular, int triquad_order,
  const Real triquad_tol, RectInfo* info)
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
  const bool hfp = dist[1] < 1e-4*L && dist[0] < L;
  if (hfp)
    integrals::calc_hfp(w, io, p, f, sigma);
  else {
    if (triquad_order <= 0)
      triquad_order = get_triquad_order(L, dist[0], triquad_tol);
    integrals::calc_integral(w, p, f, sigma, triquad_order);
  }
  if (info) {
    info->hfp = hfp;
    info->L = L;
    mv2::copy(dist, info->dist);
    info->triquad_order = triquad_order;
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
  Workspace w;
  {
#ifndef WOODLAND_ACORN_HAVE_DC3D
    printf("fs3d::unittest: Because extern/dc3d.f is not available, "
           "this unit test will purposely fail.\n");
#endif
    const Real lam = 1, mu = 1.1;
    const Real src[] = {0, 0, -0.2}, rcv[] = {0.1, 0, -0.1},
      disloc[] = {1, 1.5, 0.3};
    for (const Real n3 : {1.0, 0.99, 1.01}) {
      Real nml[] = {1, 0, n3};
      mv3::normalize(nml);
      const bool point = true;
      Real sigma_o[6], sigma_me[6];
      call_okada(true, point, lam, mu, 1, src, rcv, disloc, nml, nullptr,
                 sigma_o);
      fs3d::calc_sigma_point(lam, mu, src, nml, disloc, rcv, sigma_me);
      Real num = 0, den = 0;
      for (int i = 0; i < 6; ++i) {
        num = std::max(num, std::abs(sigma_me[i] - sigma_o[i]));
        den = std::max(den, std::abs(sigma_o[i]));
      }
      if (num > 100*mv3::eps*den) {
#ifdef WOODLAND_ACORN_HAVE_DC3D
        prc(num/den);
        prarr("sigma_o",sigma_o,6);
        prarr("sigma_m",sigma_me,6);
#endif
        ++ne;
      }
    }
  }
  {
    const Real lam = 1, mu = 1, src[] = {0,0,0}, nml[] = {0,0,1},
      tangent[] = {1,0,0}, xy_side_lens[] = {1,1}, disloc[] = {1,1,1},
      rcv[] = {0,0,0};
    Real sigma[6];
    calc_sigma_const_disloc_rect(w, lam, mu, src, nml, tangent, xy_side_lens,
                                 disloc, rcv, sigma);
  }
  return ne;
}

static void study_triquad_table () {
  Workspace w;
  const int triquad_orders[] = {1, 2, 4, 6, 8, 12};
  const int nto = sizeof(triquad_orders)/sizeof(*triquad_orders);
  const Real tols[] = {1e-6, 1e-8, 1e-10, 1e-12, 1e-14};
  const int ntol = sizeof(tols)/sizeof(*tols);

  const Real lam = 1, mu = 1;
  const Real src[3] = {0}, nml[] = {0,0,1}, tan[] = {1,0,0}, xy_side_lens[] = {1,1};
  const auto calc_sigma =
    [&] (CRPtr rcv, CRPtr disloc, const int triquad_order, RPtr sigma) {
      calc_sigma_const_disloc_rect(w, lam, mu, src, nml, tan, xy_side_lens,
                                   disloc, rcv, sigma, -1, -1, triquad_order);
    };

  const Real max_rad = 5000, max_lat = M_PI/4;
  const int nrad = 250, fac = 3, nlat = fac*10, ncirc = fac*10;
  Real disloc[3] = {0}, rcv[3], sigma_ref[6], sigma[6];
  Real rad_usable[ntol][nto] = {0};
  for (int irad = 0; irad <= nrad; ++irad) {
    const Real a = Real(irad)/nrad;
    const Real radius = std::exp((1-a)*std::log(1.5) + a*std::log(max_rad));
    Real nums[nto] = {0}, dens[nto] = {0};
    for (int ilat = 0; ilat <= nlat; ++ilat) {
      const Real lat = Real(ilat)/nlat*max_lat;
      rcv[2] = radius*std::sin(lat);
      const Real f = std::cos(lat);
      for (int icirc = 0; icirc <= ncirc; ++icirc) {
        const Real theta = Real(icirc)/ncirc*(M_PI/2);
        rcv[0] = radius*f*std::cos(theta);
        rcv[1] = radius*f*std::sin(theta);
        for (int id = 0; id < 3; ++id) {
          disloc[id] = 1;
          calc_sigma(rcv, disloc, 20, sigma_ref);
          for (int ito = 0; ito < nto; ++ito) {
            calc_sigma(rcv, disloc, triquad_orders[ito], sigma);
            for (int i = 0; i < 6; ++i) {
              dens[ito] = std::max(dens[ito], std::abs(sigma_ref[i]));
              nums[ito] = std::max(nums[ito], std::abs(sigma[i] - sigma_ref[i]));
            }
          }
          disloc[id] = 0;
        }
      }
    }
    printf("%1.5e:", radius);
    for (int ito = 0; ito < nto; ++ito) {
      const Real err = nums[ito]/dens[ito];
      for (int ie = 0; ie < ntol; ++ie)
        if (err < tols[ie] && rad_usable[ie][ito] == 0)
          rad_usable[ie][ito] = radius;
      printf(" %1.3e", err);
    }
    printf("\n");
  }
  for (int ie = ntol-1; ie >= 0; --ie) {
    if (ie > 0) {
      printf(ie == ntol-1 ? "  " : " else ");
      printf("if (tol <= %1.1e) {\n", tols[ie]);
    } else {
      printf(" else {\n");
    }
    bool done = false;
    for (int ito = nto-1; ito >= 0; --ito) {
      printf(ito == nto-1 ? "    return (" : "            ");
      if (rad_usable[ie][ito] != 0) {
        printf("dist < %6.1f*L ? %2d :\n",
               std::ceil(10*(rad_usable[ie][ito] - std::sqrt(0.5)))/10.0,
               ito == nto-1 ? 20 : triquad_orders[ito+1]);
      } else {
        printf("%d);\n", triquad_orders[ito+1]);
        done = true;
        break;
      }
    }
    if (not done) printf("            %d);\n", triquad_orders[0]);
    printf("  }");
  }
  printf("\n");
}

void study_triquad () {
  study_triquad_table();
  //todo binary search
}

struct FigureIntegrands : public acorn::CallerIntegrands {
  enum : int { ndisloc = 3, nsigma = 6 };

  FigureIntegrands ()
    : vtxs{-0.75,-0.75, 1.1,-1, 1,0.9, -0.6,0.4},
      rcv{0,0}
  {}

  int nintegrands () const override { return nsigma; }

  Real permitted_r_min (const Real r_max) const override {
    return 1e-3*r_max;
  }

  bool singular_pt (Real p[2]) {
    mv2::copy(rcv, p);
    return true;
  }

  plane::Polygon get_polygon () const {
    return plane::Polygon(vtxs, sizeof(vtxs)/sizeof(*vtxs)/2);
  }

  void eval (const int n, CRPtr xys_lcs, RPtr integrands) const override {
    const Real lam = 1, mu = 1;
    const Real disloc_lcs[] = {1,0,0};

    Real rcv_gcs[3], lcs_rcv[9];
    mv2::copy(rcv, rcv_gcs);
    eval_shape(rcv[0], rcv[1], rcv_gcs[2]);
    make_lcs(rcv[0], rcv[1], lcs_rcv);

    for (int k = 0; k < n; ++k) {
      Real src_gcs[3], src_grad[2], lcs_src[9];
      mv2::copy(&xys_lcs[2*k], src_gcs);
      eval_shape(src_gcs[0], src_gcs[1], src_gcs[2], src_grad);
      make_lcs(src_gcs[0], src_gcs[1], lcs_src);
      
      Real disloc_gcs[3];
      mv3::tmatvec(lcs_src, disloc_lcs, disloc_gcs);

      const auto integrand = &integrands[nsigma*k];
      acorn::fs3d::calc_sigma_point(lam, mu, src_gcs, &lcs_src[6],
                                    disloc_gcs, rcv_gcs, integrand);

      const Real jacdet = std::sqrt(1 + mv2::norm22(src_grad));
      for (int i = 0; i < nsigma; ++i) integrand[i] *= jacdet;
    }
  }

private:
  const Real vtxs[8], rcv[2];

  static void
  eval_shape (const Real x, const Real y, Real& z, RPtr grad = nullptr) {
    const Real a = 0.3;
    z = a*(x*x + y*y);
    if (grad) {
      grad[0] = 2*a*x;
      grad[1] = 2*a*y;
    }
  }

  static void make_lcs (const Real x, const Real y, Real lcs[9]) {
    const Real primary[] = {1,0,0};
    Real z, grad[2];
    eval_shape(x, y, z, grad);
    Real zhat[] = {-grad[0], -grad[1], 1};
    mv3::normalize(zhat);
    Real yhat[3], xhat[3];
    mv3::cross(zhat, primary, yhat);
    mv3::normalize(yhat);
    mv3::cross(yhat, zhat, xhat);
    assert(std::abs(mv3::norm2(xhat) - 1) < 1e-13);
    mv3::copy(xhat, lcs);
    mv3::copy(yhat, lcs+3);
    mv3::copy(zhat, lcs+6);
  }
};

void make_figure_data () {
#ifdef WOODLAND_ACORN_FIGURE
  integrals::fig_init();
  integrals::Options o;
  o.np_radial = 20;
  FigureIntegrands f;
  Real hfps[6];
  Workspace w;
  calc_hfp(w, o, f.get_polygon(), f, hfps);
  integrals::fig_fin();
#else
  printf("define WOODLAND_ACORN_FIGURE to enable\n");
#endif
}

} // namespace fs3d
} // namespace acorn
} // namespace woodland
