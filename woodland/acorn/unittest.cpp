#include <cstdio>

#include "woodland/acorn/compose_triquad.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/linalg.hpp"
#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/plane_geometry.hpp"
#include "woodland/acorn/hfp.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/hs3d.hpp"
#include "woodland/acorn/elastostatics.hpp"
#include "woodland/acorn/vv.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace acorn {
namespace {

typedef Matvec<2,Real> mv2;
typedef Matvec<3,Real> mv3;

#ifdef WOODLAND_ACORN_HAVE_DC3D
int test_fullspace_rect_against_okada () {
  int ne = 0;
  const Real lam = 1.1, mu = 0.9;
  // Test scaling of the problem. dc3* has issues when scl gets too small. Our
  // code seems fine.
  for (const Real scl : {1e-2, 1e6}) {
    // Simple geometry.
    for (int idx = 0; idx < 2; ++idx) {
      const Real dx = idx == 0 ? 0 : 3.5;
      const Real src[3] = {0}, nml[] = {0,0,1}, tangent[] = {1,0,0},
        xy_side_lens[] = {scl*3, scl*2.4},
        disloc[] = {scl*0.8, scl*-0.7, scl*0.2};
      int nx = 7, ny = 5;
      // Cut down a little on runtime.
      if (idx == 0) { --nx; --ny; }
      if (scl > 1) { nx -= 2; ny -= 2; }
      Real err_me[6] = {0}, err_ok1[6] = {0}, den[6] = {0};
      for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
          const Real rcv[] = {scl*(-1 + 2*Real(i)/(nx-1) + dx),
                              scl*(-1 + 2*Real(j)/(ny-1)),
                              scl*0};
          Real sigma_me[6], sigma_ok[6], sigma_ok1[6];
          fs3d::calc_sigma_const_disloc_rect(
            lam, mu, src, nml, tangent, xy_side_lens, disloc, rcv,
            sigma_me, 40, 40, -1, 1e-12);
          call_okada(
            true, false, lam, mu, 1, src, rcv, disloc, nml, xy_side_lens,
            sigma_ok);
          fs3d::calc_sigma_const_disloc_rect_okada(
            lam, mu, src, nml, tangent, xy_side_lens, disloc, rcv,
            sigma_ok1);
          for (int k = 0; k < 6; ++k)
            err_me[k] = std::max(err_me[k],
                                 std::abs(sigma_me[k] - sigma_ok[k]));
          for (int k = 0; k < 6; ++k)
            err_ok1[k] = std::max(err_ok1[k],
                                  std::abs(sigma_ok1[k] - sigma_ok[k]));
          for (int k = 0; k < 6; ++k)
            den[k] = std::max(den[k], std::abs(sigma_ok[k]));
        }
      for (int k = 0; k < 6; ++k)
        err_me[k] /= den[k];
      for (int k = 0; k < 6; ++k)
        if (err_me[k] > (idx == 0 ? 1e-7 : 5e-9)) {
          printf("test_fullspace_rect_against_okada (me): %d %d %1.2e\n",
                 idx, k, err_me[k]);
          ++ne;
        }
      for (int k = 0; k < 6; ++k)
        if (err_ok1[k] > 1e-15) {
          printf("test_fullspace_rect_against_okada (ok impl): %d %d %1.2e\n",
                 idx, k, err_ok1[k]);
          ++ne;
        }
    }
    // General geometry.
    for (int idx = 0; idx < 2; ++idx) {
      const Real dx = idx == 0 ? 0 : 3.5;
      const Real src[] = {scl*0.4, scl*-0.2, scl*-0.5},
        xy_side_lens[] = {scl*3, scl*2.4},
        // Disloc like above.
        odisloc[] = {scl*0.8, scl*-0.7, scl*0.2};
      // Set up nml and tangent.
      Real nml[] = {-0.1, 0.2, 0.8}, yhat[] = {0.5, 1, 0}, tangent[3];
      mv3::cross(yhat, nml, tangent);
      mv3::cross(nml, tangent, yhat);
      mv3::normalize(tangent);
      mv3::normalize(yhat);
      mv3::normalize(nml);
      // Rotate disloc so that, relative to the generally oriented fault, it's
      // the same as in the first case.
      Real disloc[3];
      mv3::sum3(odisloc[0], tangent, odisloc[1], yhat, odisloc[2], nml, disloc);
      int nx = 6, ny = 5;
      // Cut down a little on runtime.
      if (idx == 0) { --nx; --ny; }
      if (scl > 1) { nx -= 2; ny -= 2; }
      Real err[6] = {0}, den[6] = {0};
      for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
          // Receiver, again set up so it's the same as in the first case.
          const Real orcv[] = {scl*(-1 + 2*Real(i)/(nx-1) + dx),
                               scl*(-1 + 2*Real(j)/(ny-1)),
                               scl*0};
          Real rcv[3];
          mv3::sum3(orcv[0], tangent, orcv[1], yhat, orcv[2], nml, rcv);
          mv3::axpy(1, src, rcv);
          // Sigmas.
          Real sigma_me[6], sigma_ok[6];
          fs3d::calc_sigma_const_disloc_rect(
            lam, mu, src, nml, tangent, xy_side_lens, disloc, rcv,
            sigma_me, 40, 40);
          fs3d::calc_sigma_const_disloc_rect_okada(
            lam, mu, src, nml, tangent, xy_side_lens, disloc, rcv,
            sigma_ok);
          // Errors.
          for (int k = 0; k < 6; ++k)
            err[k] = std::max(err[k], std::abs(sigma_me[k] - sigma_ok[k]));
          for (int k = 0; k < 6; ++k)
            den[k] = std::max(den[k], std::abs(sigma_ok[k]));
        }
      for (int k = 0; k < 6; ++k)
        err[k] /= den[k];
      for (int k = 0; k < 6; ++k)
        if (err[k] > (idx == 0 ? 1e-7 : 4e-8)) {
          printf("test_fullspace_rect_against_okada (geometry): %d %d %1.2e\n",
                 idx, k, err[k]);
          ++ne;
        }
    }
  }
  return ne;
}

int test_halfspace_rect_against_okada () {
  int ne = 0;
  return ne;
}
#endif // WOODLAND_ACORN_HAVE_DC3D
} // namespace

int unittest () {
  int nerr = 0, ne;
  rununittest(util_test);
  rununittest(fs3d::unittest);
  rununittest(hs3d::unittest);
  rununittest(linalg::unittest);
  rununittest(TriangleQuadrature::unittest);
  rununittest(mv2::unittest);
  rununittest(mv3::unittest);
  rununittest(Triangle2D::unittest);
  rununittest(plane::unittest);
  rununittest(hfp::unittest);
  rununittest(integrals::unittest);
#ifdef WOODLAND_ACORN_HAVE_DC3D
  rununittest(test_fullspace_rect_against_okada);
  rununittest(test_halfspace_rect_against_okada);
#endif
#ifdef WOODLAND_ACORN_VECTORIZE
  if ( ! fs3d::time_calc_sigma_point(1000, false)) {
    printf("fs3d::time_calc_sigma_point failed\n");
    nerr += ne;
  }
#endif
  rununittest(vv::unittest);
  printf("\n%s\n", nerr ? "FAIL" : "PASS");
  return nerr;
}

} // namespace acorn
} // namespace woodland
