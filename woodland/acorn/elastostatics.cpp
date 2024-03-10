#include "woodland/acorn/acorn.hpp"
#include "woodland/acorn/elastostatics.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace acorn {

void call_okada (const bool fullspace, const bool point, const Real lam,
                 const Real mu, const int n, const Real* srcs, const Real* rcvs,
                 const Real* dislocs, const Real* nmls, const Real* strike_dip_dims,
                 Real* sigmas) {
  typedef Matvec<3,Real> mv3;
  const auto alpha = acorn::lambda_mu_to_alpha(lam, mu);
  for (int i = 0; i < n; ++i) {
    const Real* nml = &nmls[3*i];
    assert(std::abs(mv3::norm22(nml) - 1) <= std::sqrt(mv3::eps));
    const Real base = std::sqrt(square(nml[0]) + square(nml[1]));
    Real dipdeg = (180/M_PI)*atan2(nml[2], base);
    dipdeg = 90 - dipdeg;
    // Obs z <= 0 but depth is expressed as >= 0.
    Real p[3], depth;
    mv3::subtract(&rcvs[3*i], &srcs[3*i], p);
    if (fullspace) {
      if (p[2] <= 0)
        depth = 0;
      else {
        depth = p[2];
        p[2] = 0;
      }
    } else {
      depth = -srcs[2];
      p[2] = rcvs[3*i+2];
    }
    const Real* disloc = &dislocs[3*i];
    Real pot[3]; // strike, dip, tensile;
    Real along_strike[3], along_dip[3], nml_strike[3];
    const Real z[] = {0,0,1};
    const bool up = base == 0;
    if (up)
      mv3::copy(disloc, pot);
    else {
      mv3::cross(z, nml, along_strike);
      mv3::normalize(along_strike);
      mv3::cross(nml, along_strike, along_dip);
      mv3::cross(z, along_strike, nml_strike);
      // Dot disloc into the source coords.
      pot[0] = mv3::dot(disloc, along_strike);
      pot[1] = mv3::dot(disloc, along_dip);
      pot[2] = mv3::dot(disloc, nml);
      // Rotate p around the z axis into source coords.
      Real pr[3] = {0};
      pr[0] = mv3::dot(p, along_strike);
      pr[1] = mv3::dot(p, nml_strike);
      pr[2] = p[2];
      mv3::copy(pr, p);
    }
    Real u[3] = {0}, du[9] = {0};
    if (point)
      acorn::okada::dc3d0(alpha, p[0], p[1], p[2], depth, dipdeg,
                          pot[0], pot[1], pot[2], 0,
                          u, du,
                          fullspace);
    else {
      const Real al = strike_dip_dims[2*i  ]/2;
      const Real aw = strike_dip_dims[2*i+1]/2;
      const int iret =
        acorn::okada::dc3d(alpha, p[0], p[1], p[2], depth, dipdeg,
                           -al, al, -aw, aw,  // dc3d.f
                           // al, al, aw, aw, // dc3omp.f
                           pot[0], pot[1], pot[2],
                           u, du,
                           fullspace);
      if (iret != 0) printf("okada::dc3d returned %d\n", iret);
      assert(iret == 0);
    }
    Real* s = &sigmas[6*i];
    acorn::du_to_tau(lam, mu, du, s);
    if ( ! up)
      rotate_sym_tensor_3x3_RtAR(along_strike, nml_strike, z, s);
  }
}

} // namespace acorn
} // namespace woodland
