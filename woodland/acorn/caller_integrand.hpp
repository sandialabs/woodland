#ifndef INCLUDE_WOODLAND_ACORN_CALLER_INTEGRAND
#define INCLUDE_WOODLAND_ACORN_CALLER_INTEGRAND

#include <cassert>

#include "woodland/acorn/acorn.hpp"

namespace woodland {
namespace acorn {

// Caller's integrand over a domain in the parametric coordinates (u,v) in
// Real^2.
struct CallerIntegrands {
  // Number of integrands. Must be <= max_n_integrand.
  virtual int nintegrands() const = 0;

  // Evaluate each integrands[nintegrands()*i + j], j = 0 to nintegrands()-1, at
  // the points &ps[2*i], i = 0 to n-1.
  virtual void eval(const int n, CRPtr ps, RPtr integrands) const = 0;

  // >> For the case of a singular integrand:

  // Return true if there is a singularity and set p = (u_s, v_s), the location
  // of the singularity.
  virtual bool singular_pt(Real p[2]) const { return false; }

  // > Approach 1 to r->0.
  //   The user provides Green's function values in [r_min, r_max], and the
  // algorithm extrapolates to the region [0, r_min].
  //   In practice we find r_min ~ 1e-3*r_max is sufficient for a good solution.

  // Provide the smallest r = sqrt((u-u_s)^2 + (v-v_s)^2) < r_max at which this
  // integrand can be evaluated.
  virtual Real permitted_r_min(const Real r_max) const { return 0; }

  // > Approach 2 to r->0.
  //   If supports_mult_by_R3(), then eval_mult_by_R3() and eval_shape_J() must
  // be implemented.
  //   Green's function values are obtained in [0, r_max], so only interpolation
  // is needed subsequently. This is more accurate than approach 1.

  // Does the point Green's function routine support analytically multiplying by
  // R^3 to remove the singularity? Here, R is the distance from the singularity
  // in physical coordinates, whereas r is the distance in parametric
  // coordinates. If true, then permitted_r_min is not used.
  virtual bool supports_mult_by_R3 () const { return false; }

  // If it does, then this routine will get called in some cases. pdir is a unit
  // vector in parametric space pointing from p to the singularity. J_times_pdir
  // is J from eval_shape_J(p, J) times pdir. Return the same thing as eval, but
  // multiplied by R^3.
  virtual void eval_mult_by_R3 (const int n, CRPtr p, CRPtr J_times_pdir,
                                RPtr integrand) const {
    printf("Error: CallerIntegrands::eval_mult_by_R3 is not impl'ed.\n");
    assert(false);
  }

  // Evaluate J, the Jacobian of the surface over which the integrand is being
  // evaluated,
  //     (x,y,z) = f(u,v),
  // at the point p = (u,v). J is formatted as
  //     [x_u, x_v; y_u, y_v; z_u, z_v]
  // in row-major order. Examples:
  //   1. If the surface is flat and identical to (u,v), then
  //     J = [1 0; 0 1; 0 0].
  //   2. If the surface is parametric,
  //     (x,y,z) = (u,v,z(u,v)),
  // then
  //     J = [1 0; 0 1; z_u z_v].
  // Note that eval() already internally computes something similar, the area
  // element sqrt(det(J'J)); this method just provides J explicitly.
  virtual void eval_shape_J(const Real p[2], Real J[6]) const {
    printf("Error: CallerIntegrands::eval_shape_J is not impl'ed.\n");
    assert(false);
  }
};

} // namespace acorn
} // namespace woodland

#endif
