#ifndef INCLUDE_WOODLAND_ACORN_INTERACTION_INTEGRALS
#define INCLUDE_WOODLAND_ACORN_INTERACTION_INTEGRALS

#include "woodland/acorn/caller_integrand.hpp"
#include "woodland/acorn/plane_geometry.hpp"
#include "woodland/acorn/workspace.hpp"

namespace woodland {
namespace acorn {
namespace integrals {

/* Let f be well behaved or else singular at the point x_s with singularity
   1/r^3, r = distance(x - x_s). Compute the Hadamard finite part (h.f.p.) or
   integral of the function f over a convex polygon p. x_s can be inside or
   outside of p.

   If x_s is outside of p, then the h.f.p. reduces to a proper integral. If x_s
   is far away from p (in the sense of distance relative to p's size), then
   calc_integral is more accurate than calc_hfp. But if x_s is near p, then
   calc_hfp is increasingly more accurate than calc_integral as x_s gets nearer
   to p from the outside.

   WARNING: Right now, we handle only a 1/r^3 singularity. In particular, these
   routines may give very inaccurate approximations for nearly singular
   functions g, even if closely related to f ~ 1/r^3.
 */

using plane::Pt;
using plane::Polygon;

struct Options {
  int np_radial = 40, np_angular = 40;
  // Ignore a subpolygon having area relative to that of the parent polygon
  // smaller than this.
  Real relative_area_tol = 1e-8;
};

// Calc the h.f.p. of f over the convex polygon p. cc can be outside of p. May
// return false if something goes wrong, but does not check inputs.
bool calc_hfp(
  Workspace& w,
  const Options& o,
  // Convex polygon over which to integrate.
  const Polygon& p,
  // f.singular_pos must return true.
  const CallerIntegrands& f,
  // Values are += into input values.
  RPtr hfps);

// Calc the integral of f over the convex polygon p. May return false if
// something goes wrong, but does not check inputs.
bool calc_integral(
  Workspace& w,
  // Convex polygon over which to integrate.
  const Polygon& p,
  const CallerIntegrands& f,
  // Values are += into input values.
  RPtr integrals,
  const int tq_order = 20);

bool calc_integral_tensor_quadrature(
  Workspace& w, const Options& o, const Polygon& p, const CallerIntegrands& f,
  // Subdivide p at a common point in one of two ways:
  //   nearest_bdy_pt_to_anchor true: at the projection of anchor onto p's bdy;
  //                           false: at anchor.
  const Pt anchor, const bool nearest_bdy_pt_to_anchor,
  RPtr integrals);

int unittest();

void fig_init();
void fig_fin();

} // namespace integrals
} // namespace acorn
} // namespace woodland

#endif
