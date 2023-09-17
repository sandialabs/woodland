#ifndef INCLUDE_WOODLAND_ACORN_INTERACTION_INTEGRALS
#define INCLUDE_WOODLAND_ACORN_INTERACTION_INTEGRALS

#include "woodland/acorn/caller_integrand.hpp"
#include "woodland/acorn/plane_geometry.hpp"

namespace woodland {
namespace acorn {
namespace integrals {

/* Let f be well behaved or else singular at the point cc with singularity
   1/r^3, r = distance(x - cc). Compute the Hadamard finite part (h.f.p.) or
   integral of the function f over a convex polygon p. cc can be inside or
   outside of p.

   If cc is outside of p, then the h.f.p. reduces to a proper integral. If cc is
   far away from p (in the sense of distance relative to p's size), then
   calc_integral is more accurate than calc_hfp. But if cc is near p, then
   calc_hfp is increasingly more accurate than calc_integral as cc gets nearer
   to p from the outside.

   WARNING: Right now, we handle only a 1/r^3 singularity. In particular, these
   routines may give very inaccurate approximations for nearly singular
   functions g, even if closely related to f ~ 1/r^3.
 */

using plane::Pt;
using plane::Polygon;

struct Options {
  int np_radial, np_angular;
  Options();
};

// Calc the h.f.p. of f over the convex polygon p. cc can be outside of p. May
// return false if something goes wrong, but does not check inputs.
bool calc_hfp(
  const Options& o,
  // Convex polygon over which to integrate.
  const Polygon& p,
  // Location of singularity.
  const Pt cc,
  const CallerIntegrands& f,
  // Values are += into input values.
  RPtr hfps);

// Calc the integral of f over the convex polygon p. May return false if
// something goes wrong, but does not check inputs.
bool calc_integral(
  // Convex polygon over which to integrate.
  const Polygon& p,
  const CallerIntegrands& f,
  // Values are += into input values.
  RPtr integrals,
  const int tq_order = 20);

int unittest();

} // namespace integrals
} // namespace acorn
} // namespace woodland

#endif
