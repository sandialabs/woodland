#ifndef INCLUDE_WOODLAND_ACORN_CALLER_INTEGRAND
#define INCLUDE_WOODLAND_ACORN_CALLER_INTEGRAND

#include <cassert>

#include "woodland/acorn/acorn.hpp"

namespace woodland {
namespace acorn {

// Caller's integrand over a domain in Real^2.
struct CallerIntegrands {
  static const int max_n_integrand = 6;
  // Number of integrands. Must be <= max_n_integrand.
  virtual int nintegrands() const = 0;
  // Evaluate each integrand[j], j = 0 to nintegrands()-1, at the points
  // &p[2*i], i = 0 to n-1.
  virtual void eval(const int n, CRPtr p, RPtr integrand) const = 0;
  // Smallest R = sqrt(x^2 + y^2) < R_max at which this integrand can be
  // evaluated. This is for the case where there is a singularity at (0,0). In
  // eval, the caller should be careful to avoid cancellation and truncation
  // errors related to p's proximity to the singularity.
  virtual Real permitted_R_min(const Real R_max) const { return 0; }
};

} // namespace acorn
} // namespace woodland

#endif
