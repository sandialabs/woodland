#ifndef INCLUDE_WOODLAND_SQUIRREL_GREENS_FUNCTION
#define INCLUDE_WOODLAND_SQUIRREL_GREENS_FUNCTION

#include "woodland/squirrel/squirrel.hpp"

namespace woodland {
namespace squirrel {

struct GreensFnParams {
  Real lam = -1, mu = -1;
  bool halfspace = false;

  GreensFnParams () {}
  GreensFnParams (const Real lam_, const Real mu_, const bool hs)
    : lam(lam_), mu(mu_), halfspace(hs)
  {}
};

bool operator==(const GreensFnParams& a, const GreensFnParams& b);
bool operator!=(const GreensFnParams& a, const GreensFnParams& b);

} // namespace squirrel
} // namespace woodland

#endif
