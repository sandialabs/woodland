#ifndef INCLUDE_WOODLAND_ACORN_ELASTOSTATICS_INTEGRALS
#define INCLUDE_WOODLAND_ACORN_ELASTOSTATICS_INTEGRALS

#include "woodland/acorn/acorn.hpp"

namespace woodland {
namespace acorn {

// Return the triangle quadrature appropriate for a receiver at a distance from
// a polygon of a characteristic size.
int get_triquad_order(const Real size, const Real distance,
                      const Real tol = 1e-8);

} // namespace acorn
} // namespace woodland

#endif
