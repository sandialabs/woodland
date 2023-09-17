#ifndef INCLUDE_WOODLAND_ACORN_LINALG
#define INCLUDE_WOODLAND_ACORN_LINALG

#include "woodland/acorn/acorn.hpp"

#include <memory>

namespace woodland {
namespace acorn {
namespace linalg {

// Matrices are column-major to follow LAPACK/BLAS convention.

// Factorize A(1:m,1:n), m >= n. A may be overwritten. R has size >=
// n*(n+1)/2. iwrk must have size >= n.
void qr_fac(const int m, const int n, RPtr A, RPtr R, int* iwrk);
// Solve A x ls= b. Q, R, iwrk are output from qr_fac. bx is b(1:m,1:nrhs) on
// input and x(1:n,1:nrhs) on output. rwrk has size >= n*nrhs;
void ls_slv(const int m, const int n, CRPtr Q, CRPtr R, const int* iwrk,
            const int nrhs, RPtr bx, RPtr rwrk);

int unittest();

} // namespace linalg
} // namespace acorn
} // namespace woodland

#endif
