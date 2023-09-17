#include "woodland/acorn/linalg.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/dbg.hpp"

#include <limits>

namespace woodland {
namespace acorn {
namespace linalg {

void qr_fac (const int m, const int n, RPtr A, RPtr rwrk, int* iwrk) {
  assert(m >= n);
  Real* const R = rwrk;
  for (int j = 0, ptr = 0; j < n; ++j) {
    Real* const Aj = &A[j*m];
    for (int k = 0; k < j; ++k) {
      const Real* const Ak = &A[m*k];
      Real a = 0;
      for (int i = 0; i < m; ++i)
        a += Ak[i]*Aj[i];
      R[ptr+k] = a;
      for (int i = 0; i < m; ++i)
        Aj[i] -= a*Ak[i];
    }
    Real norm = 0;
    for (int i = 0; i < m; ++i)
      norm += square(Aj[i]);
    norm = std::sqrt(norm);
    R[ptr+j] = norm;
    for (int i = 0; i < m; ++i)
      Aj[i] /= norm;
    ptr += j+1;
  }
}

void ls_slv (const int m, const int n, CRPtr Q, CRPtr R, const int* iwrk,
             const int nrhs, RPtr bx, RPtr v) {
  assert(m >= n);
  for (int i = 0; i < n*nrhs; ++i)
    v[i] = 0;
  for (int j = 0; j < n; ++j) {
    const Real* const Qj = &Q[m*j];
    for (int r = 0; r < nrhs; ++r) {
      Real* const vr = &v[n*r];
      const Real* const br = &bx[m*r];
      Real a = 0;
      for (int i = 0; i < m; ++i)
        a += Qj[i]*br[i];
      vr[j] = a;
    }
  }
  for (int i = 0; i < n*nrhs; ++i)
    bx[i] = v[i];
  int ptr = (n*(n-1))/2;
  for (int j = n-1; j >= 0; --j) {
    for (int r = 0; r < nrhs; ++r) {
      Real* const xr = &bx[n*r];
      xr[j] = xr[j]/R[ptr+j];
      for (int i = 0; i < j; ++i)
        xr[i] -= R[ptr+i]*xr[j];
    }
    ptr -= j;
  }
  assert(ptr == 0);
}

static int test_qr () {
  const auto eps = std::numeric_limits<Real>::epsilon();
  int nerr = 0;
  for (const int nrhs : {1, 3})
    for (const int m : {7})
      for (const int n : {4, 7}) {
        std::vector<Real> A(m*n), Ac(m*n), Ar(m*n), R(n*(n+1)/2), bx(m*nrhs),
          bc(m*nrhs), res(m), rwrk(n*nrhs);
        std::vector<int> iwrk(n);
        for (int i = 0; i < m*n; ++i) Ac[i] = A[i] = urand() - 0.5;
        for (int i = 0; i < m*nrhs; ++i) bc[i] = bx[i] = urand() - 0.5;
        qr_fac(m, n, A.data(), R.data(), iwrk.data());
        ls_slv(m, n, A.data(), R.data(), iwrk.data(), nrhs, bx.data(),
               rwrk.data());
        for (int j = 0, ptr = 0; j < n; ++j) {
          for (int k = 0; k <= j; ++k)
            for (int i = 0; i < m; ++i)
              Ar[j*m+i] += A[k*m+i]*R[ptr+k];
          ptr += j+1;
        }
        const auto rerr = reldif(m*n, Ac.data(), Ar.data());
        if (rerr > 10*eps) ++nerr;
        for (int r = 0; r < nrhs; ++r) {
          const Real* const xr = &bx[n*r];
          const Real* const br = &bc[m*r];
          for (int i = 0; i < m; ++i)
            res[i] = -br[i];
          for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i)
              res[i] += xr[j]*Ac[m*j+i];
          for (int j = 0; j < n; ++j) {
            Real a = 0;
            for (int i = 0; i < m; ++i)
              a += Ac[m*j+i]*res[i];
            if (std::abs(a) > 10*eps) ++nerr;
          }
        }
      }
  return nerr;
}

int unittest () {
  int nerr = 0;
  nerr += test_qr();
  return nerr;
}

} // namespace linalg
} // namespace acorn
} // namespace woodland
