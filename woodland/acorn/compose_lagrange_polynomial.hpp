// Derived from https://github.com/E3SM-Project/COMPOSE; see
//     https://github.com/E3SM-Project/COMPOSE/blob/main/LICENSE
// for details of the 3-clause BSD license.

#ifndef INCLUDE_WOODLAND_ACORN_LAGRANGE_POLYGON
#define INCLUDE_WOODLAND_ACORN_LAGRANGE_POLYGON

#include "woodland/acorn/acorn.hpp"

namespace woodland {
namespace acorn {

template <typename Scalar>
Real eval_lagrange_poly (const int& n, const Scalar* xsup, const Scalar* ysup,
                         const Scalar& x) {
  Scalar y = 0;
  for (int i = 0; i < n; ++i) {
    Scalar f = 1;
    for (int j = 0; j < n; ++j)
      f *= (i == j) ?
        1 :
        (x - xsup[j]) / (xsup[i] - xsup[j]);
    y += f*ysup[i];
  }
  return y;
}

template <typename Scalar>
Real eval_lagrange_poly_derivative (const int& np, const Scalar* const xsup,
                                    const Scalar* const ysup, const Scalar& x) {
  Scalar y = 0;
  for (int i = 0; i < np; ++i) {
    Scalar f = 0;
    for (int j = 0; j < np; ++j) {
      if (j == i) continue;
      Scalar g = 1;
      for (int k = 0; k < np; ++k)
        g *= ((k == i) ? 1 :
              ((k == j ? 1 : (x - xsup[k])) /
               (xsup[i] - xsup[k])));
      f += g;
    }
    y += f*ysup[i];
  }
  return y;
}

template <typename Scalar>
void eval_lagrange_poly_basis (const int& n, const Scalar* const xsup,
                               const Scalar& x, Scalar* const y) {
  for (int i = 0; i < n; ++i) {
    Scalar f = 1;
    for (int j = 0; j < n; ++j)
      f *= (j == i ?
            1 :
            (x - xsup[j]) / (xsup[i] - xsup[j]));
    y[i] = f;
  }
}

} // namespace acorn
} // namespace woodland

#endif
