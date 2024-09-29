#ifndef INCLUDE_WOODLAND_ACORN_VECTORIZATION
#define INCLUDE_WOODLAND_ACORN_VECTORIZATION

#include "woodland/acorn/acorn.hpp"
#include <cmath>

#define woodland_acorn_pack_for(n) ompsimd for (int i = 0; i < n; ++i)

namespace woodland {
namespace acorn {

template <typename T, int N>
struct Pack {
  typedef Pack<T,N> Me;
  typedef T Scalar;
  enum : int { n = N };

  Scalar v[n];

  const T& operator[] (const int i) const { return v[i]; }
  T& operator[] (const int i) { return v[i]; }

#define woodland_acorn_pack_unary_op(op)                                \
  Me& operator op (const Me& p) { woodland_acorn_pack_for(n) v[i] op p.v[i]; } \
  Me& operator op (const Scalar& s) { woodland_acorn_pack_for(n) v[i] op s; } \
  Me& operator op (const int s) { woodland_acorn_pack_for(n) v[i] op s; }

  woodland_acorn_pack_unary_op(+=)
  woodland_acorn_pack_unary_op(-=)
  woodland_acorn_pack_unary_op(*=)
  woodland_acorn_pack_unary_op(/=)
#undef woodland_acorn_pack_unary_op

};

#define woodland_acorn_pack_binary_op(op)                       \
  template <typename T, int N>                                  \
  Pack<T,N> operator op (const Pack<T,N>& a,                    \
                         const Pack<T,N>& b) {                  \
    Pack<T,N> c;                                                \
    woodland_acorn_pack_for(N) c.v[i] = a.v[i] op b.v[i];       \
    return c;                                                   \
  }                                                             \
  template <typename T, int N>                                  \
  Pack<T,N> operator op (const Pack<T,N>& a,                    \
                         const typename Pack<T,N>::Scalar& b) { \
    Pack<T,N> c;                                                \
    woodland_acorn_pack_for(N) c.v[i] = a.v[i] op b;            \
    return c;                                                   \
  }                                                             \
  template <typename T, int N>                                  \
  Pack<T,N> operator op (const typename Pack<T,N>::Scalar& a,   \
                         const Pack<T,N>& b) {                  \
    Pack<T,N> c;                                                \
    woodland_acorn_pack_for(N) c.v[i] = a op b.v[i];            \
    return c;                                                   \
  }                                                             \
  template <typename T, int N>                                  \
  Pack<T,N> operator op (const Pack<T,N>& a,                    \
                         const int b) {                         \
    Pack<T,N> c;                                                \
    woodland_acorn_pack_for(N) c.v[i] = a.v[i] op b;            \
    return c;                                                   \
  }                                                             \
  template <typename T, int N>                                  \
  Pack<T,N> operator op (const int a,                           \
                         const Pack<T,N>& b) {                  \
    Pack<T,N> c;                                                \
    woodland_acorn_pack_for(N) c.v[i] = a op b.v[i];            \
    return c;                                                   \
  }

woodland_acorn_pack_binary_op(+)
woodland_acorn_pack_binary_op(-)
woodland_acorn_pack_binary_op(*)
woodland_acorn_pack_binary_op(/)

#undef woodland_acorn_pack_binary_op

template <typename T, int N>
Pack<T,N> operator- (const Pack<T,N>& a) {
  Pack<T,N> c;
  woodland_acorn_pack_for(N) c.v[i] = -a.v[i];
  return c;
}

template <typename T, int N>
Pack<T,N> pow (const Pack<T,N>& a, const Real& exponent) {
  Pack<T,N> c;
  woodland_acorn_pack_for(N) c.v[i] = std::pow(a.v[i], exponent);
  return c;
}

inline Real pow (const Real& a, const Real& exponent) {
  return std::pow(a, exponent);
}

typedef Pack<Real,WOODLAND_ACORN_PACK_N> RealPack;

} // namespace acorn
} // namespace woodland

#undef woodland_acorn_pack_for

#endif
