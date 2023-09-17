// Derived from https://github.com/E3SM-Project/COMPOSE; see
//     https://github.com/E3SM-Project/COMPOSE/blob/main/LICENSE
// for details of the 3-clause BSD license.

#ifndef INCLUDE_WOODLAND_ACORN_QUADRATURE
#define INCLUDE_WOODLAND_ACORN_QUADRATURE

#include "woodland/acorn/acorn.hpp"

namespace woodland {
namespace acorn {

// Gauss-Lobatto-Legendre
bool is_gll_supported(const int np);
const Real* get_x_gll(const int np);
const Real* get_w_gll(const int np);
// Gauss-Legendre
bool is_gl_supported(const int np);
const Real* get_x_gl (const int np);
const Real* get_w_gl (const int np);

// Gauss-Lobatto-Legendre
extern Real x_gll_table[];
extern Real w_gll_table[];
// Gauss-Legendre
extern Real x_gl_table[];
extern Real w_gl_table[];

struct Quadrature {
  enum : int { max_nq = 40 };
  enum Type { gll, gl };

  static bool is_supported(const Type type, const int nq);

  int nq;
  Type type;

  Quadrature(const int nq = 16, const Type type = gll);

  void get_xw (const Real*& x, const Real*& w) const {
    if (type == gll) { x = get_x_gll(nq); w = get_w_gll(nq); }
    else { x = get_x_gl(nq); w = get_w_gl(nq); }
  }
  void get_x (const Real*& x) const {
    if (type == gll) { x = get_x_gll(nq); }
    else { x = get_x_gl(nq); }
  }
  void get_w (const Real*& w) const {
    if (type == gll) { w = get_w_gll(nq); }
    else { w = get_w_gl(nq); }
  }
};

} // namespace acorn
} // namespace woodland

#endif
