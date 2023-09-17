#ifndef INCLUDE_WOODLAND_ACORN_HFP
#define INCLUDE_WOODLAND_ACORN_HFP

#include "woodland/acorn/acorn.hpp"

namespace woodland {
namespace acorn {
namespace hfp {

/* Compute the Hadamard finite part of an integral in the following two cases,
   where p(x) is a smooth, nonsingular function.

     Case 1. With a < c < b, compute hfp int_a^b p(y)/(y-c)^2 dy.
     Case 2. Compute hfp int_0^R p(x)/x^2 dx.

   The methods are based on the following expansions.

   First, separate out the lowest-order terms of the polynomial:
       p(x+c) = p(c) + p'(c) x + (p(x+c) - p(c) - p'(c) x)
       p(r) = p(0) + p'(0) r + (p(r) - p(0) - p'(0) r)

   Second, handle each piece of the hfps:
     Case 1.
       a < c < b
       hfp int_a^b p(y)/(y-c)^2 dy
         = hfp int_(a-c)^(b-c) p(x+c)/x^2 dx       (x = y-c)
         = hfp int_(a-c)^(b-c) p (c)/x^2 dx +
           cpv int_(a-c)^(b-c) p'(c)/x   dx +
               int_(a-c)^(b-c) (p(x+c) - p(c) - p'(c) x)/x^2 dx
         = p (c) (1/(a-c) - 1/(b-c)) +
           p'(c) log((b-c)/(c-a))    +
           int_(a-c)^(b-c) (p(x+c) - p(c) - p'(c) x)/x^2 dx  (third term)
     Case 2.
       0 < R
       hfp int_0^R p(r)/r^2 dr
         = hfp int_0^R p (0)/r^2 dr +
           hfp int_0^R p'(0)/r   dr +
               int_0^R (p(r) - p(0) - p'(r) r)/r^2 dr
         = -p (0)/R     +
            p'(0) log R +
            int_0^R (p(r) - p(0) - p'(0) r)/r^2 dr           (third term)

   Third, in each case's third line of the final equality, a proper integral
   remains. Compute it using Gaussian quadrature.

   If c is outside of [a,b], it can still be useful to use this technique so
   that the quadrature's accuracy is not affected by the singularity. Thus,
   calls with c < a or c > b are permitted.
 */

struct Options {
  int gl_np;
  Options () : gl_np(20) {}
};

// Representation of p(x). It must be able to calculate p(c), p'(c), and p(x)
// for x in [a,b] accurately. Thus, if c < a, for example, the representation
// over region (c,a) need not be accurate.
struct CallerP {
  // Evaluate p(x) for x in [a,b].
  virtual Real eval(const Real x) const = 0;
  // Evaluate p(c) and p'(c).
  virtual void eval_p_c(const Real c, Real& pc, Real& pdc) const = 0;
  virtual bool in_support (const Real x) const { return true; }
};

// Representation of p(x) as a polynomial in Lagrange form.
struct CallerPolyP : public CallerP {
  int n;                // nodes
  const Real* xs, * ys; // nodes x and nodal values y(x)
  CallerPolyP();
  CallerPolyP(const int n, CRPtr xs, CRPtr ys);
  Real eval(const Real x) const override;
  void eval_p_c(const Real c, Real& pc, Real& pdc) const override;
  bool in_support(const Real x) const override;
};

Real calc_hfp(
  const Options& o, const CallerP& p,
  // a < c < b (case 1) or a = c < b (case 2)
  const Real a, const Real b, const Real c);

// Public in case it is of use.
Real integrate_third_term(
  const Options& o, const CallerP& p,
  // a <= c <= b
  const Real a, const Real b, const Real c);

int unittest();

} // namespace hfp
} // namespace acorn
} // namespace woodland

#endif
