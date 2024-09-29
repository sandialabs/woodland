#ifndef INCLUDE_WOODLAND_EXAMPLES_CONVZX_EXACT
#define INCLUDE_WOODLAND_EXAMPLES_CONVZX_EXACT

#include "woodland/examples/convzx/discretization.hpp"

#include <memory>

namespace woodland {
namespace examples {
namespace convzx {

// Exact, up to quadrature error, solution for a 2D parameterized surface on a
// rectangular reference domain. This works only if the dislocation is
// sufficiently smooth. It's for testing.
struct Exact {
  typedef std::shared_ptr<Exact> Ptr;
  typedef std::shared_ptr<const Exact> CPtr;  

  struct Description {
    typedef std::shared_ptr<Description> Ptr;
    typedef std::shared_ptr<const Description> CPtr;
    // Reference-domain rectangle.
    virtual void get_rectangle_limits(
      Real& xlo, Real& xhi, Real& ylo, Real& yhi) const = 0;
    // Compute the surface z(x,y), gradient, and LCS.
    virtual void get_surface(const Real x, const Real y,
                             Real& z, RPtr z_xy = nullptr,
                             RPtr lcs = nullptr) const = 0;
    // Dislocation in the LCS.
    virtual void get_disloc(const Real x, const Real y,
                            Real disloc[3]) const = 0;
  };

  struct Options {
    int nxr = 1, nyr = 1;
    int np_radial = 40, np_angular = 40;
    int triquad_order = 20;
    Real permitted_R_min_factor = 1e-3;
  };

  void init(const Description::CPtr& d);

  // Default is 1, 1.
  void set_lam_mu(const Real lam, const Real mu);

  void set_halfspace(const bool is_halfspace);

  void set_options(const Options& o);

  void calc_stress(const Real x, const Real y, Real sigma_lcs[6]) const;

  static int unittest();

private:
  struct Integrands;
  
  Description::CPtr d;
  Options o;
  Real lam = 1, mu = 1;
  bool halfspace = false;
  mutable Workspace w;
};

} // namespace convzx
} // namespace examples
} // namespace woodland

#endif
