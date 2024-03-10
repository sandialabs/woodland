#ifndef INCLUDE_WOODLAND_EXAMPLES_CONVZX_CALC_STRESS
#define INCLUDE_WOODLAND_EXAMPLES_CONVZX_CALC_STRESS

#include "woodland/examples/convzx/discretization.hpp"

namespace woodland {
namespace examples {
namespace convzx {

struct Stress {
  struct QuadratureParams {
    int np_radial = 20;
    int np_angular = 20;
    int triquad_order = -1;

    QuadratureParams& set_np_radial(int v);
    QuadratureParams& set_np_angular(int v);
    QuadratureParams& set_triquad_order(int v);
  };

  struct Options {
    QuadratureParams qp;
    Real hfp_dist_fac = 1.0;
    bool use_calc_integral_tensor_quadrature = false;

    Options () {}

    Options& set_hfp_dist_fac(Real v);
  };

  struct Info {
    bool hfp;
    int triquad_order;
  };

  void init(const Discretization::CPtr& d);

  void init_dislocs(CRPtr dislocs_lcs);
  // Coefficients are arranged as [3*Discretization::reconstruct_ncoef]*ntri.
  const std::vector<Real>& get_dislocs_coefs () const { return coefs; }

  // Requires init_dislocs.
  void calc_s1_r1(const Real lam, const Real mu,
                  const int isrc, const int ircv, const bool self_interaction,
                  Real sigma[6],
                  const Options o = Options(), Info* info = nullptr) const;

  // Does not require init_dislocs. (isrc_dislocs, isrc_coefs) fully describe
  // the dislocation in cell isrc. These are arrays having sizes for just one
  // cell.
  void calc_s1_r1(CRPtr isrc_dislocs, CRPtr isrc_coefs,
                  const Real lam, const Real mu,
                  const int isrc, const int ircv, const bool hfp,
                  Real sigma[6],
                  const Options o = Options(), Info* info = nullptr) const;
 
  void calc_dist_over_L(const int isrc, RPtr doL) const;

private:
  enum : int { ndisloc = 3, ncoef = ndisloc*Discretization::reconstruct_ncoef };
  
  Discretization::CPtr d;
  const Real* dislocs_lcs;
  std::vector<Real> coefs;
};

} // namespace convzx
} // namespace examples
} // namespace woodland

#endif
