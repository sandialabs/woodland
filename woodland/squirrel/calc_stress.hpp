#ifndef INCLUDE_WOODLAND_SQUIRREL_CALC_STRESS
#define INCLUDE_WOODLAND_SQUIRREL_CALC_STRESS

#include "woodland/squirrel/discretization.hpp"
#include "woodland/squirrel/greens_function.hpp"

namespace woodland {
namespace squirrel {

struct Stress {
  typedef std::shared_ptr<Stress> Ptr;
  typedef std::shared_ptr<const Stress> CPtr;

  struct QuadratureParams {
    int np_radial = 20;
    int np_angular = 20;
  };

  struct Options {
    // Quadrature parameters for h.f.p. and tensor-quadrature computations.
    QuadratureParams qp_hfp, qp_tq;
    // A rcv point is near a src cell if dist < near_dist_fac L, where L is the
    // maximum distance between any two vertices of src and dist is the distance
    // between the rcv and src.
    Real near_dist_fac = 1.0;
    // Use the h.f.p. routine if rcv is near src. Otherwise, use tensor
    // quadrature. (The h.f.p. routine is always used if the src and rcv cells
    // are the same.)
    bool use_hfp_if_near = false;
    // Force use of calc_integral_tensor_quadrature rather than
    // calc_integral. Useful for nearby interactions, particularly those of
    // neighbors.
    bool force_tensor_quadrature = false;
    // Force use of the h.f.p. routine. The two force_* options can't
    // simultaneously be true.
    bool force_hfp = false;
    // Order of triangulation quadrature in calc_integral. Set to -1 for
    // automatic estimate.
    int triquad_order = -1;
    // 'tol' argument to acorn::get_triquad_order.
    Real triquad_order_tol = 1e-8;

    Options () {}
  };

  struct Info {
    bool hfp;
    int triquad_order, triquad_order_halfspace;
  };

  enum : int { ndisloc = 3, ncoef = ndisloc*Discretization::reconstruct_ncoef };

  Stress(const Discretization::CPtr& d);

  void init_dislocs(CRPtr dislocs_lcs);
  // Coefficients are arranged as [ncoef]*ntri.
  const std::vector<Real>& get_dislocs_coefs () const { return coefs; }

  // Requires init_dislocs.
  void calc_s1_r1(const GreensFnParams& gfp,
                  const int isrc, const int ircv,
                  Real sigma[6],
                  const Options o = Options(), Info* info = nullptr) const;

  // Does not require init_dislocs. (isrc_dislocs, isrc_coefs) fully describe
  // the dislocation in cell isrc. These are arrays having sizes for just one
  // cell. isrc_coefs has size ncoef.
  void calc_s1_r1(const Real isrc_disloc[3], CRPtr isrc_coefs,
                  const GreensFnParams& gfp,
                  const int isrc, const int ircv,
                  Real sigma[6],
                  const Options o = Options(), Info* info = nullptr) const;

  // Compute entry (ircv, isrc) of the BEM matrices for dislocation component
  // disloc_comp and the six receiver-local stress components.
  void calc_matrix_entries(const GreensFnParams& gfp,
                           const Idx isrc, const Idx ircv,
                           const int disloc_comp, Real sigma_rcv_lcl[6],
                           const Options o = Options(),
                           // Use o_nbrs instead of o if rcv and src are
                           // neighboring cells.
                           const Options o_nbrs = Options()) const;

  void calc_dist_over_L(const int isrc, RPtr doL) const;

  static Stress::Ptr create (const Discretization::CPtr& d) {
    return std::make_shared<Stress>(d);
  }

private:  
  Discretization::CPtr d;
  const Real* dislocs_lcs = nullptr;
  mutable Workspace w;
  std::vector<Real> coefs;
};

} // namespace squirrel
} // namespace woodland

#endif
