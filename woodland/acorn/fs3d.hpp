#ifndef INCLUDE_WOODLAND_ACORN_FS3D_SIGMA
#define INCLUDE_WOODLAND_ACORN_FS3D_SIGMA

#include "woodland/acorn/acorn.hpp"
#include "woodland/acorn/workspace.hpp"

namespace woodland {
namespace acorn {
namespace fs3d {

template <typename RealT>
void calc_sigma_point(
  // Lame coefficients.
  const Real lam, const Real mu,
  // Source.
  const RealT src[3],
  // Surface normal vector.
  const RealT nml[3],
  // Dislocation vector.
  const RealT disloc[3],
  // Receiver.
  const RealT rcv[3],
  // [sigma_11 sigma_12 sigma_13 sigma_22 sigma_23 sigma_33]
  RealT sigma[6],
  // If true, return norm(rcv - src, 2)^3 times the Green's function.
  const bool mult_by_R3 = false,
  // If provided, direction vector from src to rcv.
  const RealT* dir = nullptr);

struct RectInfo {
  bool hfp;
  Real L, dist[2];
  int triquad_order;
};

void calc_sigma_const_disloc_rect(
  Workspace& w,
  // Lame coefficients.
  const Real lam, const Real mu,
  // Source rectangle center.
  const Real src[3],
  // Source surface unit normal vector (local z).
  const Real nml[3],
  // Source surface unit tangent vector (local x).
  const Real tangent[3],
  // Rectangle local x- and y-direction lengths, with the rectangle centered on
  // the source center.
  const Real xy_side_lens[2],
  // Dislocation vector (in global coordinates).
  const Real disloc[3],
  // Receivers stored as a list of 3-vecs.
  const Real rcv[3],
  // [sigma_11 sigma_12 sigma_13 sigma_22 sigma_23 sigma_33]
  Real sigma[6],
  // For integrals::calc_hfp.
  const int hfp_np_radial = 8, const int hfp_np_angular = 6,
  // For integrals::calc_integral. If -1, set based on distance.
  int triquad_order = -1,
  // If triquad_order = -1, then use this tol in get_triquad_order.
  const Real triquad_tol = 1e-10,
  RectInfo* info = nullptr);

// Use Okada's dc3d to compute sigma.
void calc_sigma_const_disloc_rect_okada(
  const Real lam, const Real mu, const Real src[3], const Real nml[3],
  const Real tangent[3], const Real xy_side_lens[2], const Real disloc[3],
  const Real rcv[3], Real sigma[6]);

bool time_calc_sigma_point(const int n, const bool verbose = true);

void study_triquad();

void make_figure_data();

int unittest();

} // namespace fs3d
} // namespace acorn
} // namespace woodland

#endif
