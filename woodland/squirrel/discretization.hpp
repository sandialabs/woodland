#ifndef INCLUDE_WOODLAND_SQUIRREL_DISCRETIZATION
#define INCLUDE_WOODLAND_SQUIRREL_DISCRETIZATION

#include "woodland/squirrel/mesh.hpp"

#include <memory>

namespace woodland {
namespace squirrel {

struct Surface {
  typedef std::shared_ptr<Surface> Ptr;
  typedef std::shared_ptr<const Surface> CPtr;

  virtual ~Surface () {}

  virtual const char* name () const { return "Surface"; };

  virtual Idx  get_ncell() const = 0;
  virtual int  cell_nvtx(const Idx ci) const = 0;
  virtual void cell_vtxs_uv(const Idx ci, RPtr uv    ) const = 0;
  virtual void cell_ctr_uv (const Idx ci, Real uv [2]) const = 0;
  virtual void cell_ctr_xyz(const Idx ci, Real xyz[3]) const = 0;
  // Local coordinate system (LCS).
  virtual void cell_ctr_lcs(const Idx ci, Real lcs[9]) const = 0;

  // True if the Surface can provide the Jacobian matrix output in
  // cell_position.
  virtual bool supports_J () const { return false; }
  
  // Given a polygon and a local coordinate, return position data.
  void cell_position1 (
    const Idx ci, const Real uv[2],
    // Global position.
    Real xyz[3],
    // Surface-tangent local coordinate system associated with strike-slip,
    // dip-slip, and opening, respectively, ordered as (xhat, yhat, zhat) in a
    // 9-vector. zhat is the global surface normal.
    RPtr lcs = nullptr,
    // Optionally return the Jacobian determinant of the LCS -> GCS map.
    RPtr jacdet = nullptr,
    // 3x2 Jacobian matrix, row major: J = [x_u x_v; y_u y_v; z_u z_v].
    RPtr jac = nullptr) const
  { cell_position(ci, 1, uv, xyz, lcs, jacdet, jac); }

  // n instances, packed one after the other in the arrays.
  virtual void cell_position(
    const Idx ci, const int n, CRPtr uv, RPtr xyz, RPtr lcs = nullptr,
    RPtr jacdet = nullptr, RPtr jac = nullptr) const = 0;

  // Can cell_position1 be called with uv outside of element ci?
  virtual bool can_interpolate_outside_cell () const { return false; }
  // Even if !can_interpolate_outside_cell(), a meaningful LCS point can still
  // be computed in many cases. Return true if uv can be computed, false if not.
  virtual bool calc_cell_lcs(const Idx ci, const Real xyz[3], Real uv[2]) const
    { return false; }
};

struct Discretization {
  typedef std::shared_ptr<Discretization> Ptr;
  typedef std::shared_ptr<const Discretization> CPtr;
  typedef mesh::RealArray RealArray;

  Discretization(const mesh::Mesh::CPtr& m, const Surface::CPtr& srf,
                 const bool check_consistency = true);

  const mesh::Mesh::CPtr& get_mesh () const { return m; }
  const Surface::CPtr& get_surface () const { return srf; }

  enum : int { reconstruct_ncoef = 9 };
  // Calculate the coefficients associated with the cell-center values
  // fns_ctr. fns_ctr is ntri nfn-vectors. coefs is nfn
  // reconstruct_ncoef-vectors.
  void init_reconstruction(const Idx ci, const int nfn, CRPtr fns_ctr,
                           RPtr coefs) const;
  // Variant used in forming a matrix for the operator. nfn = 1. fns_ctr is 0
  // everywhere except in index fn_one_idx, at which its value is 1.
  void init_reconstruction(const Idx ci, const Idx fn_one_idx, RPtr coefs) const;

  // Evaluate the reconstruction at LCS (x,y).
  void reconstruct(const Idx ci, const int nfn, CRPtr coefs, CRPtr fns_ctr,
                   const Real uv[2], RPtr fns) const
  { reconstruct(ci, nfn, coefs, 1, fns_ctr, uv, fns); }
  void reconstruct(const Idx ci, const int nfn, CRPtr coefs, const int n,
                   // fns_ctr is for the whole grid; fns is for just cell ci.
                   CRPtr fns_ctr, CRPtr uv, RPtr fns) const;
  // Variant in which fns_ctr is just an array for cell ci, the same as
  // fns. Note that n does not apply to fns_ctr_ci.
  void reconstruct_ci(const Idx ci, const int nfn, CRPtr coefs, const int n,
                      CRPtr fns_ctr_ci, CRPtr uv, RPtr fns) const;

  // coefs arrays must be sized to hold reconstruct_ncoef, but the
  // reconstruction might use fewer coefficients, which is the following number.
  int get_disloc_ncoef() const;

  static int unittest();

  // For analysis. Values are 0 (constant), 1 (linear), 2 (quadratic).
  void set_disloc_order(const int order);

private:
  mesh::Mesh::CPtr m;
  Surface::CPtr srf;
  RealArray reconstruct_coefs;
  int disloc_order = 2;

  static int test_convergence(const bool verbose);
};

inline Discretization::Ptr
make_discretization (const mesh::Mesh::CPtr& m, const Surface::CPtr& srf,
                     const bool check_consistency = true) {
  return std::make_shared<Discretization>(m, srf, check_consistency);
}

// Curved-surface area in the domain bounded by cell ci.
Real cell_surface_area(Workspace& w, const Surface& s, const Idx ci);

} // namespace squirrel
} // namespace woodland

#endif
