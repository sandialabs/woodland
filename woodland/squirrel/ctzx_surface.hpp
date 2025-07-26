#ifndef INCLUDE_WOODLAND_SQUIRREL_CTZX_SURFACE
#define INCLUDE_WOODLAND_SQUIRREL_CTZX_SURFACE

#include "woodland/squirrel/discretization.hpp"
#include "woodland/squirrel/mesh_based_surface.hpp"
#include "woodland/squirrel/gallery.hpp"
#include "woodland/squirrel/ctzx_triangulation.hpp"

namespace woodland {
namespace squirrel {
namespace ctzx {

// Tensor product of a cubic Bezier spline in the x direction and 1 in the y
// direction.
struct ExtrudedCubicSplineSurface : public Surface {
  enum class Recon { c1_tan_exact,
                     c1_tan_2,
                     c1_tan_4,
                     c2 };
  
  ExtrudedCubicSplineSurface(const gallery::ZxFn::Shape zshape,
                             const Triangulation::CPtr& t,
                             const Recon recon = Recon::c1_tan_4);

  void init(const gallery::ZxFn::Shape zshape, const Triangulation::CPtr& t,
            const Recon recon = Recon::c1_tan_4);

  Idx get_ncell() const override;
  int cell_nvtx(const Idx ci) const override { return 3; }
  void cell_vtxs_uv(const Idx ti, Real tuv[6]) const override;
  void cell_ctr_uv (const Idx ti, Real uv [2]) const override;
  void cell_ctr_xyz(const Idx ti, Real xyz[3]) const override;
  void cell_ctr_lcs(const Idx ti, Real lcs[9]) const override;

  bool supports_J() const override;
  
  void cell_position(const Idx ti, const int n, CRPtr uv, RPtr xyz,
                     RPtr lcs = nullptr, RPtr jacdet = nullptr,
                     RPtr jac = nullptr)
    const override;

  bool calc_cell_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const override;

  static int unittest();

private:
  typedef std::vector<Real> RealArray;

  struct LclTri {
    Real xmin, xmax, v0[2], v1[2], v2[2], ctr[2];
    int seg;

    Real* vtxs () { return v0; }
    const Real* vtxs () const { return v0; }
  };
  typedef std::vector<LclTri> LclTriArray;

  gallery::ZxFn::Shape zshape;
  Triangulation::CPtr t;
  RealArray tri_ctrs_gcs, tri_lcs_at_ctrs, cs;
  LclTriArray lcl_tris;
};

// For clarity, even though other Surface classes could accommodate a 'flat'
// option, make a separate class for it.
struct FlatElementSurface : public Surface {
  FlatElementSurface(const Triangulation::CPtr& t, const Real primary[3]);

  void init(const Triangulation::CPtr& t, const Real primary[3]);

  Idx get_ncell() const override;
  int cell_nvtx(const Idx ci) const override { return 3; }
  void cell_vtxs_uv(const Idx ti, Real tuv[6]) const override;
  void cell_ctr_uv (const Idx ti, Real uv [2]) const override;
  void cell_ctr_xyz(const Idx ti, Real xyz[3]) const override;
  void cell_ctr_lcs(const Idx ti, Real lcs[9]) const override;

  bool supports_J() const override;
  
  void cell_position(const Idx ti, const int n, CRPtr uv, RPtr xyz,
                     RPtr lcs = nullptr, RPtr jacdet = nullptr,
                     RPtr jac = nullptr)
    const override;

  bool can_interpolate_outside_cell () const override { return false; }
  bool calc_cell_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const override;

  static int unittest();

private:
  typedef std::vector<Real> RealArray;

  // Local cooordinate system for each triangle.
  struct LCS {
    // xhat is the LS projection of the primary vector onto the triangle.
    // zhat is the triangle normal.
    // yhat completes the LCS.
    Real xhat[3], yhat[3], zhat[3];

    Real* data () { return xhat; }
    const Real* data () const { return xhat; }
  };
  typedef std::vector<LCS> LCSArray;

  // Triangle vertices in LCS.
  struct VtxLCS {
    // v0,1,2 are the (x,y) coordinates in the LCS. z is the LCS z coordinate
    // common to all three vertices. b is the barycentric matrix.
    Real v0[2], v1[2], v2[2], z, ctr[2], b[4];

    Real* vtxs () { return v0; }
    const Real* vtxs () const { return v0; }
  };
  typedef std::vector<VtxLCS> VtxLCSArray;

  Triangulation::CPtr t;
  RealArray tri_ctrs_gcs;
  LCSArray lcss;
  VtxLCSArray vtx_lcss;
};

} // namespace ctzx
} // namespace squirrel
} // namespace woodland

#endif
