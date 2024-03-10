#ifndef INCLUDE_WOODLAND_EXAMPLES_CONVZX_DISCRETIZATION
#define INCLUDE_WOODLAND_EXAMPLES_CONVZX_DISCRETIZATION

#include "woodland/examples/convzx/convzx.hpp"
#include "woodland/examples/convzx/triangulation.hpp"
#include "woodland/examples/convzx/gallery.hpp"

#include <memory>

namespace woodland {
namespace examples {
namespace convzx {

struct TriangulationRelations {
  typedef std::shared_ptr<TriangulationRelations> Ptr;
  typedef std::shared_ptr<const TriangulationRelations> CPtr;
  typedef std::vector<Idx> IdxArray;

  TriangulationRelations(const Triangulation::CPtr& t);

  void make_t2ts2();
  
  const IdxArray& get_v2tsi  () const { return v2tsi;  }
  const IdxArray& get_v2ts   () const { return v2ts ;  }
  const IdxArray& get_t2tsi  () const { return t2tsi;  }
  const IdxArray& get_t2ts   () const { return t2ts ;  }
  const IdxArray& get_t2etsi () const { return t2etsi; }
  const IdxArray& get_t2ets  () const { return t2ets ; }
  const IdxArray& get_t2ts2i () const { return t2ts2i; }
  const IdxArray& get_t2ts2  () const { return t2ts2 ; }

private:
  Triangulation::CPtr t;
  // xi is the pointer array into x. Thus, x[xi[k]:xi[k+1]-1] is the k'th list.
  // Names: v: vertex, e: edge, t: triangle.
  IdxArray v2tsi, v2ts;
  IdxArray t2tsi, t2ts  ; // vertex-sharing nbrs excluding self
  IdxArray t2etsi, t2ets; //   edge-sharing nbrs excluding self
  IdxArray t2ts2i, t2ts2; // optional 2-halo excluding self
};

struct Surface {
  typedef std::shared_ptr<Surface> Ptr;
  typedef std::shared_ptr<const Surface> CPtr;

  virtual Size get_ntri() const = 0;
  virtual void tri_vtxs_uv(const Idx ti, Real tuv[6]) const = 0;
  virtual void tri_ctr_uv (const Idx ti, Real uv [2]) const = 0;
  virtual void tri_ctr_xyz(const Idx ti, Real xyz[3]) const = 0;
  // Local coordinate system (LCS).
  virtual void tri_ctr_lcs(const Idx ti, Real lcs[9]) const = 0;
  
  // Given a triangle and a barycentric coordinate, return position data.
  void tri_position1 (
    const Idx ti, const Real uv[2],
    // Global position.
    Real xyz[3],
    // Surface-tangent local coordinate system associated with strike-slip,
    // dip-slip, and opening, respectively, ordered as (xhat, yhat, zhat) in a
    // 9-vector. zhat is the global surface normal.
    RPtr lcs = nullptr,
    // Optionally return the Jacobian determinant of the LCS -> GCS map.
    Real* jacdet = nullptr) const
  { tri_position(ti, 1, uv, xyz, lcs, jacdet); }
  virtual void tri_position(
    const Idx ti, const int n, CRPtr uv, RPtr xyz, RPtr lcs = nullptr,
    Real* jacdet = nullptr) const = 0;

  // Can tri_position1 be called with uv outside of element ti?
  virtual bool can_interpolate_outside_element () const { return false; }
  // Even if !can_interpolate_outside_element(), a meaningful LCS point can
  // still be computed in many cases. Return true if uv can be computed, false
  // if not.
  virtual bool calc_tri_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const
    { return false; }
};

struct Discretization {
  typedef std::shared_ptr<Discretization> Ptr;
  typedef std::shared_ptr<const Discretization> CPtr;
  typedef std::vector<Real> RealArray;

  Discretization(const Triangulation::CPtr& t,
                 const TriangulationRelations::Ptr& tr,
                 const Surface::CPtr& srf);

  const Triangulation::CPtr& get_triangulation () const { return t; }
  TriangulationRelations::CPtr get_triangulation_relations () const { return tr; }
  const Surface::CPtr& get_surface () const { return srf; }

  enum : int { reconstruct_ncoef = 9 };
  // Calculate the coefficients associated with the cell-center values
  // fns_ctr. fns_ctr is ntri nfn-vectors. coefs is nfn
  // reconstruct_ncoef-vectors.
  void tri_reconstruct_fit(const Idx ti, const int nfn, CRPtr fns_ctr,
                           RPtr coefs) const;
  // Evaluate the reconstruction at LCS (x,y).
  void tri_reconstruct(const Idx ti, const int nfn, CRPtr coefs, CRPtr fns_ctr,
                       const Real uv[2], RPtr fns) const
  { tri_reconstruct(ti, nfn, coefs, 1, fns_ctr, uv, fns); }
  void tri_reconstruct(const Idx ti, const int nfn, CRPtr coefs, const int n,
                       // fns_ctr is for the whole grid, while fns is for just
                       // triangle ti.
                       CRPtr fns_ctr, CRPtr uv, RPtr fns) const;
  // Variant in which fns_ctr is just an array for triangle ti, the same as fns.
  void tri_reconstruct_ti(const Idx ti, const int nfn, CRPtr coefs, const int n,
                          CRPtr fns_ctr_ti, CRPtr uv, RPtr fns) const;

  // coefs arrays must be sized to hold reconstruct_ncoef, but the
  // reconstruction might use fewer coefficients, which is the following number.
  int get_disloc_ncoef() const;

  static int unittest();

  // For analysis. Values are 0 (constant), 1 (linear), 2 (quadratic).
  void set_disloc_order(const int order);

private:
  Triangulation::CPtr t;
  TriangulationRelations::Ptr tr;
  Surface::CPtr srf;
  RealArray reconstruct_coefs;
  int disloc_order = 2;

  static int test_convergence(const bool verbose);
};

struct ExactParam2DSurface : public Surface {
  typedef std::vector<Real> RealArray;

  enum class Name { zero, trig1, trig2, user };

  struct UserFn {
    typedef std::shared_ptr<UserFn> Ptr;
    typedef std::shared_ptr<const UserFn> CPtr;
    // g can be nullptr.
    virtual void eval(const Real x, const Real y, Real& f, RPtr g) const = 0;
  };

  ExactParam2DSurface(const Triangulation::Ptr& t, const Name name);
  ExactParam2DSurface(const Triangulation::Ptr& t, const UserFn::CPtr& ufn);

  // Set Triangulation's z vtx values to their correct ones.
  void reset_triangulation_z();

  Size get_ntri() const override;
  void tri_vtxs_uv(const Idx ti, Real tuv[6]) const override;
  void tri_ctr_uv (const Idx ti, Real uv [2]) const override;
  void tri_ctr_xyz(const Idx ti, Real xyz[3]) const override;
  void tri_ctr_lcs(const Idx ti, Real lcs[9]) const override;
  
  void tri_position(const Idx ti, const int n, CRPtr uv, RPtr xyz,
                    RPtr lcs = nullptr, Real* jacdet = nullptr)
    const override;

  bool can_interpolate_outside_element () const override { return true; }
  bool calc_tri_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const override;

public:
  // Triangle (x,y) vertices.
  struct XyTris {
    Real v0[2], v1[2], v2[2], ctr[2];

    Real* vtxs () { return v0; }
    const Real* vtxs () const { return v0; }
  };
  typedef std::vector<XyTris> XyTrisArray;

private:
  Triangulation::Ptr t;
  Name name;
  UserFn::CPtr ufn;
  RealArray tri_ctrs_gcs, tri_lcs_at_ctrs;
  XyTrisArray xy_tris;

  void init(const Triangulation::Ptr& t, const Name name,
            const UserFn::CPtr = nullptr);
};

// Tensor product of a cubic Bezier spline in the x direction and 1 in the y
// direction.
struct ExtrudedCubicSplineSurface : public Surface {
  ExtrudedCubicSplineSurface(const gallery::ZxFn::Shape zshape,
                             const Triangulation::CPtr& t,
                             const bool use_exact_normals,
                             const int nml_recon_order = 4);

  void init(const gallery::ZxFn::Shape zshape, const Triangulation::CPtr& t,
            const bool use_exact_normals,
            const int nml_recon_order = 4);

  Size get_ntri() const override;
  void tri_vtxs_uv(const Idx ti, Real tuv[6]) const override;
  void tri_ctr_uv (const Idx ti, Real uv [2]) const override;
  void tri_ctr_xyz(const Idx ti, Real xyz[3]) const override;
  void tri_ctr_lcs(const Idx ti, Real lcs[9]) const override;
  
  void tri_position(const Idx ti, const int n, CRPtr uv, RPtr xyz,
                    RPtr lcs = nullptr, Real* jacdet = nullptr)
    const override;

  bool calc_tri_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const override;

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

  Size get_ntri() const override;
  void tri_vtxs_uv(const Idx ti, Real tuv[6]) const override;
  void tri_ctr_uv (const Idx ti, Real uv [2]) const override;
  void tri_ctr_xyz(const Idx ti, Real xyz[3]) const override;
  void tri_ctr_lcs(const Idx ti, Real lcs[9]) const override;
  
  void tri_position(const Idx ti, const int n, CRPtr uv, RPtr xyz,
                    RPtr lcs = nullptr, Real* jacdet = nullptr)
    const override;

  bool can_interpolate_outside_element () const override { return false; }
  bool calc_tri_lcs(const Idx ti, const Real xyz[3], Real uv[2]) const override;

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

// Curved-surface area in the domain bounded by triangle ti.
Real tri_surface_area(const Surface& s, const Idx ti);

void estimate_vertex_normals(const Triangulation& t, const std::vector<Idx>& v2tsi,
                             const std::vector<Idx>& v2ts, RPtr vtx_nmls,
                             // If true, estimate vtx nml using just the
                             // adjacent triangle's normals.
                             const bool from_adjacent_facets_only = false);

} // namespace convzx
} // namespace examples
} // namespace woodland

#endif
