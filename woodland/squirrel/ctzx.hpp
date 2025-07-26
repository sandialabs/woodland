#ifndef INCLUDE_WOODLAND_SQUIRREL_CTZX
#define INCLUDE_WOODLAND_SQUIRREL_CTZX

#include "woodland/squirrel/discretization.hpp"
#include "woodland/squirrel/exact.hpp"
#include "woodland/squirrel/gallery.hpp"
#include "woodland/squirrel/calc_stress.hpp"
#include "woodland/squirrel/ctzx_surface.hpp"
#ifdef WOODLAND_HAVE_HMMVP
# include "woodland/squirrel/hmatrix.hpp"
#endif

namespace woodland {
namespace squirrel {
namespace ctzx {

using gallery::ZxFn;
using gallery::Disloc;

typedef std::vector<Real> RealArray;

/* Convergence test for surfaces that are extruded 2D curves z(x,y) for (x,y) in
   [0,1]^2. z(x,y) = f(x). We can discretize such surfaces using the Okada
   routine to provide a reference solution.

   (nx,ny) is the number of rectangles (not the number of rectangle boundary
   points, (nx+1,ny+1)).
 */

struct ConvTest {
  typedef std::shared_ptr<ConvTest> Ptr;
  typedef std::shared_ptr<const ConvTest> CPtr;

  enum class EvalMethod { direct, fast, direct_hmmvp };
  static bool is_fast (const EvalMethod m) {
    return m == EvalMethod::fast;
  }
  static bool is_hmmvp (const EvalMethod m) {
    return m == EvalMethod::direct_hmmvp;
  }

  // Set the problem, without specifying the grid yet.
  void init(const ZxFn::Shape shape, const Disloc::CPtr& disloc);

  // Set grid parameters:
  //   Set nx in ZxFn.
  void set_nx(const int nx);
  //   Override ZxFn's automatically computed ny.
  void set_ny(const int ny);
  //   For looking at the effects of various things.
  void set_use_four_tris_per_rect(const bool use);
  void set_use_surface_recon(const bool use);
  void set_use_exact_tangents(const bool use);
  void set_use_c2_spline(const bool use);
  void set_use_halfspace(const bool use);
  //     2, 4
  void set_tangent_recon_order(const int order);
  //     0, 1, 2, 3
  void set_disloc_order(const int order);
  void set_use_flat_elements(const bool use);
  void set_general_lam_mu () { gfp.lam = 0.9; gfp.mu = 1.1; }
  //   Use Woodland impl of flat rect, constant disloc instead of Okada.
  void set_use_woodland_rg0c0(const bool use);
  //   Use a nonunirect discretization rather than triangles.
  void set_use_nonunirect(const bool use);
  //   Discretize uniformly in x rather than arc length.
  void set_xuniform(const bool use);

  // 0 for none.
  void set_verbosity (const int level) { verbosity = level; }

  ZxFn::CPtr get_zxfn () const { return zxfn; }
  Disloc::CPtr get_disloc () const { return disloc; }
  // Make sure general (lam,mu) work.
  Real get_lam () const { return gfp.lam; }
  Real get_mu () const { return gfp.mu; }
  bool get_use_halfspace () const { return gfp.halfspace; }
  const GreensFnParams& get_gfp () const { return gfp; }
  // nx, ny are the number of rectangles in each dimension.
  int get_nx () const { return zxfn->get_nx(); }
  int get_ny () const { return zxfn->get_ny(); }

  void print(FILE* fp = stdout) const;

  // The following routines may be called only after init and set_nx are called.

  // Evaluate sigma at the cell centers.
  void eval(RealArray& dislocs, RealArray& sigmas,
            const EvalMethod eval_method = EvalMethod::fast);

  // Valid only after calling eval().
  Discretization::CPtr get_discretization () const { return d; }
  
  // Evaluate sigma on a rectangular grid having ny x nx rectangles. Outputs are
  // ordered with x the fast direction. If ny > 0, use it; otherwise, set it to
  // ZxFn::get_ny() for the given nx.
  void eval_okada(const int nx, int& ny, RealArray& dislocs, RealArray& sigmas,
                  const EvalMethod eval_method = EvalMethod::fast) const;

  // Interpolate the Okada solution to the triangulation.
  void interp_okada(
    // src is from eval_okada
    const int nx, const int ny, const RealArray& sigmas,
    // dst is same as eval
    RealArray& isigmas) const;

  // Collect the (nx,ny)-rectangular tessellation's support points needed for
  // the current triangulation.
  struct SupportPoint {
    int ix, iy;
    bool operator< (const SupportPoint& o) const {
      if (ix < o.ix) return true;
      if (ix > o.ix) return false;
      return iy < o.iy;
    }
  };
  typedef std::vector<SupportPoint> SupportPoints;
  void collect_rect_support_points(const int nx, const int ny,
                                   SupportPoints& supports);
  // Evaluate the Okada solution at the support points.
  void eval_okada(const int nx, const int ny, const SupportPoints& supports,
                  RealArray& sigmas);
  // Interpolate to the triangulation.
  void interp_okada(const int nx, const int ny, const SupportPoints& supports,
                    const RealArray& sigmas, RealArray& isigmas);

  // Write data to a python module. All variables are key-values in 'dict'. If
  // append, open the file in append mode. One of dislocs and sigmas may be
  // empty.
  void pywrite(const std::string& python_filename, const RealArray& dislocs,
               const RealArray& sigmas, const std::string& dict,
               const bool append = true) const;
  void pywrite_okada(const std::string& python_filename, const int nx, const int ny,
                     const RealArray& dislocs, const RealArray& sigmas,
                     const std::string& dict, const bool append = true) const;

  void eval_exact_at_cell_ctrs(RealArray& sigmas, Exact::Options* o = nullptr);
  void eval_exact_at_rect_ctrs(const int nx, const int ny, RealArray& sigmas,
                               Exact::Options* o = nullptr);

#ifdef WOODLAND_HAVE_HMMVP
  const hmat::Hmatrices::Ptr& get_hmatrices () { return hmats; }
#endif

public: // for unit tests

  struct DiscretizeArgs {
    bool use_surface_recon = true;
    bool use_exact_tangents = false;
    int tan_recon_order = 4;
    bool use_flat_elements = false;
    bool use_c2_spline = false;
    bool mesh_scale = true;
    bool support_req0 = true;
    DiscretizeArgs () {}
  };

private:
  GreensFnParams gfp;
  int ntri_per_rect = 2;
  int disloc_order = 2;
  bool use_woodland_rg0c0 = false;
  bool use_nonunirect = false;
  bool xuniform = false;
  DiscretizeArgs dargs;
  ZxFn::Ptr zxfn;
  Disloc::CPtr disloc;
  Triangulation::Ptr t;
  Discretization::Ptr d;
  std::string python_fname;
  int verbosity = 0;
  mutable Exact::Ptr exact;
  mutable Workspace w;
#ifdef WOODLAND_HAVE_HMMVP
  hmat::Hmatrices::Ptr hmats;
#endif

  void eval_mesh_direct(const RealArray& dislocs, RealArray& sigmas) const;
  void eval_tri_fast(const RealArray& dislocs, RealArray& sigmas) const;
  void eval_nonunirect_fast(const RealArray& dislocs, RealArray& sigmas) const;
  void eval_mesh_hmmvp(const RealArray& dislocs, RealArray& sigmas) const;
  void eval_okada_fast(const int nx, int& ny, const SupportPoints& supports,
                       RealArray& dislocs, RealArray& sigmas) const;

public: // for unit tests
  void discretize();
  static int test_interp();
  static int unittest();

  static Triangulation::Ptr
  triangulate(const ZxFn& zxfn, const int ntri_per_rect = 2);

  static mesh::Mesh::CPtr make_nonunirect_mesh(const ZxFn& zxfn);

  void set_support_req0(const bool use);
  void set_mesh_scale(const bool use);

  static Discretization::Ptr
  discretize(const mesh::Mesh::CPtr& m,
             const Triangulation::Ptr& t,
             const ZxFn::Shape shape,
             const DiscretizeArgs args = DiscretizeArgs());

public: // for external implementations
  static void calc_rect_ctr(const ZxFn& zxfn, const int nyr, const int irect,
                            Real p[3], RPtr nml = nullptr, RPtr lengths = nullptr);
  static void fill_dislocs(const ZxFn& zxfn, const int nyr, const Disloc& disloc,
                           RealArray& dislocs);
  static void interp_okada_to_okada(
    const ZxFn& zxfn_src, const int nyr_src, const RealArray& sigmas_src,
    const ZxFn& zxfn_dst, const int nyr_dst, RealArray& sigmas_dst);
};

void run_case(const std::string& params);
void convtest_w_vs_e(const std::string& params);
void convtest_o_vs_e(const std::string& params);

// o_nbr is set to o before filling it. Thus, set fields in o that will be
// common to both.
void set_standard_options(Stress::Options& o, Stress::Options& o_nbr);

} // namespace ctzx
} // namespace squirrel
} // namespace woodland

#endif
