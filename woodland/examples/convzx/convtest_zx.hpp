#ifndef INCLUDE_WOODLAND_EXAMPLES_CONVZX_CONVTEST_ZX
#define INCLUDE_WOODLAND_EXAMPLES_CONVZX_CONVTEST_ZX

#include "woodland/examples/convzx/discretization.hpp"
#include "woodland/examples/convzx/exact.hpp"
#include "woodland/examples/convzx/gallery.hpp"

namespace woodland {
namespace examples {
namespace convzx {

using gallery::ZxFn;
using gallery::Disloc;

typedef std::vector<Real> RealArray;

/* Convergence test for surfaces that are extruded 2D curves z(x,y) for (x,y) in
   [0,1]^2. z(x,y) = const for fixed x. We can discretize such surfaces using
   the Okada routine to provide a reference solution.

   (nx,ny) is the number of rectangles (not the number of rectangle boundary
   points, (nx+1,ny+1)).
 */

struct ConvTest {
  typedef std::shared_ptr<ConvTest> Ptr;
  typedef std::shared_ptr<const ConvTest> CPtr;

  enum class EvalMethod { direct, fast, direct_hmmvp, fast_hmmvp };
  static bool is_fast (const EvalMethod m) {
    return m == EvalMethod::fast || m == EvalMethod::fast_hmmvp;
  }
  static bool is_hmmvp (const EvalMethod m) {
    return m == EvalMethod::direct_hmmvp || m == EvalMethod::fast_hmmvp;
  }

  // Set the problem, without specifying the grid yet.
  void init(const ZxFn::Shape shape, const Disloc::CPtr& disloc);

  // Set grid parameters:
  //   Set nx in ZxFn.
  void set_nx(const int nx);
  //   Override ZxFn's automatically computed ny.
  void set_ny(const int ny);
  //   4 triangles per rectangle? Otherwise, 2.
  //   Set matrix-compression relative tolerance if hmmvp is used. The default
  //   value is 1e-6.
  void set_hmmvp_tol(const Real tol);
  //   For looking at the effects of various things.
  void set_use_four_tris_per_rect(const bool use);
  void set_use_surface_recon(const bool use);
  void set_use_exact_normals(const bool use);
  //     2, 4
  void set_normal_recon_order(const int order);
  //     0, 1, 2, 3
  void set_disloc_order(const int order);
  void set_use_flat_elements(const bool use);
  void set_general_lam_mu () { lam = 0.9; mu = 1.1; }
  //   Use Woodland impl of flat rect, constant disloc instead of Okada.
  void set_use_woodland_rg0c0(const bool use);

  // 0 for none.
  void set_verbosity (const int level) { verbosity = level; }

  ZxFn::CPtr get_zxfn () const { return zxfn; }
  Disloc::CPtr get_disloc () const { return disloc; }
  // Make sure general (lam,mu) work.
  Real get_lam () const { return lam; }
  Real get_mu () const { return mu; }
  // nx, ny are the number of rectangles in each dimension.
  int get_nx () const { return zxfn->get_nx(); }
  int get_ny () const { return zxfn->get_ny(); }
  Real get_hmmvp_tol () const { return hmmvp_tol; }

  void print(FILE* fp = stdout) const;

  // The following routines may be called only after init and set_nx are called.

  // Evaluate sigma at the triangle centers.
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
                  RealArray& sigmas, const bool use_hmmvp = false);
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

  void eval_exact_at_tri_ctrs(RealArray& sigmas, Exact::Options* o = nullptr);
  void eval_exact_at_rect_ctrs(const int nx, const int ny, RealArray& sigmas,
                               Exact::Options* o = nullptr);

  static int unittest();

private:
  Real lam = 1, mu = 1;
  int ntri_per_rect = 2;
  bool use_surface_recon = true;
  bool use_exact_normals = false;
  bool use_flat_elements = false;
  int nml_recon_order = 4;
  int disloc_order = 2;
  bool use_woodland_rg0c0 = false;
  ZxFn::Ptr zxfn;
  Disloc::CPtr disloc;
  Triangulation::Ptr t;
  Discretization::Ptr d;
  std::string python_fname;
  int verbosity = 0;
  Real hmmvp_tol = 1e-6;
  mutable Exact::Ptr exact;

  void eval_direct(const RealArray& dislocs, RealArray& sigmas) const;
  void eval_fast  (const RealArray& dislocs, RealArray& sigmas) const;
  void eval_okada_fast(const int nx, int& ny, const SupportPoints& supports,
                       RealArray& dislocs, RealArray& sigmas) const;
  void init_okada_hmmvp(const EvalMethod eval_method, const ZxFn& zxfn,
                        const int nyr) const;
  void eval_okada_hmmvp(const int nx, int& ny, RealArray& dislocs,
                        RealArray& sigmas) const;
  void eval_okada_fast_hmmvp(const int nx, int& ny, const SupportPoints& supports,
                             RealArray& dislocs, RealArray& sigmas) const;

public: // for unit tests
  void discretize();
  static int test_interp();

  static Triangulation::Ptr
  triangulate(const ZxFn& zxfn, const int ntri_per_rect = 2);

  static Discretization::Ptr
  discretize(const Triangulation::Ptr& t,
             const ZxFn::Shape shape,
             const bool use_surface_recon = true,
             const bool use_exact_normals = false,
             const int nml_recon_order = 4,
             const bool use_flat_elements = false);

public: // for external implementations
  static void calc_rect_ctr(const ZxFn& zxfn, const int nyr, const int irect,
                            Real p[3], RPtr nml = nullptr, RPtr lengths = nullptr);
  static void fill_dislocs(const ZxFn& zxfn, const int nyr, const Disloc& disloc,
                           RealArray& dislocs);
  static void interp_okada_to_okada(
    const ZxFn& zxfn_src, const int nyr_src, const RealArray& sigmas_src,
    const ZxFn& zxfn_dst, const int nyr_dst, RealArray& sigmas_dst);

  struct Hmmvp;
  mutable std::shared_ptr<Hmmvp> hmmvp;
};

void run_case(const std::string& params);
void convtest_w_vs_e(const std::string& params);
void convtest_o_vs_e(const std::string& params);

} // namespace convzx
} // namespace examples
} // namespace woodland

#endif
