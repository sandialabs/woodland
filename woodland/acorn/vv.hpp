#ifndef INCLUDE_WOODLAND_ACORN_VV
#define INCLUDE_WOODLAND_ACORN_VV

#include <cmath>
#include <memory>
#include <vector>
#include <functional>

#include "woodland/acorn/acorn.hpp"

namespace woodland {
namespace acorn {
namespace vv {

// Verification and validation problems.

// (u,v) in [0,1] x [0, 2 pi].
struct Surface {
  typedef std::shared_ptr<Surface> Ptr;
  typedef std::shared_ptr<const Surface> CPtr;

  virtual bool is_v_periodic () const = 0;

  virtual void get_position (const Real uv[2], Real xyz[3]) const = 0;
  virtual void get_normal   (const Real uv[2], Real nml[3]) const = 0;
  virtual void get_vhat     (const Real uv[2], Real tan[3]) const = 0;
  virtual void get_uhat     (                 Real uhat[3]) const = 0;
  // Jacobian determinant.
  virtual Real get_jacdet   (const Real uv[2]) const = 0;

  virtual void get_uv(const Real xyz[3], Real uv[2]) const = 0;

  // Area of the surface.
  virtual Real get_area() const = 0;  
};

struct FlatRectangle : public Surface {
  typedef std::shared_ptr<FlatRectangle> Ptr;
  typedef std::shared_ptr<const FlatRectangle> CPtr;

  bool is_v_periodic () const override { return false; }

  void get_position (const Real uv[2], Real xyz[3]) const override
  { xyz[0] = uv[0]; xyz[1] = uv[1]; xyz[2] = 0; }
  void get_normal   (const Real uv[2], Real nml[3]) const override
  { nml[0] = nml[1] = 0; nml[2] = 1; }
  void get_vhat     (const Real uv[2], Real tan[3]) const override
  { tan[0] = tan[2] = 0; tan[1] = 1; }
  void get_uhat     (                 Real uhat[3]) const override
  { uhat[0] = 1; uhat[1] = uhat[2] = 0; }
  // Jacobian determinant.
  Real get_jacdet   (const Real uv[2]) const override { return 1; }

  void get_uv(const Real xyz[3], Real uv[2]) const override
  { uv[0] = xyz[0]; uv[1] = xyz[1]; }

  // Area of the surface.
  Real get_area() const override { return 2*M_PI; }
};

// Generalized cylinder centered at origin. u in [0,1] is along the GC; v in [0,
// 2pi) is around it. 'rotation' is in radians and gives a rotation about the x
// axis.
struct GeneralizedCylinder : public Surface {
  typedef std::shared_ptr<GeneralizedCylinder> Ptr;
  typedef std::shared_ptr<const GeneralizedCylinder> CPtr;
  
  enum Shape : int { circle = 0, ellipse, invalid };

  static Shape convert(const int i);
  static Shape convert(const std::string& shape);
  static std::string convert(const Shape shape);
  static bool is_valid(const Shape shape);

  typedef GeneralizedCylinder Self;

  GeneralizedCylinder();

  Self& set_shape (const Shape value) { shape  = value; return *this; }
  // Radius and ellipse parameters are independent and used only for certain
  // shapes.
  Self& set_radius (const Real value) { radius = value; return *this; }
  Self& set_ellipse(const Real a, const Real b);
  Self& set_length (const Real value) { length = value; return *this; }
  // Rotation from z axis around y in radians.
  Self& set_rotation (const Real value) { rotation = value; return *this; };

  Shape get_shape   () const { return shape; }
  Real get_radius   () const { return radius; }
  Real get_length   () const { return length; }
  Real get_rotation () const { return rotation; }

  bool is_v_periodic () const override { return true; }

  void get_position (const Real uv[2], Real xyz[3]) const override;
  void get_normal   (const Real uv[2], Real nml[3]) const override;
  void get_vhat     (const Real uv[2], Real tan[3]) const override;
  void get_uhat     (                 Real uhat[3]) const override;
  Real get_jacdet   (const Real uv[2]) const override;

  void get_uv(const Real xyz[3], Real uv[2]) const override;

  // Area of the GC excluding the ends.
  Real get_area() const override;

private:
  Shape shape;
  Real length, radius, rotation, ell_a, ell_b;
};

// Tessellation of a GC by (u,v)-aligned rectangles. The rectangle corners are
// on the GC. It might be better to make it so that each rectangle's cell center
// is tangent to the GC, but it's much harder to compute such a
// tessellation. For testing purposes, use the simple approach.
struct RectTessellation {
  typedef std::shared_ptr<RectTessellation> Ptr;
  typedef std::shared_ptr<const RectTessellation> CPtr;

  void init(const GeneralizedCylinder& gc, const int nucell, const int nvcell);

  int ncell() const;
  int get_nucell() const;
  int get_nvcell() const;
  int icell (const int iu, const int iv) const { return nvcell*iu + iv; }
  const Real* get_cellctr(const int icell) const;
  const Real* get_normal(const int icell) const;
  const Real* get_uhat(const int icell) const;
  const Real* get_vhat(const int icell) const;
  // Rectangle dimensions: (along-u length, along-v length) in (x,y,z) space.
  Real get_ulength(const int icell) const;
  Real get_vlength(const int icell) const;
  // Corners of rectangle in this order: (-u,-v), (-u,v), (u,v), (u,-v).
  void get_rect(const int icell, Real corners[12]) const;
  void get_rect_uv(const int icell, Real corners[8]) const;

private:
  int nucell, nvcell;
  Real ulength, uhat[3];
  std::vector<Real> vlengths, ctrs, nmls, vhats;
};

void global2local_sym_tensor(const RectTessellation& rt, RPtr sigmas);
void global2local_sym_tensor(const RectTessellation& rt, const int idx,
                             Real sigmas[6]);

struct CallerUvFn {
  static const int max_n_fn = 3;
  // Number of functions. Must be <= max_n_fn.
  virtual int nfunctions() const = 0;
  // Evaluate each f[i] at (u,v).
  virtual void eval(const Real uv[2], RPtr f) const = 0;
};

struct SigmaType {
  enum Enum : int { fs_okada = 0, fs, fs_exact_geometry_disloc, invalid };

  static Enum convert(const int i);
  static Enum convert(const std::string& e);
  static std::string convert(const Enum e);
  static bool is_valid(const Enum e);
};

struct QuadratureParams {
  int np_radial, np_angular, triquad_order;
  QuadratureParams () : np_radial(12), np_angular(12), triquad_order(-1) {}
  QuadratureParams& set_np_radial (int v) { np_radial = v; return *this; }
  QuadratureParams& set_np_angular (int v) { np_angular = v; return *this; }
  QuadratureParams& set_triquad_order (int v) { triquad_order = v; return *this; }
};

// Compute sigma at rect centers due to disloc. disloc.nfunctions() must be 3.
// sigmas is 6*rt.ncell(). dislocs is optional; if non-null, it's 3*rt.ncell().
void eval_sigma_at_ctrs(const GeneralizedCylinder& gc, const RectTessellation& rt,
                        const Real lam, const Real mu, const CallerUvFn& fdisloc,
                        const SigmaType::Enum stype, RPtr sigmas,
                        RPtr dislocs = nullptr,
                        const QuadratureParams qp = QuadratureParams());

void eval_sigma(const GeneralizedCylinder& gc, const RectTessellation& rt,
                const Real lam, const Real mu, const CallerUvFn& fdisloc,
                // src and rcv cell indices
                const int isrc, const int ircv,
                // nxn grid of rcv points in the rcv cell, ordered along
                // v, then u
                const int nurcv, const int nvrcv,
                const SigmaType::Enum stype,
                // src dislocs
                Real disloc[3],
                // sigmas in the rcv cell
                RPtr sigmas,
                const QuadratureParams qp = QuadratureParams());

struct SigmaExactOptions {
  QuadratureParams qp;
  int extrap_npt;
  Real hfp_dist_fac;

  SigmaExactOptions();
};

// Evaluate sigma based on exact dislocation field and geometric surface over
// the source rectangle.
void eval_sigma_exact(const GeneralizedCylinder& gc, const RectTessellation& rt,
                      const Real lam, const Real mu, const CallerUvFn& fdisloc,
                      // src cell
                      const int isrc,
                      // rcvs: points (u,v)*nuv
                      const int nuv, CRPtr uvs,
                      RPtr sigmas,
                      const SigmaExactOptions o = SigmaExactOptions(),
                      const bool omp_thread = false);

// Evaluate the source, described by rt_src and exact dislocation and geometry,
// at rt_rcv (u,v) locations corresponding to cell centers.
//   This routine can be used in a convergence test that compares
// eval_sigma_at_ctrs with rt = rt_rcv against this routine, as rt_rcv is
// refined.
void eval_sigma_exact_at_ctrs(
  const GeneralizedCylinder& gc,
  const RectTessellation& rt_src, const RectTessellation& rt_rcv,
  const Real lam, const Real mu, const CallerUvFn& fdisloc,
  // packed according to rt_rcv
  RPtr sigmas,
  const SigmaExactOptions o = SigmaExactOptions());

struct TestUvFn : public CallerUvFn {
  // Function shape, all with domain [0, 2 pi).
  //                               parameters
  enum Fn : int { constant = 0, // (value)
                  cosine_bell,  // (max value, center, width)
                  tapered,      // (max value, center, width)
                  hat,          // (max value, center, width)
                  sin,          // (a, f, p) in a sin(f x + p)
                  invalid };

  static Fn convert(const int i);
  static Fn convert(const std::string& fn);
  static std::string convert(const Fn fn);
  static bool is_valid(const Fn fn);

  TestUvFn () {}

  // The dislocation field is a tensor product of u- and v-direction functions.
  // Each function takes three params, some possibly unused. The first component
  // is along v, the second along u, the third normal.
  TestUvFn(const Fn ufns[3], const Real uparms[9],
           const Fn vfns[3], const Real vparms[9]);

  void init(const Fn ufns[3], const Real uparms[9],
            const Fn vfns[3], const Real vparms[9]);

  int nfunctions() const override;
  void eval(const Real uv[2], RPtr f) const override;

private:
  Fn ufns[3], vfns[3];
  Real uparms[9], vparms[9];

  static Real eval(const Fn fn, const Real* parms, Real x);
};

namespace convtest {

// For unit tests.
struct Results {
  // Relevant stress components.
  std::vector<int> components;
  // Errors at finest resolution.
  std::vector<std::vector<Real> > l2_errs, linf_errs;
};

namespace gencyl {

struct Config {
  GeneralizedCylinder gc;
  Real lam, mu;
  // If onedim is true, then in the u direction, everything is constant, and nu
  // = 1 always; i.e., refinement is done ony in the v direction.
  bool onedim;
  int n_u_cell_base, n_v_cell_base, n_refine;
  TestUvFn fdisloc;
  std::string python_outfn; // default empty: no output

  enum Problem {p_circle, p_cylinder, p_ellipse, p_ellipse_cylinder};

  Config(const Problem problem = p_circle); // sets default values
};

Results run(const Config& c, const bool verbose = true);

} // namespace gencyl

namespace flatstrip {

struct Config {
  Real width, length;
  Real lam, mu;
  int n_cell_base, n_refine;
  TestUvFn fdisloc;
  std::string python_outfn; // default empty: no output
  bool run_fs; // if false, run just fs_exact_geometry_disloc

  Config(); // sets default values
};

Results run(const Config& c, const bool verbose = true);

} // namespace flatstrip

namespace impl {

// Accumulate errors for multiple models. nmodel includes the reference model,
// fs_okada. Thus, two test models means nmodel = 3.
struct Errors {
  template <typename T> using A = std::vector<T>;

  A<SigmaType::Enum> stypes; // (nmodel-1)
  A<A<Real>> sigmas; // (nmodel, 6 max_ncell)
  A<A<A<Real>>> l2_num, l2_den, li_num, li_den; // (nmodel-1, 6, nrefine)

  Errors(const int nmodel, const int nrefine, const int max_ncell);

  void accum_errors(const int ncell, const std::function<Real(int)>& cell_area,
                    const int irefine);
  void collect_errors(const A<int>& sigma_components,
                      Results* results = nullptr,
                      const bool verbose = true);
};

namespace pyout {

void init(FILE* const fid, const Errors& e);
void write(FILE* const fid, const Errors& e, const Errors::A<Real>& dislocs,
           const int ncell);

} // namespace pyout
} // namespace impl
} // namespace convtest

int unittest();

} // namespace vv
} // namespace acorn
} // namespace woodland

#endif
