#ifndef INCLUDE_WOODLAND_EXAMPLES_CONVZX_GALLERY
#define INCLUDE_WOODLAND_EXAMPLES_CONVZX_GALLERY

#include "woodland/examples/convzx/convzx.hpp"

#include <memory>
#include <vector>

namespace woodland {
namespace examples {
namespace convzx {
namespace gallery {

// Function z(x) for x in [0,1].
struct ZxFn {
  typedef std::shared_ptr<ZxFn> Ptr;
  typedef std::shared_ptr<const ZxFn> CPtr;
  
  enum class Shape { zero, ramp, quadratic, trig0, trig1 };

  static Shape convert(const std::string& shape);
  static std::string convert(const Shape shape);

  ZxFn(const Shape shape_) : shape(shape_) {}

  Shape get_shape () const { return shape; }

  // Number of elements > 1.
  void set_nx(int nx);
  // Override default computed internally.
  void set_ny(int ny);

  Real calc_arclength(const Real x) const;

  // Get nx (n in init) and ny, where ny is computed to give roughly uniform
  // elements.
  int get_nx () const { return int(xs.size()) - 1; }
  int get_ny () const { return ny; }
  
  CRPtr get_xbs () const { return xs.data(); }
  CRPtr get_zbs () const { return zs.data(); }

private:
  const Shape shape;
  int ny;
  std::vector<Real> xs, zs;
};

void eval(const ZxFn::Shape shape, const Real x, Real& f, Real& g);

// Dislocation shapes as a function of (x,y).
struct Disloc {
  typedef std::shared_ptr<Disloc> Ptr;
  typedef std::shared_ptr<const Disloc> CPtr;

  enum class Shape { zero, tapered, pcosbell, stapered };

  static Shape convert(const std::string& shape);
  static std::string convert(const Shape shape);

  // Set the dislocation for dimension i (x,y,z: 0,1,2) in the local coordinate
  // system. z is opening.
  //   The dislocation is a function with reference domain (x,y) in [0,1]^2. All
  // arguments are in terms of (x,y). The dislocation is centered at ctr and has
  // axes xhat and zhat x xhat with lengths for each axis.
  //   If not set, the dislocation is 0.
  void set(const int dim, const Shape shape, const Real amplitude,
           const Real ctr_x, const Real ctr_y,
           const Real xhat_x, const Real xhat_y,
           const Real length_x, const Real length_y);

  void eval(const Real xy[2], Real disloc[3]) const;

  // Check that the dislocation is 0 on bdy [0,1]^2.
  bool is_boundary_zero(const Real threshold = 0) const;

  struct Data {
    Shape shape;
    Real amplitude, ctr[2], axis[2][2], length[2];
    Data () : shape(Shape::zero) {}
  };
  Data ds[3];
};

} // namespace gallery
} // namespace convzx
} // namespace examples
} // namespace woodland

#endif
