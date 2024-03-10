#ifndef INCLUDE_WOODLAND_EXAMPLES_CONVZX_TRIANGULATION
#define INCLUDE_WOODLAND_EXAMPLES_CONVZX_TRIANGULATION

// Triangulation of one or more cracks embedded in a 3D space.

#include "woodland/examples/convzx/convzx.hpp"

#include <cassert>
#include <memory>
#include <vector>

namespace woodland {
namespace examples {
namespace convzx {

struct Triangulation {
  typedef std::shared_ptr<Triangulation> Ptr;
  typedef std::shared_ptr<const Triangulation> CPtr;

  Triangulation () : in_mod_phase(true) {}

  // Begin modifications to the triangulation. Until end_modifications is
  // called, pointers to triangles and vertices returned from get_vtx/tri(s) are
  // invalid. A Triangulation starts in the modification phase and can re-enter
  // it multiple times.
  void begin_modifications();
  // Add a vertex. Returns the vertex's index. If no vertex vector is given,
  // space is made for the vertex, and the caller can later fill in the data
  // using get_vtx(i).
  Size append_vtx(CRPtr vtx = nullptr);
  // Add a triangle. Returns the triangle's index.
  Size append_tri(const int tri[3]);
  Size append_tri(const Idx tri[3]);
  // End modifications to the triangulation. After this call, the get_vtx/tri(s)
  // pointers are valid.
  void end_modifications();

  Size get_ntri () const { return t2vs.size()/3; }
  Size get_nvtx () const { return vtxs.size()/3; }

  // Grid geometry. Geometry can be modified at any time.
  CRPtr get_vtx (const int i) const { asserti(i, get_nvtx()); return &vtxs[3*i]; }
  RPtr get_vtx (const int i) { asserti(i, get_nvtx()); return &vtxs[3*i]; }
  CRPtr get_vtxs () const { return get_vtx(0); }

  // Grid topology. Topology data can be incrementally added to during the
  // modifications phase, but already existing topology data cannot be modified
  // or removed.
  const Idx* get_tri (const int i) const { asserti(i, get_ntri()); return &t2vs[3*i]; }
  const Idx* get_tris () const { return get_tri(0); }

  static int unittest();

private:
  bool in_mod_phase;
  std::vector<Real> vtxs;
  std::vector<Idx> t2vs;

  void asserti(const Idx i, const Size n) const {
    assert(i >= 0);
    assert(i <  n);
    assert( ! in_mod_phase);
  }
};

} // namespace convzx
} // namespace examples
} // namespace woodland

#endif
