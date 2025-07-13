#ifndef INCLUDE_WOODLAND_SQUIRREL_CTZX_TRIANGULATION
#define INCLUDE_WOODLAND_SQUIRREL_CTZX_TRIANGULATION

#include "woodland/squirrel/mesh.hpp"

#include <cassert>
#include <memory>
#include <vector>

namespace woodland {
namespace squirrel {
namespace ctzx {

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
  Idx append_vtx(CRPtr vtx = nullptr);
  // Add a triangle. Returns the triangle's index.
  template <typename Int> Idx append_tri(const Int tri[3]);
  // End modifications to the triangulation. After this call, the get_vtx/tri(s)
  // pointers are valid.
  void end_modifications();

  Idx get_ncell () const { return t2vs.size()/3; }
  Idx get_nvtx () const { return vtxs.size()/3; }

  // Grid geometry. Geometry can be modified at any time.
  CRPtr get_vtx (const int i) const { asserti(i, get_nvtx()); return &vtxs[3*i]; }
  RPtr get_vtx (const int i) { asserti(i, get_nvtx()); return &vtxs[3*i]; }
  CRPtr get_vtxs () const { return get_vtx(0); }

  // Grid topology. Topology data can be incrementally added to during the
  // modifications phase, but already existing topology data cannot be modified
  // or removed.
  const Idx* get_tri (const int i) const { asserti(i, get_ncell()); return &t2vs[3*i]; }
  const Idx* get_tris () const { return get_tri(0); }

  static int unittest();

private:
  bool in_mod_phase;
  std::vector<Real> vtxs;
  std::vector<Idx> t2vs;

  void asserti(const Idx i, const Idx n) const {
    assert(i >= 0);
    assert(i <  n);
    assert( ! in_mod_phase);
  }
};

struct TriangulationRelations {
  typedef std::shared_ptr<TriangulationRelations> Ptr;
  typedef std::shared_ptr<const TriangulationRelations> CPtr;
  typedef std::vector<Idx> IdxArray;

  TriangulationRelations(const Triangulation::CPtr& t);

  void make_t2ts2();

  bool have_edge_sharing_lists () const { return true; }
  
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
  IdxArray t2tsi, t2ts  ; // vertex-sharing nbrs excluding self
  IdxArray t2etsi, t2ets; //   edge-sharing nbrs excluding self
  IdxArray t2ts2i, t2ts2; // optional 2-halo excluding self
};

template <typename Int>
Idx Triangulation::append_tri (const Int tri[3]) {
  const auto n = t2vs.size();
  t2vs.resize(n+3);
  for (int i = 0; i < 3; ++i) {
    assert(tri[i] >= 0);
    t2vs[n+i] = Idx(tri[i]);
  }
  return get_ncell();
}

mesh::Mesh::Ptr
tri2mesh(const Triangulation& t, const TriangulationRelations::CPtr& tr);

} // namespace ctzx
} // namespace squirrel
} // namespace woodland

#endif
