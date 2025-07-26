#ifndef INCLUDE_WOODLAND_SQUIRREL_MESH
#define INCLUDE_WOODLAND_SQUIRREL_MESH

#include "woodland/squirrel/squirrel.hpp"

#include <cassert>
#include <memory>
#include <vector>

/* Mesh data structures:
     - Vertices: 2 or 3D set of vertex positions.
     - Cell: List of vertex indices corresponding to a cell.
     - Topology: List of cells.
     - Relations: Topological entity relations.
     - Mesh: Simple container for a (Vertices, Topology, Relations)
       triple. Multiple meshes can be created using the same vertices or
       topologies. The most common use case is to use the same Topology in
       multiple Meshes but different Vertices, e.g., due to multiple coordinate
       systems.
     The mesh may be conforming (no hanging nodes) or nonconforming. If it is
   conforming, all relations can be deduced. If it is not, the caller must
   provide the cell -> neighbors map and boundary edge data.
     In comments, indexing of the form a[i:j] uses python conventions; in
   particular, i:j corresonds to the set [i, i+1, ..., j-1]. If a(i:j) appears
   anywhere, the convention is F90/Matlab, and now i:j is inclusive on the right
   side.
     Many arrays come in pairs: xi and x. xi is the pointer array into x. Thus,
   x[xi[k]:xi[k+1]] is the k'th list.
 */

namespace woodland {
namespace squirrel {
namespace mesh {

typedef std::vector< Idx>  IdxArray;
typedef std::vector<char> BoolArray;
typedef std::vector<Real> RealArray;

struct Vertices {
  typedef std::shared_ptr<Vertices> Ptr;
  typedef std::shared_ptr<const Vertices> CPtr;

  // Optionally start with nvtx vertices.
  Vertices(const int ndim, const Idx nvtx = 0);

  int get_ndim () const { return ndim; }
    
  // Add a vertex. Returns the vertex's index. If no vertex vector is given,
  // space is made for the vertex, and the caller can later fill in the data
  // using get_vtx(i). Return the index of the vertex.
  Idx append_vtx(CRPtr vtx = nullptr);

  // Mesh geometry. Geometry can be modified at any time.
  Idx   get_nvtx () const { return vtxs.size()/ndim; }
  CRPtr get_vtx  (const Idx i) const { asserti(i, get_nvtx()); return &vtxs[ndim*i]; }
  RPtr  get_vtx  (const Idx i) { asserti(i, get_nvtx()); return &vtxs[ndim*i]; }
  CRPtr get_vtxs () const { return get_vtx(0); }

private:
  int ndim;
  RealArray vtxs;

  static void asserti(const Idx i, const Idx n) {
    assert(i >= 0);
    assert(i <  n);
  }
};

struct Cell {
  const int n;
  const Idx* const vtx; // vtx[0:n] describe the polygonal cell

  Idx operator[] (const int i) const { assert(i >= 0 and i < n); return vtx[i]; }

  Cell (const int n_, const Idx* const vtx_)
    : n(n_), vtx(vtx_)
  { assert(n >= 2); }
};

struct Topology {
  typedef std::shared_ptr<Topology> Ptr;
  typedef std::shared_ptr<const Topology> CPtr;

  Topology();

  // Begin modifications to the topology. Until end_modifications is called,
  // pointers to polygons and vertices returned from get_cell(s) are invalid. A
  // Mesh starts in the modification phase and can re-enter it multiple times.
  void begin_modifications();
  bool in_modifications () const { return in_mod_phase; }
  // Add a polygon. Returns the polygon's index.
  template <typename Int> Idx append_cell(const int n, const Int* vtxs);
  Idx append_cell(const Cell& cell);
  // End modifications to the topology. After this call, the get_cell(s)
  // pointers are valid.
  void end_modifications();

  Idx get_ncell () const { return static_cast<Idx>(csi.size()) - 1; }

  // This is not necessarily the number of unique vertices; rather, it's the
  // maximum vertex index + 1.
  Idx get_nvtx () const { return max_vtx+1; }

  // Mesh topology. Topology data can be incrementally added to during the
  // modifications phase, but already existing topology data cannot be modified
  // or removed.
  Cell get_cell (const Idx i) const {
    asserti(i, get_ncell());
    return Cell(csi[i+1]-csi[i], &cs[csi[i]]);
  }
  int get_cell_nvtx (const Idx i) const {
    asserti(i, get_ncell());
    return csi[i+1] - csi[i];
  }

private:
  bool in_mod_phase;
  Idx max_vtx;
  std::vector<Idx> csi, cs;

  void asserti(const Idx i, const Idx n) const {
    assert(i >= 0);
    assert(i <  n);
    assert(not in_mod_phase);
  }
};

struct Relations {
  typedef std::shared_ptr<Relations> Ptr;
  typedef std::shared_ptr<const Relations> CPtr;

  // Deduce the cell -> cell map, the cell -> edge-sharing-cell map, and the
  // boundary data. 'conforming' must be true.
  Relations(const Topology::CPtr& to,
            const bool conforming = true);

  // Provide the cell -> cell map and the boundary data. Nothing is deduced.
  Relations(const Topology::CPtr& to,
            const IdxArray& c2csi, const IdxArray& c2cs,
            const BoolArray& v2bdy, const bool check_consistency = true);

  const Topology::CPtr& get_topo () const { return to; }

  // xi and x are array pointer and data pairs as described above.
  // Names: v: vertex, e: edge, c: cell.

  // cell -> vertex-sharing nbrs excluding self
  const IdxArray& get_c2csi  () const { return c2csi ; }
  const IdxArray& get_c2cs   () const { return c2cs  ; }

  // v2bdy[i] is true if vertex i is on the boundary.
  const BoolArray& get_v2bdy () const { return v2bdy ; }

  // cell -> edge-sharing nbrs excluding self
  bool have_edge_sharing_lists() const;
  const IdxArray& get_c2ecsi () const { return c2ecsi; }
  const IdxArray& get_c2ecs  () const { return c2ecs ; }

private:
  Topology::CPtr to;
  IdxArray c2csi, c2cs, c2ecsi, c2ecs;
  BoolArray v2bdy;
};

// A Mesh holds a Vertices, a Topology, and optionally a Relations. Each of
// these can be used in multiple Meshes.
struct Mesh {
  typedef std::shared_ptr<Mesh> Ptr;
  typedef std::shared_ptr<const Mesh> CPtr;

  Mesh(const Vertices::CPtr& vtxs, const Topology::CPtr& topo,
       const Relations::CPtr& relations = nullptr,
       // If true and vtxs and topo are not consistent, raise an exception.
       const bool check_consistency = true);

  Idx get_ncell () const { return topo->get_ncell(); }

  const Vertices::CPtr& get_vtxs () const { return vtxs; }
  const Topology::CPtr& get_topo () const { return topo; }
  const Relations::CPtr& get_relations () const { return re; }

private:
  Vertices::CPtr vtxs;
  Topology::CPtr topo;
  Relations::CPtr re;
};

// Construct various relations.

void init_v2cs(const Topology& to,
               IdxArray& v2csi, IdxArray& v2cs);
void init_c2cs(const Topology& to, const IdxArray& v2csi, const IdxArray& v2cs,
               IdxArray& c2csi, IdxArray& c2cs);
void init_c2ecs(const Topology& to, const IdxArray& c2csi, const IdxArray& c2cs,
                IdxArray& c2ecsi, IdxArray& c2ecs);
// Deduce boundary edges as those having one adjacent cell rather than
// two. c2cs(i) can be c2ecs(i). c2bdy[c2bdyi[ci]: c2bdyi[ci+1]] are indices of
// edges in cell ci that are on the boundary.
void init_c2bdy(const Topology& to, const IdxArray& c2csi, const IdxArray& c2cs,
                IdxArray& c2bdyi, IdxArray& c2bdy);
// Deduce boundary vertices.
void init_v2bdy(const Topology& to, IdxArray& c2bdyi, IdxArray& c2bdy,
                BoolArray& v2bdy);

// For convenience, constructor wrappers to make Type::Ptr.

inline Vertices::Ptr make_vertices (const int ndim, const Idx nvtx = 0) {
  return std::make_shared<Vertices>(ndim, nvtx);
}

inline Topology::Ptr make_topology () { return std::make_shared<Topology>(); }

inline Relations::Ptr
make_relations (const Topology::CPtr& to, const bool conforming = true) {
  return std::make_shared<Relations>(to, conforming);
}

inline Relations::Ptr
make_relations (const Topology::CPtr& to,
                const IdxArray& c2csi, const IdxArray& c2cs,
                const BoolArray& v2bdy, const bool check_consistency = true) {
  return std::make_shared<Relations>(to, c2csi, c2cs, v2bdy, check_consistency);
}

inline Mesh::Ptr
make_mesh (const Vertices::CPtr& vtxs, const Topology::CPtr& topo,
           const Relations::CPtr& relations = nullptr,
           const bool check_consistency = true) {
  return std::make_shared<Mesh>(vtxs, topo, relations, check_consistency);
}

// Some convenience routines.

bool are_nbrs(const mesh::Relations& re, const Idx ci1, const Idx ci2);

int unittest();

template <typename Int>
Idx Topology::append_cell (const int n, const Int* vtxs_) {
  const int max_nvtx = 64;
  Idx vtxs[max_nvtx];
  assert(n > 2 and n <= max_nvtx);
  for (int i = 0; i < n; ++i) vtxs[i] = static_cast<Idx>(vtxs_[i]);
  return append_cell(Cell(n, vtxs));
}

} // namespace mesh
} // namespace squirrel
} // namespace woodland

#endif
