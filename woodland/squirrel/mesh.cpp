#include "woodland/squirrel/mesh.hpp"

#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/dbg.hpp"

#include <set>

namespace woodland {
namespace squirrel {
namespace mesh {

typedef acorn::Matvec<3,Real> mv3;

Vertices::Vertices (const int ndim_, const Idx nvtx) {
  ndim = ndim_;
  assert(ndim >= 2 and ndim <= 3);
  if (nvtx > 0) vtxs.resize(ndim*nvtx);
}

Idx Vertices::append_vtx (CRPtr vtx) {
  const auto n = vtxs.size();
  vtxs.resize(n + ndim);
  if (vtx)
    for (int i = 0; i < ndim; ++i)
      vtxs[n+i] = vtx[i];
  return get_nvtx() - 1;
}

Topology::Topology () {
  csi.push_back(0);
  begin_modifications();
}

void Topology::begin_modifications () { in_mod_phase = true; }

Idx Topology::append_cell (const Cell& p) {
  assert(in_mod_phase);

  const auto n = csi.back();
  cs.resize(n + p.n);
  for (int i = 0; i < p.n; ++i) {
    assert(p[i] >= 0);
    cs[n+i] = p[i];
  }
  csi.push_back(n + p.n);
  
  return get_ncell();
}

void Topology::end_modifications () {
  in_mod_phase = false;
  max_vtx = 0;
  for (const auto& e : cs) max_vtx = std::max(max_vtx, e);
#ifndef NDEBUG
  const auto np = get_ncell();
  for (Idx ci = 0; ci < np; ++ci) {
    const auto& c = get_cell(ci);
    for (int i = 0; i < c.n; ++i)
      assert(c[i] >= 0);
  }
#endif
}

void init_v2cs (const Topology& to, IdxArray& v2csi, IdxArray& v2cs) {
  const auto nvtx = to.get_nvtx();
  const auto nc = to.get_ncell();
  std::vector<std::set<Idx>> d(nvtx);
  for (Idx ci = 0; ci < nc; ++ci) {
    const auto c = to.get_cell(ci);
    for (int i = 0; i < c.n; ++i)
      d[c[i]].insert(ci);
  }
  v2cs.clear();
  v2csi.resize(nvtx+1);
  v2csi[0] = 0;
  for (Idx vi = 0; vi < nvtx; ++vi) {
    v2csi[vi+1] = v2csi[vi] + d[vi].size();
    for (const auto& e : d[vi])
      v2cs.push_back(e);
  }
}

void init_c2cs (const Topology& to, const IdxArray& v2csi, const IdxArray& v2cs,
                IdxArray& c2csi, IdxArray& c2cs) {
  const auto nc = to.get_ncell();
  c2cs.clear();
  c2csi.resize(nc+1);
  c2csi[0] = 0;
  std::set<Idx> d;
  for (Idx ci = 0; ci < nc; ++ci) {
    d.clear();
    const auto c = to.get_cell(ci);
    for (int i = 0; i < c.n; ++i) {
      const auto vi = c[i];
      for (Idx iidx = v2csi[vi]; iidx < v2csi[vi+1]; ++iidx) {
        if (v2cs[iidx] == ci) continue;
        d.insert(v2cs[iidx]);
      }
    }
    c2csi[ci+1] = c2csi[ci] + d.size();
    for (const auto& e : d)
      c2cs.push_back(e);
  }
}

// Assumes the mesh is conforming.
void init_c2ecs (const Topology& to, const IdxArray& c2csi, const IdxArray& c2cs,
                 IdxArray& c2ecsi, IdxArray& c2ecs_) {
  const auto nc = to.get_ncell();
  c2ecsi.resize(nc+1);
  c2ecsi[0] = 0;
  std::vector<Idx> c2ecs;
  for (Idx ci = 0; ci < nc; ++ci) {
    c2ecsi[ci+1] = c2ecsi[ci];
    const auto cei = to.get_cell(ci);
    for (Idx j = c2csi[ci]; j < c2csi[ci+1]; ++j) {
      const auto cj = c2cs[j];
      const auto cej = to.get_cell(cj);
      int n = 0;
      for (int vi = 0; vi < cei.n; ++vi)
        for (int vj = 0; vj < cej.n; ++vj)
          if (cei[vi] == cej[vj]) ++n;
      assert(n <= 2);
      if (n == 2) {
        c2ecs.push_back(cj);
        ++c2ecsi[ci+1];
      }
    }
    assert(c2ecsi[ci+1] - c2ecsi[ci] >= 0 &&
           c2ecsi[ci+1] - c2ecsi[ci] <= 3);
  }
  c2ecs_.resize(c2ecs.size());
  std::copy(c2ecs.begin(), c2ecs.end(), c2ecs_.begin());
}

void init_c2bdy (const Topology& to, const IdxArray& c2csi, const IdxArray& c2cs,
                 IdxArray& c2bdyi, IdxArray& c2bdy) {
  const auto nc = to.get_ncell();
  assert(c2csi.size() == size_t(nc) + 1);
  c2bdyi.clear(); c2bdy.clear();
  c2bdyi.push_back(0);
  std::vector<char> mark;
  for (Idx ci = 0; ci < nc; ++ci) {
    const auto cei = to.get_cell(ci);
    if (mark.size() < size_t(cei.n)) mark.resize(cei.n);
    for (int ei = 0; ei < cei.n; ++ei) mark[ei] = 0;
    for (int j = c2csi[ci]; j < c2csi[ci+1]; ++j) {
      const auto cj = c2cs[j];
      const auto cej = to.get_cell(cj);
      for (int ei = 0; ei < cei.n; ++ei) {
        const Idx edi[] = {cei[ei], cei[(ei+1)%cei.n]};
        for (int ej = 0; ej < cej.n; ++ej) {
          const Idx edj[] = {cej[ej], cej[(ej+1)%cej.n]};
          if ((edj[0] == edi[1] and edj[1] == edi[0]) or
              // We don't require consistently oriented cells.
              (edj[0] == edi[0] and edj[1] == edi[1]))
            ++mark[ei];
        }
      }
    }
    c2bdyi.push_back(c2bdyi.back());
    for (int ei = 0; ei < cei.n; ++ei) {
      assert(mark[ei] >= 0 and mark[ei] <= 1);
      if (mark[ei] == 0) {
        ++c2bdyi.back();
        c2bdy.push_back(ei);
      }
    }
  }
}

void init_v2bdy (const Topology& to, IdxArray& bi, IdxArray& b,
                 BoolArray& v2bdy) {
  const auto nc = to.get_ncell();
  v2bdy.clear();
  v2bdy.resize(to.get_nvtx(), 0);
  for (Idx ci = 0; ci < nc; ++ci) {
    const auto cei = to.get_cell(ci);
    for (int ei = bi[ci]; ei < bi[ci+1]; ++ei) {
      const auto e = b[ei];
      v2bdy[cei[e]] = 1;
      v2bdy[cei[(e+1)%cei.n]] = 1;
    }
  }
}

static void check_c2cs (const IdxArray& c2csi, const IdxArray& c2cs) {
  const Idx ncell = static_cast<Idx>(c2csi.size()) - 1;
  for (const auto& e : c2cs)
    throw_if(e < 0 or e >= ncell,
             "Relations: A c2cs entry is out of bounds.");
  ompparfor for (Idx ci = 0; ci < ncell; ++ci) {
    for (Idx j = c2csi[ci]; j < c2csi[ci+1]; ++j) {
      const auto cj = c2cs[j];
      bool fnd = false;
      for (Idx k = c2csi[cj]; k < c2csi[cj+1]; ++k)
        if (c2cs[k] == ci) {
          fnd = true;
          break;
        }
      throw_if(not fnd, "Relations:: Cell " << ci << " has cell " << cj
               << " as a neighbor, but the opposite does not hold.");
    }
  }
}

Relations::Relations (const Topology::CPtr& to_, const bool conforming)
  : to(to_)
{
  throw_if(to->in_modifications(),
           "Relations: Topology must not be in modifications.");
  throw_if(not conforming,
           "Relations: This ctor is not supported for nonconforming meshes.");
  IdxArray v2csi, v2cs;
  init_v2cs(*to, v2csi, v2cs);
  init_c2cs(*to, v2csi, v2cs, c2csi, c2cs);
#ifndef NDEBUG
  check_c2cs(c2csi, c2cs);
#endif
  init_c2ecs(*to, c2csi, c2cs, c2ecsi, c2ecs);
#ifndef NDEBUG
  check_c2cs(c2ecsi, c2ecs);
#endif
  IdxArray c2bdyi, c2bdy;
  init_c2bdy(*to, c2ecsi, c2ecs, c2bdyi, c2bdy);
  init_v2bdy(*to, c2bdyi, c2bdy, v2bdy);
}

Relations::Relations (const Topology::CPtr& to_,
                      const IdxArray& c2csi_, const IdxArray& c2cs_,
                      const BoolArray& v2bdy_, const bool check_consistency)
  : to(to_), c2csi(c2csi_), c2cs(c2cs_), v2bdy(v2bdy_)
{
  throw_if(to->in_modifications(),
           "Relations: Topology must not be in modifications.");
  throw_if(c2csi.size() != static_cast<size_t>(to->get_ncell()) + 1,
           "Relations: Topology and c2csi size are not consistent.");
  throw_if(v2bdy.size() != static_cast<size_t>(to->get_nvtx()),
           "Relations: Topology and v2bdy size are not consistent.");
  if (check_consistency) check_c2cs(c2csi, c2cs);
}

bool Relations::have_edge_sharing_lists () const { return not c2ecsi.empty(); }

Mesh::Mesh (const Vertices::CPtr& vtxs_, const Topology::CPtr& topo_,
            const Relations::CPtr& re_, const bool check_consistency)
  : vtxs(vtxs_), topo(topo_), re(re_)
{
  if (check_consistency) {
    bool ok = true;
    const auto nvtx = vtxs->get_nvtx();
    const auto nc = topo->get_ncell();
    for (Idx ci = 0; ci < nc; ++ci) {
      const auto& c = topo->get_cell(ci);
      for (int i = 0; i < c.n; ++i)
        if (not (c[i] >= 0 and c[i] < nvtx))
          ok = false;
    }
    throw_if(not ok, "Mesh: vtxs and topo are not consistent.");
    throw_if(re and re->get_topo().get() != topo.get(),
             "Mesh: Topology and Relations are not consistent.");
  }
}

static void fill (const Idx* cells, const int n, Topology& to) {
  int i0 = 0;
  while (i0 < n) {
    int i1 = i0;
    while (cells[i1] != -1) ++i1;
    to.append_cell(i1 - i0, &cells[i0]);
    i0 = i1 + 1;
  }
  to.end_modifications();
}

static void fill_arrays (const Idx* a, const int n, IdxArray& xi, IdxArray& x) {
  int i0 = 0;
  xi.push_back(0);
  while (i0 < n) {
    int i1 = i0;
    while (a[i1] != -1) ++i1;
    xi.push_back(xi.back() + (i1 - i0));
    for (int i = i0; i < i1; ++i) x.push_back(a[i]);
    i0 = i1 + 1;
  }
}

static int
check (const Topology& to, const Idx* bdys, const int n,
       const IdxArray& c2bdyi, const IdxArray& c2bdy) {
  IdxArray edgesi, edges;
  fill_arrays(bdys, n, edgesi, edges);
  const auto nc = to.get_ncell();
  if (edgesi.size() != size_t(nc) + 1) return 1;
  if (c2bdyi.size() != size_t(nc) + 1) return 1;
  int nerr = 0;
  for (Idx ci = 0; ci < nc; ++ci) {
    const auto c = to.get_cell(ci);
    for (int i = edgesi[ci]; i < edgesi[ci+1]; i += 2) {
      const Idx ei[] = {edges[i], edges[i+1]};
      bool fnd = false;
      for (int j = c2bdyi[ci]; j < c2bdyi[ci+1]; ++j) {
        const Idx ej[] = {c[c2bdy[j]], c[(c2bdy[j]+1)%c.n]};
        if (ei[0] == ej[0] and ei[1] == ej[1]) fnd = true;
      }
      if (not fnd) ++nerr;
    }
    for (int j = c2bdyi[ci]; j < c2bdyi[ci+1]; ++j) {
      const Idx ej[] = {c[c2bdy[j]], c[(c2bdy[j]+1)%c.n]};
      bool fnd = false;
      for (int i = edgesi[ci]; i < edgesi[ci+1]; i += 2) {
        const Idx ei[] = {edges[i], edges[i+1]};
        if (ei[0] == ej[0] and ei[1] == ej[1]) fnd = true;
      }
      if (not fnd) ++nerr;
    }
  }
  return nerr;
}

static int test_Relations () {
  int nerr = 0;
  { // conforming
    const Idx cells[] = {0, 1, 2, 3, -1,
                         1, 4, 2, -1,
                         2, 5, 3, -1,
                         2, 4, 5, -1,
                         3, 5, 6, -1};
    const Idx bdys[] = {0, 1, 3, 0, -1, // edges (i,j)
                        1, 4, -1,
                        -1,
                        4, 5, -1,
                        5, 6, 6, 3, -1};
    auto to = make_topology();
    fill(cells, sizeof(cells)/sizeof(*cells), *to);
    auto re = make_relations(to);
    IdxArray c2bdyi, c2bdy;
    init_c2bdy(*to, re->get_c2csi(), re->get_c2cs(), c2bdyi, c2bdy);
    nerr += check(*to, bdys, sizeof(bdys)/sizeof(*bdys), c2bdyi, c2bdy);
    const auto& v2bdy = re->get_v2bdy();
    // Only vtx 2 isn't on the boundary.
    if (v2bdy[2]) ++nerr;
    for (size_t i = 0; i < v2bdy.size(); ++i)
      if (i != 2 and not v2bdy[i]) ++nerr;
    try {
      const auto re1 = make_relations(to, re->get_c2csi(), re->get_c2cs(),
                                      re->get_v2bdy());
    } catch (...) { ++nerr; }
  }
  return nerr;
}

int unittest () {
  int nerr = 0, ne;
  rununittest(test_Relations);
  return nerr;
}

bool are_nbrs (const mesh::Relations& re, const Idx ci1, const Idx ci2) {
  const auto& c2csi = re.get_c2csi();
  const auto& c2cs = re.get_c2cs();
  for (int i = c2csi[ci1]; i < c2csi[ci1+1]; ++i)
    if (c2cs[i] == ci2)
      return true;
  return false;
}

} // namespace mesh
} // namespace squirrel
} // namespace woodland
