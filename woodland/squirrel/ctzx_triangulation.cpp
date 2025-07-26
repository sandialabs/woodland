#include "woodland/squirrel/ctzx_triangulation.hpp"
#include "woodland/acorn/matvec.hpp"

#include <set>

namespace woodland {
namespace squirrel {
namespace ctzx {

typedef acorn::Matvec<3,Real> mv3;

void Triangulation::begin_modifications () { in_mod_phase = true; }

Idx Triangulation::append_vtx (CRPtr vtx) {
  const auto n = vtxs.size();
  vtxs.resize(n+3);
  if (vtx) mv3::copy(vtx, &vtxs[n]);
  return get_nvtx();
}

void Triangulation::end_modifications () {
  in_mod_phase = false;
#ifndef NDEBUG
  const auto ntri = get_ncell();
  const auto nvtx = get_nvtx();
  for (Idx ti = 0; ti < ntri; ++ti) {
    const auto tri = get_tri(ti);
    for (int i = 0; i < 3; ++i) {
      assert(tri[i] >= 0);
      assert(tri[i] < nvtx);
    }
  }
#endif
}

int Triangulation::unittest () {
  return 0;
}

typedef TriangulationRelations::IdxArray IdxArray;

static void init_v2ts (const Triangulation& t, IdxArray& v2tsi, IdxArray& v2ts) {
  const auto nvtx = t.get_nvtx();
  const auto ntri = t.get_ncell();

  std::vector<std::set<Idx>> d(nvtx);
  for (Idx ti = 0; ti < ntri; ++ti) {
    const auto tri = t.get_tri(ti);
    for (int i = 0; i < 3; ++i)
      d[tri[i]].insert(ti);
  }

  v2tsi.resize(nvtx+1);
  v2tsi[0] = 0;
  for (Idx vi = 0; vi < nvtx; ++vi) {
    v2tsi[vi+1] = v2tsi[vi] + d[vi].size();
    for (const auto& e : d[vi])
      v2ts.push_back(e);
  }
}

static void
init_t2ts (const Triangulation& t, const IdxArray& v2tsi, const IdxArray& v2ts,
           IdxArray& t2tsi, IdxArray& t2ts, const int halo = 1) {
  assert(halo == 1 || halo == 2);
  const auto ntri = t.get_ncell();
  t2tsi.resize(ntri+1);
  t2tsi[0] = 0;
  std::set<Idx> d;
  for (Idx ti = 0; ti < ntri; ++ti) {
    d.clear();
    const auto tri = t.get_tri(ti);
    for (int i = 0; i < 3; ++i) {
      const auto vi = tri[i];
      for (Idx iidx = v2tsi[vi]; iidx < v2tsi[vi+1]; ++iidx) {
        if (v2ts[iidx] == ti) continue;
        d.insert(v2ts[iidx]);
        if (halo > 1) {
          const auto trj = t.get_tri(v2ts[iidx]);
          for (int j = 0; j < 3; ++j) {
            const auto vj = trj[j];
            for (Idx jidx = v2tsi[vj]; jidx < v2tsi[vj+1]; ++jidx) {
              if ((v2ts[jidx] == ti)) continue;
              d.insert(v2ts[jidx]);
            }
          }
        }
      }
    }
    t2tsi[ti+1] = t2tsi[ti] + d.size();
    for (const auto& e : d)
      t2ts.push_back(e);
  }
}

static void
init_t2ets (const Triangulation& t, const IdxArray& t2tsi, const IdxArray& t2ts,
            IdxArray& t2etsi, IdxArray& t2ets_) {
  const auto ntri = t.get_ncell();
  t2etsi.resize(ntri+1);
  t2etsi[0] = 0;
  std::vector<Idx> t2ets;
  for (Idx ti = 0; ti < ntri; ++ti) {
    t2etsi[ti+1] = t2etsi[ti];
    const auto tri = t.get_tri(ti);
    for (Idx j = t2tsi[ti]; j < t2tsi[ti+1]; ++j) {
      const auto tj = t2ts[j];
      const auto trj = t.get_tri(tj);
      int n = 0;
      for (int vi = 0; vi < 3; ++vi)
        for (int vj = 0; vj < 3; ++vj)
          if (tri[vi] == trj[vj]) ++n;
      assert(n <= 2);
      if (n == 2) {
        t2ets.push_back(tj);
        ++t2etsi[ti+1];
      }
    }
    assert(t2etsi[ti+1] - t2etsi[ti] >= 0 &&
           t2etsi[ti+1] - t2etsi[ti] <= 3);
  }
  t2ets_.resize(t2ets.size());
  std::copy(t2ets.begin(), t2ets.end(), t2ets_.begin());
}

TriangulationRelations::TriangulationRelations (const Triangulation::CPtr& t_) {
  t = t_;
  IdxArray v2tsi, v2ts;
  init_v2ts(*t, v2tsi, v2ts);
  init_t2ts(*t, v2tsi, v2ts, t2tsi, t2ts);
  init_t2ets(*t, t2tsi, t2ts, t2etsi, t2ets);
}

void TriangulationRelations::make_t2ts2 () {
  if (not t2ts2.empty()) return;
  IdxArray v2tsi, v2ts;
  init_v2ts(*t, v2tsi, v2ts);
  init_t2ts(*t, v2tsi, v2ts, t2ts2i, t2ts2, 2);
}

static bool same (const IdxArray& a, const IdxArray& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i] != b[i])
      return false;
  return true;
}

Triangulation::Ptr mesh2tri (const mesh::Mesh& mesh) {
  const auto v = mesh.get_vtxs();
  const auto to = mesh.get_topo();
  auto t = std::make_shared<Triangulation>();
  const auto nvtx = v->get_nvtx();
  const auto ndim = v->get_ndim();
  for (Idx i = 0; i < nvtx; ++i) {
    const auto* src = v->get_vtx(i);
    Real dst[3];
    dst[2] = 0;
    for (int d = 0; d < ndim; ++d)
      dst[d] = src[d];
    t->append_vtx(dst);
  }
  const auto nc = to->get_ncell();
  for (Idx i = 0; i < nc; ++i) {
    const auto& c = to->get_cell(i);
    throw_if(c.n != 3, "mesh2tri: topo must have only triangles");
    t->append_tri(c.vtx);
  }
  t->end_modifications();
  return t;
}

static void check (const Triangulation& t1, const TriangulationRelations& tr,
                   const mesh::Mesh& m) {
  const auto v = m.get_vtxs();
  const auto re = m.get_relations();
  const auto to = m.get_topo();
  const auto t2 = mesh2tri(m);
  throw_if_nomsg(v->get_ndim() != 3);
  throw_if_nomsg(v->get_nvtx() != t1.get_nvtx());
  throw_if_nomsg(to->get_ncell() != t1.get_ncell());
  {
    bool ok = true;
    const auto n = t1.get_nvtx();
    const auto v1 = t1.get_vtxs();
    const auto v2 = v->get_vtxs();
    for (Idx i = 0; i < n; ++i)
      if (v1[i] != v2[i]) ok = false;
    throw_if_nomsg(not ok);
  }
  {
    bool ok = true;
    const auto n = t1.get_ncell();
    for (Idx i = 0; i < n; ++i) {
      const auto c1 = t1.get_tri(i);
      const auto c2 = to->get_cell(i);
      if (c2.n != 3) ok = false;
      for (int j = 0; j < 3; ++j)
        if (c1[j] != c2[j]) ok = false;
    }
    throw_if_nomsg(not ok);
  }
  {
    throw_if_nomsg(not same(re->get_c2csi(), tr.get_t2tsi()));
    throw_if_nomsg(not same(re->get_c2cs(), tr.get_t2ts()));
    throw_if_nomsg(not same(re->get_c2ecsi(), tr.get_t2etsi()));
    throw_if_nomsg(not same(re->get_c2ecs(), tr.get_t2ets()));
    // Test other Relations ctor.
    make_relations(to, re->get_c2csi(), re->get_c2cs(),
                   re->get_v2bdy());
  }  
}

mesh::Mesh::Ptr
tri2mesh (const Triangulation& t, const TriangulationRelations::CPtr& tr) {
  const auto nvtx = t.get_nvtx();
  const auto v = mesh::make_vertices(3, nvtx);
  for (Idx i = 0; i < nvtx; ++i)
    mv3::copy(t.get_vtx(i), v->get_vtx(i));
  const auto to = mesh::make_topology();
  const auto ntri = t.get_ncell();
  for (Idx i = 0; i < ntri; ++i)
    to->append_cell(3, t.get_tri(i));
  to->end_modifications();
  const auto re = mesh::make_relations(to);
  const auto m = mesh::make_mesh(v, to, re);
  if (tr) check(t, *tr, *m);
  return m;
}

} // namespace ctzx
} // namespace squirrel
} // namespace woodland
