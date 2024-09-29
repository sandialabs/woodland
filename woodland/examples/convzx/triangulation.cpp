#include "woodland/examples/convzx/triangulation.hpp"
#include "woodland/acorn/matvec.hpp"

namespace woodland {
namespace examples {
namespace convzx {

typedef acorn::Matvec<3,Real> mv3;

void Triangulation::begin_modifications () { in_mod_phase = true; }

Size Triangulation::append_vtx (CRPtr vtx) {
  const auto n = vtxs.size();
  vtxs.resize(n+3);
  if (vtx) mv3::copy(vtx, &vtxs[n]);
  return get_nvtx();
}

Size Triangulation::append_tri (const int tri[3]) {
  const int n = t2vs.size();
  t2vs.resize(n+3);
  for (size_t i = 0; i < 3; ++i) {
    assert(tri[i] >= 0);
    t2vs[n+i] = Idx(tri[i]);
  }
  return get_ntri();
}

Size Triangulation::append_tri (const Idx tri[3]) {
  const int n = t2vs.size();
  t2vs.resize(n+3);
  for (size_t i = 0; i < 3; ++i) t2vs[n+i] = tri[i];
  return get_ntri();
}

void Triangulation::end_modifications () {
  in_mod_phase = false;
#ifndef NDEBUG
  const auto ntri = get_ntri();
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

} // namespace convzx
} // namespace examples
} // namespace woodland
