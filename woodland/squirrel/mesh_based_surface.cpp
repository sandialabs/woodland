#include "woodland/squirrel/mesh_based_surface.hpp"

#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/plane_geometry.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace squirrel {

static const int max_nvtx = 8;

typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;

static void asserti (const mesh::Mesh& m, const Idx ci) {
  assert(ci >= 0);
  assert(ci < m.get_ncell());
}

void init_xhat_from_primary (const Real zhat[3], const Real primary[3],
                             Real xhat[3]) {
  const auto alpha = mv3::dot(zhat, primary);
  mv3::axpbyz(1, primary, -alpha, zhat, xhat);
  mv3::normalize(xhat);
}

MeshBasedSurface
::MeshBasedSurface (const mesh::Mesh::CPtr& m, const bool scale)
{ init(m, scale); }

static void svd (CRPtr ps, const int n, Real sig[2], Real Ut[4]) {
  Real coef[4];
  for (int r = 0; r < 2; ++r)
    for (int c = 0; c < 2; ++c) {
      if (r == 1 and c == 0) continue;
      Real acc = 0;
      for (int i = 0; i < n; ++i)
        acc += ps[2*i+r]*ps[2*i+c];
      coef[2*r+c] = acc;
    }
  acorn::sym_2x2_eig(coef[0], coef[3], coef[1], sig, Ut, Ut+2);
  assert(sig[0] > 0);
}

static void init_scale (const Real ctr[2], const acorn::plane::Polygon& p,
                        Real lambda[2], Real Vt[4]) {
  Real ps[2*max_nvtx];
  for (int i = 0; i < p.n; ++i)
    for (int d = 0; d < 2; ++d)
      ps[2*i+d] = p.xys[2*i+d] - ctr[d];
  svd(ps, p.n, lambda, Vt);
  lambda[1] = std::sqrt(lambda[1]/lambda[0]);
  lambda[0] = 1;
}

void MeshBasedSurface::init (const mesh::Mesh::CPtr& m_, const bool scale_) {
  m = m_;
  this->scale_ = scale_;
  const auto& v = *m->get_vtxs();
  const auto& to = *m->get_topo();
  const auto nc = to.get_ncell();
  uv_ctrs.resize(2*nc);
  RealArray tlams, tVts;
  if (scale_) {
    tlams.resize(2*nc);
    tVts.resize(4*nc);
  }
  ompparfor for (Idx ci = 0; ci < nc; ++ci) {
    const auto c = to.get_cell(ci);
    Real vtxs[2*max_nvtx];
    assert(c.n <= max_nvtx);
    for (int i = 0; i < c.n; ++i)
      mv2::copy(v.get_vtx(c[i]), &vtxs[2*i]);
    acorn::plane::Polygon p(vtxs, c.n);
    acorn::plane::calc_centroid(p, &uv_ctrs[2*ci]);
    assert(not std::isnan(uv_ctrs[2*ci]));
    if (scale_) init_scale(&uv_ctrs[2*ci], p, &tlams[2*ci], &tVts[4*ci]);
  }
  if (scale_) {
    sclidxs.resize(nc);
    Idx slot = 0;
    for (Idx ci = 0; ci < nc; ++ci) {
      sclidxs[ci] = -1;
      const auto lam = &tlams[2*ci];
      if (lam[1] < 2*lam[0]) continue;
      sclidxs[ci] = slot++;
    }
    lambdas.resize(2*slot);
    Vts.resize(4*slot);
    ompparfor for (Idx ci = 0; ci < nc; ++ci) {
      const auto slot = sclidxs[ci];
      if (slot < 0) continue;
      acorn::copy(2, &tlams[2*ci], &lambdas[2*slot]);
      acorn::copy(4, &tVts[4*ci], &Vts[4*slot]);
    }
  }
}

Idx MeshBasedSurface::get_ncell () const { return m->get_ncell(); }

int MeshBasedSurface::cell_nvtx (const Idx ci) const {
  return m->get_topo()->get_cell_nvtx(ci);
}

void MeshBasedSurface
::apply_scale (const Idx ci, const bool fwd, Real p[2]) const {
  if (not scale(ci)) return;
  const auto slot = sclidxs[ci];
  assert(slot >= 0);
  const auto Vt = &Vts[4*slot];
  const auto lam = &lambdas[2*slot];
  Real tmp[2];
  for (int d = 0; d < 2; ++d)
    tmp[d] = p[d] - uv_ctrs[2*ci + d];
  if (fwd) {
    mv2::matvec(Vt, tmp, p);
    p[1] /= lam[1]; // lam[0] is 1
    mv2::tmatvec(Vt, p, tmp);
  } else {
    mv2::tmatvec(Vt, tmp, p); 
    p[1] *= lam[1];
    mv2::matvec(Vt, p, tmp);
  }
  for (int d = 0; d < 2; ++d)
    p[d] = tmp[d] + uv_ctrs[2*ci + d];
}

void MeshBasedSurface::get_jacobian_scale (const Idx ci, Real J[4]) const {
  J[0] = 1; J[1] = 0; J[2] = 0; J[3] = 1;
  if (not scale(ci)) return;
  const auto slot = sclidxs[ci];
  const auto Vt = &Vts[4*slot];
  const auto lam = &lambdas[2*slot];
  for (int i = 0; i < 2; ++i) {
    Real tmp[2];
    Real* p = &J[2*i];
    mv2::tmatvec(Vt, p, tmp);
    tmp[1] *= lam[1];
    mv2::matvec(Vt, tmp, p);
  }
}

Real MeshBasedSurface::get_jacdet_scale (const Idx ci) const {
  if (not scale(ci)) return 1;
  const auto slot = sclidxs[ci];
  return lambdas[2*slot]*lambdas[2*slot+1];
}

void MeshBasedSurface::cell_vtxs_uv (const Idx ci, RPtr uv) const {
  const auto& v = *m->get_vtxs();
  const auto& to = *m->get_topo();
  const auto c = to.get_cell(ci);
  for (int i = 0; i < c.n; ++i) {
    const auto p = v.get_vtx(c[i]);
    for (int d = 0; d < 2; ++d)
      uv[2*i + d] = p[d];
  }
  if (scale(ci))
    for (int i = 0; i < c.n; ++i)
      apply_scale(ci, true, &uv[2*i]);
}

void MeshBasedSurface::cell_ctr_uv (const Idx ci, Real p[2]) const {
  asserti(*m, ci);
  mv2::copy(&uv_ctrs[2*ci], p);
}

} // namespace squirrel
} // namespace woodland
