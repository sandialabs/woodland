#include "woodland/squirrel/global_z_surface.hpp"

#include "woodland/acorn/util.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace squirrel {

typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;
using mesh::Mesh;

static void asserti (const mesh::Mesh& m, const Idx ci) {
  assert(ci >= 0);
  assert(ci < m.get_ncell());
}

GlobalZSurface
::GlobalZSurface (const Mesh::CPtr& m, const UserFn::CPtr& ufn_,
                  const bool support_req0_, const bool mesh_scale)
  : MeshBasedSurface(m, mesh_scale), ufn(ufn_), support_req0(support_req0_)
{ init(); }

void GlobalZSurface::init () {
  const auto& to = *m->get_topo();
  const auto nc = to.get_ncell();
  z_ctrs.resize(nc);
  lcs_ctrs.resize(9*nc);
  ompparfor for (Idx ci = 0; ci < nc; ++ci) {
    Real xyz[3];
    const auto lcs = &lcs_ctrs[9*ci];
    lcs[0] = 1; lcs[1] = 0; lcs[2] = 0;
    cell_position1(ci, &uv_ctrs[2*ci], xyz, lcs);
    z_ctrs[ci] = xyz[2];
    assert(not std::isnan(z_ctrs[ci]));
  }
}

void GlobalZSurface::cell_ctr_xyz (const Idx ci, Real xyz[3]) const {
  asserti(*m, ci);
  xyz[0] = uv_ctrs[2*ci  ];
  xyz[1] = uv_ctrs[2*ci+1];
  xyz[2] = z_ctrs[ci];
}

void GlobalZSurface::cell_ctr_lcs (const Idx ci, Real lcs[9]) const {
  asserti(*m, ci);
  acorn::copy(9, &lcs_ctrs[9*ci], lcs);
}

bool GlobalZSurface
::calc_cell_lcs (const Idx ci, const Real xyz[3], Real uv[2]) const {
  mv2::copy(xyz, uv);
  apply_scale(ci, true, uv);
  return true;
}

bool GlobalZSurface::supports_J () const { return support_req0; }

void GlobalZSurface
::cell_position (const Idx ci, const int n, CRPtr uv, RPtr p_gcs_, RPtr lcs_,
                 RPtr jacdet, RPtr jac) const {
  const bool scale = this->scale(ci);
  for (int i = 0; i < n; ++i) {
    Real p_gcs[3];
    mv2::copy(&uv[2*i], p_gcs);
    if (scale) apply_scale(ci, false, p_gcs);
    ufn->eval(p_gcs[0], p_gcs[1], p_gcs[2], nullptr);
    mv3::copy(p_gcs, &p_gcs_[3*i]);
    Real grad[2];
    if (jacdet or jac or lcs_) {
      Real f;
      ufn->eval(p_gcs[0], p_gcs[1], f, grad);
    }
    // J = [1 0; 0 1; z_u z_v]
    // sqrt(det(J'J))
    //   = sqrt((1 + z_u^2) (1 + z_v^2) - (z_u z_v)^2)
    //   = sqrt(1 + z_u^2 + z_v^2)
    if (jacdet) {
      jacdet[i] = std::sqrt(1 + mv2::norm22(grad));
      if (scale) jacdet[i] *= get_jacdet_scale(ci);
    }
    if (jac) {
      auto J = &jac[6*i];
      for (int k = 0; k < 6; ++k) J[k] = 0;
      if (scale) {
        // J times (Jacobian of map from scaled to unscaled)
        Real A[4];
        get_jacobian_scale(ci, A);
        J[0] = A[0];
        J[1] = A[2];
        J[2] = A[1];
        J[3] = A[3];
        J[4] = A[0]*grad[0] + A[1]*grad[1];
        J[5] = A[2]*grad[0] + A[3]*grad[1];
      } else {
        J[0] = 1; J[3] = 1; J[4] = grad[0]; J[5] = grad[1];
      }
    }
    if (lcs_) {
      Real zhat[] = {-grad[0], -grad[1], 1};
      mv3::normalize(zhat);
      Real xhat[3], yhat[3];
      init_xhat_from_primary(zhat, &lcs_ctrs[9*ci], xhat);
      mv3::cross(zhat, xhat, yhat);
      mv3::normalize(yhat);
      const auto lcs = &lcs_[9*i];
      mv3::copy(xhat,  lcs   );
      mv3::copy(yhat, &lcs[3]);
      mv3::copy(zhat, &lcs[6]);
    }
  }
}

} // namespace squirrel
} // namespace woodland
