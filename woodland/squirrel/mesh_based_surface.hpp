#ifndef INCLUDE_WOODLAND_SQUIRREL_MESH_BASED_SURFACE
#define INCLUDE_WOODLAND_SQUIRREL_MESH_BASED_SURFACE

#include "woodland/squirrel/discretization.hpp"

namespace woodland {
namespace squirrel {

// Partial implementation of Surface for the case of a parameterization
// represented by a Mesh. The Mesh's Vertices object stores the global
// parametric surface vertices (u,v). Vertices::get_ndim() may be 3 as long as
// the first two dimensions provide (u,v).
struct MeshBasedSurface : public Surface {
  typedef mesh::RealArray RealArray;
  typedef mesh::IdxArray IdxArray;

  MeshBasedSurface(const mesh::Mesh::CPtr& m,
                   // Optionally scale the (u,v) polygon to make a roughly unit
                   // aspect ratio. This improves the interaction integral
                   // calculation. If used, then cell_position must apply the
                   // inverse to the input (u,v) values.
                   const bool scale = false);

  bool scale () const { return scale_; }
  bool scale (const Idx ci) const { return scale_ and sclidxs[ci] >= 0; }

  Idx  get_ncell   () const override;
  int  cell_nvtx   (const Idx ci) const override;
  // If scale(ci), then uv is scaled.
  void cell_vtxs_uv(const Idx ci, RPtr uv   ) const override;
  void cell_ctr_uv (const Idx ci, Real uv[2]) const override;
  
protected:
  mesh::Mesh::CPtr m;
  RealArray uv_ctrs;

  void apply_scale(const Idx ci, const bool fwd, Real p[2]) const;
  // Row-major fwd=true or col-major fwd=false.
  void get_jacobian_scale(const Idx ci, Real J[4]) const;
  Real get_jacdet_scale(const Idx ci) const;

private:
  bool scale_;
  IdxArray sclidxs;
  RealArray lambdas, Vts;

  void init(const mesh::Mesh::CPtr& m, const bool scale);
};

void init_xhat_from_primary(const Real zhat[3], const Real primary[3],
                            Real xhat[3]);

} // namespace squirrel
} // namespace woodland

#endif
