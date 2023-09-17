#include "woodland/acorn/vv.hpp"
#include "woodland/acorn/macros.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/caller_integrand.hpp"
#include "woodland/acorn/plane_geometry.hpp"
#include "woodland/acorn/interaction_integrals.hpp"
#include "woodland/acorn/elastostatics_integrals.hpp"
#include "woodland/acorn/compose_quadrature.hpp"
#include "woodland/acorn/compose_lagrange_polynomial.hpp"
#include "woodland/acorn/dbg.hpp"

namespace woodland {
namespace acorn {
namespace vv {

using plane::Pt;
typedef Matvec<3,Real> mv3;
typedef Matvec<2,Real> mv2;

GeneralizedCylinder::GeneralizedCylinder ()
  : shape(circle), length(1), radius(1), rotation(0), ell_a(2), ell_b(1)
{}

GeneralizedCylinder& GeneralizedCylinder
::set_ellipse (const Real a, const Real b) {
  ell_a = a; ell_b = b;
  return *this;
}

GeneralizedCylinder::Shape GeneralizedCylinder::convert (const int i) {
  auto shape = static_cast<Shape>(i);
  if ( ! is_valid(shape)) shape = invalid;
  return shape;
}

GeneralizedCylinder::Shape GeneralizedCylinder
::convert (const std::string& shape) {
  if (shape == "circle") return circle;
  if (shape == "ellipse") return ellipse;
  return invalid;
}

std::string GeneralizedCylinder::convert (const Shape shape) {
  switch (shape) {
  case circle: return "circle";
  case ellipse: return "ellipse";
  case invalid:
  default: return "invalid";
  }
}

bool GeneralizedCylinder::is_valid (const Shape shape) {
  return shape >= circle && shape < invalid;
}

void GeneralizedCylinder::get_position (const Real u[2], Real x[3]) const {
  x[2] = length*(u[0] - 0.5);
  const auto c = std::cos(u[1]), s = std::sin(u[1]);
  switch (shape) {
  case circle: {
    x[0] = radius*c; x[1] = radius*s;
  } break;
  case ellipse: {
    x[0] = ell_a*c; x[1] = ell_b*s;
  } break;
  default: assert(0); break;
  }
  if (rotation != 0) {
    const Real xhat[] = {1, 0, 0};
    rotate_vector_3(xhat, rotation, x);
  }
}

void GeneralizedCylinder::get_uv (const Real xyz[3], Real uv[2]) const {
  Real p[3];
  mv3::copy(xyz, p);
  if (rotation != 0) {
    const Real xhat[] = {1,0,0};
    rotate_vector_3(xhat, -rotation, p);
  }
  uv[0] = p[2]/length + 0.5;
  switch (shape) {
  case circle: {
    uv[1] = std::atan2(p[1], p[0]);
    if (uv[1] < 0) uv[1] += 2*M_PI;
  } break;
  case ellipse: {
    uv[1] = std::atan2(p[1]/ell_b, p[0]/ell_a);
    if (uv[1] < 0) uv[1] += 2*M_PI;
  } break;
  default: assert(0); break;
  }
}

void GeneralizedCylinder::get_uhat (Real uhat[3]) const {
  uhat[0] = uhat[1] = 0;
  uhat[2] = 1;
  if (rotation != 0) {
    const Real xhat[] = {1, 0, 0};
    rotate_vector_3(xhat, rotation, uhat);
  }
}

void GeneralizedCylinder::get_normal (const Real u[2], Real nml[3]) const {
  nml[2] = 0;
  const auto c = std::cos(u[1]), s = std::sin(u[1]);
  switch (shape) {
  case circle: {
    nml[0] = c; nml[1] = s;
  } break;
  case ellipse: {
    nml[0] = ell_b*c;
    nml[1] = ell_a*s;
    mv2::normalize(nml);
  } break;
  default: assert(0); break;
  }
  mv3::normalize(nml);
  if (rotation != 0) {
    const Real xhat[] = {1, 0, 0};
    rotate_vector_3(xhat, rotation, nml);
  }
}

void GeneralizedCylinder::get_vhat (const Real u[2], Real tan[3]) const {
  tan[2] = 0;
  const auto c = std::cos(u[1]), s = std::sin(u[1]);
  switch (shape) {
  case circle: {
    tan[0] = -s; tan[1] = c;
  } break;
  case ellipse: {
    tan[0] = -ell_a*s;
    tan[1] =  ell_b*c;
    mv2::normalize(tan);
  } break;
  default: assert(0); break;
  }
  mv3::normalize(tan);
  if (rotation != 0) {
    const Real xhat[] = {1, 0, 0};
    rotate_vector_3(xhat, rotation, tan);
  }
}

Real GeneralizedCylinder::get_jacdet (const Real u[2]) const {
  switch (shape) {
  case circle: {
    return length*radius;
  } break;
  case ellipse: {
    const auto c = std::cos(u[1]), s = std::sin(u[1]);
    return length*std::sqrt(square(ell_a*s) + square(ell_b*c));
  } break;
  default: return 0;
  }
}

Real GeneralizedCylinder::get_area () const {
  switch (shape) {
  case circle: {
    return length*2*M_PI*radius;
  } break;
  case ellipse: {
    Real area = 0;
    acorn::Quadrature q(40, acorn::Quadrature::gl);
    const Real* x, * w;
    q.get_xw(x, w);
    Real uv[2] = {0};
    for (int i = 0; i < q.nq; ++i) {
      uv[1] = M_PI*(1 + x[i]); // 2 pi [(1 + x)/2]
      area += 0.5*w[i]*get_jacdet(uv);
    }
    return 2*M_PI*area;
  } break;
  default: return 0;
  }
}

void RectTessellation
::init (const GeneralizedCylinder& gc, const int nucell_, const int nvcell_) {
  nucell = nucell_;
  nvcell = nvcell_;
  ulength = gc.get_length()/nucell;
  const auto ncell = nucell*nvcell;
  vlengths.resize(nvcell);
  ctrs.resize(3*ncell);
  nmls.resize(3*nvcell);
  vhats.resize(3*nvcell);
  for (int iu = 0; iu < nucell; ++iu) {
    const Real us[] = {Real(iu)/nucell, Real(iu+1)/nucell};
    for (int iv = 0; iv < nvcell; ++iv) {
      const auto icell = this->icell(iu, iv);
      const Real vs[] = {2*M_PI*Real(iv)/nvcell, 2*M_PI*Real(iv+1)/nvcell};
      Real ctr[3] = {0};
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
          const Real uv[] = {us[i], vs[j]};
          Real p[3];
          gc.get_position(uv, p);
          mv3::axpy(1, p, ctr);
        }
      mv3::scale(0.25, ctr);
      mv3::copy(ctr, &ctrs[3*icell]);
      if (iu == 0) {
        Real base[3], uhat[3], vhat[3], nml[3];
        { const Real uv[] = {us[0], vs[0]}; gc.get_position(uv, base); }
        { const Real uv[] = {us[1], vs[0]}; gc.get_position(uv, uhat); }
        { const Real uv[] = {us[0], vs[1]}; gc.get_position(uv, vhat); }
        mv3::axpy(-1, base, uhat);
        mv3::axpy(-1, base, vhat);
        mv3::cross(vhat, uhat, nml);
        mv3::normalize(nml);
        mv3::copy(nml, &nmls[3*iv]);
        vlengths[iv] = mv3::norm2(vhat);
        mv3::normalize(vhat);
        mv3::copy(vhat, &vhats[3*iv]);
        if (iv == 0) {
          mv3::normalize(uhat);
          mv3::copy(uhat, this->uhat);
        }
      }
    }
  }
}

int RectTessellation::ncell () const { return nucell*nvcell; }
int RectTessellation::get_nucell () const { return nucell; }
int RectTessellation::get_nvcell () const { return nvcell; }

static void ic2iuiv (const int nvcell, const int icell, int& iu, int& iv) {
  iu = icell / nvcell;
  iv = icell % nvcell;
}

const Real* RectTessellation::get_cellctr (const int icell) const {
  return &ctrs[3*icell];
}

const Real* RectTessellation::get_normal (const int icell) const {
  int iu, iv;
  ic2iuiv(nvcell, icell, iu, iv);
  return &nmls[3*iv];
}

const Real* RectTessellation::get_uhat (const int icell) const { return uhat; }

const Real* RectTessellation::get_vhat (const int icell) const {
  int iu, iv;
  ic2iuiv(nvcell, icell, iu, iv);
  return &vhats[3*iv];
}

Real RectTessellation::get_ulength (const int icell) const { return ulength; }

Real RectTessellation::get_vlength (const int icell) const {
  int iu, iv;
  ic2iuiv(nvcell, icell, iu, iv);
  return vlengths[iv];
}

void RectTessellation::get_rect (const int icell, Real corners[12]) const {
  const auto ctr = get_cellctr(icell);
  const auto vhat = get_vhat(icell);
  const auto uhat = get_uhat(icell);
  const auto ulen = get_ulength(icell);
  const auto vlen = get_vlength(icell);
  mv3::sum3(1, ctr, -ulen/2, uhat, -vlen/2, vhat, corners  );
  mv3::sum3(1, ctr, -ulen/2, uhat,  vlen/2, vhat, corners+3);
  mv3::sum3(1, ctr,  ulen/2, uhat,  vlen/2, vhat, corners+6);
  mv3::sum3(1, ctr,  ulen/2, uhat, -vlen/2, vhat, corners+9);
}

void RectTessellation::get_rect_uv (const int icell, Real corners[8]) const {
  int iu, iv;
  ic2iuiv(nvcell, icell, iu, iv);
  const Real du = 1.0/nucell, dv = 2*M_PI/nvcell;
  corners[0] =  iu   *du; corners[1] =  iv   *dv;
  corners[2] =  iu   *du; corners[3] = (iv+1)*dv;
  corners[4] = (iu+1)*du; corners[5] = (iv+1)*dv;
  corners[6] = (iu+1)*du; corners[7] =  iv   *dv;
}

void global2local_sym_tensor (const RectTessellation& rt, RPtr sigmas) {
  const auto ncell = rt.ncell();
  ompparfor for (int i = 0; i < ncell; ++i)
    global2local_sym_tensor(rt, i, &sigmas[6*i]);
}

void global2local_sym_tensor (const RectTessellation& rt, const int idx,
                              Real sigma[6]) {
  rotate_sym_tensor_3x3_RARt(rt.get_vhat(idx), rt.get_uhat(idx),
                             rt.get_normal(idx), sigma);
}


SigmaType::Enum SigmaType::convert (const int i) {
  auto s = static_cast<SigmaType::Enum>(i);
  if ( ! is_valid(s)) s = invalid;
  return s;
}

SigmaType::Enum SigmaType::convert (const std::string& e) {
  if (e == "fs_okada") return fs_okada;
  if (e == "fs") return fs;
  if (e == "fs_exact_geometry_disloc") return fs_exact_geometry_disloc;
  return invalid;
}

std::string SigmaType::convert (const SigmaType::Enum e) {
  switch (e) {
  case fs_okada: return "fs_okada";
  case fs: return "fs";
  case fs_exact_geometry_disloc: return "fs_exact_geometry_disloc";
  case invalid:
  default: return "invalid";
  }
}

bool SigmaType::is_valid (const Enum e) {
  return e >= fs_okada && e < invalid;
}

TestUvFn::Fn TestUvFn::convert (const int i) {
  auto f = static_cast<Fn>(i);
  if ( ! is_valid(f)) f = invalid;
  return f;
}

TestUvFn::Fn TestUvFn::convert (const std::string& fn) {
  if (fn == "constant") return Fn::constant;
  if (fn == "cosine_bell") return Fn::cosine_bell;
  if (fn == "tapered") return Fn::tapered;
  if (fn == "hat") return Fn::hat;
  if (fn == "sin") return Fn::sin;
  return Fn::invalid;
}

std::string TestUvFn::convert (const TestUvFn::Fn fn) {
  switch (fn) {
  case Fn::constant: return "constant";
  case Fn::cosine_bell: return "cosine_bell";
  case Fn::tapered: return "tapered";
  case Fn::hat: return "hat";
  case Fn::sin: return "sin";
  case Fn::invalid:
  default: return "invalid";
  }
}

bool TestUvFn::is_valid (const Fn fn) {
  return fn >= constant && fn < invalid;
}

TestUvFn::TestUvFn (const Fn ufns[3], const Real uparms[9],
                    const Fn vfns[3], const Real vparms[9]) {
  init(ufns, uparms, vfns, vparms);
}

void TestUvFn::init (const Fn ufns_[3], const Real uparms_[9],
                     const Fn vfns_[3], const Real vparms_[9]) {
  copy(3, ufns_, ufns);
  copy(3, vfns_, vfns);
  copy(9, uparms_, uparms);
  copy(9, vparms_, vparms);
}

int TestUvFn::nfunctions() const { return 3; }

void TestUvFn::eval (const Real uv[2], RPtr f) const {
  for (int i = 0; i < 3; ++i)
    f[i] = (eval(ufns[i], &uparms[3*i], 2*M_PI*uv[0])*
            eval(vfns[i], &vparms[3*i],        uv[1]));
}

Real TestUvFn::eval (const Fn fn, const Real p[3], Real x) {
  while (x < 0) x += 2*M_PI;
  switch (fn) {
  case constant: return p[0];
  case cosine_bell:
  case tapered:
  case hat: {
    const Real ctr = p[1], w = p[2];
    assert(w > 0);
    if (x >= ctr + w/2 || x <= ctr - w/2) return 0;
    return (fn == cosine_bell ?
            0.5*p[0]*(1 + std::cos(2*M_PI*(x - ctr)/w)) :
            fn == tapered ?
            p[0]*std::pow(1 - square(2*(x - ctr)/w), 1.5):
            p[0]);
  }
  case sin: {
    assert(p[1] != 0);
    return p[0]*std::sin(p[1]*x + p[2]);
  }
  case invalid:
  default: assert(0); return 0;
  }
}

void eval_sigma_at_ctrs (
  const GeneralizedCylinder& gc, const RectTessellation& rt, const Real lam,
  const Real mu, const CallerUvFn& dislocf, const SigmaType::Enum stype,
  RPtr sigmas, RPtr dislocs, const QuadratureParams qp)
{
  assert(SigmaType::is_valid(stype));
  const auto ncell = rt.ncell();
  ompparfor for (int ir = 0; ir < ncell; ++ir) {
    auto* const sigma = &sigmas[6*ir];
    for (int i = 0; i < 6; ++i) sigma[i] = 0;
    const auto* rcv = rt.get_cellctr(ir);
    for (int is = 0; is < ncell; ++is) {
      Real s[6];
      const auto* const src = rt.get_cellctr(is);
      const auto* const nml = rt.get_normal(is);
      const auto* const tan = rt.get_vhat(is);
      const auto* const uhat = rt.get_uhat(is);
      Real suv[2], disloc_uv[3];
      gc.get_uv(src, suv);
      dislocf.eval(suv, disloc_uv);
      if (dislocs && is == ir) copy(3, disloc_uv, &dislocs[3*is]);
      Real disloc[3];
      mv3::sum3(disloc_uv[0], tan, disloc_uv[1], uhat, disloc_uv[2], nml,
                disloc);
      const Real xy_side_lens[] = {rt.get_vlength(is), rt.get_ulength(is)};
      switch (stype) {
      case SigmaType::fs_okada:
#ifdef WOODLAND_ACORN_HAVE_DC3D
        {
          fs3d::calc_sigma_const_disloc_rect_okada(
            lam, mu, src, nml, tan, xy_side_lens, disloc, rcv, s);
        } break;
#endif
      case SigmaType::fs: {
        fs3d::calc_sigma_const_disloc_rect(
          lam, mu, src, nml, tan, xy_side_lens, disloc, rcv, s,
          qp.np_radial, qp.np_angular, qp.triquad_order);
      } break;
      case SigmaType::fs_exact_geometry_disloc: {
        Real cc_uv[2];
        gc.get_uv(rcv, cc_uv);
        SigmaExactOptions o;
        o.qp = qp;
        eval_sigma_exact(gc, rt, lam, mu, dislocf, is, 1, cc_uv, s, o);
      } break;
      case SigmaType::invalid:
      default: assert(0);
      }
      for (int i = 0; i < 6; ++i) sigma[i] += s[i];
    }
  }
}

void eval_sigma (
  const GeneralizedCylinder& gc, const RectTessellation& rt,
  const Real lam, const Real mu, const CallerUvFn& dislocf,
  const int is, const int ircv, const int nurcv, const int nvrcv,
  const SigmaType::Enum stype, Real disloc_uv[3], RPtr sigmas,
  const QuadratureParams qp)
{
  assert(SigmaType::is_valid(stype));
  const auto ncell = rt.ncell();
  assert(is >= 0 && is < ncell);
  assert(ircv >= 0 && ircv < ncell);

  const auto* const src = rt.get_cellctr(is);
  const auto* const nml = rt.get_normal(is);
  const auto* const tan = rt.get_vhat(is);
  const auto* const uhat = rt.get_uhat(is);
  Real suv[2];
  gc.get_uv(src, suv);
  dislocf.eval(suv, disloc_uv);
  Real disloc[3];
  mv3::sum3(disloc_uv[0], tan, disloc_uv[1], uhat, disloc_uv[2], nml, disloc);

  Real rcrnrs[12];
  rt.get_rect(ircv, rcrnrs);
  
  const Real xy_side_lens[] = {rt.get_vlength(is), rt.get_ulength(is)};

  for (int iu = 0, kr = 0; iu < nurcv; ++iu) {
    const Real a = Real(iu + 1)/(nurcv + 1);
    for (int iv = 0; iv < nvrcv; ++iv, ++kr) {
      const Real b = Real(iv + 1)/(nvrcv + 1);
      Real rcv[3] = {0};
      mv3::axpy((1-a)*(1-b), rcrnrs  , rcv);
      mv3::axpy((1-a)*   b , rcrnrs+3, rcv);
      mv3::axpy(   a *   b , rcrnrs+6, rcv);
      mv3::axpy(   a *(1-b), rcrnrs+9, rcv);
      Real s[6];
      switch (stype) {
      case SigmaType::fs_okada:
#ifdef WOODLAND_ACORN_HAVE_DC3D
        fs3d::calc_sigma_const_disloc_rect_okada(
          lam, mu, src, nml, tan, xy_side_lens, disloc, rcv, s);
        break;
#endif
      case SigmaType::fs:
        fs3d::calc_sigma_const_disloc_rect(
          lam, mu, src, nml, tan, xy_side_lens, disloc, rcv, s,
          qp.np_radial, qp.np_angular, qp.triquad_order);
        break;
      case SigmaType::invalid:
      default: assert(0);
      }
      auto* const sigma = &sigmas[6*kr];
      for (int i = 0; i < 6; ++i) sigma[i] += s[i];
    }
  }
}

SigmaExactOptions::SigmaExactOptions ()
  : qp(), extrap_npt(5), hfp_dist_fac(0.5)
{}

static void
rcv_distance (const Real ur[2], const Real vr[2], Real uv[2], Real& dist) {
  Real u_dist = std::max(ur[0] - uv[0], uv[0] - ur[1]);
  Real v_dist = 2*M_PI, v_min = uv[1];
  for (const Real v_os : {-2*M_PI, 0.0, 2*M_PI}) {
    const Real v = uv[1] + v_os;
    const Real dist = std::max(vr[0] - v, v - vr[1]);
    if (dist < v_dist) {
      v_dist = dist;
      v_min = v;
    }
  }
  uv[1] = v_min;
  dist = std::max(u_dist, v_dist);
}

struct SigmaExactIntegrands : public CallerIntegrands {
  SigmaExactIntegrands (
    const GeneralizedCylinder& gc_, const RectTessellation& rt_,
    const Real lam_, const Real mu_, const CallerUvFn& fdisloc_,
    const int isrc_, const Pt rcv_uv_, const SigmaExactOptions& o)
    : npt(o.extrap_npt), gc(gc_), rt(rt_), fdisloc(fdisloc_), lam(lam_),
      mu(mu_), isrc(isrc_)
  {
    assert(npt >= 2);
    assert(npt <= max_npt);
    gc.get_position(rcv_uv_, rcv);
    rt.get_rect_uv(isrc, rect_uv);
    rcv_uv[0] = rcv_uv_[0];
    rcv_uv[1] = rcv_uv_[1];
    const Real u_range[] = {rect_uv[0], rect_uv[4]};
    const Real v_range[] = {rect_uv[1], rect_uv[3]};
    rcv_distance(u_range, v_range, rcv_uv, dist);
    L = std::max(u_range[1] - u_range[0],
                 v_range[1] - v_range[0]);

    if (dist > 0) {
      // For self-other interactions, establish a smooth dislocation
      // approximation valid both inside and outside the src cell.
      Quadrature q(npt, Quadrature::gl);
      const Real* qx;
      q.get_x(qx);
      for (int iu = 0; iu < npt; ++iu) {
        const auto a = (1 + qx[iu])/2;
        ip.u[iu] = (1-a)*rect_uv[0] + a*rect_uv[4];
        for (int iv = 0; iv < npt; ++iv) {
          const auto b = (1 + qx[iv])/2;
          ip.v[iv] = (1-b)*rect_uv[1] + b*rect_uv[3];
          const Pt p_uv{ip.u[iu], ip.v[iv]};
          Real disloc_lcl[3];
          fdisloc.eval(p_uv, disloc_lcl);
          const int i = npt*iu + iv;
          for (int j = 0; j < 3; ++j)
            ip.f[j][i] = disloc_lcl[j];
        }
      }
    }
  }

  virtual int nintegrands () const override { return 6; }

  virtual Real permitted_R_min (const Real R_max) const override {
    return 1e-3*R_max;
  }

  virtual void eval (const int n, CRPtr p_uvs, RPtr integrand) const override {
    for (int i = 0; i < n; ++i) {
      Real p[3], vhat[3], nml[3], uhat[3], disloc_lcl[3], disloc[3];
      const auto p_uv = &p_uvs[2*i];
      gc.get_position(p_uv, p);
      gc.get_vhat(p_uv, vhat);
      gc.get_normal(p_uv, nml);
      mv3::cross(nml, vhat, uhat);
      if (dist > 0) {
        switch (npt) {
        case 2: extrap<2>(p_uv, disloc_lcl); break;
        case 3: extrap<3>(p_uv, disloc_lcl); break;
        case 4: extrap<4>(p_uv, disloc_lcl); break;
        case 5: extrap<5>(p_uv, disloc_lcl); break;
        case 6: extrap<6>(p_uv, disloc_lcl); break;
        case 7: extrap<7>(p_uv, disloc_lcl); break;
        case 8: extrap<8>(p_uv, disloc_lcl); break;
        default: assert(0);
        }
      } else {
        fdisloc.eval(p_uv, disloc_lcl);
      }
      mv3::sum3(disloc_lcl[0], vhat, disloc_lcl[1], uhat, disloc_lcl[2], nml,
                disloc);
      const auto ig = &integrand[6*i];
      fs3d::calc_sigma_point(lam, mu, p, nml, disloc, rcv, ig);
      const Real jacdet = gc.get_jacdet(p_uv);
      for (int j = 0; j < 6; ++j) ig[j] *= jacdet;
    }
  }

  plane::Polygon get_src_polygon () const {
    return plane::Polygon(rect_uv, 4);
  }

  int npt;
  const GeneralizedCylinder& gc;
  const RectTessellation& rt;
  const CallerUvFn& fdisloc;
  const Real lam, mu;
  const int isrc;
  Pt rcv_uv;
  Real rcv[3];
  Real rect_uv[8];
  Real dist, L;

private:
  enum : int { max_npt = 8 };
  struct {
    Real u[max_npt], v[max_npt], f[3][max_npt*max_npt];
  } ip;

  template<int N>
  void extrap (const Pt p_uv, Real disloc_lcl[3]) const {
    for (int j = 0; j < 3; ++j) {
      Real f[N];
      for (int i = 0; i < N; ++i)
        f[i] = eval_lagrange_poly(N, ip.v, &ip.f[j][N*i], p_uv[1]);
      disloc_lcl[j] = eval_lagrange_poly(N, ip.u, f, p_uv[0]);
    }    
  }
};

static void eval_sigma_exact (const SigmaExactIntegrands& s, RPtr sigma_out,
                              const SigmaExactOptions o) {
  Real sigma[6] = {0};
  const plane::Polygon p = s.get_src_polygon();
  if (s.dist < o.hfp_dist_fac*s.L) {
    integrals::Options io;
    io.np_radial = o.qp.np_radial;
    io.np_angular = o.qp.np_angular;
    integrals::calc_hfp(io, p, s.rcv_uv, s, sigma);
  } else {
    const auto triquad_order = o.qp.triquad_order <= 0 ?
      get_triquad_order(s.L, s.dist) : o.qp.triquad_order;
    integrals::calc_integral(p, s, sigma, triquad_order);
  }
  for (int i = 0; i < 6; ++i) sigma_out[i] = sigma[i];
}

void eval_sigma_exact (
  const GeneralizedCylinder& gc, const RectTessellation& rt,
  const Real lam, const Real mu, const CallerUvFn& fdisloc,
  const int isrc, const int nuv, CRPtr uvs, RPtr sigmas,
  const SigmaExactOptions o, const bool thread)
{
  const auto ncell = rt.ncell();
  assert(isrc >= 0 && isrc < ncell);

  if (thread) {
    ompparfor for (int ir = 0; ir < nuv; ++ir) {
      SigmaExactIntegrands sei(gc, rt, lam, mu, fdisloc, isrc, &uvs[2*ir], o);
      eval_sigma_exact(sei, &sigmas[6*ir], o);
    }
  } else {
    for (int ir = 0; ir < nuv; ++ir) {
      SigmaExactIntegrands sei(gc, rt, lam, mu, fdisloc, isrc, &uvs[2*ir], o);
      eval_sigma_exact(sei, &sigmas[6*ir], o);
    }
  }
}

void eval_sigma_exact_at_ctrs (
  const GeneralizedCylinder& gc,
  const RectTessellation& rt_src, const RectTessellation& rt_rcv,
  const Real lam, const Real mu, const CallerUvFn& fdisloc,
  RPtr sigmas, const SigmaExactOptions o)
{
  const auto ncell_rcv = rt_rcv.ncell();
  const auto ncell_src = rt_src.ncell();
  ompparfor for (int ir = 0; ir < ncell_rcv; ++ir) {
    const auto* rcv = rt_rcv.get_cellctr(ir);
    Real cc_rcv_uv[2];
    gc.get_uv(rcv, cc_rcv_uv);
    auto* const sigma = &sigmas[6*ir];
    for (int i = 0; i < 6; ++i) sigma[i] = 0;
    for (int is = 0; is < ncell_src; ++is) {
      Real s[6];
      eval_sigma_exact(gc, rt_src, lam, mu, fdisloc, is, 1, cc_rcv_uv, s, o);
      for (int i = 0; i < 6; ++i) sigma[i] += s[i];
    }
  }
}

namespace convtest {
namespace impl {

Errors::Errors (const int nmodel, const int nrefine, const int max_ncell) {
  stypes.resize(nmodel-1);
  sigmas.resize(nmodel);
  for (int i = 0; i < nmodel; ++i) sigmas[i].resize(6*max_ncell, 0);
  const auto resize = [&] (A<A<A<Real>>>& a) {
    a.resize(nmodel-1);
    for (int i = 0; i < nmodel-1; ++i) {
      a[i].resize(6);
      for (int j = 0; j < 6; ++j)
        a[i][j].resize(nrefine, 0);
    }
  };
  resize(l2_num);
  resize(l2_den);
  resize(li_num);
  resize(li_den);
}

void Errors
::accum_errors (const int ncell, const std::function<Real(int)>& cell_area,
                const int iref) {
  for (size_t sti = 0; sti < stypes.size(); ++sti) {
    const auto rs = sigmas[0];
    const auto ms = sigmas[sti+1];
    for (int si = 0; si < 6; ++si) {
      for (int i = 0; i < ncell; ++i) {
        const auto dx = cell_area(i);
        const auto d = ms[6*i+si] - rs[6*i+si];
        l2_num[sti][si][iref] += square(d)*dx;
        l2_den[sti][si][iref] += square(rs[6*i+si])*dx;
        li_num[sti][si][iref] = std::max(li_num[sti][si][iref], std::abs(d));
        li_den[sti][si][iref] = std::max(li_den[sti][si][iref], std::abs(rs[6*i+si]));
      }
      l2_num[sti][si][iref] = std::sqrt(l2_num[sti][si][iref]);
      l2_den[sti][si][iref] = std::sqrt(l2_den[sti][si][iref]);
    }
  }
}

void Errors::collect_errors (const A<int>& sigma_components, Results* r,
                             const bool verbose) {
  if (r) {
    const auto ns = stypes.size();
    const auto nc = sigma_components.size();
    r->components = sigma_components;
    r->l2_errs.resize(ns, std::vector<Real>(nc));
    r->linf_errs.resize(ns, std::vector<Real>(nc));
  }
  for (size_t sti = 0; sti < stypes.size(); ++sti) {
    if (verbose) printf("%s\n", SigmaType::convert(stypes[sti]).c_str());
    for (size_t ci = 0; ci < sigma_components.size(); ++ci) {
      const auto si = sigma_components[ci];
      if (verbose) printf("  %d\n", si);
      const auto nref = l2_num[0][0].size();
      for (size_t iref = 0; iref < nref; ++iref) {
        const auto l2n = l2_num[sti][si][iref];
        const auto l2d = l2_den[sti][si][iref];
        const auto lin = li_num[sti][si][iref];
        const auto lid = li_den[sti][si][iref];
        const auto l2_rel = l2d == 0 ? l2n : l2n/l2d;
        const auto li_rel = lid == 0 ? lin : lin/lid;
        if (verbose)
          printf(" %9.3e (%9.3e) | %9.3e (%9.3e) \n", l2_rel, l2d, li_rel, lid);
        if (r && iref+1 == nref) {
          r->l2_errs[sti][ci] = l2_rel;
          r->linf_errs[sti][ci] = li_rel;
        }
      }
    }
  }
}

namespace pyout {

void init (FILE* const fid, const Errors& e) {
  fprintf(fid, "stypes = ['fs_okada']\n");
  for (const auto& s : e.stypes)
    fprintf(fid, "stypes.append('%s')\n", SigmaType::convert(s).c_str());
  fprintf(fid, "dislocs = []\n");
  fprintf(fid, "sigmas = []\n");
}

void write (FILE* const fid, const Errors& e, const Errors::A<Real>& dislocs,
            const int ncell) {
  assert(dislocs.size() >= (size_t) 3*ncell);
  fprintf(fid, "dislocs.append([");
  for (int i = 0; i < 3*ncell; ++i)
    fprintf(fid, "%1.15e,", dislocs[i]);
  fprintf(fid, "])\n");
  fprintf(fid, "tmp = []\n");
  for (const auto& s : e.sigmas) {
    assert(s.size() >= (size_t) 6*ncell);
    fprintf(fid, "tmp.append([");
    for (int i = 0; i < 6*ncell; ++i)
      fprintf(fid, "%1.15e,", s[i]);
    fprintf(fid, "])\n");
  }
  fprintf(fid, "sigmas.append(tmp)\n");
}

} // namespace pyout
} // namespace impl

namespace gencyl {

Config::Config (const Problem problem) {
  lam = 1;
  mu = 1;
  onedim = problem == p_circle || problem == p_ellipse;
  gc.set_shape(problem == p_circle || problem == p_cylinder ?
               gc.circle : gc.ellipse);
  if (onedim) {
    n_u_cell_base = 1;
    n_v_cell_base = 16;
    n_refine = 7;
    const auto c = TestUvFn::Fn::constant;
    const TestUvFn::Fn ufns[] = {c, c, c};
    const Real uparms[] = {1, 0, 0, 0, 0, 0, 0, 0, 0};
#if 1
    const TestUvFn::Fn vfns[] = {TestUvFn::Fn::tapered, c, c};
    const Real vparms[] = {1, M_PI/2, 1*M_PI, 0, 0, 0, 0, 0, 0};
#else
    const TestUvFn::Fn vfns[] = {TestUvFn::Fn::sin, c, c};
    const Real vparms[] = {1, 1, 0, 0, 0, 0, 0, 0, 0};
#endif
    fdisloc.init(ufns, uparms, vfns, vparms);
  } else {
    n_u_cell_base = 8;
    n_v_cell_base = 32;
    n_refine = 3;
    const auto cb = TestUvFn::Fn::cosine_bell;
    const auto sn = TestUvFn::Fn::sin;
    const TestUvFn::Fn ufns[] = {cb, cb, cb};
    const Real uparms[] = {1, M_PI, 2*M_PI, 1, M_PI, 2*M_PI, 1, M_PI, 2*M_PI};
    const TestUvFn::Fn vfns[] = {sn, sn, sn};
    const Real vparms[] = {1, 1, 0, 0.2, 1, 0, 0.05, 1, 0};
    fdisloc.init(ufns, uparms, vfns, vparms);
  }
  python_outfn = (std::string("pyout_vv_") + (onedim ? "onedim" : "gencyl") +
                  ".py");
}

Results run (const Config& c, const bool verbose) {
  const int max_ncell = ((c.onedim ?
                          1 :
                          (c.n_u_cell_base*(1 << (c.n_refine-1))))*
                         (c.n_v_cell_base*(1 << (c.n_refine-1))));
  impl::Errors e(2, c.n_refine, max_ncell);
  std::vector<Real> dislocs(3*max_ncell, 0);

  e.stypes[0] = SigmaType::fs_exact_geometry_disloc;

  QuadratureParams qp;
  qp.set_np_radial(10).set_np_angular(6);

  FILE* const fid = c.python_outfn.empty() ? nullptr :
    fopen(c.python_outfn.c_str(), "w");
  if (fid) impl::pyout::init(fid, e);
  
  for (int iref = 0; iref < c.n_refine; ++iref) {
    const int nucell = c.onedim ? 1 : c.n_u_cell_base*(1 << iref);
    const int nvcell = c.n_v_cell_base*(1 << iref);
    RectTessellation rt;
    rt.init(c.gc, nucell, nvcell);
    {
      const auto stype = SigmaType::fs_okada;
      eval_sigma_at_ctrs(c.gc, rt, c.lam, c.mu, c.fdisloc, stype,
                         e.sigmas[0].data(), dislocs.data());
      global2local_sym_tensor(rt, e.sigmas[0].data());
    }
    {
      const int fac = 8;
      RectTessellation rt_src;
      assert(c.onedim || nucell % fac == 0);
      assert(nvcell % fac == 0);
      rt_src.init(c.gc, c.onedim ? 1 : nucell/fac, nvcell/fac);
      SigmaExactOptions o;
      o.qp = qp;
      // The circle problem creates skinny rectangle cells as iref
      // increases. Need to be careful with extrap. This shouldn't happen in
      // real problems.
      if (c.onedim) o.extrap_npt = 3;
      eval_sigma_exact_at_ctrs(c.gc, rt_src, rt, c.lam, c.mu, c.fdisloc,
                               e.sigmas[1].data(), o);
      global2local_sym_tensor(rt, e.sigmas[1].data());
    }
    const Real area = ((c.gc.get_length()/nucell)*
                       (2*M_PI*c.gc.get_radius()/nvcell));
    const auto ncell = nucell*nvcell;
    e.accum_errors(ncell, [&] (int) { return area; }, iref);
    if (fid) impl::pyout::write(fid, e, dislocs, ncell);
  }

  if (fid) fclose(fid);

  if (verbose)
    printf("%s %s\nStress components are in element-local coordinates.\n",
           c.onedim ? "1D" : "gencyl", c.gc.convert(c.gc.get_shape()).c_str());
  Results r;
  const auto comps = (c.onedim ?
                      std::vector<int>({0,2,3,5}) :
                      std::vector<int>({0,1,2,3,4,5}));
  e.collect_errors(comps, &r, verbose);

  return r;
}

} // namespace gencyl
} // namespace convtest

int unittest () {
#ifndef WOODLAND_ACORN_HAVE_DC3D
  printf("vv::unittest: Because extern/dc3d.f is not available, "
         "this unit test will purposely fail.\n");
#endif
  int ne = 0;
  {
    GeneralizedCylinder gc;
    for (const auto shape : {gc.circle, gc.ellipse}) {
      gc.set_shape(shape).set_length(1.5).set_radius(0.8)
        .set_ellipse(1.7, 0.9).set_rotation(-0.1);
      for (int i = 0; i < 100; ++i) {
        const Real uv[] = {0.3, 2*M_PI*Real(i)/100};
        Real xyz[3];
        gc.get_position(uv, xyz);
        Real uv1[2];
        gc.get_uv(xyz, uv1);
        if (reldif(2, uv, uv1) > 2*mv3::eps) ++ne;
      }
      {
        RectTessellation rt;
        rt.init(gc, 3, 15);
        const auto nc = rt.ncell();
        if (nc != 45) ++ne;
        for (int ic = 0; ic < nc; ++ic) {
          Real rect[12];
          rt.get_rect(ic, rect);
        }
      }
      {
        RectTessellation rt;
        rt.init(gc, 1, 3);
        Real sigmas[3*6];
        using Fn = TestUvFn::Fn;
        const Fn ufns[] = {Fn::constant, Fn::constant, Fn::constant};
        const Fn vfns[] = {Fn::cosine_bell, Fn::constant, Fn::constant};
        Real uparms[9] = {0}; uparms[0] = 1;
        Real vparms[9] = {0}; vparms[0] = 1; vparms[1] = M_PI; vparms[2] = M_PI/5;
        TestUvFn disloc(ufns, uparms, vfns, vparms);
        eval_sigma_at_ctrs(gc, rt, 1, 1, disloc, SigmaType::fs, sigmas);
      }
    }
  }
  const Real lo = 0.99, hi = 1.01;
  const bool verbose = false;
  const auto cmp_l2s = [&] (CRPtr l2s, const convtest::Results& r, const int n,
                            int& ne) {
    for (int i = 0; i < n; ++i) {
      const auto v = r.l2_errs[0][i];
      if (v > hi*l2s[i] || v < lo*l2s[i])
        ++ne;
    }                         
  };
  {
    convtest::gencyl::Config c(convtest::gencyl::Config::p_circle);
    c.python_outfn = "";
    c.n_v_cell_base = 128;
    c.n_refine = 1;
    const auto r = run(c, verbose);
    const Real l2s[] = {6.558e-02, 2.115e-02, 8.357e-02, 3.435e-02};
    cmp_l2s(l2s, r, 4, ne);
  }
  {
    convtest::gencyl::Config c(convtest::gencyl::Config::p_ellipse);
    c.python_outfn = "";
    c.n_v_cell_base = 128;
    c.n_refine = 1;
    const auto r = run(c, verbose);
    const Real l2s[] = {1.609e-01, 1.580e-02, 1.414e-01, 3.157e-02};
    cmp_l2s(l2s, r, 4, ne);
  }
  {
    convtest::flatstrip::Config c;
    c.python_outfn = "";
    c.n_cell_base = 32;
    c.n_refine = 1;
    {
      const auto r = run(c, verbose);
      const Real l2s[] = {1.163e-09, 9.063e-03};
      for (int i = 0; i < 2; ++i) {
        const auto v = r.l2_errs[i][0];
        if (v > hi*l2s[i] || v < lo*l2s[i])
          ++ne;
      }
    }
    {
      c.n_cell_base = 64;
      c.run_fs = false;
      const auto r = run(c, verbose);
      const Real l2s[] = {2.933e-03};
      for (int i = 0; i < 1; ++i) {
        const auto v = r.l2_errs[i][0];
        if (v > hi*l2s[i] || v < lo*l2s[i])
          ++ne;
      }
    }
  }
  {
    convtest::gencyl::Config c(convtest::gencyl::Config::p_cylinder);
    c.python_outfn = "";
    c.n_refine = 1;
    const auto r = run(c, verbose);
    const Real l2s[] = {1.884e-01, 5.113e-02, 3.871e-02, 1.275e-01, 6.760e-02,
                        1.393e-01};
    cmp_l2s(l2s, r, 6, ne);
  }
  {
    convtest::gencyl::Config c(convtest::gencyl::Config::p_ellipse_cylinder);
    c.python_outfn = "";
    c.n_refine = 1;
    const auto r = run(c, verbose);
    const Real l2s[] = {3.766e-01, 8.234e-02, 5.171e-02, 2.463e-01, 1.086e-01,
                        2.883e-01};
    cmp_l2s(l2s, r, 6, ne);
  }
  return ne;
}

} // namespace vv
} // namespace acorn
} // namespace woodland
