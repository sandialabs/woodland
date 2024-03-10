#include "woodland/acorn/compose_quadrature.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/dbg.hpp"

#include "woodland/examples/convzx/convtest_zx.hpp"
#include "woodland/examples/convzx/convtest_zx_hmmvp.hpp"
#include "woodland/examples/convzx/pywrite.hpp"

#include <cmath>

namespace woodland {
namespace examples {
namespace convzx {

typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;

namespace {

void pywrite (FILE* fp, const std::string& dict, const Triangulation& t) {
  fprintf(fp, "%s['tri'] = npy.array([\n", dict.c_str());
  const auto ntri = t.get_ntri();
  for (int ti = 0; ti < ntri; ++ti) {
    const auto tri = t.get_tri(ti);
    fprintf(fp, "[%d, %d, %d],\n", int(tri[0]), int(tri[1]), int(tri[2]));
  }
  fprintf(fp, "], dtype=npy.int32)\n");
  fprintf(fp, "%s['vtx'] = npy.array([\n", dict.c_str());
  const auto nvtx = t.get_nvtx();
  for (int vi = 0; vi < nvtx; ++vi) {
    const auto vtx = t.get_vtx(vi);
    fprintf(fp, "[%15.8e, %15.8e, %15.8e],\n", vtx[0], vtx[1], vtx[2]);
  }
  fprintf(fp, "])\n");
}

} // namespace

void ConvTest::init (const ZxFn::Shape shape, const Disloc::CPtr& disloc_) {
  zxfn = std::make_shared<ZxFn>(shape);
  disloc = disloc_;
#ifndef WOODLAND_ACORN_HAVE_DC3D
  use_woodland_rg0c0 = true;
#endif
}

void ConvTest::set_nx (const int nx) {
  assert(nx >= 1);
  if (zxfn->get_nx() != nx) { t = nullptr; d = nullptr; }
  zxfn->set_nx(nx);
}

void ConvTest::set_ny (const int ny) {
  assert(ny >= 1);
  if (zxfn->get_ny() != ny) { t = nullptr; d = nullptr; }
  zxfn->set_ny(ny);
}

void ConvTest::set_use_four_tris_per_rect(const bool use) {
  if ((use ? 4 : 2) != ntri_per_rect) { t = nullptr; d = nullptr; }
  ntri_per_rect = use ? 4 : 2;
}

void ConvTest::set_use_surface_recon (const bool use) {
  if (use != use_surface_recon) { t = nullptr; d = nullptr; }
  use_surface_recon = use;
}

void ConvTest::set_use_exact_normals (const bool use) {
  if (use != use_exact_normals) { t = nullptr; d = nullptr; }
  use_exact_normals = use;
}

void ConvTest::set_hmmvp_tol (const Real tol) {
  hmmvp_tol = tol;
}

void ConvTest::set_normal_recon_order (const int order) {
  if (order != nml_recon_order) { t = nullptr; d = nullptr; }
  nml_recon_order = order;
}

void ConvTest::set_disloc_order (const int order) {
  if (order != disloc_order) { t = nullptr; d = nullptr; }
  disloc_order = order;
}

void ConvTest::set_use_flat_elements (const bool use) {
  if (use != use_flat_elements) { t = nullptr; d = nullptr; }
  use_flat_elements = use;
}

void ConvTest::set_use_woodland_rg0c0 (const bool use) {
#ifndef WOODLAND_ACORN_HAVE_DC3D
  if (not use) {
    printf("WARNING: Using Woodland rg0c0 impl b/c dc3d impl is missing.");
    use_woodland_rg0c0 = true;
    return;
  }
#endif
  use_woodland_rg0c0 = use;
}

void ConvTest::print (FILE* fp) const {
  fprintf(fp, "ct> lam %1.3e mu %1.3e\n", lam, mu);
  fprintf(fp, "ct> ntriperrect %d dislocorder %d exactsrf %d flatelem %d"
          " exactnml %d nmlorder %d\n",
          ntri_per_rect, disloc_order, int(not use_surface_recon),
          int(use_flat_elements), int(use_exact_normals), nml_recon_order);
  fprintf(fp, "ct> zshape %s\n", ZxFn::convert(zxfn->get_shape()).c_str());
  if (use_woodland_rg0c0) printf("ct> using woodland impl of rg0c0\n");
  if (disloc)
    for (int i = 0; i < 3; ++i) {
      const auto& d = disloc->ds[i];
      const Real* const a = &d.amplitude;
      bool all0 = true;
      for (int k = 0; k < 9; ++k)
        if (a[k] != 0)
          all0 = false;
      if (all0) continue;
      fprintf(fp, "ct> disloc %d %s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f "
              "%6.3f %6.3f\n",
              i, Disloc::convert(d.shape).c_str(),
              a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8]);
    }
}

Triangulation::Ptr ConvTest
::triangulate (const ZxFn& zxfn, const int ntri_per_rect) {
  assert(zxfn.get_nx() >= 1 && zxfn.get_ny() >= 1);

  const auto t = std::make_shared<Triangulation>();
  t->begin_modifications();

  const auto nx = zxfn.get_nx(), ny = zxfn.get_ny(), pnx = nx+1,
    os = pnx*(ny+1);

  for (int yi = 0, pat = 0; yi < ny; ++yi) {
    for (int xi = 0; xi < nx; ++xi, ++pat) {
      if (ntri_per_rect == 4) {
        const int tris[][3] =
          {{pnx* yi    + xi+1, pnx*(yi+1) + xi+1, os + nx*yi + xi},
           {pnx*(yi+1) + xi  , pnx* yi    + xi  , os + nx*yi + xi},
           {pnx* yi    + xi  , pnx* yi    + xi+1, os + nx*yi + xi},
           {pnx*(yi+1) + xi+1, pnx*(yi+1) + xi  , os + nx*yi + xi}};
        for (int i = 0; i < 4; ++i) t->append_tri(tris[i]);
      } else {
        // Don't want the alternation: breaks symmetry and then can't use
        // fast-eval.
        if (1) { //(pat % 2 == 0) {
          const int tris[][3] =
            {{pnx*yi + xi, pnx*(yi+1) + xi+1, pnx*(yi+1) + xi  },
             {pnx*yi + xi, pnx* yi    + xi+1, pnx*(yi+1) + xi+1}};
          for (int i = 0; i < 2; ++i) t->append_tri(tris[i]);
        } else {
          const int tris[][3] =
            {{pnx*yi + xi,   pnx* yi    + xi+1, pnx*(yi+1) + xi},
             {pnx*yi + xi+1, pnx*(yi+1) + xi+1, pnx*(yi+1) + xi}};
          for (int i = 0; i < 2; ++i) t->append_tri(tris[i]);
        }
      }
    }
  }

  for (int yi = 0; yi <= ny; ++yi)
    for (int xi = 0; xi <= nx; ++xi)
      t->append_vtx();
  if (ntri_per_rect == 4)
    for (int yi = 0; yi < ny; ++yi)
      for (int xi = 0; xi < nx; ++xi)
        t->append_vtx();

  t->end_modifications();

  const auto xs = zxfn.get_xbs();
  const auto zs = zxfn.get_zbs();
  for (int xi = 0; xi <= nx; ++xi) {
    const auto x = xs[xi];
    const auto z = zs[xi];
    for (int yi = 0; yi <= ny; ++yi) {
      const Real y = Real(yi)/ny;
      Real* vtx = t->get_vtx(pnx*yi + xi);
      vtx[0] = x;
      vtx[1] = y;
      vtx[2] = z;
    }
  }
  if (ntri_per_rect == 4) {
    const auto zshape = zxfn.get_shape();
    for (int xi = 0; xi < nx; ++xi) {
      const auto x = (xs[xi] + xs[xi+1])/2;
      Real z, unused;
      gallery::eval(zshape, x, z, unused);
      for (int yi = 0; yi < ny; ++yi) {
        const Real y = (Real(yi)/ny + Real(yi+1)/ny)/2;
        Real* vtx = t->get_vtx(os + nx*yi + xi);
        vtx[0] = x;
        vtx[1] = y;
        vtx[2] = z;
      }
    }
  }

  return t;
}

namespace {
struct Param2DUserFn : public ExactParam2DSurface::UserFn {
  Param2DUserFn (const ZxFn::Shape shape_) : shape(shape_) {}
  
  void eval (const Real x, const Real y, Real& f, RPtr g) const override {
    Real gx;
    gallery::eval(shape, x, f, gx);
    if (g) {
      g[0] = gx;
      g[1] = 0;
    }
  }

  ZxFn::Shape shape;
};
} // namespace

Discretization::Ptr ConvTest
::discretize (const Triangulation::Ptr& t, const ZxFn::Shape shape,
              const bool use_surface_recon, const bool use_exact_normals,
              const int nml_recon_order, const bool use_flat_elements) {
  const auto tr = std::make_shared<TriangulationRelations>(t);
  Surface::CPtr s;
  if (use_surface_recon) {
    if (use_flat_elements) {
      const Real primary[] = {1,0,0};
      s = std::make_shared<FlatElementSurface>(t, primary);
    } else {
      s = std::make_shared<ExtrudedCubicSplineSurface>(
        shape, t, use_exact_normals, nml_recon_order);
    }
  } else {
    const auto ufn = std::make_shared<Param2DUserFn>(shape);
    const auto srf = std::make_shared<ExactParam2DSurface>(t, ufn);
    srf->reset_triangulation_z();
    s = srf;
  }
  return std::make_shared<Discretization>(t, tr, s);
}

void ConvTest
::pywrite (const std::string& python_filename, const RealArray& dislocs,
           const RealArray& sigmas, const std::string& dict,
           const bool append) const {
  assert(dislocs.empty() || 3*t->get_ntri() == int(dislocs.size()));
  assert(sigmas .empty() || 6*t->get_ntri() == int(sigmas .size()));
  assert(not (dislocs.empty() && sigmas.empty()));
  FILE* fp = fopen(python_filename.c_str(), append ? "a" : "w");
  if (not append) pywrite_header(fp);
  fprintf(fp, "%s = {}\n", dict.c_str());
  fprintf(fp, "%s['type'] = 'tri'\n", dict.c_str());
  const auto nx = zxfn->get_nx();
  fprintf(fp, "%s['nx'] = %d\n", dict.c_str(), nx);
  std::string tdi = dict + "['t']";
  fprintf(fp, "%s = {}\n", tdi.c_str());
  convzx::pywrite(fp, tdi, *t);
  std::string ddi = dict + "['disloc']";
  if (dislocs.empty())
    fprintf(fp, "%s = None\n", ddi.c_str());
  else
    pywrite_double_array(fp, ddi, dislocs.size()/3, 3, dislocs.data());
  std::string sdi = dict + "['sigma']";
  if (sigmas.empty())
    fprintf(fp, "%s = None\n", sdi.c_str());
  else
    pywrite_double_array(fp, sdi, sigmas.size()/6, 6, sigmas.data());
  fclose(fp);
}

void ConvTest
::pywrite_okada (const std::string& python_filename, const int nxr, const int nyr,
                 const RealArray& dislocs, const RealArray& sigmas,
                 const std::string& dict, const bool append) const {
  const int nrect = nxr*nyr;
  assert(dislocs.empty() || int(dislocs.size()) == 3*nrect);
  assert(sigmas .empty() || int(sigmas .size()) == 6*nrect);
  assert(not (dislocs.empty() && sigmas.empty()));
  FILE* fp = fopen(python_filename.c_str(), append ? "a" : "w");
  if (not append) pywrite_header(fp);
  fprintf(fp, "%s = {}\n", dict.c_str());
  fprintf(fp, "%s['type'] = 'rect'\n", dict.c_str());
  fprintf(fp, "%s['nx'] = %d\n", dict.c_str(), nxr);
  fprintf(fp, "%s['ny'] = %d\n", dict.c_str(), nyr);
  std::string ddi = dict + "['disloc']";
  if (dislocs.empty())
    fprintf(fp, "%s = None\n", ddi.c_str());
  else
    pywrite_double_array(fp, ddi, nrect, 3, dislocs.data());
  std::string sdi = dict + "['sigma']";
  if (sigmas.empty())
    fprintf(fp, "%s = None\n", sdi.c_str());
  else
    pywrite_double_array(fp, sdi, nrect, 6, sigmas.data());
  fclose(fp);
}

static void
pywrite_zxfn (const std::string& python_filename, const ZxFn& zxfn,
              const std::string& dict, const bool append) {
  FILE* fp = fopen(python_filename.c_str(), append ? "a" : "w");
  if (not append) pywrite_header(fp);
  fprintf(fp, "%s = {}\n", dict.c_str());
  fprintf(fp, "%s['type'] = 'zxfn'\n", dict.c_str());
  fprintf(fp, "%s['nx'] = %d\n", dict.c_str(), zxfn.get_nx());
  std::string ddi = dict + "['x']";
  pywrite_double_array(fp, ddi, 1, zxfn.get_nx(), zxfn.get_xbs());
  ddi = dict + "['z']";
  pywrite_double_array(fp, ddi, 1, zxfn.get_nx(), zxfn.get_zbs());
  fclose(fp);
}

void ConvTest::discretize () {
  t = triangulate(*zxfn, ntri_per_rect);
  d = discretize(t, zxfn->get_shape(), use_surface_recon,
                 use_exact_normals, nml_recon_order, use_flat_elements);
  d->set_disloc_order(disloc_order);
}

static Idx search (const ZxFn& zxfn, const Triangulation& t, const Real p[2]) {
  assert(p[0] > 0 && p[0] < 1 && p[1] > 0 && p[1] < 1);
  const auto nxr = zxfn.get_nx(), nyr = zxfn.get_ny();
  const auto xs = zxfn.get_xbs();
  const int ix = std::lower_bound(xs, xs + nxr + 1, p[0]) - xs - 1;
  assert(ix >= 0); assert(ix < nxr);
  const int iy = std::floor(p[1]*nyr);
  const int ntri_per_rect = t.get_ntri() / (nxr*nyr);
  assert(ntri_per_rect == 2 || ntri_per_rect == 4);
  const int ti_lo = ntri_per_rect*(nxr*iy + ix);
  Real lam[3], lam_lo = 1;
  int i_lo = 0;
  for (int i = 0; i < ntri_per_rect; ++i) {
    Real tvs[6], b[4];
    const auto tri = t.get_tri(ti_lo + i);
    for (int j = 0; j < 3; ++j)
      mv2::copy(t.get_vtx(tri[j]), &tvs[2*j]);
    acorn::Triangle2D::calc_barycentric_matrix(tvs, b);
    acorn::Triangle2D::xy_to_barycentric(&tvs[4], b, p, lam);
    Real out = 0;
    for (int j = 0; j < 3; ++j)
      if (lam[j] > 1 || lam[j] < 0)
        out = std::max(out, std::max(-lam[j], lam[j]-1));
    if (out < lam_lo) {
      lam_lo = out;
      i_lo = i;
    }
  }
  return ti_lo + i_lo;
}

static void
setup (const int testcase, ZxFn::Shape& zshape, Disloc::Ptr& disloc,
       const Real r = 1) {
  disloc = std::make_shared<Disloc>();
  const auto dshape = Disloc::Shape::stapered;
  switch (testcase) {
  case -1:
    zshape = ZxFn::Shape::trig1;
    disloc->set(0, dshape, 1, 0.5, 0.5, 1, 0, r, r);
    break;
  case 10: case 11: case 12:
    zshape = ZxFn::Shape::zero;
    disloc->set(testcase - 10, dshape, 1, 0.5, 0.5, 1, 0, r, r);
    break;
  case 20: case 21: case 22:
    zshape = ZxFn::Shape::trig1;
    disloc->set(testcase - 20, dshape, 1, 0.5, 0.5, 1, 0, r, r);
    break;
  case 30: case 31: case 32:
    zshape = ZxFn::Shape::trig0;
    disloc->set(testcase - 30, dshape, 1, 0.5, 0.5, 1, 0, r, r);
    break;
  case 0: case 1: case 2: case 3:
  default:
    zshape = (testcase == 0 ? ZxFn::Shape::trig1 :
              testcase == 1 ? ZxFn::Shape::zero :
              testcase == 2 ? ZxFn::Shape::trig0 :
              ZxFn::Shape::ramp);
    disloc->set(0, dshape,  0.7, 0.5, 0.5, 1, -0.4, 0.75, 0.65);
    disloc->set(1, dshape, -0.3, 0.4, 0.6, 1,    1, 0.7 , 0.4 );
    disloc->set(2, dshape,  0.1, 0.5, 0.5, 1,    0, 0.75, 0.75);
    break;
  }
  if (not disloc->is_boundary_zero(1e-12)) {
    printf("Disloc is not 0 on boundary; exiting.\n");
    exit(-1);
  }
}

static void
calc_errors (const RealArray& s_true, const RealArray& s, Real l2_err[6],
             Real li_err[6], Surface::CPtr srf = nullptr) {
  const int ncell = int(s_true.size()/6);
  assert(s.size() == s_true.size());
  const auto nthr = acorn::get_max_threads();
  const int naccum = 24;
  std::vector<Real> a(nthr*naccum, 0);
  ompparfor for (int ti = 0; ti < ncell; ++ti) {
    const auto tid = acorn::get_thread_num();
    auto ai = &a[tid*naccum];
    const auto area = srf ? tri_surface_area(*srf, ti) : 1;
    for (int i = 0; i < 6; ++i) {
      const auto den = s_true[6*ti+i];
      const auto diff = s[6*ti+i] - den;
      ai[4*i+0] += acorn::square(diff)*area;
      ai[4*i+1] += acorn::square(den)*area;
      ai[4*i+2] = std::max(ai[4*i+2], std::abs(diff));
      ai[4*i+3] = std::max(ai[4*i+3], std::abs(den));
    }
  }
  for (int i = 0; i < 6; ++i) {
    Real l2n = 0, l2d = 0, lin = 0, lid = 0;
    for (int tid = 0; tid < nthr; ++tid) {
      auto ai = &a[tid*naccum + 4*i];
      l2n += ai[0];
      l2d += ai[1];
      lin = std::max(lin, ai[2]);
      lid = std::max(lid, ai[3]);
    }
    if (l2d == 0) l2d = 1;
    if (lid == 0) lid = 1;
    l2_err[i] = std::sqrt(l2n/l2d);
    li_err[i] = lin/lid;
  }
}

static void print_errors (const Real l2_err[6], const Real l2ep[6],
                          const Real li_err[6], const Real liep[6],
                          const Real fac = 0) {
  assert(not (l2ep || liep) || fac > 0);
  printf("l2:");
  for (int i = 0; i < 6; ++i) {
    Real ooa = 0;
    if (l2ep && l2_err[i] != 0)
      ooa = -std::log(l2_err[i]/l2ep[i])/std::log(fac);
    printf(" %8.2e (%4.2f)", l2_err[i], ooa);
  }
  printf("\nli:");
  for (int i = 0; i < 6; ++i) {
    Real ooa = 0;
    if (liep && l2_err[i] != 0)
      ooa = -std::log(li_err[i]/liep[i])/std::log(fac);
    printf(" %8.2e (%4.2f)", li_err[i], ooa);
  }
  printf("\n");
}

struct ConvData {
  Discretization::CPtr d;
  int ny;
  RealArray dislocs, sigmas;
  ZxFn::CPtr zxfn;
};

[[maybe_unused]] static void
interp_tt (const ConvData& dref, const ConvData& dmeas, RealArray& si) {
  const int nfn = 6;
  const auto nc = nfn*Discretization::reconstruct_ncoef;
  const auto ntrid = dref.d->get_triangulation()->get_ntri();
  std::vector<Real> coefs(nc*ntrid);
  ompparfor for (int ti = 0; ti < ntrid; ++ti)
    dref.d->tri_reconstruct_fit(ti, nfn, dref.sigmas.data(), &coefs[ti*nc]);
  const auto ntrin = dmeas.d->get_triangulation()->get_ntri();
  si.resize(nfn*ntrin);
  ompparfor for (int ti = 0; ti < ntrin; ++ti) {
    Real p[3];
    dmeas.d->get_surface()->tri_ctr_xyz(ti, p);
    const auto ti_ref = search(*dref.zxfn, *dref.d->get_triangulation(), p);
    dref.d->tri_reconstruct(ti_ref, nfn, &coefs[ti_ref*nc], dref.sigmas.data(),
                            p, &si[nfn*ti]);
  }
}

[[maybe_unused]] static void
interp_to (const ConvTest& ct, const RealArray& tsigmas,
           const ZxFn& ozxfn, const int nyr, RealArray& isigmas) {
  const auto& d = *ct.get_discretization();
  const auto& tzxfn = *ct.get_zxfn();
  const int nfn = 6;
  const auto nc = nfn*Discretization::reconstruct_ncoef;
  const auto ntrid = d.get_triangulation()->get_ntri();
  std::vector<Real> coefs(nc*ntrid);
  ompparfor for (int ti = 0; ti < ntrid; ++ti)
    d.tri_reconstruct_fit(ti, nfn, tsigmas.data(), &coefs[ti*nc]);
  const auto nxr = ozxfn.get_nx();
  const auto nrect = nxr*nyr;
  isigmas.resize(nfn*nrect);
  ompparfor for (int ri = 0; ri < nrect; ++ri) {
    Real p[3];
    ConvTest::calc_rect_ctr(ozxfn, nyr, ri, p);
    const auto ti_ref = search(tzxfn, *d.get_triangulation(), p);
    d.tri_reconstruct(ti_ref, nfn, &coefs[ti_ref*nc], tsigmas.data(), p,
                      &isigmas[nfn*ri]);
  }
}

namespace {

struct ConvTestSettings {
  int testcase = 0;
  bool use_four_tris_per_rect = false;
  bool use_surface_recon = true;
  bool use_flat_elements = false;
  bool use_exact_normals = false;
  int nml_recon_order = 4;
  int disloc_order = 2;
  int e_max = -1;
  Real r = 0.8;
  bool use_woodland_rg0c0 = false;

  void set(const std::string& params);
};

static bool in (const std::string sub, const std::string sup) {
  return sup.find(sub) != std::string::npos;
}

void ConvTestSettings::set (const std::string& params) {
  const auto toks = acorn::split(params, ",");
  for (const auto& t : toks) {
    const auto keyval = acorn::split(t, "=");
    if (keyval.size() != 2) {
      printf("Token error: %s\n", t.c_str());
      exit(-1);
    }
    const auto& key = keyval[0];
    const auto& val = keyval[1];
    const int ival = std::stoi(val);
    if (in("testcase", key)) testcase = ival;
    else if (in("ntri", key)) use_four_tris_per_rect = ival == 4;
    else if (in("srfrecon", key)) use_surface_recon = ival;
    else if (in("flatelem", key)) use_flat_elements = ival;
    else if (in("exactnml", key)) use_exact_normals = ival;
    else if (in("nmlorder", key)) nml_recon_order = ival;
    else if (in("dislocorder", key)) disloc_order = ival;
    else if (in("nres", key)) e_max = ival - 1;
    else if (in("woodlandrect", key)) use_woodland_rg0c0 = ival;
    else {
      printf("Unrecognized key-value pair: %s %s\n",
             keyval[0].c_str(), keyval[1].c_str());
      printf("Valid pairs:\n"
             "  testcase: int [0]\n"
             "  ntri: 2, 4 [2]\n"
             "  srfrecon: bool (0,1) [1]\n"
             "  flatelem: bool [0]\n"
             "  exactnml: bool [0]\n"
             "  nmlorder: 2, 4 [4]\n"
             "  dislocorder: 0, 1, 2, 3 [2]\n"
             "  nres: int > 0 [-1]\n"
             "  woodlandrect: bool [0]\n");
      exit(-1);
    }
  }
#ifndef WOODLAND_ACORN_HAVE_DC3D
  use_woodland_rg0c0 = true;
#endif
}

void set_settings (const ConvTestSettings& s, ConvTest& ct) {
  if (s.use_four_tris_per_rect && s.use_surface_recon && not s.use_flat_elements) {
    printf("ExtrudedCubicSplineSurface cannot be used with ntri=4");
    throw_if(true, "set_settings: incompatible settings");
  }
  ct.set_verbosity(1);
  ct.set_use_four_tris_per_rect(s.use_four_tris_per_rect);
  ct.set_use_surface_recon(s.use_surface_recon);
  ct.set_use_exact_normals(s.use_exact_normals);
  ct.set_use_flat_elements(s.use_flat_elements);
  ct.set_normal_recon_order(s.nml_recon_order);
  ct.set_disloc_order(s.disloc_order);
  ct.set_general_lam_mu();
  ct.set_use_woodland_rg0c0(s.use_woodland_rg0c0);
}

} // namespace

static void conv_test_exact (const bool my_method, const ConvTestSettings& s) {
  ZxFn::Shape zshape;
  Disloc::Ptr disloc;
  setup(s.testcase, zshape, disloc, s.r);
  int e_max = s.e_max;

  ConvTest ct;
  ct.init(zshape, disloc);
  set_settings(s, ct);
  ct.print();

  const int nxfac = 10;
  int ny0 = 0;
  Real l2ep[6], liep[6];
  for (int e = 0; e <= e_max; ++e) {
    const bool last = e == e_max;
    const bool output = last;
    const auto resfac = (1 << e);
    const int nx = nxfac*resfac;

    ct.set_nx(nx);
    if (e == 0) ny0 = ct.get_ny();
    const int ny = ny0*resfac;
    printf("nx ny %d %d\n", nx, ny);
    RealArray dislocs, sigmas, esigmas;
    if (my_method) {
      ct.eval(dislocs, sigmas);
      ct.eval_exact_at_tri_ctrs(esigmas);
      if (0) {
        // Test that e vs e with different configs can get very accurate
        // agreement.
        Exact::Options o;
        o.nxr = 23; o.nyr = 23;
        o.np_radial = 40; o.np_angular = 40;
        o.triquad_order = -1;
        RealArray e1sigmas;
        ct.eval_exact_at_tri_ctrs(e1sigmas, &o);
        Real l2_err[6], li_err[6];
        calc_errors(esigmas, e1sigmas, l2_err, li_err,
                    my_method ? ct.get_discretization()->get_surface() : nullptr);
        printf("e vs e:\n");
        print_errors(l2_err, nullptr, li_err, nullptr, 2);
        exit(0);
      }
    } else {
      int iny = ny;
      ct.eval_okada(nx, iny, dislocs, sigmas);
      assert(iny == ny);
      ct.eval_exact_at_rect_ctrs(nx, ny, esigmas);
    }

    Real l2_err[6], li_err[6];
    calc_errors(esigmas, sigmas, l2_err, li_err,
                my_method ? ct.get_discretization()->get_surface() : nullptr);
    print_errors(l2_err, e == 0 ? nullptr : l2ep,
                 li_err, e == 0 ? nullptr : liep, 2);
    acorn::copy(6, l2_err, l2ep);
    acorn::copy(6, li_err, liep);

    if (output) {
      if (my_method) {
        ct.pywrite("ctzx_tri_exact.py", dislocs, esigmas, "d", false);
        ct.pywrite("ctzx_tri.py", dislocs, sigmas, "d", false);
      } else {
        ZxFn ozxfn(ct.get_zxfn()->get_shape());
        ozxfn.set_nx(nx);
        ct.pywrite_okada("ctzx_rect.py", nx, ny, dislocs, sigmas, "d", false);
        pywrite_zxfn("ctzx_rect.py", ozxfn, "z", true);
        ct.pywrite_okada("ctzx_rect_exact.py", nx, ny, dislocs, esigmas, "d",
                         false);
        pywrite_zxfn("ctzx_rect_exact.py", ozxfn, "z", true);
      }
    }
  }
}

void convtest_w_vs_e (const std::string& params) {
  ConvTestSettings s;
  s.e_max = 3;
  s.set(params);
  printf("convtest_w_vs_e %d\n", s.testcase);
  conv_test_exact(true, s);
}

void convtest_o_vs_e (const std::string& params) {
  ConvTestSettings s;
  s.e_max = 5;
  s.set(params);
  printf("convtest_o_vs_e %d\n", s.testcase);
  conv_test_exact(false, s);
}

void run_case (const std::string& params) {
  //unittest();

  const bool eval_tri = 1;
  const bool eval_oka = 0;
  const bool eval_tri_exact = 1;
  const auto okada_method = ConvTest::EvalMethod::fast;

  ConvTestSettings s;
  s.set(params);
  ZxFn::Shape zshape;
  Disloc::Ptr disloc;
  setup(s.testcase, zshape, disloc, 1.0);

  ConvTest ct;
  ct.init(zshape, disloc);
  set_settings(s, ct);
  ct.set_nx(30);
  printf("nx %d\n", ct.get_nx());
  ct.set_verbosity(1);
  ct.set_hmmvp_tol(1e-4);
  ct.print();

  RealArray tsigmas, tdislocs;
  if (eval_tri) {
    const auto t0 = acorn::dbg::gettime();
    ct.eval(tdislocs, tsigmas);
    const auto t1 = acorn::dbg::gettime();
    printf("ct.eval et %1.3e\n", t1 - t0);
    ct.pywrite("ctzx_tri.py", tdislocs, tsigmas, "d", false);
  } else {
    ct.discretize();
  }
  if (eval_tri_exact) {
    RealArray esigmas;
    const auto t0 = acorn::dbg::gettime();
    ct.eval_exact_at_tri_ctrs(esigmas);
    const auto t1 = acorn::dbg::gettime();
    printf("ct.eval_exact_at_tri_ctrs et %1.3e\n", t1 - t0);
    ct.pywrite("ctzx_tri_exact.py", tdislocs, esigmas, "d", false);
  }
  if (eval_oka) {
#if 1
    // nx has to be large to handle the geometry well.
    const int nx = 11*ct.get_nx();
    int ny = 5*ct.get_ny();
#else
    int nx = 100, ny = 40;
#endif
    pr("rect" pu(nx) pu(ny));
    RealArray idislocs, isigmas;
    const auto t0 = acorn::dbg::gettime();
    switch (okada_method) {
    case ConvTest::EvalMethod::direct:
    case ConvTest::EvalMethod::direct_hmmvp: {
      RealArray dislocs, sigmas;
      ct.eval_okada(nx, ny, dislocs, sigmas, okada_method);
      ct.pywrite_okada("ctzx_rect.py", nx, ny, dislocs, sigmas, "d", false);
      ct.interp_okada(nx, ny, sigmas, isigmas);
    } break;
    case ConvTest::EvalMethod::fast:
    case ConvTest::EvalMethod::fast_hmmvp: {
      ConvTest::SupportPoints supports;
      ct.collect_rect_support_points(nx, ny, supports);
      RealArray sigmas;
      ct.eval_okada(nx, ny, supports, sigmas, ConvTest::is_hmmvp(okada_method));
      ct.interp_okada(nx, ny, supports, sigmas, isigmas);
    } break;
    default: assert(0);
    }
    const auto t1 = acorn::dbg::gettime();
    printf("ct.eval/interp_okada et %1.3e\n", t1 - t0);
    ct.pywrite("ctzx_tri_interp.py", idislocs, isigmas, "d", false);
    if (not (tsigmas.empty() || isigmas.empty())) {
      Real l2_err[6], li_err[6];
      calc_errors(isigmas, tsigmas, l2_err, li_err,
                  ct.get_discretization()->get_surface());
      print_errors(l2_err, nullptr, li_err, nullptr);
    }
  }
}

static int test_arclength () {
  int nerr = 0;
  ZxFn zxfn(ZxFn::Shape::trig1);
  const auto shape = zxfn.get_shape();
  const Real xend = 0.7;
  const auto L = zxfn.calc_arclength(xend);
  const int n = 100;
  Real Lsum = 0, xprev = 0, fprev = 0;
  for (int i = 0; i <= n; ++i) {
    const Real x = xend*Real(i)/n;
    Real f, g;
    eval(shape, x, f, g);
    if (i > 0)
      Lsum += std::sqrt(acorn::square(x - xprev) + acorn::square(f - fprev));
    xprev = x;
    fprev = f;
  }
  if (std::abs(L - Lsum) > 1e-4*Lsum) ++nerr;
  if (nerr) printf("convtest_zx: test_arclength failed\n");
  return nerr;
}

static int test_zxfn () {
  int nerr = 0;
  ZxFn zxfn(ZxFn::Shape::trig1);
  const int nx = 274;
  zxfn.set_nx(nx);
  const auto xs = zxfn.get_xbs();
  const auto zs = zxfn.get_zbs();
  Real segmin = 1, segmax = 0;
  for (int i = 0; i < nx; ++i) {
    const Real seg = std::sqrt(acorn::square(xs[i+1] - xs[i]) +
                               acorn::square(zs[i+1] - zs[i]));
    segmin = std::min(segmin, seg);
    segmax = std::max(segmax, seg);
  }
  if (segmax - segmin > 1e-3*segmax) ++nerr;
  if (nerr) printf("convtest_zx: test_zxfn failed\n");
  return nerr;  
}

static Real testfn(const int fno, const Real x, const Real y) {
  switch (fno) {
  case 0:
  default:
    return std::cos(2*M_PI*x)*std::sin(2*M_PI*y);
  }
}

int ConvTest::test_interp () {
  const bool verbose = false;
  int nerr = 0;

  ZxFn zxfn(ZxFn::Shape::trig1), zxfnr(zxfn.get_shape());
  const auto disloc = std::make_shared<Disloc>();
  ConvTest ct;
  ct.init(zxfn.get_shape(), disloc);
  const int fno = 0;

  Real pl2e = 0, plie = 0;
  for (int e = 0; e < 6; ++e) {
    const int nx = 5*(1 << e);
    zxfn.set_nx(nx);
    const int ny = zxfn.get_ny();
    const int onx = std::ceil(3*nx), ony = std::ceil(1.6*ny);
    ct.set_nx(nx);
    ct.discretize();

    RealArray tfs; {
      RealArray rfs(6*ony*onx, 0);
      zxfnr.set_nx(onx);
      const auto rxs = zxfnr.get_xbs();
      for (int iy = 0, k = 0; iy < ony; ++iy) {
        const auto y = (iy + 0.5)/ony;
        for (int ix = 0; ix < onx; ++ix, ++k) {
          const auto x = (rxs[ix] + rxs[ix+1])/2;
          rfs[6*k] = testfn(fno, x, y);
        }
      }
      ct.interp_okada(onx, ony, rfs, tfs);
    }

    const auto& t = *ct.d->get_triangulation();
    const auto& srf = *ct.d->get_surface();
    const auto ntri = t.get_ntri();
    RealArray tfs_true(ntri);
    for (int ti = 0; ti < ntri; ++ti) {
      Real p[3];
      srf.tri_ctr_xyz(ti, p);
      tfs_true[ti] = testfn(fno, p[0], p[1]);
    }

    Real l2n = 0, l2d = 0, lin = 0, lid = 0;
    for (int ti = 0; ti < ntri; ++ti) {
      const auto d = tfs[6*ti] - tfs_true[ti];
      l2n += acorn::square(d);
      l2d += acorn::square(tfs_true[ti]);
      lin = std::max(lin, std::abs(d));
      lid = std::max(lid, std::abs(tfs_true[ti]));
    }
    const auto l2e = std::sqrt(l2n/l2d), lie = lin/lid;
    const auto l2_ooa = e == 0 ? 0 : std::log(pl2e/l2e)/std::log(2);
    const auto li_ooa = e == 0 ? 0 : std::log(plie/lie)/std::log(2);

    if (e > 4) {
      if (l2_ooa < 1.98) ++nerr;
      if (li_ooa < 1.95) ++nerr;
    }
    if (verbose)
      printf("%2d nx %3d ny %3d l2 %9.3e (%5.3f) linf %9.3e (%5.3f)\n",
             e, nx, ny, l2e, l2_ooa, lie, li_ooa);

    pl2e = l2e;
    plie = lie;
  }
  if (nerr) printf("convtest_zx: test_interp failed\n");
  return nerr;
}

static int test_fast_okada_methods (const bool use_own_impl) {
  int nerr = 0;

  ZxFn::Shape zshape;
  Disloc::Ptr disloc;
  setup(0, zshape, disloc);

  ConvTest ct;
  ct.init(zshape, disloc);
  ct.set_verbosity(0);
  const int nx = use_own_impl ? 7 : 11;
  ct.set_nx(nx);
  ct.discretize();
  const int ny = ct.get_zxfn()->get_ny();
  ct.set_use_woodland_rg0c0(use_own_impl);

  const int onx = 4*nx;
  int ony = ny + 9;

  RealArray s2, is2; {
    RealArray dislocs, s1;
    ct.eval_okada(onx, ony, dislocs, s1, ConvTest::EvalMethod::direct);
    ct.eval_okada(onx, ony, dislocs, s2);
    const int n = int(s1.size());
    for (int i = 0; i < n; ++i) {
      const auto e = std::abs(s2[i] - s1[i]);
      if (e > (use_own_impl ? 1e-7 : 1e-9) ||
          e > (use_own_impl ? 1e-5 : 1e-6)*std::abs(s1[i])) {
        if (nerr < 10) pr(puf(i) pu(s1[i]) pu(s2[i]) pu(s2[i] - s1[i]));
        ++nerr;
      }
    }
    ct.interp_okada(onx, ony, s2, is2);
  }

  {
    ConvTest::SupportPoints supports;
    ct.collect_rect_support_points(onx, ony, supports);
    RealArray s3, is3;
    ct.eval_okada(onx, ony, supports, s3);
    const int ns = int(supports.size());
    for (int i = 0; i < ns; ++i) {
      const auto& sp = supports[i];
      const int idx = sp.iy*onx + sp.ix;
      for (int k = 0; k < 6; ++k) {
        const int i3 = 6*i + k, i2 = 6*idx + k;
        if (s3[i3] != s2[i2]) {
          if (nerr < 10)
            pr(puf(i) pu(idx) pu(k) pu(s2[i2]) pu(s3[i3]) pu(s3[i3] - s2[i2]));
          ++nerr;
        }
      }
    }
    ct.interp_okada(onx, ony, supports, s3, is3);
    const int n = int(is2.size());
    for (int i = 0; i < n; ++i)
      if (is3[i] != is2[i]) {
        if (nerr < 10) pr(puf(i) pu(is2[i]) pu(is3[i]) pu(is3[i] - is2[i]));
        ++nerr;
      }
  }

  if (nerr) printf("convtest_zx: test_fast_okada_methods failed\n");
  return nerr;
}

static int test_fast_woodland_methods () {
  int nerr = 0;

  ZxFn::Shape zshape;
  Disloc::Ptr disloc;
  setup(0, zshape, disloc);

  ConvTest ct;
  ct.init(zshape, disloc);
  ct.set_verbosity(0);
  ct.set_use_surface_recon(0);
  ct.set_use_exact_normals(1);
  ct.set_disloc_order(2);
  printf("test_fast_woodland_methods: skipping order 3 b/c too expensive\n");

  for (const bool use_four : {false, true}) {
    ct.set_use_four_tris_per_rect(use_four);
    ct.set_nx(use_four ? 5 : 6);

    RealArray dd, sd, df, sf, se;
    const auto t0 = acorn::dbg::gettime();
    ct.eval(dd, sd, ConvTest::EvalMethod::direct);
    const auto t1 = acorn::dbg::gettime();
    ct.eval(df, sf, ConvTest::EvalMethod::fast);
    const auto t2 = acorn::dbg::gettime();
    if (0) {
      // We don't expect speedup for a test-problem size. Each disloc component
      // requires 2*6 GF evaluations (2 for each side of the source, 6 for the 6
      // components of the quadratic fit), so a fully general problem requires
      // 2*6*3 = 36 GF evaluations. Thus, breakeven is roughly nyrect = 36. For a
      // single disloc component, breakeven is 12.
      printf("eval_direct/fast et %1.3e %1.3e speedup %1.2f\n",
             t1 - t0, t2 - t1, (t1 - t0)/(t2 - t1));
    }
    ct.pywrite("ctzx_tri.py", df, sf, "d", false);
    const auto ntri = ct.get_discretization()->get_triangulation()->get_ntri();
    for (int i = 0; i < 6*ntri; ++i) {
      Real e = std::abs(sf[i] - sd[i]);
      if (e > 1e-8 || e > 1e-6*std::abs(sd[i])) {
        if (nerr < 10) pr(puf(i) pu(sd[i]) pu(sf[i]) pu(sf[i] - sd[i]));
        ++nerr;
      }
    }
  }

  if (nerr) printf("convtest_zx: test_fast_woodland_methods failed\n");
  return nerr;
}

int ConvTest::unittest () {
  const auto t0 = acorn::dbg::gettime();
  int nerr = 0;
  nerr += test_arclength();
  nerr += test_zxfn();
  nerr += ExtrudedCubicSplineSurface::unittest();
  nerr += FlatElementSurface::unittest();
  nerr += ConvTest::test_interp();
  nerr += test_fast_okada_methods(false);
  nerr += test_fast_okada_methods(true);
  nerr += test_fast_woodland_methods();
  const auto t1 = acorn::dbg::gettime();
  printf("convzx::unittest et %1.2e\n", t1 - t0);
  return nerr;
}

} // namespace convzx
} // namespace examples
} // namespace woodland
