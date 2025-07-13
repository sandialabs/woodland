#include "woodland/acorn/compose_quadrature.hpp"
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/triangle.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/dbg.hpp"

#include "woodland/squirrel/mesh.hpp"
#include "woodland/squirrel/global_z_surface.hpp"
#include "woodland/squirrel/ctzx.hpp"
#include "woodland/squirrel/pywrite.hpp"
#include "woodland/squirrel/mesh.hpp"

#include <cmath>
#include <set>
#include <map>

namespace woodland {
namespace squirrel {
namespace ctzx {

typedef acorn::Matvec<2,Real> mv2;
typedef acorn::Matvec<3,Real> mv3;

namespace {

static const bool use_calc_integral_tensor_quadrature = true;

void pywrite (FILE* fp, const std::string& dict, const Triangulation& t) {
  fprintf(fp, "%s['topo'] = npy.array([\n", dict.c_str());
  const auto ntri = t.get_ncell();
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

void pywrite (FILE* fp, const std::string& dict, const mesh::Mesh& m) {
  fprintf(fp, "%s['topo'] = npy.array([\n", dict.c_str());
  const auto nc = m.get_ncell();
  const auto to = *m.get_topo();
  const auto v = *m.get_vtxs();
  int nv = 0;
  for (int ci = 0; ci < nc; ++ci) nv = std::max(nv, to.get_cell_nvtx(ci));
  for (int ci = 0; ci < nc; ++ci) {
    const auto c = to.get_cell(ci);
    fprintf(fp, "[");
    for (int i = 0; i < c.n; ++i) fprintf(fp, "%d,", int(c[i]));
    for (int i = c.n; i < nv; ++i) fprintf(fp, "%d,", int(c[c.n-1]));
    fprintf(fp, "],\n");
  }
  fprintf(fp, "], dtype=npy.int32)\n");
  fprintf(fp, "%s['vtx'] = npy.array([\n", dict.c_str());
  const auto nvtx = v.get_nvtx();
  for (int vi = 0; vi < nvtx; ++vi) {
    const auto vtx = v.get_vtx(vi);
    fprintf(fp, "[%15.8e, %15.8e, %15.8e],\n", vtx[0], vtx[1], vtx[2]);
  }
  fprintf(fp, "])\n");
}

} // namespace

void ConvTest::init (const ZxFn::Shape shape, const Disloc::CPtr& disloc_) {
  zxfn = std::make_shared<ZxFn>(shape);
  disloc = disloc_;
  gfp.lam = gfp.mu = 1;
  gfp.halfspace = false;
#ifndef WOODLAND_ACORN_HAVE_DC3D
  use_woodland_rg0c0 = true;
#endif
}

void ConvTest::set_nx (const int nx) {
  assert(nx >= 1);
  if (zxfn->get_nx() != nx or zxfn->is_uniform_in_x() != xuniform) {
    t = nullptr; d = nullptr;
  }
  zxfn->set_nx(nx, xuniform);
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
  if (use != dargs.use_surface_recon) { t = nullptr; d = nullptr; }
  dargs.use_surface_recon = use;
}

void ConvTest::set_use_exact_tangents (const bool use) {
  if (use != dargs.use_exact_tangents) { t = nullptr; d = nullptr; }
  dargs.use_exact_tangents = use;
}

void ConvTest::set_use_c2_spline (const bool use) {
  if (use != dargs.use_c2_spline) { t = nullptr; d = nullptr; }
  dargs.use_c2_spline = use;
}

void ConvTest::set_use_halfspace (const bool use) {
  gfp.halfspace = use;
}

void ConvTest::set_tangent_recon_order (const int order) {
  if (order != dargs.tan_recon_order) { t = nullptr; d = nullptr; }
  dargs.tan_recon_order = order;
}

void ConvTest::set_disloc_order (const int order) {
  if (order != disloc_order) { t = nullptr; d = nullptr; }
  disloc_order = order;
}

void ConvTest::set_use_flat_elements (const bool use) {
  if (use != dargs.use_flat_elements) { t = nullptr; d = nullptr; }
  dargs.use_flat_elements = use;
}

void ConvTest::set_use_woodland_rg0c0 (const bool use) {
#ifndef WOODLAND_ACORN_HAVE_DC3D
  if (not use) {
    printf("WARNING: Using Woodland rg0c0 impl b/c dc3d impl is missing.\n");
    use_woodland_rg0c0 = true;
    return;
  }
#endif
  use_woodland_rg0c0 = use;
}

void ConvTest::set_use_nonunirect (const bool use) {
  if (use != use_nonunirect) { t = nullptr; d = nullptr; }
  use_nonunirect = use;
}

void ConvTest::set_xuniform (const bool use) {
  if (use != xuniform) { t = nullptr; d = nullptr; }
  xuniform = use;
  if (zxfn->get_nx() > 0) set_nx(zxfn->get_nx());
}

void ConvTest::set_support_req0 (const bool use) {
  if (use != dargs.support_req0) { t = nullptr; d = nullptr; }
  dargs.support_req0 = use;
}

void ConvTest::set_mesh_scale (const bool use) {
  if (use != dargs.mesh_scale) { t = nullptr; d = nullptr; }
  dargs.mesh_scale = use;
}

void ConvTest::print (FILE* fp) const {
  fprintf(fp, "ct> lam %1.3e mu %1.3e halfspace %d\n",
          gfp.lam, gfp.mu, int(gfp.halfspace));
  fprintf(fp, "ct> xuniform %d nonunirect %d ntriperrect %d dislocorder %d\n"
          "ct> exactsrf %d flatelem %d exacttan %d tanorder %d c2spline %d\n",
          int(xuniform), use_nonunirect, ntri_per_rect, disloc_order,
          int(not dargs.use_surface_recon), int(dargs.use_flat_elements),
          int(dargs.use_exact_tangents), dargs.tan_recon_order,
          int(dargs.use_c2_spline));
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
  assert(zxfn.get_nx() >= 1 and zxfn.get_ny() >= 1);

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
struct Param2DUserFn : public GlobalZSurface::UserFn {
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
::discretize (const mesh::Mesh::CPtr& m, const Triangulation::Ptr& t,
              const ZxFn::Shape shape, const DiscretizeArgs a) {
  Surface::CPtr s;
  if (a.use_surface_recon) {
    assert(t);
    if (a.use_flat_elements) {
      const Real primary[] = {1,0,0};
      s = std::make_shared<FlatElementSurface>(t, primary);
    } else {
      using Recon = ExtrudedCubicSplineSurface::Recon;
      const auto recon = (a.use_c2_spline ? Recon::c2 :
                          a.use_exact_tangents ? Recon::c1_tan_exact :
                          a.tan_recon_order == 2 ? Recon::c1_tan_2 :
                          Recon::c1_tan_4);
      assert(recon != Recon::c1_tan_4 || a.tan_recon_order == 4);
      s = std::make_shared<ExtrudedCubicSplineSurface>(shape, t, recon);
    }
  } else {
    const auto ufn = std::make_shared<Param2DUserFn>(shape);
    const auto srf = std::make_shared<GlobalZSurface>(
      m, ufn, a.support_req0, a.mesh_scale);
    s = srf;
  }
  return make_discretization(m, s);
}

void ConvTest
::pywrite (const std::string& python_filename, const RealArray& dislocs,
           const RealArray& sigmas, const std::string& dict,
           const bool append) const {
  const auto& mesh = *d->get_mesh();
  assert(dislocs.empty() || 3*mesh.get_ncell() == int(dislocs.size()));
  assert(sigmas .empty() || 6*mesh.get_ncell() == int(sigmas .size()));
  assert(not (dislocs.empty() && sigmas.empty()));
  FILE* fp = fopen(python_filename.c_str(), append ? "a" : "w");
  if (not append) pywrite_header(fp);
  fprintf(fp, "%s = {}\n", dict.c_str());
  fprintf(fp, "%s['nx'] = %d\n", dict.c_str(), zxfn->get_nx());
  if (t) {
    fprintf(fp, "%s['type'] = 'mesh'\n", dict.c_str());
    std::string tdi = dict + "['m']";
    fprintf(fp, "%s = {}\n", tdi.c_str());
    ctzx::pywrite(fp, tdi, *t);
  } else {
    fprintf(fp, "%s['type'] = 'mesh'\n", dict.c_str());
    std::string tdi = dict + "['m']";
    fprintf(fp, "%s = {}\n", tdi.c_str());
    ctzx::pywrite(fp, tdi, mesh);
  }
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
  pywrite_double_array(fp, ddi, 1, zxfn.get_nx()+1, zxfn.get_xbs());
  ddi = dict + "['z']";
  pywrite_double_array(fp, ddi, 1, zxfn.get_nx()+1, zxfn.get_zbs());
  fclose(fp);
}

// Refine alternating y-strips by one level to test convergence when using a
// nonuniform qmesh.
//   Cells are grouped by y-strip.
mesh::Mesh::CPtr ConvTest::make_nonunirect_mesh (const ZxFn& zxfn) {
  assert(zxfn.get_nx() >= 1 && zxfn.get_ny() >= 1);

  const bool uniform = true; // optionally uniform
  if (uniform)
    printf("NOTE, make_nonunirect_mesh: logically rectangular grid.\n");
  const auto nx = zxfn.get_nx(), ny = zxfn.get_ny();

  const auto v = mesh::make_vertices(3);
  const auto to = mesh::make_topology();
  mesh::BoolArray v2bdy;

  // Build topo and vtxs.
  {
    // Given (x,y), compute z and add the vertex (x,y,z) to the vertices.
    const auto zshape = zxfn.get_shape();
    const auto xs = zxfn.get_xbs();
    const auto addvtx =
      [&] (const Real x, const Real y, const bool on_bdy) -> Idx {
        Real z, unused;
        gallery::eval(zshape, x, z, unused);
        const Real p[] = {x, y, z};
        v2bdy.push_back(on_bdy);
        return v->append_vtx(p);
      };
    // Given grid integer coordinates (i,j), return the associated vertex.
    using std::make_pair;
    std::map<std::pair<Idx,Idx>,Idx> c2vi;
    const auto getvi = [&] (const Idx i, const Idx j) -> Idx {
      const auto p = make_pair(i,j);
      const auto it = c2vi.find(p);
      assert(it != c2vi.end());
      return it->second;
    };
    // Build grid.
    Idx cell[4];
    for (Idx xi = 0; xi < nx; ++xi) {
      for (Idx yi = 0; yi < ny; ++yi) {
        if (xi == 0) {
          if (yi == 0) c2vi[make_pair(0,0)] = addvtx(xs[0], 0, true);
          c2vi[make_pair(0,2*(yi+1))] = addvtx(xs[0], Real(yi+1)/ny, true);
        }
        if (yi == 0) c2vi[make_pair(2*(xi+1),0)] = addvtx(xs[xi+1], 0, true);
        c2vi[make_pair(2*(xi+1), 2*(yi+1))] =
          addvtx(xs[xi+1], Real(yi+1)/ny, xi+1 == nx or yi+1 == ny);
      }
      if (xi % 2 == 0 or uniform) {
        for (Idx yi = 0; yi < ny; ++yi) {
          cell[0] = getvi(2*xi, 2*yi);
          cell[1] = getvi(2*(xi+1), 2*yi);
          cell[2] = getvi(2*(xi+1), 2*(yi+1));
          cell[3] = getvi(2*xi, 2*(yi+1));
          to->append_cell(mesh::Cell(4, cell));
        }
      } else {
        const Real xl = xs[xi], xr = xs[xi+1], xm = (xl + xr)/2;
        for (Idx yi = 0; yi < ny; ++yi) {
          const auto ym = Real(2*yi+1)/(2*ny);
          c2vi[make_pair(2*xi+1, 2*yi)] = addvtx(xm, Real(yi)/ny, yi == 0);
          c2vi[make_pair(2*xi, 2*yi+1)] = addvtx(xl, ym, false);
          c2vi[make_pair(2*xi+1, 2*yi+1)] = addvtx(xm, ym, false);
          c2vi[make_pair(2*(xi+1), 2*yi+1)] = addvtx(xr, ym, false);
        }
        c2vi[make_pair(2*xi+1, 2*ny)] = addvtx(xm, 1, true);
        for (int k = 0; k < 2; ++k) {
          for (Idx yi = 0; yi < 2*ny; ++yi) {
            cell[0] = getvi(2*xi+k, yi);
            cell[1] = getvi(2*xi+k+1, yi);
            cell[2] = getvi(2*xi+k+1, yi+1);
            cell[3] = getvi(2*xi+k, yi+1);
            to->append_cell(mesh::Cell(4, cell));
          }
        }
      }
    }
    to->end_modifications();
  }

  const auto ncell = to->get_ncell(), nvtx = v->get_nvtx();

  // Assert various entity counts of the grid.
  assert(to->get_nvtx() == nvtx);
  assert(v2bdy.size() == static_cast<size_t>(nvtx));
#ifndef NDEBUG
  int nbdy = 0;
  for (const auto& e : v2bdy)
    if (e) ++nbdy;
#endif
  if (uniform) {
    assert(ncell == nx*ny);
    assert(nvtx == (nx+1)*(ny+1));
    assert(nbdy == 2*(nx + ny));
  } else {
    assert(ncell == (nx+1)/2*ny + 4*(nx/2)*ny);
    assert(nvtx == (nx+1)*(ny+1) + 4*(nx/2)*ny + nx/2);
    assert(nbdy == 2*(nx + ny) + 2*(nx/2));
  }

  // Make neighbor sets. In this grid, every nbr is a vtx nbr.
  std::vector<std::set<Idx>> nbrss(ncell); {
    std::vector<std::set<Idx>> vi2ci(nvtx);
    for (Idx ci = 0; ci < ncell; ++ci) {
      const auto c = to->get_cell(ci);
      for (int i = 0; i < c.n; ++i)
        vi2ci[c[i]].insert(ci);
    }
    for (Idx ci = 0; ci < ncell; ++ci) {
      const auto c = to->get_cell(ci);
      for (int i = 0; i < c.n; ++i) {
        const auto& cjs = vi2ci[c[i]];
        for (const auto& cj : cjs)
          if (cj != ci) nbrss[ci].insert(cj);
      }
    }
  }
  // Convert nbr sets to arrays.
  mesh::IdxArray c2csi(ncell+1), c2cs;
  c2csi[0] = 0;
  for (Idx ci = 0; ci < ncell; ++ci) {
    const auto& nbrs = nbrss[ci];
    const Idx nnbr = static_cast<Idx>(nbrs.size());
    c2csi[ci+1] = c2csi[ci] + nnbr;
    for (const auto& cj : nbrs)
      c2cs.push_back(cj);
  }

  const auto re = mesh::make_relations(to, c2csi, c2cs, v2bdy);
  return mesh::make_mesh(v, to, re);
}

void ConvTest::discretize () {
  mesh::Mesh::CPtr m;
  if (not use_nonunirect) {
    t = triangulate(*zxfn, ntri_per_rect);
    const auto tr = std::make_shared<TriangulationRelations>(t);
    m = tri2mesh(*t, tr);
  } else {
    m = make_nonunirect_mesh(*zxfn);
    throw_if(dargs.use_surface_recon,
             "nonunirect=1 does not support srfrecon=1");
  }
  d = discretize(m, t, zxfn->get_shape(), dargs);
  d->set_disloc_order(disloc_order);
#ifdef WOODLAND_HAVE_HMMVP
  hmats = hmat::make_hmatrices(d);
#endif
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
  case 40: case 41: case 42:
    zshape = ZxFn::Shape::steep;
    disloc->set(testcase - 40, dshape, 1, 0.5, 0.5, 1, 0, r, r);
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
  Workspace w;
  ompparfor for (int ti = 0; ti < ncell; ++ti) {
    const auto tid = acorn::get_thread_num();
    auto ai = &a[tid*naccum];
    const auto area = srf ? cell_surface_area(w, *srf, ti) : 1;
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

namespace {

struct ConvTestSettings {
  int testcase = 0;
  bool use_four_tris_per_rect = false;
  bool use_surface_recon = true;
  bool use_flat_elements = false;
  bool use_exact_tangents = false;
  bool use_c2_spline = false;
  bool use_woodland_rg0c0 = false;
  bool use_halfspace = false;
  bool use_nonunirect = false;
  bool use_hmmvp = false;
  bool xuniform = false;
  int tan_recon_order = 4;
  int disloc_order = 2;
  int e_max = -1;
  Real r = 0.8;
  int ny_override = -1;

  void set(const std::string& params);
};

static bool in (const std::string sub, const std::string sup) {
  return sup.find(sub) != std::string::npos;
}

static void print_valid_pairs () {
  printf("Valid pairs:\n"
         "  testcase: int [0]\n"
         "  ntri: 2, 4 [2]\n"
         "  srfrecon: bool (0,1) [1]\n"
         "  flatelem: bool [0]\n"
         "  exacttan: bool [0]\n"
         "  tanorder: 2, 4 [4]\n"
         "  c2spline: bool [0]\n"
         "  dislocorder: 0, 1, 2, 3 [2]\n"
         "  nres: int > 0 [-1]\n"
         "  woodlandrect: bool [0]\n"
         "  halfspace: bool [0]\n"
         "  nonunirect: bool [0]\n"
         "  hmmvp: bool [0]\n"
         "  xuniform: bool [0]\n");
}

void ConvTestSettings::set (const std::string& params) {
  const auto toks = acorn::split(params, ",");
  for (const auto& t : toks) {
    const auto keyval = acorn::split(t, "=");
    if (keyval.size() != 2) {
      printf("Token error: %s\n", t.c_str());
      print_valid_pairs();
      exit(-1);
    }
    const auto& key = keyval[0];
    const auto& val = keyval[1];
    const int ival = std::stoi(val);
    if (in("testcase", key)) testcase = ival;
    else if (in("ntri", key)) use_four_tris_per_rect = ival == 4;
    else if (in("srfrecon", key)) use_surface_recon = ival;
    else if (in("flatelem", key)) use_flat_elements = ival;
    else if (in("exacttan", key)) use_exact_tangents = ival;
    else if (in("tanorder", key)) tan_recon_order = ival;
    else if (in("c2spline", key)) use_c2_spline = ival;
    else if (in("dislocorder", key)) disloc_order = ival;
    else if (in("nres", key)) e_max = ival - 1;
    else if (in("woodlandrect", key)) use_woodland_rg0c0 = ival;
    else if (in("halfspace", key)) use_halfspace = ival;
    else if (in("nonunirect", key)) use_nonunirect = ival;
    else if (in("hmmvp", key)) use_hmmvp = ival;
    else if (in("xuniform", key)) xuniform = ival;
    // Unadvertised debug option with unspecified behavior.
    else if (in("nyoverride", key)) ny_override = ival;
    else {
      printf("Unrecognized key-value pair: %s %s\n",
             keyval[0].c_str(), keyval[1].c_str());
      print_valid_pairs();
      exit(-1);
    }
  }
#ifndef WOODLAND_ACORN_HAVE_DC3D
  use_woodland_rg0c0 = true;
#endif
}

void set_settings (const ConvTestSettings& s, ConvTest& ct) {
  if (s.use_four_tris_per_rect and s.use_surface_recon
      and not s.use_flat_elements) {
    printf("ExtrudedCubicSplineSurface cannot be used with ntri=4");
    throw_if(true, "set_settings: incompatible settings");
  }
  ct.set_verbosity(1);
  ct.set_use_four_tris_per_rect(s.use_four_tris_per_rect);
  ct.set_use_surface_recon(s.use_surface_recon);
  ct.set_use_exact_tangents(s.use_exact_tangents);
  ct.set_use_c2_spline(s.use_c2_spline);
  ct.set_use_flat_elements(s.use_flat_elements);
  ct.set_use_woodland_rg0c0(s.use_woodland_rg0c0);
  ct.set_use_halfspace(s.use_halfspace);
  ct.set_use_nonunirect(s.use_nonunirect);
  ct.set_xuniform(s.xuniform);
  ct.set_tangent_recon_order(s.tan_recon_order);
  ct.set_disloc_order(s.disloc_order);
  ct.set_general_lam_mu();
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

  const auto method = (s.use_hmmvp ?
                       ConvTest::EvalMethod::direct_hmmvp :
                       ConvTest::EvalMethod::fast);
  if (s.use_hmmvp) printf("ct> use_hmmvp\n");

  const int nxfac = 10;
  int ny0 = 0;
  Real l2ep[6], liep[6];
  Real t_accum[2] = {0};
  for (int e = 0; e <= e_max; ++e) {
    const bool last = e == e_max;
    const bool output = last;
    const auto resfac = (1 << e);
    const int nx = nxfac*resfac;

    ct.set_nx(nx);
    if (e == 0) ny0 = ct.get_ny();
    const int ny = s.ny_override > 0 ? s.ny_override : ny0*resfac;
    printf("nx ny %d %d\n", nx, ny);
    RealArray dislocs, sigmas, esigmas;
    Real t0, t1, t2;
    if (my_method) {
      t0 = acorn::dbg::gettime();
      ct.eval(dislocs, sigmas, method);
      t1 = acorn::dbg::gettime();
      ct.eval_exact_at_cell_ctrs(esigmas);
      t2 = acorn::dbg::gettime();
    } else {
      int iny = ny;
      t0 = acorn::dbg::gettime();
      ct.eval_okada(nx, iny, dislocs, sigmas);
      t1 = acorn::dbg::gettime();
      assert(iny == ny);
      ct.eval_exact_at_rect_ctrs(nx, ny, esigmas);
      t2 = acorn::dbg::gettime();
    }

    t_accum[0] += t1 - t0;
    t_accum[1] += t2 - t1;
    printf("t: method %1.2e (%1.2e) exact %1.2e (%1.2e)\n",
           t1 - t0, t_accum[0], t2 - t1, t_accum[1]);

    Real l2_err[6], li_err[6];
    calc_errors(esigmas, sigmas, l2_err, li_err,
                my_method ? ct.get_discretization()->get_surface() : nullptr);
    print_errors(l2_err, e == 0 ? nullptr : l2ep,
                 li_err, e == 0 ? nullptr : liep, 2);
    acorn::copy(6, l2_err, l2ep);
    acorn::copy(6, li_err, liep);

    if (output) {
      if (my_method) {
        ct.pywrite("ctzx_wl_exact.py", dislocs, esigmas, "d", false);
        ct.pywrite("ctzx_wl.py", dislocs, sigmas, "d", false);
      } else {
        ZxFn ozxfn(ct.get_zxfn()->get_shape());
        ozxfn.set_nx(nx);
        ct.pywrite_okada("ctzx_rect.py", nx, ny, dislocs, sigmas, "d", false);
        pywrite_zxfn("ctzx_rect.py", ozxfn, "z", true);
        ct.pywrite_okada("ctzx_rect_exact.py", nx, ny, dislocs, esigmas, "d", false);
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
    ct.eval_exact_at_cell_ctrs(esigmas);
    const auto t1 = acorn::dbg::gettime();
    printf("ct.eval_exact_at_cell_ctrs et %1.3e\n", t1 - t0);
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
    case ConvTest::EvalMethod::direct: {
      RealArray dislocs, sigmas;
      ct.eval_okada(nx, ny, dislocs, sigmas, okada_method);
      ct.pywrite_okada("ctzx_rect.py", nx, ny, dislocs, sigmas, "d", false);
      ct.interp_okada(nx, ny, sigmas, isigmas);
    } break;
    case ConvTest::EvalMethod::fast: {
      ConvTest::SupportPoints supports;
      ct.collect_rect_support_points(nx, ny, supports);
      RealArray sigmas;
      ct.eval_okada(nx, ny, supports, sigmas);
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

static Real testfn(const int fno, const Real x, const Real y) {
  switch (fno) {
  case 0:
  default:
    return std::cos(2*M_PI*x)*std::sin(2*M_PI*y);
  }
}

void set_standard_options (Stress::Options& o, Stress::Options& o_nbr) {
  // Use nbr relations rather than geometry to determine when to call calc_hfp
  // on other-interactions.
  o.near_dist_fac = 0;
  o_nbr = o;
  if (use_calc_integral_tensor_quadrature) {
    o_nbr.qp_tq.np_radial = 40;
    o_nbr.qp_tq.np_angular = 40;
    o_nbr.force_tensor_quadrature = true;
  } else {
    o_nbr.force_hfp = true;
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

    const auto& m = *ct.d->get_mesh();
    const auto& srf = *ct.d->get_surface();
    const auto ntri = m.get_ncell();
    RealArray tfs_true(ntri);
    for (int ti = 0; ti < ntri; ++ti) {
      Real p[3];
      srf.cell_ctr_xyz(ti, p);
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
  if (nerr) printf("ctzx: test_interp failed\n");
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

  if (nerr) printf("ctzx: test_fast_okada_methods failed\n");
  return nerr;
}

static int
check_errors (const RealArray& sd, const RealArray& sf,
              const Real abs_tol = 2e-8, const Real rel_tol = 1e-6,
              const int display_cnt = 100, const bool verbose = false) {
  int nerr = 0;
  const auto nc = sd.size() / 6;
  Real max_abs = 0;
  for (size_t i = 0; i < 6*nc; ++i) {
    const Real e = std::abs(sf[i] - sd[i]);
    max_abs = std::max(max_abs, e);
    if (std::isnan(sf[i]) or std::isnan(sd[i]) or e > abs_tol or
        e > rel_tol*std::abs(sd[i])) {
      if (nerr < display_cnt)
        pr(puf(i/6) pu(i % 6) pu(sd[i]) pu(sf[i]) pu(sf[i] - sd[i]));
      ++nerr;
    }
  }
  if (nerr) printf("check_errors max_abs %1.3e\n", max_abs);
  return nerr;
}

static int test_fast_woodland_methods (const bool use_nonunirect = false) {
  int nerr = 0;

  ZxFn::Shape zshape;
  Disloc::Ptr disloc;
  setup(0, zshape, disloc);

  ConvTest ct;
  ct.init(zshape, disloc);
  ct.set_verbosity(0);
  ct.set_use_surface_recon(0);
  ct.set_use_exact_tangents(1);
  ct.set_disloc_order(2);
  ct.set_use_nonunirect(use_nonunirect);

  for (const bool use_four : {false, true}) {
    if (use_nonunirect and use_four) continue;
    ct.set_use_four_tris_per_rect(use_four);
    ct.set_nx(use_four ? 5 : 6);
    ct.set_use_halfspace(not use_four);

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
    ct.pywrite(use_nonunirect ? "ctzx_qmesh.py" : "ctzx_tri.py",
               df, sf, "d", false);
    nerr += check_errors(sd, sf);
  }

  if (nerr) printf("ctzx: test_fast_woodland_methods failed\n");
  return nerr;
}

static int test_fast_nonunirect () {
  return test_fast_woodland_methods(true);
}

static int test_matrix_form () {
  int nerr = 0;

  ZxFn::Shape zshape;
  Disloc::Ptr disloc;
  setup(20, zshape, disloc);
  const int disloc_comp = 0;

  ConvTest ct;
  ct.init(zshape, disloc);
  ct.set_verbosity(0);
  ct.set_use_surface_recon(0);
  ct.set_disloc_order(2);
  ct.set_nx(4);
  ct.set_use_halfspace(false);
  ct.discretize();

  const auto d = ct.get_discretization();
  const auto& gfp = ct.get_gfp();

  Stress stress(d);

  Stress::Options o, onbr;
  set_standard_options(o, onbr);

  // Form 6 nc x nc matrices.
  const auto& m = *d->get_mesh();
  const auto nc = m.get_ncell();
  RealArray A(6*nc*nc);
  ompparfor for (Idx si = 0; si < nc; ++si) {
    for (Idx ri = 0; ri < nc; ++ri) {
      Real s[6];
      stress.calc_matrix_entries(gfp, si, ri, disloc_comp, s, o, onbr);
      for (int i = 0; i < 6; ++i)
        A[i*nc*nc + ri*nc + si] = s[i];
    }
  }

  RealArray dd, sd;
  ct.eval(dd, sd, ConvTest::EvalMethod::direct);

  // 6 matvecs.
  RealArray sm(6*nc, 0);
  ompparfor for (Idx ri = 0; ri < nc; ++ri) {
    for (int comp = 0; comp < 6; ++comp) {
      Real accum = 0;
      for (Idx si = 0; si < nc; ++si)
        accum += A[comp*nc*nc + ri*nc + si] * dd[3*si + disloc_comp];
      sm[6*ri + comp] = accum;
    }
  }

  nerr += check_errors(sd, sm, 1e5*mv3::eps, 1e7*mv3::eps);

  return nerr;
}

static int test_hmatrix () {
  int nerr = 0;

#ifdef WOODLAND_HAVE_HMMVP
  ZxFn::Shape zshape;
  Disloc::Ptr disloc;
  setup(0, zshape, disloc);

  ConvTest ct;
  ct.init(zshape, disloc);
  ct.set_verbosity(0);
  ct.set_use_surface_recon(0);
  ct.set_nx(4);
  ct.set_use_halfspace(false);
  ct.set_use_nonunirect(true);
  ct.discretize();

  RealArray dd, sd;
  ct.eval(dd, sd, ConvTest::EvalMethod::direct);

  const auto d = ct.get_discretization();
  RealArray sh;
  nerr += hmat::Hmatrices::unittest(d, &ct.get_gfp(), &dd, &sh);

  nerr += check_errors(sd, sh, 1e5*mv3::eps, 1e-6);
#endif

  return nerr;
}

static int test_req0_and_scale () {
  int nerr = 0;

  ZxFn::Shape zshape;
  Disloc::Ptr disloc;
  setup(0, zshape, disloc);

  ConvTest ct;
  ct.init(zshape, disloc);
  ct.set_verbosity(0);
  ct.set_use_surface_recon(0);
  ct.set_use_exact_tangents(1);
  ct.set_disloc_order(2);
  ct.set_use_nonunirect(false);

  for (const bool use_halfspace : {false, true}) {
    ct.set_use_halfspace(use_halfspace);
    ct.set_nx(use_halfspace ? 5 : 6);
    ct.set_ny(3);
    RealArray dr, sr; // ref
    ct.set_mesh_scale(false);
    ct.set_support_req0(false);
    ct.eval(dr, sr, ConvTest::EvalMethod::direct);
    for (const bool scale : {false, true})
      for (const bool req0 : {false, true}) {
        if (not (scale or req0)) continue; // same as ref
        ct.set_mesh_scale(scale);
        ct.set_support_req0(req0);
        RealArray dt, st;
        ct.eval(dt, st, ConvTest::EvalMethod::direct);
        const int ne =
          scale ?
          check_errors(sr, st, 5e-5, 3e-3, 10) :
          check_errors(sr, st, 1e-8, 5e-7, 10);
        if (ne)
          printf("use_halfspace %d scale %d req0 %d\n",
                 int(use_halfspace), int(scale), int(req0));
        nerr += ne;
      }
  }

  return nerr;
}

int ConvTest::unittest () {
  int nerr = 0, ne;
  rununittest(ZxFn::unittest);
  rununittest(ExtrudedCubicSplineSurface::unittest);
  rununittest(FlatElementSurface::unittest);
  rununittest(ConvTest::test_interp);
  rununittest(test_req0_and_scale);
#ifdef WOODLAND_ACORN_HAVE_DC3D
  nerr += test_fast_okada_methods(false);
#endif
  nerr += test_fast_okada_methods(true);
  rununittest(test_fast_woodland_methods);
  rununittest(test_fast_nonunirect);
  rununittest(test_matrix_form);
  rununittest(test_hmatrix);
  return nerr;
}

} // namespace ctzx
} // namespace squirrel
} // namespace woodland
