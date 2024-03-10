#include "woodland/examples/convzx/convtest_zx_hmmvp.hpp"

#ifdef WOODLAND_HAVE_HMMVP
#include "woodland/acorn/matvec.hpp"
#include "woodland/acorn/util.hpp"
#include "woodland/acorn/fs3d.hpp"
#include "woodland/acorn/dbg.hpp"

#include "hmmvp/include/Compress.hpp"
#include "hmmvp/include/Hmat.hpp"
#endif

namespace woodland {
namespace examples {
namespace convzx {

#ifdef WOODLAND_HAVE_HMMVP

using acorn::zero;
typedef acorn::Matvec<3,Real> mv3;
typedef hmmvp::UInt UInt;

struct ConvTest::Hmmvp::OkadaData {
  EvalMethod eval_method = EvalMethod::direct_hmmvp;
  Real lam = -1, mu = -1, tol = -1;
  ZxFn::Ptr zxfn;
  int nyr = -1;
  hmmvp::Compressor::TolMethod tol_method = hmmvp::Compressor::tm_mrem_fro;
  enum : int { nhmat = 18 };
  hmmvp::Hmat* hmats[nhmat] = {0};

  void reset_hmats () {
    for (int i = 0; i < nhmat; ++i)
      if (hmats[i]) {
        hmmvp::DeleteHmat(hmats[i]);
        hmats[i] = nullptr;
      }
  }

  ~OkadaData () { reset_hmats(); }
};

struct OkadaGreensFn : public hmmvp::GreensFn {
  OkadaGreensFn (const ConvTest::Hmmvp::OkadaData& d_, const int di, const int si)
    : d(d_), disloc_idx(di), sigma_idx(si), nrect(d.zxfn->get_nx()*d.nyr)
  {}

  bool Call (const int nr, const int nc, const UInt* const rs,
             const UInt* const cs, double* B) const {
    neval += nr*nc;
    const Real yhat[] = {0, 1, 0};
    for (int ridx = 0; ridx < nr; ++ridx) {
      assert(rs[ridx] > 0 && rs[ridx] <= nrect);
      const auto ir = rs[ridx] - 1;
      Real rcv[3], rnml[3], rxhat[3];
      ConvTest::calc_rect_ctr(*d.zxfn, d.nyr, ir, rcv, rnml);
      mv3::cross(yhat, rnml, rxhat);
      for (int cidx = 0; cidx < nc; ++cidx) {
        assert(cs[cidx] > 0 && cs[cidx] <= nrect);
        const auto is = cs[cidx] - 1;
        const auto Bidx = nr*cidx + ridx;
        if (is == ir && sigma_diag) {
          assert(not all_sigma_mode);
          B[Bidx] = sigma_diag[6*is + sigma_idx];
          continue;
        }
        Real src[3], nml[3], lengths[2], xhat[3], dlcl[3] = {0}, dglbl[3];
        ConvTest::calc_rect_ctr(*d.zxfn, d.nyr, is, src, nml, lengths);
        mv3::cross(yhat, nml, xhat);
        dlcl[disloc_idx] = 1;
        acorn::tmatvec(xhat, yhat, nml, dlcl, dglbl);
        Real sigma[6];
        acorn::fs3d::calc_sigma_const_disloc_rect_okada(
          d.lam, d.mu, src, nml, xhat, lengths, dglbl, rcv, sigma);
        acorn::rotate_sym_tensor_3x3_RARt(rxhat, yhat, rnml, sigma);
        if (all_sigma_mode)
          acorn::copy(6, sigma, &B[6*Bidx]);
        else
          B[Bidx] = sigma[sigma_idx];
      }
    }
    return true;
  }
  
  virtual bool
  Call (const hmmvp::CompressBlockInfo& cbi, const std::vector<UInt>& rs,
        const std::vector<UInt>& cs, double* B) const {
    assert(not all_sigma_mode || (rs.size() == size_t(1) && cs.size() == size_t(1)));
    return Call(rs.size(), cs.size(), rs.data(), cs.data(), B);
  }

  // To get diagonal elements efficiently:
  void set_calc_all_sigma_mode (const bool on) { all_sigma_mode = on; }

  void set_sigma_diag (const Real* const sd) { sigma_diag = sd; }

  UInt get_neval () const { return neval; }

private:
  const ConvTest::Hmmvp::OkadaData& d;
  const int disloc_idx, sigma_idx;
  bool all_sigma_mode = false;
  const Real* sigma_diag = nullptr;
  UInt nrect;
  mutable UInt neval = 0;
};

static std::string
make_file_name (const std::string method, const ZxFn& zxfn, const int nyr,
                const int disloc_idx, const int sigma_idx,
                hmmvp::Compressor::TolMethod tol_method, const Real tol) {
  const int nbuf = 16;
  char tol_str[nbuf];
  snprintf(tol_str, nbuf, "%1.3e", tol);
  return std::string("hmats/") + method + "_" + ZxFn::convert(zxfn.get_shape()) +
    "_" + std::to_string(zxfn.get_nx()) + "_" + std::to_string(nyr) + "_" +
    std::to_string(disloc_idx) + "_" + std::to_string(sigma_idx) + "_" +
    std::to_string(int(tol_method)) + "_" + tol_str + ".hmat";
}

void ConvTest::Hmmvp
::compress_okada_if_not (const EvalMethod eval_method,
                         const ZxFn& zxfn, const int nyr) {
  assert(ConvTest::is_hmmvp(eval_method));
  const bool force_new_hmats = false;

  if (not od) od = std::make_shared<OkadaData>();
  auto& d = *od;
  if (d.eval_method == eval_method && d.lam == lam && d.mu == mu &&
      d.tol == tol && d.nyr == nyr && d.zxfn &&
      d.zxfn->get_nx() == zxfn.get_nx() && not force_new_hmats)
    return;
  d.reset_hmats();
  d.eval_method = eval_method;
  d.lam = lam; d.mu = mu; d.tol = tol;
  d.zxfn = std::make_shared<ZxFn>(zxfn);
  const auto nxr = d.zxfn->get_nx();
  d.nyr = nyr;

  const bool direct = eval_method == EvalMethod::direct_hmmvp;
  const std::string method_str = std::string("okada_") +
    (direct ? "direct" : "fast");

  bool all = not force_new_hmats;
  if (all) {
    for (int di = 0; di < 3; ++di) {
      for (int si = 0; si < 6; ++si) {
        const auto filename =
          make_file_name(method_str, *d.zxfn, nyr, di, si, d.tol_method, d.tol);
        FILE* fp = fopen(filename.c_str(), "r");
        if (fp)
          fclose(fp);
        else {
          all = false;
          break;
        }
        try {
          d.hmats[6*di + si] = hmmvp::NewHmat(filename);
        } catch (...) {
          all = false;
          printf("compress_okada_if_not: Error reading %s\n", filename.c_str());
          break;
        }
      }
      if (not all) break;
    }
    if (all && not force_new_hmats) return;
    d.reset_hmats();
  }

  const UInt nrect = nxr*nyr;
  hmmvp::Hd* hd; {
    hmmvp::Matrix<Real> C(3, nrect);
    ompparfor for (auto ir = zero(nrect); ir < nrect; ++ir) {
      Real p[3];
      ConvTest::calc_rect_ctr(*d.zxfn, nyr, ir, p);
      for (int d = 0; d < 3; ++d)
        C(d+1,ir+1) = p[d];
    }
    if (direct)
      hd = hmmvp::NewHd(C);
    else {
      hmmvp::Matrix<Real> R(3, nxr);
      ompparfor for (auto ir = zero(nxr); ir < nxr; ++ir) {
        Real p[3];
        ConvTest::calc_rect_ctr(*d.zxfn, nyr, ir, p);
        for (int d = 0; d < 3; ++d)
          R(d+1,ir+1) = p[d];
      }
      hd = hmmvp::NewHd(R, C);
    }
  }

  ompparfor for (int k = 0; k < 18; ++k) {
    const int di = k / 6, si = k % 6;
    OkadaGreensFn ogf(d, di, si);
    hmmvp::Compressor* c = hmmvp::NewCompressor(hd, &ogf);
    c->SetTolMethod(d.tol_method);
    c->SetNumberOfRejectedZeroRowsBeforeStopping(4);
    c->SetOutputLevel(0);
    c->AvoidRedundantGfCalls(true);
    c->SetOmpNthreads(1);
    const double norm_estimate = c->EstimateBfro();
    c->SetBfroEstimate(norm_estimate);
    c->SetTol(d.tol);
    const auto filename = make_file_name(method_str, *d.zxfn, nyr, di, si,
                                         c->GetTolMethod(), d.tol);
    FILE* fp = fopen(filename.c_str(), "r");
    const bool new_hmat = not fp || force_new_hmats;
    if (new_hmat)
      c->CompressToFile(filename);
    else
      fclose(fp);
    DeleteCompressor(c);
    d.hmats[k] = hmmvp::NewHmat(filename);
    if (new_hmat) {
      const auto Bnnz = acorn::square(UInt(nxr*nyr));
      printf("%d %d: B nnz %lu H nnz %lu (%1.3e, %1.3e) neval %lu (%1.3e)\n"
             "fro norm %d %d %1.3e (est %1.3e)\n",
             di, si, Bnnz, d.hmats[k]->GetNnz(),
             Real(d.hmats[k]->GetNnz()) / Bnnz,
             Real(d.hmats[k]->GetNnz()) / (nxr*nyr),
             ogf.get_neval(),
             Real(ogf.get_neval()) / Bnnz,
             di, si, std::sqrt(d.hmats[k]->NormFrobenius2()),
             norm_estimate);
    }
  }

  if (0) {
    prc(nrect);
    bool hdm = false;
    if (acorn::dbg::getenv("hmmvp_dump_matrices", hdm) && hdm && nrect < 2000) {
      printf("Writing full matrices.\n");
      std::vector<UInt> cs(nrect);
      std::vector<hmmvp::Blint> ics(nrect);
      for (auto i = zero(nrect); i < nrect; ++i)
        cs[i] = i + 1;
      for (auto i = zero(nrect); i < nrect; ++i)
        ics[i] = i;
      const auto nthr = omp_get_max_threads();
      std::vector<std::vector<Real>> data(nthr), Bs(6);
      for (int i = 0; i < 6; ++i) Bs[i].resize(nrect*nrect);
      FILE* fp = fopen("Bmats.py", "w");
      fprintf(fp, "import numpy as npy\nB = {}\n");
      for (int di = 0; di < 3; ++di) {
        if (di != 0) continue;
        printf("Writing di %d\n", di);
        for (int imat = 0; imat < 2; ++imat) {
          if (imat == 0) {
            OkadaGreensFn ogf(d, di, -1);
            ogf.set_calc_all_sigma_mode(true);
            ompparfor for (auto ir = zero(nrect); ir < nrect; ++ir) {
              const auto tid = omp_get_thread_num();
              data[tid].resize(6*nrect);
              Real* const sigma = data[tid].data();
              const auto irp1 = ir + 1;
              ogf.Call(1, nrect, &irp1, cs.data(), sigma);
              for (int j = 0; j < 6; ++j)
                for (auto i = zero(nrect); i < nrect; ++i)
                  Bs[j][nrect*ir + i] = sigma[6*i + j];
            }
          } else {
            for (int j = 0; j < 6; ++j) {
              const auto h = d.hmats[6*di + j];
              //h->TurnOffPermute();
              h->Extract(ics, ics, Bs[j].data());
              //h->TurnOnPermute();
            }
          }
          for (int j = 0; j < 6; ++j) {
            if (j != 4 && j != 2) continue;
            fprintf(fp, "B[(%d,%d,%d)] = npy.array([", imat, di, j);
            for (auto r = zero(nrect); r < nrect; ++r) {
              fprintf(fp, "[");
              for (auto c = zero(nrect); c < nrect; ++c) {
                fprintf(fp, " %1.3e,", Bs[j][r*nrect+c]);
                if (c % 20 == 0) fprintf(fp, "\n");
              }
              fprintf(fp, "],\n");
            }
            fprintf(fp, "])\n");
          }
        }
      }
      fclose(fp);
      printf("Done.\n");
    }
  }
}

void ConvTest::Hmmvp
::eval_okada(const Disloc& disloc, RealArray& dislocs, RealArray& sigmas) const {
  assert(od);
  const auto& d = *od;
  const int nxr = d.zxfn->get_nx();
  const int nrect = nxr*d.nyr;

  fill_dislocs(*d.zxfn, d.nyr, disloc, dislocs);

  sigmas.clear();
  sigmas.resize(6*nrect, 0);
  std::vector<Real> x(nrect), y(nrect);
  for (int di = 0; di < 3; ++di) {
    for (int i = 0; i < nrect; ++i)
      x[i] = dislocs[3*i + di];
    for (int si = 0; si < 6; ++si) {
      d.hmats[6*di + si]->Mvp(x.data(), y.data(), 1);
      for (int i = 0; i < nrect; ++i)
        sigmas[6*i + si] += y[i];
    }
  }
}

void ConvTest::Hmmvp
::eval_okada_fast (const SupportPoints& supports, const Disloc& disloc,
                   RealArray& dislocs, RealArray& sigmas) const {
  throw_if(true, "not impl'ed");
}

struct ConvTest::Hmmvp::WoodlandData {
  Real lam = -1, mu = -1, tol = -1;
  Discretization::CPtr d;
  Stress::Options cso;
};

void ConvTest::Hmmvp
::compress_if_not (const Discretization::CPtr& d, Stress::Options& o) {
  
}

#else // WOODLAND_HAVE_HMMVP

void ConvTest::Hmmvp
::compress_okada_if_not (const EvalMethod eval_method, const ZxFn& zxfn,
                         const int nyr)
{}

void ConvTest::Hmmvp
::compress_if_not (const Discretization::CPtr& d, Stress::Options& o) {}

void ConvTest::Hmmvp
::eval_okada (const Disloc& disloc, RealArray& dislocs,
              RealArray& sigmas) const
{}

void ConvTest::Hmmvp
::eval_okada_fast (const SupportPoints& supports, const Disloc& disloc,
                   RealArray& dislocs, RealArray& sigmas) const
{}

#endif // WOODLAND_HAVE_HMMVP

} // namespace convzx
} // namespace examples
} // namespace woodland
