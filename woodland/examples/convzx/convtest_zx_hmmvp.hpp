#ifndef INCLUDE_WOODLAND_EXAMPLES_CONVZX_CONVTEST_ZX_HMMVP
#define INCLUDE_WOODLAND_EXAMPLES_CONVZX_CONVTEST_ZX_HMMVP

#include "woodland/examples/convzx/convtest_zx.hpp"
#include "woodland/examples/convzx/calc_stress.hpp"

namespace woodland {
namespace examples {
namespace convzx {

struct ConvTest::Hmmvp {
  typedef std::shared_ptr<Hmmvp> Ptr;
  typedef std::shared_ptr<const Hmmvp> CPtr;

  void set_lam_mu (const Real lam_, const Real mu_) { lam = lam_; mu = mu_; }
  void set_tol (const Real tol_) { tol = tol_; }

  Real get_lam () const { return lam; }
  Real get_mu  () const { return mu ; }
  Real get_tol () const { return tol; }

  // Compress the Okada operator.
  void compress_okada_if_not(const EvalMethod eval_method,
                             const ZxFn& zxfn, const int nyr);
  // Compress the Woodland operator.
  void compress_if_not(const Discretization::CPtr& d, Stress::Options& o);

  void eval_okada(const Disloc& disloc, RealArray& dislocs,
                  RealArray& sigmas) const;
  void eval_okada_fast(const SupportPoints& supports, const Disloc& disloc,
                       RealArray& dislocs, RealArray& sigmas) const;

public: // for external implementations
  struct OkadaData;
  struct WoodlandData;

private:
  Real lam = -1, mu = -1, tol = -1;
  std::shared_ptr<OkadaData> od;
  std::shared_ptr<WoodlandData> wd;
};

} // namespace convzx
} // namespace examples
} // namespace woodland

#endif
