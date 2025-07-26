#ifndef INCLUDE_WOODLAND_SQUIRREL_SOLVER1D
#define INCLUDE_WOODLAND_SQUIRREL_SOLVER1D

namespace woodland {
namespace squirrel {

template <typename Real>
struct Solver1d {
  struct Fn {
    virtual void eval(const Real& x, Real& f, Real& fp) const = 0;
  };

  struct Tols {
    // Absolute tolerances on x bracket width and |f| deviation from 0.
    Real xtol, ftol;
    Tols();
    Tols(const Real xtol, const Real ftol);
  };

  struct Info {
    enum Result : int {
      not_run = -3,
      input_f_same_sign = -2, tol_not_reached = -1,
      success_on_xtol = 0, success_on_ftol = 1
    };

    Real x, f, fp;
    int neval, nbisect;
    Result result;

    Info () : x(-1), f(-1), fp(-1), neval(0), nbisect(0), result(not_run) {}

    bool success () const {
      return (result == success_on_xtol or
              result == success_on_ftol);
    }
  };

  static void fzero(const Tols& tols, const Fn& fn,
                    const Real& xlo, const Real& xhi,
                    const Real& flo, const Real& fhi,
                    Info& info);

  static int unittest(const bool verbose = false);
};

} // namespace squirrel
} // namespace woodland

#endif
