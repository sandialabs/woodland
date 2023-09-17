#include "woodland/acorn/matvec.hpp"

namespace woodland {
namespace acorn {

template <int dim, typename Real>
int Matvec<dim,Real>::unittest () {
  return 0;
}

template struct Matvec<2,double>;
template struct Matvec<3,double>;

} // namespace acorn
} // namespace woodland
