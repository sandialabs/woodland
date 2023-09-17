#include "woodland/acorn/elastostatics_integrals.hpp"

namespace woodland {
namespace acorn {

int get_triquad_order (const Real L, const Real dist) {
  return (dist <  3*L ? 20 :
          dist <  4*L ? 12 :
          dist <  6*L ?  8 :
          dist < 16*L ?  4 :
          dist < 24*L ?  2 :
          1);
}

} // namespace acorn
} // namespace woodland
