#include "woodland/examples/convzx/unittest.hpp"
#include "woodland/examples/convzx/discretization.hpp"
#include "woodland/examples/convzx/exact.hpp"
#include "woodland/examples/convzx/convtest_zx.hpp"

namespace woodland {
namespace examples {
namespace convzx {

#define rununittest(f) do {                     \
    ne = f();                                   \
    if (ne) printf(#f " ne %d\n", ne);          \
    nerr += ne;                                 \
  } while (0)

int unittest () {
  int nerr = 0, ne;
  rununittest(Exact::unittest);
  rununittest(Discretization::unittest);
  rununittest(ConvTest::unittest);
  printf("\n%s\n", nerr ? "FAIL" : "PASS");
  return nerr;
}

} // namespace convzx
} // namespace examples
} // namespace woodland
