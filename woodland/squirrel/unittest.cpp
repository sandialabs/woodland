#include "woodland/squirrel/unittest.hpp"
#include "woodland/squirrel/mesh.hpp"
#include "woodland/squirrel/discretization.hpp"
#include "woodland/squirrel/exact.hpp"
#include "woodland/squirrel/ctzx.hpp"
#include "woodland/squirrel/solver1d.hpp"
#include "woodland/squirrel/bezier_cubic.hpp"

#include <sys/stat.h>

namespace woodland {
namespace squirrel {

int unittest () {
  int nerr = 0, ne;
  // Several tests write files in tmp/.
  mkdir("tmp", 0777);
  rununittest(Solver1d<Real>::unittest);
  rununittest(BezierCubic<Real>::unittest);
  rununittest(mesh::unittest);
  rununittest(Exact::unittest);
  rununittest(Discretization::unittest);
  rununittest(ctzx::ConvTest::unittest);
  printf("\n%s\n", nerr ? "FAIL" : "PASS");
  return nerr;
}

} // namespace squirrel
} // namespace woodland
