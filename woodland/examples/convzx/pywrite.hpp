#ifndef INCLUDE_WOODLAND_EXAMPLES_CONVZX_PYWRITE
#define INCLUDE_WOODLAND_EXAMPLES_CONVZX_PYWRITE

#include "woodland/examples/convzx/convzx.hpp"

#include <cstdio>

namespace woodland {
namespace examples {
namespace convzx {

void pywrite_header(FILE* fp);

void pywrite_double_array(FILE* fp, const std::string& var,
                          const int n, CRPtr a);

void pywrite_double_array(FILE* fp, const std::string& var,
                          const int m, const int n, CRPtr a);

} // namespace convzx
} // namespace examples
} // namespace woodland

#endif
