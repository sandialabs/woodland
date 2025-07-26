#ifndef INCLUDE_WOODLAND_SQUIRREL_PYWRITE
#define INCLUDE_WOODLAND_SQUIRREL_PYWRITE

#include "woodland/squirrel/squirrel.hpp"

#include <cstdio>

namespace woodland {
namespace squirrel {

void pywrite_header(FILE* fp);

void pywrite_double_array(FILE* fp, const std::string& var,
                          const int n, CRPtr a);

void pywrite_double_array(FILE* fp, const std::string& var,
                          const int m, const int n, CRPtr a);

} // namespace squirrel
} // namespace woodland

#endif
