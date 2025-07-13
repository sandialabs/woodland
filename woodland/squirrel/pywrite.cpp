#include "woodland/squirrel/pywrite.hpp"

namespace woodland {
namespace squirrel {

void pywrite_header (FILE* fp) {
  fprintf(fp, "import numpy as npy\n");
}

void pywrite_double_array (FILE* fp, const std::string& var,
                           const int n, CRPtr a) {
  fprintf(fp, "%s = npy.array([", var.c_str());
  for (int j = 0; j < n; ++j) {
    fprintf(fp, "%12.5e,", a[j]);
    if ((j+1) % 8 == 0) fprintf(fp, "\n");
  }
  fprintf(fp, "])\n");
}

void pywrite_double_array (FILE* fp, const std::string& var,
                           const int m, const int n, CRPtr a) {
  fprintf(fp, "%s = npy.array([", var.c_str());
  for (int i = 0, k = 0; i < m; ++i) {
    fprintf(fp, "[");
    for (int j = 0; j < n; ++j, ++k) {
      fprintf(fp, "%16.9e,", a[k]);
      if ((j+1) % 8 == 0) fprintf(fp, "\n");
    }
    fprintf(fp, "],\n");
  }
  fprintf(fp, "])\n");
}

} // namespace squirrel
} // namespace woodland
