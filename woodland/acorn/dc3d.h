#ifndef INCLUDE_ACORN_DC3D
#define INCLUDE_ACORN_DC3D

#include "woodland/acorn/types.hpp"

// Interface to extern/dc3d.f

#ifdef WOODLAND_ACORN_HAVE_DC3D

extern "C" void dc3d0_(
  // in
  const char* SPACE, const double* ALPHA,
  const double* X, const double* Y, const double* Z,
  const double* DEPTH, const double* DIP, const double* POT1,
  const double* POT2, const double* POT3, const double* POT4,
  // out
  double* UX, double* UY, double* UZ,
  double* UXX, double* UYX, double* UZX,
  double* UXY, double* UYY, double* UZY,
  double* UXZ, double* UYZ, double* UZZ,
  int* iret);

extern "C" void dc3d_(
  // in
  const char* SPACE, const double* ALPHA,
  const double* X, const double* Y, const double* Z,
  const double* DEPTH, const double* DIP,
  const double* AL1, const double* AL2, const double* AW1, const double* AW2,
  const double* DISL1, const double* DISL2, const double* DISL3,
  // out
  double* UX, double* UY, double* UZ,
  double* UXX, double* UYX, double* UZX,
  double* UXY, double* UYY, double* UZY,
  double* UXZ, double* UYZ, double* UZZ,
  int* iret);

#endif // WOODLAND_ACORN_HAVE_DC3D

#endif
