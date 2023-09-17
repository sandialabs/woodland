// Derived from https://github.com/E3SM-Project/COMPOSE; see
//     https://github.com/E3SM-Project/COMPOSE/blob/main/LICENSE
// for details of the 3-clause BSD license.

#ifndef INCLUDE_WOODLAND_ACORN_TRIQUAD
#define INCLUDE_WOODLAND_ACORN_TRIQUAD

#include "woodland/acorn/acorn.hpp"

namespace woodland {
namespace acorn {

struct TriangleQuadrature {
  static bool is_order_supported(const int order);
  
  static void get_coef(const int order,
                       // quadrature barycentric coordinates; list of triples
                       const Real*& coord,
                       // quadrature point weights; sum to 1
                       const Real*& weight,
                       // number of quadrature points
                       int& n);

  static int unittest();

private:
  static const Real trisym_order1_coord_  [  3];
  static const Real trisym_order1_weight_ [  1];
  static const Real trisym_order2_coord_  [  9];
  static const Real trisym_order2_weight_ [  3];
  static const Real trisym_order4_coord_  [ 18];
  static const Real trisym_order4_weight_ [  6];
  static const Real tritay_order6_coord_  [ 33];
  static const Real tritay_order6_weight_ [ 11];
  static const Real trisym_order8_coord_  [ 48];
  static const Real trisym_order8_weight_ [ 16];
  static const Real tritay_order12_coord_ [ 96];
  static const Real tritay_order12_weight_[ 32];
  static const Real trisym_order20_coord_ [264];
  static const Real trisym_order20_weight_[ 88];
};

} // namespace acorn
} // namespace woodland

#endif
