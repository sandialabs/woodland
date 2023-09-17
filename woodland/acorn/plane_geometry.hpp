#ifndef INCLUDE_WOODLAND_ACORN_PLANE_GEOMETRY
#define INCLUDE_WOODLAND_ACORN_PLANE_GEOMETRY

#include <cassert>

#include "woodland/acorn/acorn.hpp"

namespace woodland {
namespace acorn {
namespace plane {

typedef Real Pt[2];  // A point is an (x,y) 2-vector.
typedef Real Tri[6]; // A triangle is a list of three (x,y) 2-vectors.

// A 2D polygon described by a periodic list of vertex coordinates with edges
// implied between each pair of vertices.
struct Polygon {
  CRPtr xys;       // [x1,y1,...,xn,yn]
  const int n;
  Polygon (CRPtr xys_, const int n_) : xys(xys_), n(n_) { assert(n >= 3); }
};

// Does the triangle (a,b,c) have counterclockwise (CCW) orientation?
bool is_ccw(const Pt a, const Pt b, const Pt c);

// p must be convex.
bool is_ccw(const Polygon& p,
            // Assume there are no collinear vertices.
            const bool assume_noncollinear=true);

// Is p to the right of the ray v1->v2?
bool is_point_to_right_of_ray(const Pt p, const Pt v1, const Pt v2);

bool is_point_inside_convex_polygon(const Polygon& p, const Pt c);

// Project p onto the line (v1 + v2)/2 + a (v2 - v1). The projection point is on
// the segment if a in [-0.5, 0.5].
void project_p_onto_line(const Pt p, const Pt v1, const Pt v2,
                         Real& a, Pt projection);

// Return the radius^2 of the largest circle centered at s that fits inside
// convex polygon p. It is assumed, but not checked, that p is convex and s is
// inside p.
//   If support_sides is provided, a side to which the circle is tangent is
// provided.
Real inscribe_largest_circle(const Polygon& p, const Pt s,
                             int* support_side = nullptr);

// Shortest distance between the point p and the line
//     v1 + a (v2 - v1), -inf < a < inf.
Real calc_shortest_distance2_between_point_and_line(
  const Pt p, const Pt v1, const Pt v2);

/* calc_tangent_points
     Let c = cc - a, s = norm(c), L = sqrt(s^2 - r^2). The tangent condition is
       dot(L [cos t, sin t], L [cos t, sin t] - c) = 0.
   =>
       L^2 - L c'[cos t, sin t] = 0.
   Let
       d = L c(1)
       e = L c(2)
       f = L^2
   in
       d cos t + e sin t = f.                         (trig eq)
   Then
       v = -1 or 1
       d v sqrt(1 - sin^2 t) + e sin t = f
       x = sin t
       d^2 (1 - x^2) = (f - e x)^2 = f^2 - 2 e f x + e^2 x^2
       (e^2 + d^2) x^2 - 2 e f x + f^2 - d^2 = 0 => (x1, x2)
   This gives four possible solutions:
       (cos t, sin t) in {(v sqrt(1 - xi^2), xi)}_{i=1,2, v=-1,1}.
   Choose the two satisfying (trig eq). Return false if a is inside the circle.
 */
bool calc_tangent_points(const Pt a, const Pt cc, const Real radius2,
                         Pt pts[2]);

// Return false if the lines are collinear, otherwise true.
bool intersect_line_line(
  const Pt xa, const Pt xb, // line xa + a (xb - xa)  
  const Pt ya, const Pt yb, // line ya + b (yb - ya)
  Pt p,                     // intersection
  // If requested, compute p_xb with row-major order.
  Real* p_xb = nullptr);

// Intersection of xa + a (xb - xa) with the nearer of the two points to xa on
// the circle. Return false if the line does not intersect the circle, otherwise
// true. If assume_intersect, then compute the roots assuming the discriminant
// is >= 0.
bool intersect_line_circle_nearer(
  const Pt xa, const Pt xb,        // line xa + a (xb - xa)
  const Pt cc, const Real radius2, // circle centered at cc; provide radius^2
  Pt p,                            // intersection point nearer xa
  const bool assume_intersect = false,
  // If requested, compute p_xb with row-major order.
  Real* p_xb = nullptr);

int unittest();

} // namespace plane
} // namespace acorn
} // namespace woodland

#endif
