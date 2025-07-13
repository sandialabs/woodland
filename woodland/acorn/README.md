## Acorn

`woodland::acorn` is a library containing methods to integrate the elastostatic
Displacement Discontinuity Method Green's function over convex polygonal
elements on a curved fault.

It is largely complete.

### Interface

Acorn provides routines to compute the self-interaction (`calc_hfp` in
`interaction_integrals.hpp`) and other-interaction (`calc_integral` in
`interaction_integrals.hpp`) integrals in the DDM between curved polygons with
high-order dislocation representation.

A curved polygon is a polygon on a two-dimensional manifold. It is the
_physical_ polygon. Each curved polygon has a _reference_ polygon that is a
flat, straight-edged polygon and a map between reference and physical
coordinates. `calc_hfp` and `calc_integral` are given reference polygons.

The caller's integrand object, `CallerIntegrands` in `caller_integrand.hpp`,
needs to know how to map reference coordinate to physical, to evaluate the
Jacobian determinant associated with this map, and to evaluate dislocation at
the given reference coordinate.

`vv.hpp` is a set of verification and validation problems and shows how to use
these routines and construct a problem-specific `CallerIntegrands` structure.

### Methods

`calc_integral` computes a proper integral using triangular quadrature.

The self-interaction integral is the crux of a DDM. The integrand is
hypersingular, and thus the integral must be evaluated as a Hadamard finite part
(f.p.). The integrand has singularity $r^{-3}$ for radial distance from the
receiver point $r$. The f.p. is evaluated as follows.

Decompose each reference polygon into triangles rooted at the receiver
point. Since this is a self-interaction integral, this point is inside the
polygon. For a DDM, it is well separated from the polygon edges. Evalute the
f.p. over each triangle and then sum the results.

Now consider a triangle. Inscribe the largest circular sector centered at the
receiver point. This leaves either one or two regions, each having a specific
structure: one edge that is a radial segment, one that is a triangle edge, and
one that is an arc of the circle. In the code, we refer to such a shape as a
_radtricirc_. The integral is proper over each such region and is evaluated
using a quadrature grid suited to the region's specific structure.

This leaves an f.p. calculation over the circular sector. First, integrate in
the angular direction, reducing the two-dimensional $r^{-3}$ singularity to a
one-dimensional $r^{-2}$ one. This integral is proper and smooth, so standard
quadrature can be used. Next, evaluate the f.p. of the one-dimensional
integral. This uses the methods implemented in `hfp.*pp`. See the header
comments in `hfp.hpp` for details. Finally, sum the values from the circular
sector and the one or two radtricircs.
