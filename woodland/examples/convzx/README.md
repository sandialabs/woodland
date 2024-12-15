## examples/convzx

This directory contains a convergence test that demonstrates use of
Woodland::Acorn.

The surface is described by $z(x,y) = f(x)$ for $x,y \in [0,1]$. $f(x)$ is
usually $f(x) = 0.3 \ x \ \sin(2 \pi x) - 0.4$, which is interesting but simple
enough to permit efficient convergence testing. -0.4 makes the halfspace
Green's function relevant.

This surface can be tessellated by smoothly varying rectangles, permitting use
of the Okada Green's function (constant dislocation over a flat rectangle).

Each rectangle is then divided into triangles to test Woodland::Acorn. This is
not an optimal tessellation of the surface; rather, it's meant as a test of
handling triangles.

The shape of the surface and the regularity of the tessellations permits reuse
of Green's functions, substantially speeding up the computations and permitting
the convergence curves to be extended to high accuracy.

### Examples

Woodland::Acorn method; slip in the $x$ direction; full-space solution; use the exact surface when computing integrals; run three resolutions:
```
$ OMP_NUM_THREADS=4 ./bin/woodland_examples_convzx ct_we_zx 'testcase=20, halfspace=0, srfrecon=0, nres=3'

#threads 4
convtest_w_vs_e 20
ct> lam 9.000e-01 mu 1.100e+00 halfspace 0
ct> ntriperrect 2 dislocorder 2 exactsrf 1 flatelem 0 exacttan 0 tanorder 4 c2spline 0
ct> zshape trig1
ct> disloc 0 stapered  1.000  0.500  0.500  1.000  0.000 -0.000  1.000  0.800  0.800
nx ny 10 8
t: method 3.76e-01 (3.76e-01) exact 4.34e-01 (4.34e-01)
l2: 4.92e-02 (0.00) 8.16e-02 (0.00) 3.72e-02 (0.00) 5.52e-02 (0.00) 8.14e-02 (0.00) 4.54e-02 (0.00)
li: 7.10e-02 (0.00) 9.42e-02 (0.00) 2.88e-02 (0.00) 8.96e-02 (0.00) 9.83e-02 (0.00) 6.46e-02 (0.00)
nx ny 20 16
t: method 1.41e+00 (1.78e+00) exact 1.74e+00 (2.18e+00)
l2: 8.39e-03 (2.55) 1.33e-02 (2.62) 1.03e-02 (1.86) 1.00e-02 (2.46) 1.94e-02 (2.07) 7.96e-03 (2.51)
li: 1.39e-02 (2.35) 2.40e-02 (1.97) 8.21e-03 (1.81) 1.54e-02 (2.54) 2.24e-02 (2.14) 1.25e-02 (2.36)
nx ny 40 32
t: method 6.86e+00 (8.64e+00) exact 6.93e+00 (9.11e+00)
l2: 1.02e-03 (3.04) 1.88e-03 (2.82) 2.54e-03 (2.02) 1.26e-03 (2.99) 4.53e-03 (2.10) 9.95e-04 (3.00)
li: 2.52e-03 (2.46) 4.59e-03 (2.38) 2.52e-03 (1.71) 2.81e-03 (2.45) 6.17e-03 (1.86) 2.15e-03 (2.54)
```

The `ct>` lines show information about the specific problem. For details, see
`convtest_zx.cpp`, particularly `ConvTestSettings` and `ConvTest::print`. For
details about the test slip functions, see `setup` in that file.

The `nx ny` lines give the number of rectangles; `ntriperrect 2` then means each
rectangle is divided into two triangles.

The `t:` line displays computation time:
```
t: method this-resolution-in-sec (cumulative) exact this-resolution-in-sec (cumulative)
```
`exact` refers to the reference solution.

The `l2` ($l_2$-norm) and `li` ($l_{\infty}$-norm) lines show errors in the six
stress components: $xx$, $xy$, $xz$, $yy$, $yz$, $zz$. In each pair of numbers,
the first is the relative error. The second is the empirical order of accuracy,
computed with respect to the previous resolution; hence OOA is 0 for the
coarsest resolution.

Flat rectangles with Okada Green's function; slip in the $x$ direction;
full-space solution; run four resolutions:
```
$ OMP_NUM_THREADS=4 ./bin/woodland_examples_convzx ct_oe_zx 'testcase=20, halfspace=0, nres=4'

#threads 4
convtest_o_vs_e 20
ct> lam 9.000e-01 mu 1.100e+00 halfspace 0
ct> ntriperrect 2 dislocorder 2 exactsrf 0 flatelem 0 exacttan 0 tanorder 4 c2spline 0
ct> zshape trig1
ct> disloc 0 stapered  1.000  0.500  0.500  1.000  0.000 -0.000  1.000  0.800  0.800
nx ny 10 8
t: method 7.36e-04 (7.36e-04) exact 2.17e-01 (2.17e-01)
l2: 5.02e-01 (0.00) 3.10e-01 (0.00) 9.40e-02 (0.00) 2.61e-01 (0.00) 2.73e-01 (0.00) 3.82e-01 (0.00)
li: 4.98e-01 (0.00) 4.44e-01 (0.00) 8.60e-02 (0.00) 2.78e-01 (0.00) 3.57e-01 (0.00) 4.43e-01 (0.00)
nx ny 20 16
t: method 2.33e-03 (3.07e-03) exact 8.73e-01 (1.09e+00)
l2: 3.51e-01 (0.52) 1.25e-01 (1.31) 3.78e-02 (1.31) 1.49e-01 (0.81) 1.15e-01 (1.24) 1.94e-01 (0.98)
li: 4.31e-01 (0.21) 1.83e-01 (1.28) 4.03e-02 (1.09) 1.77e-01 (0.65) 1.44e-01 (1.31) 2.52e-01 (0.81)
nx ny 40 32
t: method 1.82e-02 (2.13e-02) exact 3.47e+00 (4.56e+00)
l2: 2.05e-01 (0.78) 5.16e-02 (1.27) 2.13e-02 (0.83) 1.00e-01 (0.57) 5.24e-02 (1.14) 9.53e-02 (1.03)
li: 2.73e-01 (0.66) 7.59e-02 (1.27) 2.84e-02 (0.50) 1.32e-01 (0.43) 6.26e-02 (1.20) 1.21e-01 (1.06)
nx ny 80 64
t: method 1.84e-01 (2.05e-01) exact 1.37e+01 (1.83e+01)
l2: 1.10e-01 (0.90) 2.27e-02 (1.18) 1.17e-02 (0.86) 5.78e-02 (0.79) 2.50e-02 (1.07) 4.71e-02 (1.02)
li: 1.51e-01 (0.85) 3.37e-02 (1.17) 1.72e-02 (0.73) 8.00e-02 (0.72) 2.89e-02 (1.11) 6.22e-02 (0.96)
```

Options:
```
option: type [default]

testcase: int [0]
ntri: 2, 4 [2]
srfrecon: bool (0,1) [1]
flatelem: bool [0]
exacttan: bool [0]
tanorder: 2, 4 [4]
c2spline: bool [0]  dislocorder: 0, 1, 2, 3 [2]
nres: int > 0 [-1]
woodlandrect: bool [0]
halfspace: bool [0]
```