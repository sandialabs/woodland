Woodland will provide methods to construct the Displacement Discontinuity Method
(DDM) elastostatic operator for curved faults. This operator can then be used in
quasidynamic rate-state friction earthquake simulators.

Woodland is a work in progress and is not yet ready for general use. Currently,
low-level calculations are implemented in the directory `woodland/acorn` for the
self-interaction Hadamard finite part and the other-interaction proper integral
for convex polygons on curved fractures with general dislocations, where the
fault shape and the dislocation must be smooth within a polygon.

In the future, the directory `woodland/squirrel` will contain a DDM
discretization based on these tools. Currently, it is a work in progress and
contains only a discretization to support a convergence test.

Finally, the directory `woodland/oak` will contain a simple quasidynamic
rate-state friction simulator to demonstrate the use of the operator.

## Build and run unit test

At the command line, run the following commands, where `~/tmp/woodland_install`
is an example of the location to install the library files.
```
cmake PATH_TO_WOODLAND_DIRECTORY \
      -D CMAKE_BUILD_TYPE=RelWithDebInfo \
      -D CMAKE_INSTALL_PREFIX=~/tmp/woodland_install;
make -j4 install
OMP_NUM_THREADS=4 ctest -VV
```
The tests will purposely fail if they do not have access to Y. Okada's `dc3d.f`
code. See `extern/README.md` for instructions to obtain and patch this file.
Once it is available and patched, the above commands will detect the file, and
the tests should pass.

## References

If you use Woodland, please cite

```
@misc{woodland-software,
  title={{Woodland: Methods to construct the elasticity operator for curved fractures in the Displacement Discontinuity Method}},
  author={Andrew M. Bradley},
  howpublished={[Computer Software] \url{https://github.com/ambrad/woodland}},
  year={2023}
}
```
