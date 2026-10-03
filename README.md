Woodland provides methods to construct the Displacement Discontinuity Method
(DDM) elastostatic operator for curved faults. This operator can then be used in
quasidynamic rate-state friction earthquake simulators.

Low-level calculations are implemented in the directory `woodland/acorn`. These
include the self-interaction Hadamard finite part and the other-interaction
proper integral for convex polygons on curved fractures with general
dislocations in elastic whole- and half-spaces.

`woodland/squirrel` provides a basic DDM discretization API built on these
tools. A convergence test demonstrates and exercises this API.

The directory `woodland/oak` may one day contain a simple quasidynamic
rate-state friction simulator to demonstrate the use of the operator. Currently,
it is empty.

Finally, `doc/theory.pdf` describes the methods used in this software.

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
