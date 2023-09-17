This directory is for external code. If a code is publicly available but does
not have a license, we provide instructions to obtain it.

## Okada rectangular dislocation

This code is described in

> Y. Okada, Internal deformation due to shear and tensile faults in a half-space, Bull. Seism. Soc. Am., 1992.

The original version is on the website of Japan's National Research Institute for Earth Science and Disaster Resilience (NIED) [here](https://www.bosai.go.jp/information/dc3d_e.html).

The patch file 0001-dc3d.f-Modify-original-source-in-four-ways.patch makes small
modifications to the original code for thread safety and to support a full-space
option. Apply it like this:

1. Download the original code from [here](https://www.bosai.go.jp/information/pdf/copy_of_DC3Dfortran.txt); call it `dc3d.f`.
2. You will likely need to convert the newline convention from DOS to UNIX. On a typical Linux setup, you can run `dos2unix dc3d.f`
3. `git add dc3d.f`
4. `git commit dc3d.f -m "original dc3d.f file"`
5. `git am 0001-dc3d.f-Modify-original-source-in-four-ways.patch`
6. If you want to keep the repo unmodified, `git reset HEAD~2` to restore the original state while still keeping the updated `dc3d.f` file.

To enable all unit tests, `dc3d.f` must be added to this directory and then
updated with the patch before building Woodland.
