# Building DHSVM on DCC

This records the recipe that builds the DHSVM binary on the Duke Compute
Cluster. The build is reproducible. The resulting binary crashes at runtime;
see Status. The build is documented here so the binary can be rebuilt from
persistent source at any time.

## Source

Source tree: `/hpc/group/abmurraylab/ys451/dhsvm/DHSVM-26`

The source lives in group space, which is snapshotted and persistent. The build
output goes to `/work`, which has no backup and a 75-day auto-purge. Binary
safety rests on the source in group space plus this recipe, not on the built
artifact in `/work`.

## Toolchain

No GCC module is needed. The system compiler is sufficient.

- gcc 11.5.0 at `/usr/bin/gcc`
- cmake via `module load cmake/3.28.3`

## Build

From the DHSVM source directory:

```bash
module load cmake/3.28.3
mkdir -p build && cd build
cmake -DDHSVM_USE_NETCDF=OFF -DDHSVM_USE_X11=OFF ..
make
```

NetCDF and X11 are both off. This is justified. The CA config uses Format=BIN,
so NetCDF output is not needed. DCC compute nodes are headless, so X11 is not
available. Turning both off removes the only optional dependencies and keeps the
build self-contained.

## Status

The binary builds clean. It crashes at runtime about two seconds in, at the same
point each time, just after this log line:

```
No spatial LAI provided, generating from vegetation table
```

The crash is in the DCC build, not in the inputs or config. The same inputs that
crash this binary run clean on the Mac binary. Diagnosis is parked. The next
step is to attach gdb on a compute node and read the backtrace. The login node
denies ptrace, so gdb is unusable there; `srun` to a compute node first.
