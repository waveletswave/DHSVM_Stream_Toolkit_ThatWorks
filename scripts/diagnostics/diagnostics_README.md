# Diagnostics

Read-only analysis and verification tools. None of these are part of the DHSVM input pipeline; they support correctness audits and parameter studies. All read rasters and print results (some also write a CSV or a `*_filled.tif`); none run DHSVM.

- `diagnose_slope_gaps.py` — locate and characterize the missing-slope cells: the clipped-DEM boundary ring vs cells ringing former interior DEM voids, with patch count. The diagnosis that identified the source of the soil-depth boundary artifact.
- `fix_slope_nodata.py` — standalone slope NoData fill (neighbour interpolation) plus a soil-depth before/after report. The fill logic is now in `qgis_CA/slope_fill.py` (used by the pipeline); this is kept for the before/after diagnostic and audit trail.
- `compare_soildepth.py` — compare two `soildepth.tif` (old vs new) cell by cell, with a cross-check of the changed cells against the old slope NoData mask.
- `check_slope_consumers.py` — whether any stream cell overlaps the conditioned cells, i.e. whether the channel files (`stream.class.dat`, `stream.network.dat`) could change.
- `soildepth_sensitivity.py` — B-Q2 soil-depth parameter sweep over MAX_SLOPE, POW_SLOPE, WT_SOURCE. `--verify` checks the formula against an existing `soildepth.tif`; `--fixgap` shows the distribution with the boundary artifact removed.

Most accept `--help`. They expect the CA-style `Intermediate_GIS` layout, or explicit `--slope` / `--dem` / `--fac` paths.
