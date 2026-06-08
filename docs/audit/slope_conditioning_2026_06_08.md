# Slope NoData conditioning (gap-fill) fix — 2026-06-08

## Context

This follows the Tier B unit fix (2026-05-12). During the B-Q2 soil-depth sensitivity analysis on the CA (Camp Branch) watershed, the sensitivity sweep surfaced a residual correctness problem in the soil-depth field that is the same character as the Tier B bugs: silent, no runtime error, and it distorts a model input. It is a slope-raster completeness issue, not a parameter-tuning issue, so it is recorded here as a correctness fix rather than as part of the deferred B-Q1/B-Q2 parameter work.

## Issue

The CA `soildepth.tif` carried a ring of spuriously deep cells. 318 cells (7.34% of the 4334-cell basin) reached 4.95 m mean, 5.23 m max, while the surrounding terrain-driven depth sat near 2.3 m. All of the 20 deepest cells in the basin were these. The artifact inflated the reported basin-mean depth (2.55 m) above the artifact-free central value (median 2.36 m), and the entire upper tail of the depth distribution (p95 4.87 m, max 5.23 m) was artifact rather than real cove soil. Because the affected cells have slope identically 0, they are insensitive to all soil-depth parameters, so any parameter tuning done on the un-conditioned field would chase a distribution whose top 7% is fixed and unphysical.

## Root cause

`stream_slope.tif` is produced by GRASS `r.slope.aspect`. Two things leave it NoData on cells that are valid in the DEM:

1. The outermost ring of the clipped watershed has no full 3x3 neighbourhood, so `r.slope.aspect` returns NoData there (291 cells).
2. A few small interior DEM voids existed when slope was computed; `r.slope.aspect` returned NoData on the cells ringing them, and the voids were later filled in the DEM, so those NoData cells now appear interior rather than on a boundary (27 cells, in 4 small patches).

`soildepthscript.py` then filled the missing slope with 0: `slope_clean = np.where(valid & slope_valid & slope>=0, slope, 0.0)`. In the PNNL depth formula `term_slope = WT_SLOPE * (1 - (slope/MAX_SLOPE)**POW_SLOPE)`, a slope of 0 yields the maximum slope bonus (WT_SLOPE = 0.7), i.e. the deepest soil. The gap-fill conflated "slope unknown" with "slope zero", so edge cells received the maximum depth instead of a value consistent with their surroundings.

## Fix

`slope_fill.py` conditions the slope raster by neighbour interpolation: each missing-slope cell (DEM-valid, slope-NoData) is filled with the mean of its valid 8-neighbours, iterated inward for the few multi-cell patches. This recovers a physically reasonable slope (the boundary cells are real terrain whose slope simply was not computable at the edge) and leaves the depth formula and its parameters untouched. It is pure NumPy + GDAL with no QGIS dependency, so the same module serves the QGIS pipeline and the planned standalone/cluster pipeline.

It is folded into `prep_dhsvm_inputs.py` as a conditioning step immediately after the `r.slope.aspect` block and before any consumer reads the slope, so `stream_slope.tif` is complete on disk for every downstream step (channel slope sampling, stream network write, soil depth). The step reports the number of cells filled.

## Validation

On the CA data the conditioning filled all 318 missing-slope cells (291 boundary, 27 interior), 0 unfilled, with filled-cell slope mean 24.64 deg against an interior mean of 25.91 deg. The reported slope mean of 24.01 deg before conditioning was the 4334-cell average pulled down by the 318 zeros; excluding them gives 25.91 deg, which matches the Tier B post-fix figure.

Regression on the soil-depth output (`compare_soildepth.py`, old vs new `soildepth.tif`): 318 cells changed, 4016 cells byte-identical, max change confined to the changed cells (old depth mean 4.95 m / max 5.23 m, new mean 2.32 m / max 3.02 m). Cross-check against the old slope's NoData mask: every changed cell was a conditioned cell and only those changed (0 changed cells outside the conditioned set). So the change is exactly the conditioning, with no side effects elsewhere.

## Scope (what changed, what did not)

Slope feeds soil depth and the channel files. Everything else is independent of slope.

- `soildepth.bin` / `soildepth.tif`: changed at the 318 cells (confirmed above).
- `stream.class.dat`: unchanged. All 41 segments are class 13 (steep band, smallest area bin) both before and after; the classification is robust to the slope change at this scale.
- `stream.network.dat`: the column-3 tan(slope) differs against the pre-Tier-B backup, but that difference is the Tier B percent-to-degrees unit fix, not this conditioning. The old/new ratio is about 1.9 across all segments, which is tan(percent-value-misread-as-degrees) / tan(true-degrees), the signature of the unit fix. The conditioning itself touches only the channel sampling for at most 7 stream-adjacent cells, and Tier A showed channel routing has near-zero effect on CA streamflow, so the conditioning's effect on the channel files is negligible.
- `Channel.State`: downstream of the channel files (class widths x length); unchanged because `stream.class.dat` is unchanged.
- `dem.bin`, `mask.bin`, `soil.bin`, `veg.bin`, `soildepth_uniform_*.bin`, `veg_bs*.bin`, `stream.map.dat`, the grid state files, and the `[AREA]` config: no slope dependency, not affected.

## Relation to the audit and to B-Q2

This completes the slope-correctness part of Tier B alongside the `format:1 -> 0` unit fix. The B-Q2 soil-depth decision (uniform vs variable depth) is a separate, still-open question and is independent of this fix: the conditioning makes the slope raster correct regardless of which depth approach is chosen, and the variable-depth field (if used) is now free of the edge artifact.

## Files

- Fix: `slope_fill.py` (folded into `prep_dhsvm_inputs.py`).
- Diagnostics used (kept under `scripts/diagnostics/`): `fix_slope_nodata.py` (standalone fill plus before/after), `soildepth_sensitivity.py` (B-Q2 parameter sweep), `compare_soildepth.py` (old/new soil-depth regression), `check_slope_consumers.py` (stream-cell overlap with conditioned cells).
- Companion report: `docs/audit/bq2_soildepth_sensitivity_2026_06_08.md`.
