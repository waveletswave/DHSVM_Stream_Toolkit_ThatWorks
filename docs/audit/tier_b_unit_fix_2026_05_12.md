# Tier B Unit-Fix Audit Log

**Date**: 2026-05-12
**Scope**: QGIS DHSVM preprocessing pipeline (`qgis_CA/`)
**Mainline commits**:
- `71d667e` — GRASS slope format 1→0 (percent→degrees)
- `4b70f9e` — Consumer A sample() unpacking fix

## Context

Tier B was originally scoped as a pure numerical-reasonableness review
of channel hydraulic geometry (Q1) and soildepthscript weights (Q2),
not intended to modify code. Two silent unit-correctness bugs were
uncovered during the audit and required correction before any
parameter tuning could be meaningfully discussed. Same character as
Tier A: silent structural errors, runtime no warning. Q1/Q2 parameter
tuning is deferred to a separate follow-up commit.

## Bug 1: GRASS r.slope.aspect format parameter

The pipeline section header read "Slope raster (DEGREES)" but the
algorithm was invoked with `format=1`, which produces percent. The
RAT label on the resulting `stream_slope.tif` confirmed this with
explicit per-class labels "N percent" / "zero slope".

Verification: `algorithmHelp("grass7:r.slope.aspect")` documents
`0=degrees, 1=percent`. An independent reference TIF generated with
`format=0` had Min 4.34, Mean 25.91, Max 44.91 degrees — an exact
arctan match to the pre-fix percent values (Min 7.58, Mean 49.64,
Max 99.70).

Three downstream consumers depended on the slope raster, each with
its own assumption about the unit:
- `_sample_mean_slope_deg` (line 316): wrote raw values into a field
  named `slope_deg`. Misinformed `channelclass.py` heuristic.
- `_sample_mean_slope_tan_along_line` (line 543): applied
  `tan(radians(val))` treating percent as degrees, producing
  tan values 2-4x too large in `stream.network.dat` column 3.
- `soildepthscript.generate_soildepth`: applied `MAX_SLOPE=30` clip
  with degrees intent but received percent; ~80% of cells exceeded
  the clip threshold, collapsing simulated catena to mean 2.43m,
  std 0.72m (vs Coweeta saprolite ~3-6m, Swank & Douglass 1975).

Fix: changed `'format':1` to `'format':0` in `prep_dhsvm_inputs.py`
line 250 (commit `71d667e`).

## Bug 2: Consumer A sample() unpacking order

`_sample_mean_slope_deg` unpacked the result of
`QgsRasterDataProvider.sample()` as `ok,val`, but the QGIS API
documents the return order as `(value, ok_bool)` — see PyQGIS
Cookbook "Using Raster Layers" section. Consumer B
(`_sample_mean_slope_tan_along_line`, line 547) had the correct
`val,ok` order, confirming this was a copy-paste typo rather than
API ambiguity.

Pre-fix consequence: every segment received `float(True)=1.0` as
its `slope_deg` field, completely bypassing real terrain sampling.
The `channelclass.py` heuristic (which selects degrees/percent/tan
mode based on the sample max) then mistakenly selected "tan" mode
because the max was 1.0. On CA — a genuinely all-steep basin — this
bug was empirically silent because the heuristic still routed every
segment to class 13 via a different logical path (tan_s=1.0 → steep
band → class 13). For any non-all-steep basin, this bug would have
caused silent misclassification.

Fix: swapped `ok,val` → `val,ok` in `prep_dhsvm_inputs.py` line 324
(commit `4b70f9e`).

## Post-fix file-level changes (CA basin, full pipeline rerun)

Three input files modified:

| File | Pre-fix | Post-fix |
|---|---|---|
| `Intermediate_GIS/stream_slope.tif` | percent, Mean 49.64, Max 99.70 | degrees, Mean 25.91, Max 44.91 |
| `Intermediate_GIS/soildepth.tif` | Mean 2.425m, StdDev 0.724 | Mean 2.549m, StdDev 0.714 (+5%) |
| `DHSVM_input_streams/stream.network.dat` col 3 | unphysical max 2.68 (theta=69.5 deg) | physical max ~0.70 |

Files byte-identical post-fix (md5 verified across 28-file manifest):
`stream.class.dat`, `stream.map.dat`, all 13 `DHSVM_input_binaries`,
all 6 `modelstate`, `elev_clipped.tif`, `stream_raster.tif`,
`flow_dir.tif`, `flow_acc.tif`. Total changed: 3 files, all on the
slope-derivation path.

## DHSVM regression (LAI20 baseline, three-year simulation)

Three-way comparison via `compare_three_runs.py`:

| Run | Tag | Stage | Total Q (mm) | Mean Q (mm/d) | Peak Q (mm/d) |
|---|---|---|---|---|---|
| CA | CA_LAI20 | pre-Tier-A | 3879.26 | 3.5395 | 40.14 |
| 508 | 508_LAI20 | post-Tier-A | 3879.26 | 3.5395 | 40.14 |
| 512 | 512_LAI20 | post-Tier-B | 3855.22 | 3.5175 | 39.93 |

**Tier A re-verified**: r(CA, 508) = 1.000000, vol delta = +0.000%.
Second independent confirmation of Tier A byte-identical streamflow
on CA.

**Tier B effect**: vol delta = -0.620%, Pearson r(508, 512) = 0.99951.
A 24 mm streamflow reduction over 3 years (3879 → 3855 mm), absorbed
by increased soil storage (soildepth +5% mean) rather than ET.
Aggregated.Values P / ET / T / I unchanged at 0.01 mm precision;
the storage absorption is verifiable via Mass.Balance, not parsed
in compare_three_runs.py.

**Spatial-temporal pattern**: hydrograph difference concentrated in
year 1 (visible spin-up artifact as deeper soil column equilibrates,
peak negative ΔQ ~-2 mm/d in 2016 spring); years 2-3 near-steady.
Implication for future runs: first year should be treated as
spin-up and excluded from calibration metrics.

**Additivity**: cumulative (-0.620%) = Tier A (+0.000%) + Tier B
(-0.620%) to round-off precision. Two fixes are independent in
effect, no interaction.

## Pending follow-up (Tier B Q1 / Q2, separate commit)

Parameter tuning questions deferred because they require judgment
calls now unblocked by unit correctness:

1. **`MAX_SLOPE = 30 deg` in soildepthscript.py** still clips ~25%
   of CA cells (mean slope 25.9°, range up to 45°), suppressing
   catena even post-fix. Candidate revision to 45° or 60°.
2. **`WT_SOURCE = 0.0`** disables flow-accumulation contribution to
   soil depth. Now that flow_acc sign convention and MAX_SOURCE
   magnitude mismatch (~30x) are documented, an informed decision
   on whether to enable this term is possible.
3. **`channelclass.py CLASS_TABLE area_edges`** at `[1e6, 1e7, 2e7,
   3e7, 4e7]` m² designed for >10 km² basins; CA at 3.4 km² only
   exercises bins 0-1 of 6. Candidate to rescale for CA-class
   headwater.

## Cross-cutting flags noted but not addressed

- `DHSVM_input_binaries/{veg,mask,soil}.bin` share identical md5
  (same byte content), suggesting single-class basin assumption.
  Confirm whether spatially distributed soil/veg classification is
  intended for future LAI scaling experiments.
- `flow_acc.tif` GRASS convention produces all-negative values
  (drainage from outside DEM); soildepthscript already handles via
  `np.abs()` but `MAX_SOURCE=100000` is ~30x larger than actual
  |max|=3388.
- soildepth.tif Max=5.23m is a boundary artifact from
  `np.where(slope_arr is valid, slope_arr, 0.0)` gap-fill, not a
  real cove signal. Real cove cells suppressed by MAX_SLOPE clip.
