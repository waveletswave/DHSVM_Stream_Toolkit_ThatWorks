# Drop analysis: Strahler order cached across the threshold sweep (2026-09-23)

## Finding

`diagnostics/drop_analysis.py` sweeps support-area thresholds and, for each,
assigns Strahler order to the extracted streams and tests first-order drops
against higher-order drops (Broscoe 1959; Tarboton et al. 1991). It called
`pyflwdir.FlwdirRaster.stream_order(type="strahler", mask=...)` once per
threshold on one flow object. That method (pyflwdir 0.5.5, 2022, through
0.5.12, 2026-07; `flwdir.py`) caches the order map on the object and returns
the cached map on every later call, whatever mask is passed. Every threshold
after the first therefore used the Strahler orders of the first, densest
network (10 cells). The t-test compared groups defined on the wrong network.

Found on 2026-09-23 while porting the analysis into WW-DHSVM. Signature in
every published sweep table: maximum Strahler order 4 at all thresholds,
including 300 cells, where only 7 streams exist (order 4 needs at least 8
first-order streams).

## Verification

Cloud environment, pyflwdir 0.5.12, on the CA 28 m fixture
(`standalone_CA/tests/fixtures/CA_28m/elev_clipped.tif`, the same DEM as the
June analysis). The unmodified script reproduces the documented table
exactly (objective 50 cells, band 50 to 100, max order 4 throughout). With
the order recomputed per threshold the same script gives:

| cells | km2 | streams | order 1 | higher | max order | drainage density | abs t | pass |
|---|---|---|---|---|---|---|---|---|
| 10 | 0.0079 | 191 | 101 | 90 | 4 | 6.97 | 8.01 | no |
| 20 | 0.0159 | 80 | 42 | 38 | 4 | 4.05 | 6.11 | no |
| 30 | 0.0238 | 59 | 31 | 28 | 4 | 3.29 | 4.82 | no |
| 40 | 0.0317 | 48 | 25 | 23 | 4 | 2.71 | 4.14 | no |
| 50 | 0.0396 | 36 | 19 | 17 | 4 | 2.40 | 3.43 | no |
| 60 | 0.0476 | 30 | 16 | 14 | 3 | 2.18 | 2.80 | no |
| 70 | 0.0555 | 25 | 13 | 12 | 3 | 1.98 | 2.66 | no |
| 80 | 0.0634 | 21 | 11 | 10 | 3 | 1.84 | 2.47 | no |
| 90 | 0.0714 | 19 | 10 | 9 | 3 | 1.75 | 2.02 | no |
| 100 | 0.0793 | 17 | 9 | 8 | 3 | 1.66 | 2.20 | no |
| 110 | 0.0872 | 17 | 9 | 8 | 3 | 1.61 | 2.07 | no |
| 120 | 0.0951 | 17 | 9 | 8 | 3 | 1.57 | 1.92 | yes |
| 150 | 0.1189 | 15 | 8 | 7 | 3 | 1.39 | 1.33 | yes |
| 200 | 0.1586 | 9 | 5 | 4 | 3 | 1.18 | 1.10 | yes |
| 300 | 0.2379 | 7 | 4 | 3 | 2 | 0.95 | 1.10 | yes |
| 600 | 0.4757 | 5 | 3 | 2 | 2 | 0.66 | 0.46 | yes |

Objective threshold (start of the first band of three consecutive passes):
120 cells = 95143 m2 = 0.0951 km2; every threshold from 120 to 600 cells
passes. abs t falls monotonically from 8 at 10 cells to 1.9 at 120 cells,
where the network has 17 streams (9 first order, 8 higher), so the test has
little power there and the pass is in part a small-sample effect: on a
4334-cell basin the constant-drop criterion can say where the law is
violated (below about 100 cells) more firmly than where it holds.

## Consequences

- `standalone_CA/docs/stream_threshold.md` stated an objective of 50 cells
  with a passing band of 50 to 100 and placed the 60-cell default at 0.83 x
  the objective. Corrected: the objective is 120 cells; the 60-cell default
  is 0.5 x the objective and does not satisfy the law on CA at 28 m
  (abs t 2.80).
- `standalone_CA/docs/validation_log.md`: the CA 10 m objective (500 cells,
  band 0.050 to 0.10 km2), the AR 10 m objective (220 cells) and the AR
  30 m objective (100 cells), and the resolution-portability reading built
  on them, were computed with the cached orders. They are withdrawn, not
  redone: the 10 m DEMs were in scratch space and are gone, and the check
  is not needed for the manuscript basins (decision 2026-09-23).
- The June AR portability runs used those objectives as their A_c. They
  are validation exercises of the general entry point, not manuscript
  inputs, and their pipeline outputs remain valid input sets at a
  threshold that is simply no longer called objective.
- Not affected: the manuscript networks for CA and AR (both use the 60-cell
  default, 47571.5 m2, a modelling choice) and the Tier E audit (network
  topology at a fixed threshold).

## Fix

`drop_analysis.py` now computes the order with
`pyflwdir.streams.strahler_order(flw.idxs_ds, flw.idxs_seq, mask=...)` for
every threshold (`strahler_order()` helper); nothing else in the method
changed. `standalone_CA/tests/test_drop_analysis.py` pins the corrected
behaviour on the CA fixture: the maximum order falls along the sweep (4, 3,
3, 2 at 10, 60, 120, 300 cells) and the objective is 120 cells. The
[dev] extras now include pyflwdir and scipy so the test runs in CI.

## Default support area

The default stays at 60 cells (47571.5 m2 at 28.158 m) for CA and AR. It is
the value the manuscript runs and the Tier E audit used, and it is now
described as what it always was, a visual choice, with the corrected
objective and the small-sample caveat beside it. Moving A_c is a separate
decision: it changes the channel cell set and therefore the channel
interception term, and belongs with the manuscript's input-set decision.

## Upstream

Filed as https://github.com/Deltares/pyflwdir/issues/123 (2026-09-23)
with a two-mask reproducer (call `stream_order` with a dense mask, then a
sparse one; the second call returns the first map) and the one-line fix.
WW-DHSVM is not affected: it computes Strahler order itself.
