# B-Q2 Soil-Depth Sensitivity Report — CA (Camp Branch)

Date: 2026-06-08. Toolkit repo HEAD `72c837c`. Target reference: `soildepthscript.py` (`qgis_CA/`), which produced the current `Intermediate_GIS/soildepth.tif`.

Inputs (grid-matched, from `DEM_CA_0406/Intermediate_GIS/`): `elev_clipped.tif`, `stream_slope.tif`, `flow_acc.tif`. Sweep run with `soildepth_sensitivity.py`, a read-only reimplementation of the soil-depth formula (no `.bin`/`.tif` written, DHSVM not run).

## Verification

The baseline reconstruction (MS30 / POW0.25 / WS0, the current pipeline values) was compared cell-by-cell against the existing `soildepth.tif`: 4334 cells, max|diff| 4.07e-07, mean|diff| 7.91e-08, correlation 1.000000, 0 cells differing by more than 1e-3. The reimplementation reproduces the pipeline output, so every sweep row below is faithful to `soildepthscript.py`.

Reproducibility note (migration item): the repo `qgis_CA/soildepthscript.py` resolves its inputs from root-level candidates (`<WS>/stream_slope.tif`, `<WS>/flow_acc.tif`), but the actual rasters live under `Intermediate_GIS/`. A fresh clone-and-run would not find them. Worth reconciling when the toolkit moves to the standalone CLI.

## Domain

4334 valid cells. Slope mean 24.01 deg, max 44.91 deg, 29.4% above 30 deg, 0% above 45 deg. The reported slope mean of 24.01 deg includes the 318 gap-filled cells whose slope was set to 0; excluding them gives 25.91 deg, which matches the prior handoff figure. Flow accumulation: max 3388, p99 1181, p95 125, p90 35, median 4. The fac distribution is heavy-tailed, so the channel network is a small fraction of the basin.

## Finding 1 — the reported Mean 2.55 m overstates real depth; the deep tail is an artifact

318 cells (7.34% of the domain) are mask-edge cells that are valid in the DEM but have no valid slope. The current gap-fill sets their slope to 0, and since `term_slope = WT_SLOPE * (1 - (slope/MAX_SLOPE)^POW)`, a slope of 0 yields the maximum slope bonus (0.7). These cells therefore receive 4.95–5.23 m depth. All 20 of the deepest cells in the basin are these edge cells, not coves.

Consequences. The artifact inflates the basin mean: the artifact-free central depth (the median) is 2.36 m, against the reported mean of 2.55 m, a 0.19 m overstatement. The p95 (4.87 m) and max (5.23 m) are the artifact ceiling, not real deep soil. Running the sweep with `--fixgap` collapses p95 toward the real terrain maximum (in the demo it drops from ~5.06 to ~2.89, with the median essentially unchanged), confirming the artifact lives entirely in the tail.

These cells are also untunable: a slope of 0 makes `term_slope` insensitive to MAX_SLOPE and POW_SLOPE, so they keep their spurious depth no matter how the parameters are set. Any parameter tuning done on the current field is chasing a distribution whose top 7% is fixed and fake. The gap-fill should be fixed before parameters are chosen.

## Finding 2 — the catena is real but compressed and shallow

Baseline gentle (bottom slope quartile, valley-like) 3.37 m, steep (top slope quartile, ridge-like) 2.15 m, contrast 1.22 m. Valleys are deeper than ridges, so a catena exists, but the absolute depths sit below the Coweeta 3–6 m band and the artifact-free central depth is 2.36 m.

Correction to an earlier framing. The clip at MAX_SLOPE=30 does not pin steep cells at the 2 m floor: only 2.0% of cells reach MIN_DEPTH. Clipping sets `term_slope` to 0 for the 29.4% of cells above 30 deg, but the elevation term still lifts them above 2 m. MAX_SLOPE compresses the slope contribution rather than flooring depth.

## Finding 3 — lever behavior (faithful sweep, WT_SOURCE=0)

Depth in metres. gentle = bottom slope quartile, steep = top slope quartile, contrast = gentle − steep.

| config | mean | gentle | steep | contrast |
|---|---|---|---|---|
| MS30 POW0.25 (current) | 2.55 | 3.37 | 2.15 | 1.22 |
| MS30 POW0.5 | 2.68 | 3.64 | 2.15 | 1.49 |
| MS30 POW1.0 | 2.88 | 4.05 | 2.15 | 1.90 |
| MS45 POW0.25 | 2.76 | 3.53 | 2.32 | 1.21 |
| MS45 POW0.5 | 3.05 | 3.90 | 2.48 | 1.42 |
| MS45 POW1.0 | 3.51 | 4.37 | 2.76 | 1.62 |
| MS60 POW0.25 | 2.91 | 3.63 | 2.50 | 1.13 |
| MS60 POW0.5 | 3.31 | 4.05 | 2.81 | 1.24 |
| MS60 POW1.0 | 3.88 | 4.54 | 3.31 | 1.23 |

POW_SLOPE is the dominant lever for both basin mean and catena contrast. At MS30, raising POW from 0.25 to 1.0 lifts the mean from 2.55 to 2.88 and the contrast from 1.22 to 1.90, by pulling the slope response from near-constant back toward near-linear.

MAX_SLOPE lifts absolute depths, most visibly on ridges: at POW1.0, steep-cell depth rises from 2.15 (MS30) to 2.76 (MS45) to 3.31 (MS60). But it reduces contrast, because widening the clip compresses the normalized slope range that gentle and steep cells map onto. There is a tradeoff: MS30 maximizes contrast but leaves steep cells shallow and flat, while MS60 raises the floor but flattens the gradient.

The deepest slope-only setting is MS60 POW1.0 (mean 3.88, gentle 4.54, steep 3.31, median 3.84). This places the bulk distribution in the lower part of the Coweeta band. No slope-only combination reaches 5–6 m, and `%@6m` is 0 across the entire sweep.

## Finding 4 — the source term is a channel-targeted micro-lever, dead at MAX_SOURCE=100000

WT_SOURCE barely moves the basin mean: at MS45 POW1.0 it changes the mean from 3.51 (WS0) to 3.52 (WS0.3). It only touches high-fac channel-proximal cells, which are a small fraction of the basin. At MS45 POW1.0 it shifts the top fac-decile mean from 3.97 (WS0) to 4.12 (WS0.3) and contrast from 1.62 to 1.66.

MAX_SOURCE controls whether the source term does anything at all. At MS45 POW1.0 WS0.2, the top fac-decile mean moves with the normalization:

| MAX_SOURCE | topfac mean | p95 | %@6m | contrast |
|---|---|---|---|---|
| 100000 (current) | 3.98 | 4.87 | 0.0 | 1.62 |
| max (3388) | 4.07 | 4.88 | 0.0 | 1.65 |
| p99 (1181) | 4.18 | 4.90 | 0.0 | 1.68 |
| p95 (125) | 4.59 | 5.02 | 0.1 | 1.77 |

With MAX_SOURCE at 100000, the source term adds essentially nothing (the largest valley cell, fac 3388, reaches only `0.2 * 3388/100000` of the source range). It becomes meaningful only when MAX_SOURCE is set near the real fac distribution (p95 to max). Even at p95, `%@6m` is only 0.1%. The current sweep range (WT_SOURCE up to 0.3, POW_SOURCE 1.0) nudges the channel cells but does not populate the deep coves. If the target requires strong cove deepening to 5–6 m, a second sweep with higher source weight or lower POW_SOURCE would be needed.

## Finding 5 — de-artifacted results (gap-fill = zero)

Re-running with the gap-fill fix (the 318 cells given no slope bonus) changes the picture in two ways. The deep tail disappears: at the current parameters, p95 drops from 4.87 to 2.78 m and max from 5.23 to 3.31 m. And the catena contrast collapses from 1.22 to 0.40 m, because the 318 slope-zero cells all fall in the bottom slope quartile and were inflating the gentle (valley) mean by 0.82 m. The steep (ridge) mean is unchanged at 2.15 m. Most of the apparent catena in the current field was the edge artifact, not terrain.

The real, de-artifacted field at the current parameters is shallow and nearly flat: median 2.32 m, gentle 2.55 m, steep 2.15 m, contrast 0.40 m. Both the Coweeta depth band and a meaningful ridge-to-valley gradient are absent.

De-artifacted slope-only sweep (WT_SOURCE=0), depth in metres:

| config | mean | gentle | steep | contrast | p50 |
|---|---|---|---|---|---|
| MS30 POW0.25 (current) | 2.34 | 2.55 | 2.15 | 0.40 | 2.32 |
| MS30 POW0.5 | 2.47 | 2.82 | 2.15 | 0.67 | 2.37 |
| MS30 POW1.0 | 2.68 | 3.23 | 2.15 | 1.08 | 2.48 |
| MS45 POW0.25 | 2.55 | 2.71 | 2.32 | 0.39 | 2.56 |
| MS45 POW0.5 | 2.85 | 3.07 | 2.48 | 0.60 | 2.85 |
| MS45 POW1.0 | 3.30 | 3.55 | 2.76 | 0.79 | 3.33 |
| MS60 POW0.25 | 2.71 | 2.81 | 2.50 | 0.31 | 2.74 |
| MS60 POW0.5 | 3.11 | 3.23 | 2.81 | 0.42 | 3.14 |
| MS60 POW1.0 | 3.67 | 3.71 | 3.31 | 0.41 | 3.76 |

The depth-versus-contrast tradeoff is now clean. At POW1.0, raising MAX_SLOPE lifts the mean (2.68 to 3.30 to 3.67) but flattens the catena (contrast 1.08 to 0.79 to 0.41), because the steep cells deepen faster than the gentle ones as the clip widens. The slope and elevation terms alone cannot produce both deep absolute soil and a strong catena.

The source term is the lever that can break that tradeoff, since it deepens high-fac valley cells without touching ridges. De-artifacted, at MS45 POW1.0 WS0.2, lowering MAX_SOURCE from 100000 to p95 raises contrast from 0.80 to 0.95 and the channel-cell mean from 3.95 to 4.56. Within the current sweep range (WT_SOURCE up to 0.3) the catena stays moderate; reaching the deep coves would need a wider source sweep.

## Recommended ordering

1. Fix the gap-fill first. Decide between excluding the 318 edge cells from the domain or keeping them with no slope bonus, then re-run the sweep with `--fixgap` so parameters are chosen on the de-artifacted field.
2. Set the slope parameters (MAX_SLOPE, POW_SLOPE) to the depth target.
3. If coves need to reach the upper Coweeta band, tune the source term, with MAX_SOURCE set to the real fac distribution rather than 100000.

## Decision points

1. Gap-fill fix: `exclude` (drop the edge cells) vs `zero` (keep them, no slope bonus). This depends on whether the 318 cells are genuine terrain that happens to lack a slope value or are mask padding at the boundary. Worth checking their spatial location before deciding.
2. Target framing: match a basin-mean depth in the middle of 3–6 m, or maximize the catena gradient (thin ridges, thick valleys). These point to different parameter sets, since contrast and absolute depth trade off against each other under MAX_SLOPE.
3. Source-term normalization: raw max, p99, or p95 for MAX_SOURCE, and whether to widen the WT_SOURCE / POW_SOURCE sweep, if deep coves are required.
4. Weight sum when WT_SOURCE > 0: keep the additive 0.7 / 0.3 weights (sum exceeds 1, valleys can clip at 6 m) or renormalize so the three weights sum to 1 (`--renorm`, default off, matches current code).

## Candidate parameter sets (de-artifacted field)

These assume the gap-fill is fixed. Choose by target.

- Match basin-mean depth, catena secondary: MS60 POW1.0, mean 3.67 (gentle 3.71, steep 3.31, contrast 0.41). The deepest single-knob option, but nearly flat.
- Balance depth and catena: MS45 POW1.0, mean 3.30 (gentle 3.55, steep 2.76, contrast 0.79). Adding WT_SOURCE 0.2 with MAX_SOURCE at p95 raises it to mean 3.39, contrast 0.95, channel cells 4.56.
- Maximum catena from slope alone: MS30 POW1.0, contrast 1.08 (gentle 3.23, steep 2.15), but mean only 2.68.
- Strong catena into the upper Coweeta band (5–6 m coves): none of the above reach it. A second sweep with WT_SOURCE up to roughly 0.5–0.6 and MAX_SOURCE at p95–p99 is the next experiment.

## Next step

Pick the gap-fill fix and a parameter direction. The code change then covers the gap-fill fix and the parameter update in `qgis_CA/soildepthscript.py`, ported identically to `qgis/` and `standalone/` (the formula is identical across all three), followed by a byte-level regression of the new `soildepth.tif` against the rerun.
