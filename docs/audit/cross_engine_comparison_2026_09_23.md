# Cross-engine comparison: the Toolkit's stream network against WW-DHSVM's (2026-09-23)

**Date**: 2026-09-23
**Scope**: the stream network stage, compared with an independent implementation on the same grid, DEM and support area, and the effect of the two networks on the manuscript's DHSVM runs
**Branch**: `feat-cross-engine`
**Companion**: WW-DHSVM branch `feat-channel-initiation` (github.com/waveletswave/WW-DHSVM, fork of UConn-EFC/WW-DHSVM at `19a6987`), which adds the constant-drop analysis, the outlet invariants and a bring-your-own-DEM entry to Zhi Li's package

## Question

Tier E rebuilt the Toolkit's network stage from the raster (`tier_e_network_orientation_2026_09_22.md`). Its invariants say the network is well formed; they do not say it is the network another method would draw. WW-DHSVM builds channel networks for DHSVM by a different route (pyflwdir priority-flood fill and D8 on the basin mask, its own raster-native segment builder, a basin-scaled 18-class table), so running both engines on identical inputs measures how much of the network is method and how much is terrain, and the DHSVM reruns measure whether that difference reaches the manuscript's numbers.

## Method

`scripts/diagnostics/cross_engine/compare_engines.py` runs WW-DHSVM's terrain and network stages on the Toolkit's grid (the `elev_clipped.tif` of the Toolkit run defines rows, columns, cell size, CRS and the basin mask; the elevations come from the same reprojected 3DEP tile the clip was cut from) at the Toolkit's support area, A_c = 47571.5 m2 = 60 cells at 28.158 m, with no stream burning in either engine. It writes WW-DHSVM's three stream files and a `Channel.State` in the Toolkit's convention (width x length x 0.05 m), rasters in the Toolkit's conventions so `check_network_orientation.py` can read the WW-DHSVM network, and a comparison at the level DHSVM sees: channel cells, segments, outlets, lengths, ranks, classes and the outlet invariants. It also runs the corrected drop analysis (`drop_analysis_strahler_cache_2026_09_23.md`) as ported to WW-DHSVM.

Inputs (sha256, first 16 hex): CA, the regression fixture `standalone_CA/tests/fixtures/CA_28m/` (elev_clipped 181539c95da679c9, stream files 26983e5342993a03 / 6a01bd2f6c1d6404 / a4f0f02ea78d8ab5, segments.csv d283299e3bb00372); AR, the fixed AR 28 m run (elev_clipped a0e570a8e6ee6ca2, stream files 26983e5342993a03 / 9d5fdca8f1222bbf / 94cf497545aa1fad, segments.csv 8a2d8025e0c170f0); tile `USGS_1_n36w084_20220725_UTM17.tif` ad34883c05cb0c3d for both. Outputs are under `docs/audit/cross_engine_2026_09_23/`.

The DHSVM side reuses the Tier E rerun scripts (`scripts/diagnostics/tierE_dhsvm/`), which gained a third network kind, `ww`, prefix `wA_`: the April manuscript inputs with only the three stream files and the `Channel.State` replaced by the WW-DHSVM ones, the grid states copied unchanged from the April set (sha256-checked). Runs: CA `wA_S4h`, `wA_LAI70`, `wA_LAI20` (base `apr`) and AR `wA_S4h`, DHSVM3.2 of 2025-01-13 on the Mac, 2026-09-23.

## Networks

CA, 74 x 82, 4334 basin cells (3.436 km2):

| | Toolkit (GRASS MFD, r.stream.extract, Tier E) | WW-DHSVM (pyflwdir D8) |
|---|---|---|
| channel cells / segments | 231 / 22 | 266 / 30 |
| map records, cells under two segments | 231, 0 | 266, 0 |
| total length (m), mean segment (m) | 7600.8, 345.5 | 8901.2, 296.7 |
| drainage density (km/km2) | 2.212 | 2.590 |
| outlet segments (tail cell, z) | 12 at (67, 68) 829.61 m and 22 at (65, 68) 838.0 m | 29 at (68, 68) 826.81 m |
| max routing rank, max in-degree | 8, 2 | 10, 3 |
| classes used | 13 x 17, 14 x 5 | 6 x 1, 13 x 11, 14 x 5, 15 x 3, 16 x 4, 17 x 2, 18 x 4 |
| max-|acc| cell (cells), in an outlet segment | (63, 67), 3387.5, yes | (68, 68), 4334, yes |
| lowest channel cell (z), in an outlet segment | (67, 68), 829.61 m, yes | (68, 68), 826.81 m, yes |

AR, 55 x 72, 2870 basin cells (2.276 km2):

| | Toolkit | WW-DHSVM |
|---|---|---|
| channel cells / segments | 170 / 19 | 203 / 29 |
| map records, cells under two segments | 170, 0 | 203, 0 |
| total length (m), mean segment (m) | 5498.3, 289.4 | 6602.4, 227.7 |
| drainage density (km/km2) | 2.416 | 2.902 |
| outlet segments (tail cell, z) | 19 at (53, 41) 796.61 m | 29 at (53, 41) 796.61 m |
| max routing rank, max in-degree | 7, 2 | 8, 2 |
| classes used | 13 x 16, 14 x 3 | 13 x 12, 14 x 4, 15 x 4, 16 x 5, 17 x 1, 18 x 3 |
| max-|acc| cell (cells), in an outlet segment | (53, 41), 2602.6, yes | (53, 41), 2870, yes |
| lowest channel cell (z), in an outlet segment | (53, 41), 796.61 m, yes | same |

Channel-cell overlap, by the Toolkit's MFD |acc| band (common of Toolkit / WW-DHSVM cells in the band):

| |acc| band (cells) | CA | AR |
|---|---|---|
| 0 to 60 | 1 of 5 / 77 | 2 of 3 / 36 |
| 60 to 120 | 31 of 56 / 52 | 45 of 56 / 57 |
| 120 to 300 | 37 of 73 / 46 | 49 of 55 / 55 |
| 300 to 1000 | 36 of 44 / 43 | 27 of 28 / 27 |
| 1000 and above | 48 of 53 / 48 | 28 of 28 / 28 |
| all | 153 common, 78 Toolkit only, 113 WW only, Jaccard 0.445 | 151 common, 19 Toolkit only, 52 WW only, Jaccard 0.68 |

Reading:

- The two engines agree on the trunk and disagree in the headwaters. Above 300 cells of MFD accumulation 84 of 97 Toolkit cells on CA and 55 of 56 on AR are WW-DHSVM channel cells; the WW-only cells are mostly below 60 cells of MFD accumulation (76 of 113 on CA, 34 of 52 on AR), which is where D8 concentrates flow that MFD spreads. At the same nominal A_c, D8 initiates channels higher on the hillslope, so the WW-DHSVM networks have 15 to 19 percent more channel cells, 17 to 20 percent more length, and more, shorter segments.
- The CA mouth. The Toolkit ends in two outlet segments at (67, 68) and (65, 68) because MFD accumulation is not monotone along the extracted stream (the earlier record); WW-DHSVM ends in one outlet at (68, 68), the cell that holds the whole basin's D8 accumulation (4334) and the lowest cell in the mask. The Toolkit's max-|acc| cell (63, 67) is a WW-DHSVM channel cell; WW-DHSVM's outlet cell (68, 68) is not a Toolkit channel cell. On AR the two engines drain to the same outlet cell.
- Both WW-DHSVM networks pass WW-DHSVM's own topology and outlet checks (no errors, no warnings, no depression-fill change to the raw DEM) and the Toolkit's read-only diagnostic (`docs/audit/cross_engine_2026_09_23/*/orientation_*_ww.txt`: I1, I2 and I5 PASS, 0 hard invariants failed).
- The class tables differ by design: the Toolkit assigns classes 13 and 14 (width 0.5 and 1.0 m, cut depth 0.10 m, Manning n 0.045); WW-DHSVM scales its 18-class table to the basin (widths 0.72 to 3.57 m on CA, 0.74 to 3.01 m on AR, cut depths 0.15 to 0.35 m, n 0.10; one CA segment in class 6, n 0.03). The DHSVM reruns therefore carry both the network and the class table.
- Drop analysis, as ported (Strahler order recomputed per threshold): CA objective 120 cells (0.0951 km2), band 120 to 600, the same as the Toolkit's corrected objective on the same DEM; AR has isolated passes at 60 to 70 and 160 to 170 cells and its first sustained band starts at 320 cells with 8 streams left, which is not a usable objective, as with the Toolkit on AR. One arithmetic note: the sweep thresholds the upstream cell count with `>=`, while `extractNetwork` thresholds the upstream area in m2, so 47571.5 m2 excludes the three CA cells whose count is exactly 60 (269 stream cells in the sweep row, 266 in the network).

## DHSVM: the WW-DHSVM network on the manuscript inputs

Mass.Final.Balance (mm; `oA` control on the April inputs, `nA` Tier E network, `wA` WW-DHSVM network):

| run | ET | ChannelInt | Initial Storage | Final Storage | Mass Error |
|---|---|---|---|---|---|
| CA S4h oA / nA / wA | 2808.941 / 2808.850 / 2807.320 | 2210.152 / 2210.246 / 2211.787 | 606.150 / 606.140 / 605.944 | 1021.259 / 1021.235 / 1021.054 | -0.071 / -0.081 / -0.056 |
| CA LAI70 oA / nA / wA | 2325.027 / 2324.950 / 2323.671 | 2686.421 / 2686.508 / 2687.797 | 606.150 / 606.140 / 605.944 | 1029.104 / 1029.081 / 1028.880 | -0.052 / -0.056 / -0.051 |
| CA LAI20 oA / nA / wA | 1119.144 / 1119.112 / 1118.575 | 3878.913 / 3878.942 / 3879.535 | 606.150 / 606.140 / 605.944 | 1042.750 / 1042.733 / 1042.497 | -0.032 / -0.041 / -0.025 |
| AR S4h oA / nA / wA | 2875.769 / 2875.726 / 2874.328 | 2124.688 / 2124.803 / 2126.427 | 651.039 / 651.034 / 650.867 | 1079.481 / 1079.409 / 1079.018 | 0.136 / 0.140 / 0.143 |

Channel lateral inflow (the manuscript's simulated Q) and the routed outflow, Tier E versus WW-DHSVM:

| run | Q, three years | Q, daily | routed outflow, hourly | R12 on wA |
|---|---|---|---|---|
| CA S4h | +1.532 mm (+0.069%) | max 0.284 mm/d, RMSE 0.032 mm/d, r 0.99991 | +1.599 mm; max step 434.7 m3 (6.2% of the peak), r 0.99889 | one SAVE outlet, sum minus total 2e-35 m3/step |
| CA LAI70 | +1.296 mm (+0.048%) | max 0.320, RMSE 0.033, r 0.99992 | +1.358 mm; max step 410.3 m3 (5.6%), r 0.99909 | 1e-34 |
| CA LAI20 | +0.583 mm (+0.015%) | max 0.325, RMSE 0.035, r 0.99993 | +0.642 mm; max step 450.5 m3 (5.5%), r 0.99923 | 1e-34 |
| AR S4h | +1.627 mm (+0.077%) | max 0.443, RMSE 0.046, r 0.99975 | +1.712 mm; max step 349.2 m3 (8.1%), r 0.99785 | 1e-34 |

Against the April control the numbers are the same to the second decimal (CA S4h +1.632 mm, LAI70 +1.377, LAI20 +0.616, AR +1.735). The cumulative channel lateral inflow and the routed outflow agree to 0.01% in every wA run (ratio 1.0000, AR 1.0001).

Metrics with the manuscript's own evaluation code:

| | NSE | r | RMSE (mm/d) | PBIAS (%) |
|---|---|---|---|---|
| CA stitched all, oA / nA / wA | 0.681 / 0.682 / 0.681 | 0.841 / 0.841 / 0.841 | 2.310 / 2.309 / 2.311 | -0.3 / -0.3 / -0.3 |
| CA stitched 2017 | 0.507 / 0.509 / 0.504 | 0.750 / 0.751 / 0.750 | 1.094 / 1.091 / 1.097 | -0.4 / -0.4 / -0.4 |
| CA stitched 2018 | 0.691 / 0.691 / 0.691 | 0.858 / 0.858 / 0.858 | 3.021 / 3.020 / 3.021 | -0.2 / -0.2 / -0.2 |
| AR overall, S4h / nA / wA | 0.600 / 0.601 / 0.596 | 0.777 / 0.777 / 0.774 | 1.964 / 1.962 / 1.975 | -3.5 / -3.5 / -3.5 |
| AR 2017 | 0.298 / 0.300 / 0.294 | 0.627 / 0.627 / 0.623 | 0.916 / 0.915 / 0.919 | +14.1 / +14.1 / +14.2 |
| AR 2018 | 0.521 / 0.522 / 0.515 | 0.728 / 0.729 / 0.724 | 2.573 / 2.570 / 2.588 | -8.5 / -8.5 / -8.5 |

The per-year optimal leaf-area multipliers of the CA stitch (LAI20 for 2017, LAI70 for 2018) do not move.

Reading:

- Swapping the whole network engine, 15 to 19 percent more channel cells with wider and deeper channels, moves the manuscript's three-year Q by 0.015 to 0.077 percent and its metrics by at most 0.007 in NSE, 0.018 mm/d in RMSE and 0.1 point in PBIAS. The extra 0.6 to 1.6 mm of channel inflow is matched by a fall in ET of nearly the same amount. The network is not a source of uncertainty at the level the manuscript reports.
- The Tier E effect (control versus nA, up to 0.005% of Q) is an order of magnitude smaller than the engine effect (nA versus wA), which puts the size of the Tier E correction in context: the orientation defect it removed had no effect on Q, and even a different method of drawing the channels barely has one.
- The hourly routed hydrograph differs more (single steps up to 5 to 8 percent of the peak, r 0.998 to 0.999): one outlet instead of two on CA, more channel storage (end-of-run 442 to 501 m3 against 254 to 287 m3 on CA, 234 against 148 m3 on AR), and Manning n 0.10 instead of 0.045 slow and merge the routed response. The manuscript's Q is the lateral inflow, which this does not touch; the routed hydrograph is now usable in both engines (R12 holds to floating-point precision on the single WW-DHSVM outlet).
- Initial Storage differs by 0.2 mm on both basins; the channel state and the channel cells' initialization are the only inputs that changed (grid states sha256-identical).

## What this settles

- The Toolkit's Tier E network is not an artefact of the GRASS route: an independent engine draws the same trunk from the same DEM, and the two agree wherever the accumulation exceeds a few hundred cells. The headwater disagreement is the MFD versus D8 difference in where a given support area is reached, which is a modelling choice, not an error in either engine.
- For the manuscript, the network question is closed in both forms: neither the Tier E correction nor a change of engine reaches the reported numbers beyond the third decimal. The open input-set decision (`tier_e_network_orientation_2026_09_22.md`, Consequences) is unchanged by this record.
- For the joint tool: the WW-DHSVM fork branch carries the drop analysis, the outlet invariants (`checkOutlets`, merged into `checkTopology`) and `demFromRaster`, with tests that pin the Camp Branch network (266 cells, 30 segments, outlet (68, 68)) and the corrected drop objective (120 cells). This is the material for the pull request to UConn-EFC/WW-DHSVM.

## Follow-ups

- The CA double outlet remains a Toolkit property (MFD accumulation is not monotone along the stream). A single-outlet option was recorded and not adopted in the Tier E record; the cross-engine result gives it a concrete target, the WW-DHSVM outlet at (68, 68).
- A_c in cells versus m2: both engines take the support area as an area; the Toolkit rounds it to a cell count for GRASS (`hydrology.py`, `round(area / cell_area)`), WW-DHSVM compares the upstream area in m2 directly, and its sweep compares the cell count. A shared convention (an integer cell count, or an area with the rounding stated) belongs in the joint tool.
- pyflwdir issue for the Strahler cache: https://github.com/Deltares/pyflwdir/issues/123.

## Files

- Script: `scripts/diagnostics/cross_engine/compare_engines.py` (needs the WW-DHSVM fork branch and its import-time dependencies).
- Evidence: `docs/audit/cross_engine_2026_09_23/{CA,AR}/` (comparison `.md` and `.json`, the WW-DHSVM drop sweep, the WW-DHSVM stream files and `Channel.State`, the Toolkit diagnostic on the WW-DHSVM network) and `dhsvm_eval_{CA,AR}.csv` (the metrics tables above). The WW-DHSVM `stream.map.dat` header carries a generation timestamp, so its bytes differ between runs while its body does not.
- DHSVM reruns: `scripts/diagnostics/tierE_dhsvm/` (the `ww` kind; the `wA_` outputs live in `TestCase/CA/output/` and `TestCase/AR/output/` on the Mac, with the `tierE_compare_*_vs_wA_*.csv` column tables).
- Working records: the project documents STEP2_CROSS_ENGINE_PLAN and DROP_ANALYSIS_CACHE_DEFECT of 2026-09-23.
