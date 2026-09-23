# Tier E: stream network orientation audit (2026-09-22)

**Date**: 2026-09-22
**Scope**: the stream network stage of `standalone_CA/` (and, by descent, the same stage of `qgis_CA/`, which is not changed)
**Branch**: `feat-tierE-raster-network`
**Commits**:
- `0fc5fda` Export the r.stream.extract direction raster
- `71b8a4c` Build the stream network from the raster
- this commit: the audit record, the validation log entry, the README updates, the DHSVM rerun scripts, and the retirement of `iso_check_network.py`

## Context

While comparing this toolkit with WW-DHSVM (Zhi Li, UConn) in September 2026, a synthetic test of `stream_network.py` (`test_orientation_stream_network`) showed that the network it writes points every segment the wrong way. A read-only diagnostic (`standalone_CA/diagnostics/check_network_orientation.py`) then confirmed the defect on the real CA files: on the QGIS reference tree (`qgis_CA_ref`, the June 8 outputs) and on a fresh standalone baseline regenerated from the group-space inputs at commit `1aaee3d`, both of which are the same 41-segment network.

This is the same character as Tier A, Tier B and the slope conditioning fix: silent, no runtime error, a structurally wrong model input. It sat behind the earlier validations because every check compared the two pipelines with each other (both carry the defect), and because the manuscript's simulated streamflow is the channel lateral inflow, which does not depend on routing (below).

## Issue

On a basin-clipped DEM `r.watershed` (without `-a`) writes negative accumulation on every cell: the sign marks "possibly runoff from outside region", the magnitude is the accumulation. `vector_attrs.py` took the endpoint of each `r.to.vect` line with the larger raw value as downstream, that is, the smaller magnitude, so every segment was reversed. `stream_network.py` then re-linked the reversed lines by geometry and merged the resulting sinks, which made a headwater the written outlet.

What the CA 28 m network looked like (diagnostic, D1 to D7):

- 41 segments, 19 of them one cell long; the written outlet (segment 41, `down 0`) is a headwater at row 19 whose "downstream" end carries |acc| 60.5 cells, the initiation threshold. The true mouth cells (rows 65 to 67, columns 67 to 68, |acc| up to 3387.5) sat inside segments pointing away from the mouth.
- in-degree histogram {0: 18, 1: 20, 2: 2, 16: 1}: sixteen merged sinks feeding one node.
- 275 map records over 235 distinct cells: 29 cells under two segments, 4 records off the stream raster.
- one channel class (13) on all 41 segments: the `meanmsq` sampler kept only positive accumulation, found none, and fell back to the cell area.
- a second finding closed the first repair plan: MFD accumulation is not monotone along the extracted stream (8 of 41 segments; the main stem reads 3387.5 at (63, 67), 1953 at (65, 68) and 1559 at (67, 68)), so orienting the `r.to.vect` lines by |acc| would still misplace segments. Six segments were also ambiguous under the `r.watershed` drainage raster.

The earlier "Tier A" note that the merged-sink segment was the true basin outlet is withdrawn: it is a headwater.

## Root cause

Two decisions in the old stage: reading the sign of `r.watershed` accumulation as a magnitude, and rebuilding the topology from line geometry instead of following the flow direction that `r.stream.extract` had already computed for every stream cell.

## Fix

The network is now built from the rasters, the method WW-DHSVM uses.

- `run_hydrology_grass.sh` exports the `r.stream.extract` `direction` raster as `stream_dir.tif` (D8 codes 1 to 8 counter-clockwise from north-east; negative where the stream leaves the raster). `hydrology.py` stamps its CRS like the other exports; `paths.py` gains `STREAM_DIR`, `SEGMENTS_CSV`, `STREAM_CELLS`.
- `segments_from_raster.py` (new stage 6): each stream cell's successor comes from `stream_dir.tif`; a segment starts at every head (no inflow) and below every confluence (two or more inflows); a segment ends where its successor leaves the raster, is NoData, or is not a stream cell (outlet); a Kahn topological order gives the propagated rank, dense from 1. Attributes keep the earlier conventions: slope is the mean tan(slope_filled) over the segment's cells (0.01 where undefined), `meanmsq` is the largest |acc| in the segment times the cell area, per-cell length is the D8 step (the outlet cell carries one cell size), aspect is the D8 azimuth. Outputs: `streamfile_attr.shp`, `segments.csv`, `stream_cells.csv`. Cycles and unassigned cells raise.
- `channelclass_standalone.py` is unchanged and now sees real drainage areas (5.2e4 to 2.69e6 m2 on CA), so the PNNL bins engage.
- `network_files.py` (new): checks the invariants before writing (every `down` exists, `down 0` agrees with `is_outlet`, positive slope and length, class ids exist, ranks dense and increasing downstream, no duplicate cells, the lowest stream cell lies in an outlet segment, map records equal stream cells) and writes `stream.network.dat` with `SAVE "outlet_k"` on every outlet row and `stream.map.dat` in the fixed-width format used before.
- Outlets are written as the raster has them. CA 28 m has two mouth cells (segment 22, the main stem, tail (65, 68) at 838.0 m, rank 8; segment 12, a first-order tributary along the south edge, tail (67, 68) at 829.6 m, rank 1). DHSVM takes both (verified below). Forcing a single outlet stays an option, not a need.
- `vector_attrs.py`, `stream_network.py` and `iso_check_network.py` are retired to `scripts/legacy/standalone_CA_pre_tierE/`. `iso_check_network.py` matched the standalone network to the reference by (length, slope, class) signature; with the reference itself reversed and the corrected network a different segmentation, the check has no meaning.
- `standalone_CA/tests/test_segments_from_raster.py` covers the stage without GRASS (six synthetic tests, including the file format end to end).
- `standalone_CA/diagnostics/check_network_orientation.py` is the read-only diagnostic (invariants I1 to I9; I2, I5, I8 and I9 are hard, the rest informational).

## Validation, pipeline side (DCC, CA 28 m, `/work/ys451/dhsvm_ca/tierE/`)

Byte gate against the baseline: dem, mask, soil, veg, soildepth, the five uniform soil depths and the three grid state files are sha256-identical. Changed, as designed: `stream.class.dat` (26983e53..), `stream.map.dat` (a4f0f02e..), `stream.network.dat` (6a01bd2f..), `Channel.State` (819fa1f7..). The four pre-existing rasters are unchanged; `stream_dir.tif` is new (7e1c7acb..).

| | before (r.to.vect + geometric re-linking) | after (raster-native) |
|---|---|---|
| segments | 41 (19 one cell long) | 22 |
| map records / distinct cells / cells under two segments / cells off the raster | 275 / 235 / 29 / 4 | 231 / 231 / 0 / 0 |
| total channel length | 7526.0 m | 7600.8 m (exact D8 steps; each outlet cell carries a full cell) |
| outlets | 1, a headwater at row 19 | 2 mouth cells (segments 22 and 12) |
| in-degree histogram | {0: 18, 1: 20, 2: 2, 16: 1} | {0: 12, 2: 10} |
| rank histogram | {1: 18, 2: 12, 3: 4, 4: 2, 5: 2, 6: 2, 7: 1} | {1: 12, 2: 2, 3: 2, 4: 2, 5: 1, 6: 1, 7: 1, 8: 1} |
| channel classes | 13 on all 41 | 13 on 17, 14 on 5 |
| SAVE rows | none | outlet_1 (segment 12), outlet_2 (segment 22) |

Diagnostic v3 on the fixed outputs: 0 hard invariants failed; all 22 segments descend in elevation along `stream_dir.tif` (20 of 22 also in |acc|, the two exceptions being the MFD non-monotonicity above); the two outlet rows are the two outlet cells; I7 (single outlet) is informational and fails by design. The old outlet segment is now segment 5, a rank-1 head draining to segment 13.

## Validation, DHSVM side (Mac, DHSVM3.2 built 2025-01-13, `TestCase/CA`)

Scripts: `scripts/diagnostics/tierE_dhsvm/` (`tierE_rerun_CA.py`, `tierE_compare_CA.py`, `tierE_eval_CA.py`; the working copies live in `DHSVM-PNNL-2025/TestCase/CA/`).

**What the manuscript's Q is.** `DHSVM_CA_Stitch_Optimal_LAI_V2.py` reads `parts[2]` of the Stream.Flow "Totals" rows, which `channel.c` writes as `date 0 total_lateral_inflow total_outflow total_storage total_storage_change total_error`. The manuscript's daily Q is therefore the water entering the channel network, summed to daily mm; it equals the water-balance term ChannelInt (ratio 1.00009 on all three runs). Routing only redistributes water within the day: on a May run whose network routed (`512_LAI20`), daily routed outflow against daily lateral inflow over 2017-02 to 2018-12 gives NSE 1.0000, r 0.99999, PBIAS -0.001%, max 0.09 mm/d. In the April manuscript runs the routed outflow was zero at every step and the channel error absorbed the inflow (the pre-Tier-A network never routed a segment); nothing in the manuscript used it.

**Input provenance.** The April manuscript inputs survive as `DEM_CA_apr` (a copy of "DEM_CA_0406 copy", files dated 2026-04-06 16:17). They differ from today's `DEM_CA_0406` (the June 8 outputs, equal to `qgis_CA_ref` and to the standalone baseline) in two files only: `soildepth.bin` (7443b833.. versus ff76ec11..; the slope unit fix of Tier B and the slope conditioning of June 8) and `stream.network.dat` (d0973965.. versus e6ab190f..; Tier A and B). Initial Storage in Mass.Final.Balance dates each set: 606.150 mm (April), 637.306 (May 12, Tier B), 589.065 (June 8).

**Reproduction.** The April base control runs (`oA_S4h`, `oA_LAI70`, `oA_LAI20`: the templates with all nine input files and the state directory pointed at `DEM_CA_apr`) reproduce the April outputs byte for byte: all fifteen output files carry the April sha256.

**Network effect, the manuscript inputs with only the network swapped** (`nA_*` = `DEM_CA_apr` binaries, Tier E stream files and Channel.State; the grid states are the same bytes). Mass.Final.Balance, mm, three years:

| run | ET | ChannelInt | Initial Storage | Final Storage | Mass Error |
|---|---|---|---|---|---|
| CA_S4h = oA_S4h | 2808.941 | 2210.152 | 606.150 | 1021.259 | -0.071 |
| nA_S4h | 2808.850 | 2210.246 | 606.140 | 1021.235 | -0.081 |
| CA_LAI70 = oA_LAI70 | 2325.027 | 2686.421 | 606.150 | 1029.104 | -0.052 |
| nA_LAI70 | 2324.950 | 2686.508 | 606.140 | 1029.081 | -0.056 |
| CA_LAI20 = oA_LAI20 | 1119.144 | 3878.913 | 606.150 | 1042.750 | -0.032 |
| nA_LAI20 | 1119.112 | 3878.942 | 606.140 | 1042.733 | -0.041 |

Daily lateral inflow (the manuscript's Q): S4h +0.100 mm over three years (+0.0045%), daily RMSE 0.0054 mm/d, max daily difference 0.057 mm/d; LAI70 +0.082 mm, RMSE 0.0052, max 0.060; LAI20 +0.032 mm, RMSE 0.0044, max 0.050; r 0.999998 or better. The same pair on today's inputs (`old_*` versus `new_*`) gives +0.095, +0.075 and +0.033 mm with the same daily statistics. The hourly routed hydrograph does change (max single-step difference about 600 m3 per step, 8% of the peak; RMSE 18 m3 per step; r 0.9991) because the direction is now right and there are two mouth outlets; the routed outflow is usable for the first time (ratio to ChannelInt 1.0000).

R12, both SAVE outlets against the routed total: |outlet_1 + outlet_2 - total_outflow| is at most 0.1 m3 per step over 26280 steps, 1.2e-5 to 1.4e-5 of the peak, the rounding of the five significant digits DHSVM prints. Outlet shares: main stem 84.4%, south-edge tributary 15.6%.

**Manuscript metrics**, computed with the manuscript's own code (`tierE_eval_CA.py` imports `DHSVM_CA_Stitch_Optimal_LAI_V2.py` and uses its parsers and `calc_metrics`; the manuscript base reproduces `CA_Stitched_Summary.csv` and the Evaluation_Summary tables exactly). Stitched series, 2017 from the LAI20 run and 2018 from the LAI70 run, window 2017-02-01 to 2018-12-31:

| inputs, network | NSE | r | RMSE (mm/d) | PBIAS (%) |
|---|---|---|---|---|
| April, April (manuscript = oA) | 0.681499 | 0.841120 | 2.31034 | -0.322 |
| April, Tier E (nA) | 0.681898 | 0.841491 | 2.30889 | -0.323 |
| June, post-Tier-A+B (old) | 0.682523 | 0.841351 | 2.30662 | +0.374 |
| June, Tier E (new) | 0.682915 | 0.841716 | 2.30520 | +0.373 |

At the manuscript's precision: 0.681 / 0.841 / 2.310 / -0.3 becomes 0.682 / 0.841 / 2.309 / -0.3 with the network alone, and 0.683 / 0.842 / 2.305 / +0.4 with the June inputs and the network. The per-year optimum cannot move: the candidate neighbours sit at PBIAS -9.1 (LAI30) and +9.2 (LAI10) in 2017 and -4.3 (LAI80) and +4.2 (LAI60) in 2018, against a network effect of 0.001 and a soil-depth effect of 0.7.

**A DHSVM limit found on the way.** This DHSVM3.2 build crashes (SIGTRAP from the fortified `sprintf`) when the Output Directory value exceeds 78 characters: `InitDump.c` writes it into `char sumoutfile[100]` with `failure_summary.txt` appended, and `RouteSubSurface.c` into `char satoutfile[100]` with `saturation_extent.txt` (21 characters) appended. The rerun script enforces the limit.

## Scope

- Changed: `stream.class.dat`, `stream.network.dat`, `stream.map.dat`, `Channel.State` (segment count and widths), and the new `stream_dir.tif`, `segments.csv`, `stream_cells.csv`.
- Unchanged: every grid binary, every grid state file, `flow_acc.tif`, `flow_dir.tif`, `stream_raster.tif`, `streamfile.shp`, the slope and soil-depth stages, `channelclass_standalone.py`, `states.py` (it reads the new files and writes a 22-segment Channel.State).
- The manuscript's streamflow, ET partition and calibrated leaf-area multipliers are unaffected beyond the third decimal; the quantities above are the record.
- `qgis_CA/` still carries the defect in its step 7 and is not changed; it is the reference implementation, kept for the record.

## Consequences for the manuscript (decision pending)

Two ways to carry this into the manuscript are on the table; Song and Brad decide.

1. Keep the April runs as published and cite this audit: the network could not have changed any reported number beyond the third decimal, and the simulated Q never depended on routing.
2. Rerun the CA trio on the current pipeline outputs (June soil depth, Tier E network) and refresh Fig 2 (0.683 / 0.842 / 2.305 / +0.4). Recommended: the published numbers then come from inputs the released toolkit reproduces.

## Follow-ups

- AR: regenerate the AR network with the fixed pipeline on the grid the manuscript AR run used (the AR reference network was built by the same stage) and rerun AR S4h. The AR source DEM and polygon have to be restaged on DCC first.
- The joint-tool plan: port the drop analysis to WW-DHSVM and cross-compare both engines on CA and AR.
- Recorded alternatives, not adopted: segment slope as drop over length instead of the mean tan(slope raster); a forced single outlet.

## Files

- Fix: `standalone_CA/pipeline/run_hydrology_grass.sh`, `hydrology.py`, `paths.py`, `segments_from_raster.py`, `network_files.py`, `run_pipeline.sh`.
- Tests and diagnostics: `standalone_CA/tests/test_segments_from_raster.py`, `standalone_CA/diagnostics/check_network_orientation.py`.
- DHSVM reruns: `scripts/diagnostics/tierE_dhsvm/`.
- Retired: `scripts/legacy/standalone_CA_pre_tierE/` (`vector_attrs.py`, `stream_network.py`, `iso_check_network.py`).
- Working records: the project documents TIER_E_NETWORK_ORIENTATION_PLAN, TIER_E_DIAGNOSTIC_RESULTS, TIER_E_FIXED_RUN and TIER_E_MAC_RERUN of 2026-09-22.
