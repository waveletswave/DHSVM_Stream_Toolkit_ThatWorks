# QGIS preprocessing pipeline (`qgis_CA/`)

The original implementation of the DHSVM preprocessing pipeline, run from the QGIS Python Console. It clips the DEM, builds the hydrologically routed stream network with GRASS through the QGIS Processing framework, and coordinates a set of Python sub-modules to write the gridded binaries, the topological stream files, and the initial model states.

This is the implementation the headless pipeline was built from and validated against. For batch or HPC use, prefer the no-QGIS port in [`../standalone_CA/`](../standalone_CA/), which reproduces these outputs (byte-identical binaries, data-equivalent rasters) without a desktop GUI. Use `qgis_CA/` when you want the interactive, visual-debugging path or a reference for what the standalone reproduces. The narrative writeup of the workflow is in the repository root (`DHSVM_Workflow_merged.tex` and the accompanying PDF); note that the writeup predates some of the audited fixes listed below, so this README is the current description.

## Quick start

Run the master script inside the QGIS Python Console, with the Processing plugin and the GRASS provider enabled:

```python
exec(open("prep_dhsvm_inputs.py").read())
```

## Inputs (auto-discovered)

Inputs are discovered relative to the parent directory of the scripts folder.

- `Reprojected_DEM/`: the projected DEM (a UTM grid). The script tries `Cropped_DEM.tif`, `dem.tif`, `elev.tif`, and `elev_clipped.tif`, then any other `.tif` in the folder.
- `Reprojected_Watersheds/`: the watershed mask polygon (optional). Shapefiles whose name contains "watershed" are preferred. If no polygon is found, the full DEM extent is used.

## Outputs

- `DHSVM_input_binaries/`: flat binary grids (`dem.bin`, `mask.bin`, `soil.bin`, `veg.bin`, `soildepth.bin`, plus uniform soil-depth baselines).
- `DHSVM_input_streams/`: topological routing files (`stream.network.dat`, `stream.map.dat`, `stream.class.dat`).
- `modelstate/`: initial Interception, Snow, Soil, and Channel states.
- `Intermediate_GIS/`: diagnostic rasters and shapefiles (`elev_clipped.tif`, `flow_acc.tif`, `flow_dir.tif`, `stream_raster.tif`, `stream_slope.tif`, `streamfile.shp`, `soildepth.tif`).

## Pipeline stages

`prep_dhsvm_inputs.py` runs these stages in order, reporting progress to the console:

1. Clip the DEM to the watershed polygon (`gdal:cliprasterbymasklayer`).
2. GRASS `r.watershed` for flow direction and accumulation.
3. GRASS `r.stream.extract` for the stream raster, at a support threshold of `MIN_SRC_CELLS` cells (default 60).
4. GRASS `r.to.vect` to vectorize the streams.
5. GRASS `r.slope.aspect` for slope, in degrees.
6. Condition the slope raster (`slope_fill.py`): fill the NoData ring `r.slope.aspect` leaves on the DEM edge and around former interior voids, by the mean of valid neighbours, before any consumer reads it.
7. Build the directed network and assign Row/Col, inline in the master script (`_ensure_fields_rowcol_len`, `_write_stream_map`); write `stream.network.dat` and `stream.map.dat`.
8. Classify channels (`channelclass.py`); write `stream.class.dat`.
9. Compute the topographic soil depth (`soildepthscript.py`); write `soildepth.bin` and `soildepth.tif`.
10. Serialize the base maps and uniform soil-depth baselines (`dem_to_dhsvm_bins.py`).
11. Generate the initial states (`generate_dhsvm_states.py`).

A few stage parameters are set at the top of the master script: `MIN_SRC_CELLS` (stream support threshold, default 60), `NIN_MODE` (stream-order mode, default `propagated` for the dense order DHSVM requires, see Tier A below), and `BASE_SLOPE_SAMPLES` (samples per segment when computing slope, default 12).

## Sub-modules

`prep_dhsvm_inputs.py` loads these through `importlib`, so each can be updated independently.

| Sub-module | Purpose | Output |
|---|---|---|
| `slope_fill.py` | fill missing-slope cells by neighbour mean (slope conditioning) | conditioned `stream_slope.tif` |
| `channelclass.py` | slope and area channel classification | `stream.class.dat` |
| `soildepthscript.py` | weighted slope, area, and elevation soil depth | `soildepth.bin` |
| `dem_to_dhsvm_bins.py` | DEM, mask, soil, veg base maps and uniform depths | `*.bin` |
| `generate_dhsvm_states.py` | physically consistent initial conditions | `*.State.*` |

`slope_fill.py` is pure GDAL and NumPy (no QGIS), so the same module is called by the standalone pipeline.

## Standalone utilities (not part of the pipeline)

- `dhsvm_area_config_from_dem.py`: prints the `[AREA]` block for the DHSVM configuration file (a second copy lives under `DHSVM_Input_Workflow/`).
- `plot_mainstem.py`: traces the physical main stem from the outlet and exports a longitudinal profile for geomorphic QA.
- `make_2class_veg_map.py`, `make_4class_veg_map.py`, `rdnbr_to_veg_bs.py`: project-specific burn-severity vegetation maps.

## Audited corrections

The pipeline has been through several correctness audits, recorded in [`../docs/audit/`](../docs/audit/). The fixes are silent-failure cases: no runtime error, but a distorted model input.

- **Tier A** (stream topology): `r.stream.extract` outlet sentinel and stream-order semantics in `stream.network.dat`. `NIN_MODE` is `propagated` so the stream order is dense; DHSVM's `channel_route_network` breaks its outer loop on the first empty order, which left the outlet unrouted (basin outflow stuck at zero) under indegree or Shreve ordering.
- **Tier B** (slope units): `r.slope.aspect` output switched from percent (`format=1`) to degrees (`format=0`), with the consumer's raster `sample()` unpacking fixed to match. See `tier_b_unit_fix_2026_05_12.md`.
- **Slope conditioning** (2026-06-08): the slope-completeness fix in step 6 above. `soildepthscript` previously filled the missing-slope ring with slope 0, which receives the maximum slope bonus in the depth formula and produced a spurious deep ring at the basin edge. `slope_fill.py` fills those cells with the neighbour mean instead. See `slope_conditioning_2026_06_08.md`.
- **Tier D** (Row/Col origin): Row/Col uses a top-left origin throughout, and the unused `rowcolmap.py` and `roadaspect.py` were archived (see `_archive/` below).

## Row/Col origin convention

DHSVM uses a top-left origin for grid indexing: row 0 is the northern edge, and the row index increases downward. Two places compute Row/Col, both top-left:

1. `_ensure_fields_rowcol_len`: Row/Col on the stream shapefile's attribute table (advisory only).
2. `_write_stream_map`: Row/Col in `stream.map.dat`, which DHSVM reads.

If you modify one, modify the other.

## `_archive/`

`_archive/rowcolmap.py` and `_archive/roadaspect.py` are not imported by the pipeline; the master script implements Row/Col inline instead. They are kept as design history and as a cautionary example: `rowcolmap.py` computes the row index from the bottom-left origin, which would silently flip stream rows and break `stream.map.dat` if reused. See [`_archive/README.md`](./_archive/README.md). Do not import from `_archive/`.

## On previously generated outputs

`stream.map.dat` produced by earlier versions of this pipeline is correct: its Row/Col values come from `_write_stream_map`, which has always used the top-left convention. The Row/Col attribute fields on the intermediate shapefile were bottom-left until Tier D2 (May 2026), but those fields are not consumed by DHSVM. No re-run is required on account of that fix alone.
