# DHSVM Spatial & Stream Network Toolkit (That Actually Works)

A toolkit for generating the spatial inputs, stream networks, and initial model states that the Distributed Hydrology Soil Vegetation Model (DHSVM) requires, starting from a watershed polygon and a DEM. It produces the gridded terrain and soil binaries, the topological stream-routing files, and the initial state files, plus the diagnostic GeoTIFFs for inspection.

The repository holds two implementations of the same preprocessing logic:

| Implementation | Location | Use |
| --- | --- | --- |
| **Standalone CLI** | [`standalone_CA/`](./standalone_CA/) | Current and validated. Runs from any terminal, including HPC compute nodes with no GUI. Recommended for new work. |
| **QGIS Python Console** | [`qgis_CA/`](./qgis_CA/) | The original implementation. Interactive desktop work and visual debugging. The reference the standalone was validated against. |

The standalone version uses `rasterio` + `geopandas` + `shapely` + `pyflwdir` and calls GRASS GIS for the slope and hydrology stages, instead of the QGIS Python API, so it runs headless and is scriptable into batch jobs.

## Repository layout

- `standalone_CA/`: the current pipeline (CLI / HPC). See below.
- `qgis_CA/`: the original QGIS Python Console scripts and its own `README.md`.
- `scripts/`: archived material. `legacy/` holds earlier QGIS and standalone versions of the scripts, and `diagnostics/` holds one-off analysis utilities. `README_legacy_v0.md` documents the older layout.
- `docs/audit/`: the audit trail. The standalone rebuild inventory and the per-stage audits (slope conditioning, soil-depth sensitivity, the Tier B unit fix).
- `example_outputs/`: a sample stream-profile figure and its data.
- `DHSVM_Workflow_merged.tex` and `Workflow_for_Preparing_DHSVM_Input_Files_...pdf`: the written narrative of the full workflow.

## The standalone pipeline (`standalone_CA/`)

One command turns a source DEM and a watershed polygon into a complete DHSVM input set. The stages run in dependency order, each checked before the next, driven by environment variables read in `pipeline/paths.py` so a new basin or machine needs no code edits.

```
standalone_CA/
  prep/         fetch_dem.py     fetch a DEM over the basin from USGS 3DEP (the one networked step)
                prep_dem.py      reproject to the pipeline CRS at the target resolution, clip to the polygon
  pipeline/     clip.py          CA 28 m byte-match reproducer for the regression test
                slope.py         GRASS r.slope.aspect, then CRS stamp and fill conditioning (slope_fill.py)
                bins.py          dem / mask / soil / veg binaries plus uniform soil-depth baselines
                hydrology.py     GRASS r.watershed -> r.stream.extract -> r.to.vect (flow, stream raster, vector)
                soildepth.py     dynamic soil depth from the PNNL weighting
                vector_attrs.py  stream attributes; channelclass_standalone.py -> stream.class.dat;
                                 stream_network.py -> stream.network.dat + stream.map.dat
                states.py        initial Interception / Snow / Soil / Channel states
                paths.py         single source of truth for paths and config (env-overridable)
                run_pipeline.sh  the one-command orchestrator (fail-fast); the GRASS stages are run_*_grass.sh
                compare.py, compare_bin.py   raster and binary comparison helpers for validation
  diagnostics/  quicklook.py     six-panel diagnostic figure (DEM + network, soil depth, slope, flow, location)
                drop_analysis.py constant stream drop analysis to pick the stream support area A_c
                plot_drop.py     plot the drop sweep; iso_check_*.py and check_slope_units.py are isolation checks
  run/          make_dhsvm_config.py  render the DHSVM .dhs config (DEM-derived [AREA], declared meridian);
                                      area.py and CA.dhs.template support it
  docs/         validation_log.md, NEW_WATERSHED_GUIDE.md, DCC_SETUP.md, BUILD_DCC.md, stream_threshold.md
```

### Running it

For a new watershed, follow `standalone_CA/docs/NEW_WATERSHED_GUIDE.md`. The short version: fetch the DEM, set the stream support area from the drop analysis, run the pipeline, and eyeball the quick-look.

```bash
# set the basin (see paths.py for the full list of variables)
export DHSVM_EPSG=32617
export DHSVM_WATERSHED=/path/your_watershed.shp
export DHSVM_OUT=/path/outputs_yourbasin_10m
export DHSVM_DEM_SOURCE=/path/src_3dep_10m.tif
export DHSVM_DEM_RES=10

python3 standalone_CA/prep/fetch_dem.py "$DHSVM_DEM_SOURCE" --res 10 --watershed "$DHSVM_WATERSHED"
python3 standalone_CA/prep/prep_dem.py "$DHSVM_DEM_SOURCE"
export DHSVM_STREAM_SOURCE_AREA_M2=$(python3 standalone_CA/diagnostics/drop_analysis.py "$DHSVM_OUT/elev_clipped.tif" --emit-area)
bash standalone_CA/pipeline/run_pipeline.sh
python3 standalone_CA/diagnostics/quicklook.py --run-dir "$DHSVM_OUT" --boundary "$DHSVM_WATERSHED" --out "$DHSVM_OUT/quicklook.png" --title "yourbasin"
```

To reproduce the development basin's calibrated 28 m grid for the regression test, use `clip.py` as the entry instead of `prep_dem.py`, then run the pipeline.

## Validation

The standalone pipeline was validated against the audited QGIS reference on the CA (Camp Branch) test case, stage by stage and end to end. Grid binaries are byte-identical to the reference; rasters are checked for data equivalence, since GRASS and rasterio write different TIFF containers; the stream network is verified against DHSVM's topological requirements rather than the QGIS tie-break. End to end, hourly streamflow over 2016 to 2018 is numerically equivalent to the QGIS pipeline (NSE 0.99999510, PBIAS -0.000059 percent). The two input sets differ only in two documented numerical respects, a 1-ULP float difference in soil depth and a tie-break in stream extraction, which leave every aggregate metric unchanged. Portability was confirmed by building a full input set on a second basin, AR (Arrowwood), at 10 m and 30 m from its polygon alone. The full record is in `standalone_CA/docs/validation_log.md`.

## Output directory structure

- `DHSVM_input_binaries/`: flat binary grids (`dem.bin`, `mask.bin`, `soil.bin`, `veg.bin`, `soildepth.bin`, plus uniform soil-depth variants) for the `[TERRAIN]` and `[SOILS]` blocks.
- `DHSVM_input_streams/`: topological routing files (`stream.class.dat`, `stream.network.dat`, `stream.map.dat`).
- `modelstate/`: initial-condition files for the `[STATE]` block. `Interception.State.[Date].bin`, `Snow.State.[Date].bin`, `Soil.State.[Date].bin` (Float32), and `Channel.State.[Date]` (ASCII).
- intermediate GeoTIFFs (`elev_clipped.tif`, `slope_filled.tif`, `flow_acc.tif`, and so on) for inspection and for the quick-look.

## Key features and fixes

- **Direct binary serialization.** Float32 and Int8 flat binaries are written natively from Python, with no ESRI ASCII export and no legacy `myconvert` C utility.
- **Physically consistent initial states.** Endianness is handled and channel widths are parsed to produce physically accurate initial channel water volumes, avoiding numerical shock at spin-up.
- **Physical-area stream threshold.** The channel-initiation threshold is set as a support area A_c rather than a fixed cell count, so drainage density is held fixed as resolution changes. `drop_analysis` picks A_c as the start of the first sustained passing band of the constant stream drop law.
- **Tight clip by default.** `prep_dem` clips the grid to the basin with a zero-cell margin, matching the calibrated CA grid; `--pad-cells` adds a buffer if wanted.
- **Single source of truth for paths.** Every stage reads `paths.py`, with env-overridable roots, so a case or machine change needs only the `DHSVM_*` variables. Both the Python and the GRASS shell layers are parameterized.
- **One-command orchestration.** `run_pipeline.sh` runs the full pipeline fail-fast and was verified to match a step-by-step run byte for byte.
- **Top-left origin convention.** Row/Col assignment uses the top-left origin `stream.map.dat` requires. The earlier bottom-left version is preserved under `scripts/legacy/` and is deprecated.

## The QGIS implementation (`qgis_CA/`)

The original scripts, run from the QGIS Python Console with the Processing plugin and the GRASS provider enabled. This is the implementation the standalone was validated against; see [`qgis_CA/README.md`](./qgis_CA/README.md) for its own notes. The QGIS workflow is interactive and useful for visual debugging, but it needs a desktop GUI and cannot run on a headless compute node.

---

**Author:** Y. Song
**Organization:** Duke University
**Reference physics:** Wigmosta et al. (1994), *Water Resour. Res.* 30(6), 1665-1679
