# QGIS Scripts

QGIS Python implementation of the DHSVM preprocessing pipeline.
Run from the QGIS Python Console; entry point is `prep_dhsvm_inputs.py`.

## Quick start

```python
exec(open("prep_dhsvm_inputs.py").read())
```

Inputs are auto-discovered relative to the parent directory of the scripts folder:

- `Reprojected_DEM/` — projected DEM (e.g., UTM)
- `Reprojected_Watersheds/` — watershed mask polygon (optional)

Outputs land in:

- `DHSVM_input_binaries/` — flat binary grids (DEM, Mask, Soil, Veg, soildepth)
- `DHSVM_input_streams/` — `stream.network.dat`, `stream.map.dat`, `stream.class.dat`
- `modelstate/` — Interception / Snow / Soil / Channel initial states
- `Intermediate_GIS/` — diagnostic rasters and shapefiles

## Pipeline

`prep_dhsvm_inputs.py` orchestrates these submodules:

| Submodule | Purpose | Output |
|---|---|---|
| `channelclass.py` | slope × area channel classification | `stream.class.dat` |
| `soildepthscript.py` | weighted slope + elevation soil depth | `soildepth.bin` |
| `dem_to_dhsvm_bins.py` | DEM / Mask / Soil / Veg base maps | `*.bin` |
| `generate_dhsvm_states.py` | initial conditions | `*.State.*` |

Standalone utilities (not part of the pipeline):

- `dhsvm_area_config_from_dem.py` — prints the `[AREA]` block for the DHSVM INPUT file
- `plot_mainstem.py` — longitudinal profile diagnostic
- `make_2class_veg_map.py`, `make_4class_veg_map.py`, `rdnbr_to_veg_bs.py` — project-specific burn-severity vegetation maps

## Row/Col origin convention

DHSVM uses a **top-left** origin for grid indexing: row 0 is the northern edge, row index increases downward. Two places in `prep_dhsvm_inputs.py` compute Row/Col:

1. `_ensure_fields_rowcol_len()` — Row/Col on the stream shapefile's attribute table (advisory only).
2. `_write_stream_map()` — Row/Col in `stream.map.dat` (read by DHSVM).

Both use the top-left origin. If you ever modify one, modify the other. See the docstring on `_ensure_fields_rowcol_len` for historical context.

## `_archive/`

`_archive/rowcolmap.py` and `_archive/roadaspect.py` are not imported by the pipeline. They are retained as historical references and as a documented cautionary example (`rowcolmap.py` has a bottom-left origin bug that would silently flip stream rows if reused). See `_archive/README.md`. Do not import from `_archive/`.

## On previously generated outputs

`stream.map.dat` produced by earlier versions of this pipeline is correct: its Row/Col values come from `_write_stream_map()`, which has always used the top-left convention. The Row/Col **attribute fields** on the intermediate shapefile were bottom-left until Tier D2 (May 2026), but those attribute fields are not consumed by DHSVM. No re-run is required on account of this fix alone.
