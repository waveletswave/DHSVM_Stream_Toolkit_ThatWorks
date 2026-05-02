# DHSVM Spatial & Stream Network Toolkit (That Actually Works)

A QGIS + Python toolkit for generating all spatial inputs, stream networks, and initial model states required by the Distributed Hydrology Soil Vegetation Model (DHSVM).

This repository now contains **two parallel implementations** of the same preprocessing pipeline:

| Implementation                | Location                      | Use Case                                                  |
| ----------------------------- | ----------------------------- | --------------------------------------------------------- |
| **QGIS Python Console** | [`qgis/`](./qgis/)             | Interactive desktop work, visual debugging                |
| **Standalone CLI**      | [`standalone/`](./standalone/) | HPC / Duke DCC, batch processing, reproducible Slurm jobs |

The standalone version uses `rasterio` + `geopandas` + `shapely` instead of the QGIS Python API, so it can be run from any terminal — including HPC compute nodes that have no GUI / X server.

## Standalone migration progress

| Module                      | QGIS file                                                    | Standalone file                    | Status                                             |
| --------------------------- | ------------------------------------------------------------ | ---------------------------------- | -------------------------------------------------- |
| Channel classification      | `qgis/channelclass.py`                                     | `standalone/channelclass.py`     | ✅ Verified byte-identical (DEM_CA_0406)           |
| ROW/COL & SegID assignment  | `qgis/rowcolmap_robust.py`                                 | `standalone/rowcolmap.py`        | ✅ Verified byte-identical (DEM_CA_0406)           |
| Main stem profile           | `qgis/plot_mainstem.py`                                    | `standalone/plot_mainstem.py`    | 🟡 Draft — pending verification                   |
| Direct binary serialization | `qgis/dem_to_dhsvm_bins.py`                                | (already standalone, no QGIS deps) | ✅ Already QGIS-free                               |
| State initialization        | `qgis/generate_dhsvm_states.py`                            | (already standalone, no QGIS deps) | ✅ Already QGIS-free                               |
| Soil depth                  | `qgis/soildepthscript.py`                                  | (already standalone, no QGIS deps) | ✅ Already QGIS-free                               |
| Pipeline orchestration      | `qgis/prep_dhsvm_inputs.py`                                | —                                 | 📋 Not yet ported                                  |
| Vegetation maps             | `qgis/rdnbr_to_veg_bs.py`, `qgis/make_2class_veg_map.py` | —                                 | 📋 Not yet ported                                  |
| Road/stream aspect          | (none)                                                       | (none)                             | ❌ Not needed for this toolkit's target watersheds |

## Key Improvements and Bug Fixes (both implementations)

* **Direct Binary Serialization (No More `myconvert`):** Replaces the cumbersome ESRI ASCII export workflow. The toolkit natively serializes Float32 and Int8 flat binaries directly from Python, completely bypassing the legacy `myconvert` C-utility.
* **Physically Consistent Initial States:** Replaces `MakeModelStateBin` and legacy bash scripts. Safely handles Endianness and parses dynamically computed channel widths to produce physically accurate initial channel water volumes, preventing severe numerical shocks during model spin-up.
* **Uniform Soil Depth Baselines:** Generates a suite of pixel-aligned uniform soil depth grids for subsurface stormflow sensitivity testing.
* **Top-left Origin Convention (DHSVM-correct):** The robust rowcolmap uses top-left origin as required by `stream.map.dat`. An earlier bottom-left version is preserved (`qgis/rowcolmap_legacy_bottom_left.py`) but **deprecated** — see [`qgis/README.md`](./qgis/README.md) for details.
* **Floating-point Parity Across Implementations:** The standalone `rowcolmap` uses an explicit `(y_top - y) / dy` division (matching the QGIS robust version) rather than `rasterio.transform.rowcol()`, ensuring byte-identical Row/Col output even at exact cell-boundary centroids.

## Output Directory Structure

* `/DHSVM_input_binaries/`: Flat binary grids (`dem.bin`, `mask.bin`, `soildepth.bin`, etc.) for the `[TERRAIN]` and `[SOILS]` blocks.
* `/DHSVM_input_streams/`: Topological routing files (`stream.map.dat`, `stream.network.dat`, `stream.class.dat`).
* `/modelstate/`: Initial condition files for the `[STATE]` block:
  * `Interception.State.[Date].bin`, `Snow.State.[Date].bin`, `Soil.State.[Date].bin` (Float32 binary)
  * `Channel.State.[Date]` (ASCII, no extension)
* `/Intermediate_GIS/`: Diagnostic `.tif` and `.shp` files for QGIS visualization.

## Usage

### QGIS workflow (interactive)

Run the scripts in `qgis/` from the QGIS Python Console. Ensure the Processing plugin and GRASS provider are enabled. Auto-discovery resolves paths from `Reprojected_DEM/elev_clipped.tif` upward.

### Standalone workflow (CLI / HPC)

```bash
# Set up the conda environment once
conda env create -f standalone/environment.yml
conda activate dhsvm-prep

# Run individual modules
python standalone/channelclass.py <streamfile.shp> <output_dir>
python standalone/rowcolmap.py annotate <streamfile.shp> <dem.tif>
python standalone/plot_mainstem.py <streamfile.shp> <dem.tif> --output-dir <dir>
```

---

**Author:** Y. Song
**Organization:** Duke University
**Reference physics:** Wigmosta et al. (1994), *Water Resour. Res.* 30(6), 1665–1679
