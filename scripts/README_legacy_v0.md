# DHSVM Spatial & Stream Network Toolkit (That Actually Works)

This repository provides a comprehensive QGIS, GRASS GIS, and Python workflow for generating all spatial inputs, stream networks, and initial model states required by the Distributed Hydrology Soil Vegetation Model (DHSVM). 

This toolkit builds upon standard DHSVM preprocessing methods but introduces critical bug fixes, removes dependencies on legacy C-utilities, and enforces strict workspace organization to ensure physical accuracy and methodological reproducibility. 

## Key Improvements and Bug Fixes

* **Direct Binary Serialization (No More `myconvert`):** Replaces the cumbersome ESRI ASCII export workflow. The toolkit natively serializes Float32 and Int8 flat binaries directly from Python, completely bypassing the legacy `myconvert` C-utility.
* **Physically Consistent Initial States:** Replaces `MakeModelStateBin` and legacy bash scripts. It safely handles Endianness (byte-order) for binary grid states and parses dynamically computed channel widths to generate physically accurate initial channel water volumes, preventing severe numerical shocks during the model spin-up period.
* **Uniform Soil Depth Baselines:** Automatically generates a suite of uniform soil depth grids (e.g., 2.0m, 3.0m) that are perfectly pixel-aligned with the master mask. These serve as reliable control groups for subsurface stormflow parameter sensitivity testing.
* **Topology Sampling Correction:** Fixes a common PyQGIS raster sampling bug where boolean success flags are misparsed. This ensures the routing algorithm accurately computes the true physical gradient ($\tan(\theta)$) for every segment in `stream.network.dat`.

## Included Scripts

### 1. Main Preparation Pipeline
* `prep_dhsvm_inputs.py`: The main execution script. Coordinates DEM clipping, GRASS topological processing, and triggers all sub-modules to automatically generate ready-to-use binaries and stream files.
* `dem_to_dhsvm_bins.py`: Handles direct `.bin` generation for DEM, mask, soil type, vegetation, and uniform depth baselines.
* `generate_dhsvm_states.py`: Generates the initial condition grid binaries and the text-based channel state file.
* `soildepthscript.py`: Generates a dynamic topographic soil depth raster using slope and source area weighting.

### 2. Stream Topology & Classification
* `channelclass.py`: Defines stream classes by contributing area and slope. Writes `stream.class.dat`.
* `rowcolmap_robust.py`: Robust per-cell ROW/COL & SegID assignment using a top-left origin convention (matches DHSVM stream.map.dat format). Replaces the legacy `rowcolmap.py` and the misnamed `roadaspect.py` from earlier iterations.
* `rowcolmap_legacy_bottom_left.py`: ⚠️ Deprecated. Earlier version that used a bottom-left origin and produced vertically-flipped Row indices. Kept for historical reference only — **do not use** for new DHSVM runs.

*(Note: a true road/stream aspect calculator was never required by the watersheds this toolkit was developed for, since DHSVM road routing is not used in those configurations. The `roadaspect.py` filename in earlier iterations was a documentation/naming error — its actual content was a robust rowcolmap, now properly renamed.)*

### 3. Hydrologic Analysis Tools
* `plot_mainstem.py`: Automatically identifies the watershed outlet and traces the true physical main stem upstream using DEM elevation routing. Identifies knickpoints, plots longitudinal profiles, and exports detailed `.csv` spatial data.

*(Note: Previous iterations of the pipeline are archived in the `scripts/legacy/` directory for historical reference and methodological transparency).*

## Output Directory Structure

To keep the workspace clean, the pipeline automatically routes outputs into four distinct directories:

* `/DHSVM_input_binaries/`: Flat binary grids (`dem.bin`, `mask.bin`, `soildepth.bin`, etc.) ready for the `[TERRAIN]` and `[SOILS]` blocks in the `.dhs` file.
* `/DHSVM_input_streams/`: Topological routing files (`stream.map.dat`, `stream.network.dat`, `stream.class.dat`).
* `/modelstate/`: Initial condition files for the `[STATE]` block. The pipeline explicitly generates the following four files (appending `.bin` to grid states to satisfy DHSVM's internal I/O hardcoding):
  * `Interception.State.[Date].bin` (Float32 binary)
  * `Snow.State.[Date].bin` (Float32 binary)
  * `Soil.State.[Date].bin` (Float32 binary)
  * `Channel.State.[Date]` (ASCII text format, no extension)
* `/Intermediate_GIS/`: Diagnostic `.tif` and `.shp` files for visualization and debugging in QGIS.

## Usage

Run the scripts inside the QGIS Python Console. Ensure the Processing plugin and GRASS provider are enabled.

1. Execute `prep_dhsvm_inputs.py` to run the complete spatial and state generation pipeline.
2. Execute `plot_mainstem.py` to validate your slope calculations and visualize the main stem geomorphology. 

The scripts utilize an auto-discovery mechanism, so no manual path hardcoding is required as long as your foundational DEM sits within the project directory tree (e.g., `Reprojected_DEM/elev_clipped.tif`).

---
**Author:** Y. Song <br>
**Organization:** Duke University