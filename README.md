# DHSVM Stream Network Toolkit (That Actually Works)

This repository provides a QGIS and GRASS GIS workflow for generating DHSVM stream network input files from a DEM and watershed polygon.

This toolkit builds upon standard DHSVM preprocessing methods. It includes corrections for spatial sampling errors and an automated tool for extracting, analyzing, and visualizing main stem longitudinal profiles.

## Key Improvements and Bug Fixes

PyQGIS raster sampling functions (`provider.sample`) can return a boolean success flag alongside the data value. If parsed incorrectly, this assigns identical flat slopes to downstream routing segments. This pipeline corrects the boolean parsing order to ensure the routing algorithm accurately computes the true physical gradient ($\tan(\theta)$) for every segment in `stream.network.dat`.

## Included Scripts

### 1. Core DHSVM Pipeline

* `createstreamnetwork.py`: The main driver script. Orchestrates DEM clipping, GRASS topological processing, and safe flow accumulation routing. It writes all DHSVM stream network dat files.
* `channelclass.py`: Defines stream classes by contributing area and slope. Writes `stream.class.dat`.
* `roadaspect.py`: Computes road and stream network aspect (includes a stable fallback algorithm).
* `rowcolmap.py`: Utilities for precise mapping between DEM grid coordinates and DHSVM stream map indices.
* `wshdslope.py`: Computes accurate watershed and stream slope rasters.
* `soildepthscript.py`: Generates a soil depth raster consistent with the DEM grid using slope and source area weighting.

### 2. Hydrologic Analysis Tools

* `extract_mainstem_profile.py`: Automatically identifies the watershed outlet, traces the main stem upstream using maximum cumulative drainage area (`meanmsq`), and extracts the longitudinal profile.
  * Visualization: Plots a high quality scatter plot of the profile, applying a color gradient to represent local slope.
  * Knickpoint Detection: Automatically flags statistically significant abrupt changes in localized slope along the main stem.
  * Data Export: Exports a detailed `.csv` containing distance, elevation, slope, and knickpoint boolean flags for further statistical analysis.

## Output Artifacts

When successfully executed in the QGIS Python Console, the pipeline produces the following outputs in your workspace:

* `stream.map.dat`: Snake style, per cell, 0 to 360 degree aspect, north bearing, clockwise.
* `stream.network.dat`: Flow accumulation directed network, single stable outlet with DownSegID equal to negative one.
* `stream.class.dat`: Channel class structure.
* `stream_profile.png`: Visual validation of the main stem geomorphology.
* `stream_profile_data.csv`: Node by node spatial data for the main stem.

## Usage

Run the scripts sequentially inside the QGIS Python Console. Ensure the Processing plugin and GRASS provider are enabled.

1. Execute `createstreamnetwork.py` first to generate the topological shapefiles and `.dat` files.
2. Execute `extract_mainstem_profile.py` to validate the slope calculations and visualize your main stem. The script utilizes an auto discovery mechanism so no manual path hardcoding is required as long as it sits within the same project directory tree.

---

Author: Y. Song
Organization: Duke University
