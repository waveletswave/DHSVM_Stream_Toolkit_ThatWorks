#!/usr/bin/env bash
# Stage 2 (GRASS part): slope from clipped DEM via r.slope.aspect, degrees.
# Mirrors qgis_CA prep_dhsvm_inputs.py r.slope.aspect call (format=0 = degrees).
# Region is locked to the clipped DEM so the output lands on the same 74x82
# grid as the qgis_CA reference. Run through the py3 shim. CRS stamping and
# slope_fill conditioning happen afterward in slope.py.
# Usage:  bash run_slope_grass.sh <clipped_dem.tif> <slope_raw_out.tif> <shim.py>
set -euo pipefail

CLIPPED="$1"        # clipped DEM (EPSG:32617), drives region + grid
SLOPE_RAW="$2"      # r.slope.aspect output (Float32 degrees), export target
SHIM="$3"          # grass76_py3.py compatibility shim (was hardcoded)
GRASS="python $SHIM"
LOC="/tmp/gslope_$$"   # throwaway location, built from the CRS-bearing DEM

# Fresh location built from the clipped DEM (carries CRS as WKT; no gcs.csv needed).
rm -rf "$LOC"
$GRASS -c "$CLIPPED" "$LOC" --exec g.version

# Import the clipped DEM and LOCK THE REGION to it (critical for grid alignment).
$GRASS "$LOC/PERMANENT" --exec r.in.gdal input="$CLIPPED" output=elev_clipped --overwrite
$GRASS "$LOC/PERMANENT" --exec g.region raster=elev_clipped -p

# Slope in DEGREES (format=0), matching qgis_CA. Defaults otherwise (Horn-like,
# zscale=1) — same as the QGIS r.slope.aspect call, which set only format.
$GRASS "$LOC/PERMANENT" --exec r.slope.aspect elevation=elev_clipped slope=slope_raw format=degrees --overwrite

# Export to GeoTIFF (CRS tag will be missing due to gcs.csv; stamped in slope.py).
$GRASS "$LOC/PERMANENT" --exec r.out.gdal input=slope_raw output="$SLOPE_RAW" format=GTiff type=Float32 --overwrite

echo "GRASS slope stage done -> $SLOPE_RAW"
