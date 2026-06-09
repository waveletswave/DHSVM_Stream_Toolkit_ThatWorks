#!/usr/bin/env bash
# Stage 2 (GRASS part): slope from clipped DEM via r.slope.aspect, degrees.
# Mirrors qgis_CA prep_dhsvm_inputs.py r.slope.aspect call (format=0 = degrees).
# Region is locked to the clipped DEM so the output lands on the same 74x82
# grid as the qgis_CA reference. Run through the py3 shim. CRS stamping and
# slope_fill conditioning happen afterward in slope.py.
set -euo pipefail

GRASS="python /hpc/group/abmurraylab/ys451/bin/grass76_py3.py"
CLIPPED=/work/ys451/dhsvm_ca/standalone_dev/outputs/elev_clipped.tif
SLOPE_RAW=/work/ys451/dhsvm_ca/standalone_dev/outputs/slope_raw.tif
LOC=/work/ys451/dhsvm_ca/grass/gloc_slope

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
