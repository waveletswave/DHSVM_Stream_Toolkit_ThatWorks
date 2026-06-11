#!/usr/bin/env bash
# =====================================================================
# run_hydrology_grass.sh  -  standalone hydrology core (DCC, path a)
#
# Single GRASS location, three modules chained in one mapset (option i):
#   r.watershed  ->  r.stream.extract  ->  r.to.vect
# Mirrors the three grass7: processing.run calls in prep_dhsvm_inputs.py.
# Intermediate rasters stay inside GRASS; only the four reference outputs
# are exported (flow_acc, flow_dir, stream_raster, streamfile vector).
#
# Run through the Py3 shim, region locked to the clipped DEM, same pattern
# as the slope stage. CRS is stamped back onto the GeoTIFF exports by the
# Python wrapper (hydrology.py), because r.out.gdal cannot write the CRS tag
# on this install (gcs.csv missing).
#
# Usage:  bash run_hydrology_grass.sh <elev_clipped.tif> <out_dir> <shim.py>
# =====================================================================
set -euo pipefail

ELEV="$1"          # clipped DEM (EPSG:32617), drives region + grid
OUTDIR="$2"        # where GRASS exports land (GeoTIFF + shp)
SHIM="$3"          # grass76_py3.py compatibility shim (was hardcoded)
LOC=/tmp/ghydro_$$         # throwaway location, built from the CRS-bearing DEM
THRESH=60                  # MIN_SRC_CELLS in prep
CONV=5                     # r.watershed convergence in prep
MEM=300                    # memory in prep

mkdir -p "$OUTDIR"
rm -rf "$LOC"               # immunize against a leftover location of the same name

# Build the location from the DEM (option b: bypass the EPSG gcs.csv lookup),
# then run the whole chain inside one --exec batch so the intermediate rasters
# (accumulation, drainage) persist across modules in the same mapset.
python "$SHIM" -c "$ELEV" "$LOC" --exec bash -s <<EOF
set -euo pipefail

# Import the DEM under a known name, then lock region + resolution to it.
# Everything lands on the 74x82 grid.
r.in.gdal input="$ELEV" output=elev_in --overwrite --quiet
g.region raster=elev_in

# --- r.watershed: MFD default (no -s / -a / -m), convergence=5, memory=300 ---
# prep: {'elevation','accumulation','drainage','convergence':5,'memory':300}
# No -a: keep the negative flow-accumulation values for edge/underestimate cells.
r.watershed elevation=elev_in accumulation=flow_acc drainage=flow_dir \
    convergence=$CONV memory=$MEM --overwrite --quiet

# --- r.stream.extract: elevation + accumulation in; stream_raster + vector out ---
# prep: {'elevation','accumulation','threshold':60,'stream_raster','stream_vector','-m':True}
# NOTE 1: prep also passes 'direction': flow_dir, but in r.stream.extract 'direction'
# is an OUTPUT, not an input. The QGIS grass7: wrapper tolerated that key; the
# native module does not consume an external direction. So we do NOT feed flow_dir
# here -- elevation + accumulation is the real input set, matching what QGIS ran.
# NOTE 2: prep's '-m':True is NOT a valid r.stream.extract flag (it is an
# r.watershed flag). The QGIS wrapper silently dropped it, so r.stream.extract
# actually ran without -m; we drop it too. memory= is the only memory control.
r.stream.extract elevation=elev_in accumulation=flow_acc threshold=$THRESH \
    stream_raster=stream_rast stream_vector=stream_vec_native \
    memory=$MEM --overwrite --quiet

# --- r.to.vect fallback: stream raster -> lines ---
# In prep, r.stream.extract's stream_vector came out empty on CA (0 features),
# so _force_lines_from_raster ran r.to.vect on the stream raster. trials order is
# [point, line, line, line+smooth]; first that yields lines wins -> type=line, no -s.
# Reproduce that here. Validation against streamfile.shp (41 lines) is the check.
r.to.vect input=stream_rast output=streamfile_vec type=line --overwrite --quiet

# --- exports: only the four reference artifacts ---
# flow_acc is DCELL (negative on CA), exported as Float64 by default.
# flow_dir and stream_rast are CELL (int).
r.out.gdal input=flow_acc    output="$OUTDIR/flow_acc.tif"      format=GTiff --overwrite --quiet
r.out.gdal input=flow_dir    output="$OUTDIR/flow_dir.tif"      format=GTiff --overwrite --quiet
r.out.gdal input=stream_rast output="$OUTDIR/stream_raster.tif" format=GTiff --overwrite --quiet

# Vector export to shapefile (the fallback stream lines).
v.out.ogr input=streamfile_vec output="$OUTDIR/streamfile.shp" format=ESRI_Shapefile --overwrite --quiet
EOF

echo "[ok] GRASS hydrology chain done. Exports in $OUTDIR"
echo "     (CRS stamping of the GeoTIFFs is handled by hydrology.py)"
