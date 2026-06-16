#!/usr/bin/env bash
# =====================================================================
# run_pipeline.sh  -  end-to-end standalone DHSVM preprocessing (CA, DCC)
#
# One command: source DEM + watershed  ->  full DHSVM input set.
# Runs every stage in dependency order, fail-fast (set -e), checking each
# stage's key output exists before continuing.
#
# Order (decision B puts states AFTER vector, since channel state reads the
# standalone stream files):
#   1. DEM entry       elev_clipped.tif
#        default            clip.py            (CA 28 m byte-match reproducer)
#        DHSVM_DEM_SOURCE   prep_dem.py        (general; optional fetch_dem first)
#   2. slope  (GRASS)  slope_raw.tif       (run_slope_grass.sh)
#      slope  (py)     slope_filled.tif    (stamp CRS + slope_fill + validate)
#   3. bins            DHSVM_input_binaries/{dem,mask,soil,veg,soildepth_uniform_*}.bin
#   4. hydrology       flow_acc/flow_dir/stream_raster.tif + streamfile.shp
#                      (hydrology.py runs run_hydrology_grass.sh internally + CRS stamp)
#   5. soildepth       DHSVM_input_binaries/soildepth.bin
#   6. vector_attrs    streamfile_attr.shp
#      channelclass    DHSVM_input_streams/stream.class.dat (+chanclass write-back)
#      stream_network  DHSVM_input_streams/stream.{network,map}.dat
#   7. states          modelstate/{Interception,Snow,Soil}.State... + Channel.State...
#
# Paths come from paths.py (env-overridable: DHSVM_INPUTS/REF/OUT/...). This
# script does NOT redefine paths; it reads them from paths.py so the checks and
# the stages agree, and passes the slope stage's DEM/out/shim to the GRASS .sh.
# Both GRASS .sh now take their paths and shim as args, so the DHSVM_* env vars
# retarget the whole pipeline (Python and GRASS) end to end.
#
# DEM entry (stage 1) has two modes, switched by DHSVM_DEM_SOURCE:
#   unset          -> clip.py, the CA 28 m byte-match reproducer (default)
#   set to a path  -> general entry: prep_dem.py reprojects that source DEM to
#                     paths.EPSG at DHSVM_DEM_RES (default 10 m) and masks to the
#                     watershed; set DHSVM_FETCH=1 to fetch the source from 3DEP
#                     first (needs a networked node). prep_dem and clip.py both
#                     write paths.ELEV_CLIPPED, so stages 2..7 are unchanged.
# For a data-driven stream threshold, set A_c once before running:
#   export DHSVM_STREAM_SOURCE_AREA_M2=$(python3 ../diagnostics/drop_analysis.py \
#          "$OUT/elev_clipped.tif" --emit-area --tmin 50 --tmax 800)
#
# Usage:
#   bash run_pipeline.sh
#   DHSVM_OUT=/path/to/out bash run_pipeline.sh                       # retarget run
#   DHSVM_DEM_SOURCE=/path/src.tif DHSVM_DEM_RES=10 bash run_pipeline.sh   # general
# =====================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# Resolve OUT (and a couple of key paths) from paths.py so checks match stages.
OUT="$(python3 -c 'import paths; print(paths.OUT)')"
BIN_DIR="$(python3 -c 'import paths; print(paths.BIN_DIR)')"
STREAMS_DIR="$(python3 -c 'import paths; print(paths.STREAMS_DIR)')"
STATE_DIR="$(python3 -c 'import paths; print(paths.STATE_DIR)')"
SLOPE_SH="$HERE/run_slope_grass.sh"
# Args for the GRASS slope stage, resolved from paths.py (the same single source).
SLOPE_DEM="$(python3 -c 'import paths; print(paths.ELEV_CLIPPED)')"
SLOPE_RAW_OUT="$(python3 -c 'import paths; print(paths.SLOPE_RAW)')"
GRASS_SHIM="$(python3 -c 'import paths; print(paths.SHIM)')"
PREP_DIR="$HERE/../prep"

echo "=========================================================="
echo "  STANDALONE DHSVM PREPROCESSING PIPELINE"
echo "  OUT = $OUT"
echo "=========================================================="

need() {  # need <file/dir> <stage-name-that-should-have-made-it>
  if [ ! -e "$1" ]; then
    echo "[FAIL] expected output missing after $2: $1" >&2
    exit 1
  fi
}

step() { echo ""; echo ">>> [$1] $2"; }

# ---------------------------------------------------------------- 1. DEM entry
# Two modes, switched by DHSVM_DEM_SOURCE (see header). Both write
# paths.ELEV_CLIPPED, so everything downstream is identical.
if [ -n "${DHSVM_DEM_SOURCE:-}" ]; then
  step 1 "DEM entry (general: prep_dem)"
  SRC="$DHSVM_DEM_SOURCE"
  if [ "${DHSVM_FETCH:-0}" = "1" ]; then
    echo "    fetch DEM from 3DEP -> $SRC  (needs a networked node)"
    python3 "$PREP_DIR/fetch_dem.py" "$SRC"
    need "$SRC" "fetch_dem"
  fi
  echo "    prep_dem $SRC -> elev_clipped.tif  (res ${DHSVM_DEM_RES:-10} m)"
  python3 "$PREP_DIR/prep_dem.py" "$SRC"
else
  step 1 "DEM entry (CA 28 m reproducer: clip.py)"
  python3 clip.py
fi
need "$OUT/elev_clipped.tif" "DEM entry"

# ---------------------------------------------------------------- 2. slope
step 2 "slope (GRASS r.slope.aspect)"
bash "$SLOPE_SH" "$SLOPE_DEM" "$SLOPE_RAW_OUT" "$GRASS_SHIM"
need "$OUT/slope_raw.tif" "slope GRASS"
step 2 "slope (CRS stamp + fill + validate)"
python3 slope.py
need "$OUT/slope_filled.tif" "slope py"

# ---------------------------------------------------------------- 3. bins
step 3 "base-map binaries"
python3 bins.py
need "$BIN_DIR/dem.bin" "bins"
need "$BIN_DIR/mask.bin" "bins"

# ---------------------------------------------------------------- 4. hydrology
step 4 "hydrology (GRASS chain + CRS stamp)"
python3 hydrology.py
need "$OUT/flow_acc.tif" "hydrology"
need "$OUT/streamfile.shp" "hydrology"

# ---------------------------------------------------------------- 5. soildepth
step 5 "soil depth"
python3 soildepth.py
need "$BIN_DIR/soildepth.bin" "soildepth"

# ---------------------------------------------------------------- 6. vector chain
step 6 "vector attributes"
python3 vector_attrs.py
need "$OUT/streamfile_attr.shp" "vector_attrs"
step 6 "channel class + chanclass write-back"
python3 channelclass_standalone.py
need "$STREAMS_DIR/stream.class.dat" "channelclass"
step 6 "stream network + map"
python3 stream_network.py
need "$STREAMS_DIR/stream.network.dat" "stream_network"
need "$STREAMS_DIR/stream.map.dat" "stream_network"

# ---------------------------------------------------------------- 7. states
step 7 "initial states (grid + channel)"
python3 states.py
need "$STATE_DIR/Channel.State.01.01.2016.00.00.00" "states"

# ---------------------------------------------------------------- summary
echo ""
echo "=========================================================="
echo "  PIPELINE COMPLETE -- DHSVM input set:"
echo "=========================================================="
echo "[grid binaries]  $BIN_DIR"
ls -1 "$BIN_DIR" | sed 's/^/    /'
echo "[stream files]   $STREAMS_DIR"
ls -1 "$STREAMS_DIR" | sed 's/^/    /'
echo "[model state]    $STATE_DIR"
ls -1 "$STATE_DIR" | sed 's/^/    /'
echo "=========================================================="
