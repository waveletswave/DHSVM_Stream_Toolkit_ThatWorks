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
#   1. clip            elev_clipped.tif
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
# Usage:
#   bash run_pipeline.sh
#   DHSVM_OUT=/path/to/out bash run_pipeline.sh    # retargets the whole run
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

# ---------------------------------------------------------------- 1. clip
step 1 "clip DEM to watershed"
python3 clip.py
need "$OUT/elev_clipped.tif" "clip"

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
