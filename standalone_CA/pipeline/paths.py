"""Path constants for the standalone DHSVM preprocessing pipeline.

Single source of truth for every stage's paths. The three roots (INPUTS, REF,
OUT) plus the GRASS shim and EPSG are env-overridable, so the same code runs on
a different case or machine by setting environment variables -- no code edits:

    export DHSVM_INPUTS=/path/to/inputs
    export DHSVM_REF=/path/to/qgis_reference     # optional, validation only
    export DHSVM_OUT=/path/to/outputs
    export DHSVM_GRASS_SHIM=/path/to/grass76_py3.py
    export DHSVM_EPSG=32617
    export DHSVM_SRC_DEM=USGS_..._UTM17.tif       # input DEM filename (in INPUTS)
    export DHSVM_WATERSHED=cabr_watershed_UTM17.shp
    export DHSVM_CASE=CA                          # watershed tag; drives RUN_OUT
    export DHSVM_RUN_OUT=/work/ys451/dhsvm/CA     # DHSVM run dir (config + output)
    export DHSVM_MET_FILE=/path/to/met.txt        # met forcing, given in full
    export DHSVM_BIN=/work/ys451/dhsvm/bin/DHSVM  # compiled DHSVM binary

Defaults reproduce the CA case on DCC.
"""
import os
from pathlib import Path

# ----------------------------- roots (env-overridable) ----------------------
INPUTS = Path(os.environ.get("DHSVM_INPUTS", "/hpc/group/abmurraylab/ys451/dhsvm_ca/inputs"))
REF    = Path(os.environ.get("DHSVM_REF",    "/hpc/group/abmurraylab/ys451/dhsvm_ca/qgis_CA_ref"))
OUT    = Path(os.environ.get("DHSVM_OUT",    "/work/ys451/dhsvm_ca/standalone_dev/outputs"))
SHIM   = Path(os.environ.get("DHSVM_GRASS_SHIM", "/hpc/group/abmurraylab/ys451/bin/grass76_py3.py"))
EPSG   = int(os.environ.get("DHSVM_EPSG", "32617"))

OUT.mkdir(parents=True, exist_ok=True)

# ----------------------------- inputs ---------------------------------------
SRC_DEM   = INPUTS / os.environ.get("DHSVM_SRC_DEM", "USGS_1_n36w084_20220725_UTM17.tif")
WATERSHED = INPUTS / os.environ.get("DHSVM_WATERSHED", "cabr_watershed_UTM17.shp")

# ----------------------------- output directories ---------------------------
BIN_DIR     = OUT / "DHSVM_input_binaries"   # all DHSVM grid binaries: dem/mask/soil/veg + soildepth (decision A: unified)
STATE_DIR   = OUT / "modelstate"             # initial states
STREAMS_DIR = OUT / "DHSVM_input_streams"    # stream.class/network/map.dat

# ----------------------------- derived raster/vector outputs ----------------
# clip stage
ELEV_CLIPPED  = OUT / "elev_clipped.tif"
# slope stage
SLOPE_RAW     = OUT / "slope_raw.tif"        # r.slope.aspect output before fill
SLOPE_FILLED  = OUT / "slope_filled.tif"     # after slope_fill conditioning
# hydrology stage
FLOW_ACC      = OUT / "flow_acc.tif"
FLOW_DIR      = OUT / "flow_dir.tif"
STREAM_RASTER = OUT / "stream_raster.tif"
STREAMFILE    = OUT / "streamfile.shp"       # stream vector (geometry only)
# vector_attrs stage
STREAMFILE_ATTR = OUT / "streamfile_attr.shp"
# soildepth stage
SOILDEPTH_TIF = OUT / "soildepth.tif"        # diagnostic GeoTIFF (bin goes in BIN_DIR)

# ----------------------------- reference (validation only) ------------------
REF_ELEV_CLIPPED = REF / "elev_clipped.tif"
REF_SLOPE        = REF / "stream_slope_filled.tif"   # conditioned slope reference
REF_INTERMED     = REF / "Intermediate_GIS"          # flow_acc/flow_dir/stream_raster/streamfile + soildepth.tif
REF_BINARIES     = REF / "DHSVM_input_binaries"      # soildepth.bin etc.
REF_STREAMS      = REF / "DHSVM_input_streams"        # stream.class/network/map.dat + modelstate ref lives under REF/modelstate
REF_MODELSTATE   = REF / "modelstate"

# ----------------------------- DHSVM run paths (model run / config gen) ------
# Used by standalone_CA/run/make_dhsvm_config.py to render the .dhs. CASE drives
# RUN_OUT, so switching watershed is one env var. MET_FILE is given in full, not
# derived: some adjacent watersheds share one met file, so the name stays
# explicit rather than following a {CASE} convention.
CASE      = os.environ.get("DHSVM_CASE", "CA")               # CA / AR / AR_UP / CA_TO ...
RUN_OUT   = Path(os.environ.get("DHSVM_RUN_OUT", f"/work/ys451/dhsvm/{CASE}"))
MET_FILE  = Path(os.environ.get("DHSVM_MET_FILE",
    "/hpc/group/abmurraylab/ys451/dhsvm/met/DHSVM_met_input_CS01_hourly_0411_mean_wind_interpolation.txt"))
DHSVM_BIN = Path(os.environ.get("DHSVM_BIN", "/work/ys451/dhsvm/bin/DHSVM"))

# Stream initiation: critical support area at the channel head, as a physical
# area so the cell threshold scales with DEM resolution. Default is the CA
# 60-cell equivalent (0.04757 km2), which keeps the byte-identical validation.
# For any other watershed, set from drop analysis or override via the env var.
STREAM_SOURCE_AREA_M2 = float(os.environ.get("DHSVM_STREAM_SOURCE_AREA_M2", 47571.5))
