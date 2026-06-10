# -*- coding: utf-8 -*-
# =====================================================================
# states.py  -  standalone (no-QGIS, DCC) port of generate_dhsvm_states.py
#
# Stage: initial states (grid + channel). Pure NumPy + GDAL.
# Ports MakeModelStateBin.c / MakeChannelState.sh, same as the qgis_CA version.
#
# This is a byte-parity port for the GRID states. Only changed from the
# original: path resolution (now via paths.py) and run-on-import -> run().
# The grid-state write logic and the channel-state parsing/format are copied
# verbatim.
#
# Grid-state parity: every state array is spatially uniform, so the output
# bytes are fixed by (array count, write order, dtype, constant values) alone.
# All four are copied verbatim, so the grid binaries are byte-identical to the
# qgis_CA reference given the same grid dims.
#
# Channel state now reads the STANDALONE stream files (STREAMS_DIR), produced
# by the hydrology + vector stages (decision B), not the qgis_CA reference.
# Because the standalone stream network differs from the reference at two d=0
# tie links (see validation_log sub-step C), Channel.State is self-consistent
# with the standalone network and is NOT expected to be byte-identical to the
# qgis_CA Channel.State -- that is correct: the channel state must match the
# stream network actually used. Run order: hydrology -> vector -> states.
# =====================================================================

import numpy as np
from osgeo import gdal

from paths import ELEV_CLIPPED, STREAMS_DIR, STATE_DIR

# stream files now come from the standalone vector stage (decision B)
DEM_TIF      = ELEV_CLIPPED
CLASS_FILE   = STREAMS_DIR / "stream.class.dat"
NETWORK_FILE = STREAMS_DIR / "stream.network.dat"

# ==============================================================
# State config  (copied verbatim from generate_dhsvm_states.py)
# ==============================================================
MODEL_START   = "01/01/2016-00"     # MM/DD/YYYY-HH
INITIAL_DEPTH = 0.05                 # 5 cm, to prevent initial shocks

N_VEG_LAYERS          = 2
RAIN_INTERCEPTION     = [0.0, 0.0]
SNOW_INTERCEPTION_TOP = 0.0
TEMP_INT_STORAGE      = 0.0

SNOW_MASK             = 0.0
DAYS_SINCE_LAST_SNOW  = 0.0
SNOW_WATER_EQUIV      = 0.0
LW_BOTTOM             = 0.0
T_BOTTOM              = 0.0
LW_TOP                = 0.0
T_TOP                 = 0.0
COLD_CONTENT          = 0.0

N_SOIL_LAYERS         = 3
SOIL_MOISTURE         = [0.25, 0.25, 0.25, 0.25]   # NSoil + 1 (saturated zone)
SOIL_TEMP_SURF        = 5.0
SOIL_TEMP_LAYERS      = [5.0, 5.0, 5.0]
GROUND_HEAT           = 0.0
INITIAL_RUNOFF        = 0.0

# ==============================================================
# Helpers
# ==============================================================
def get_grid_dimensions():
    if not DEM_TIF.exists():
        raise FileNotFoundError(f"[ERROR] DEM not found: {DEM_TIF}")
    ds = gdal.Open(str(DEM_TIF))
    return ds.RasterYSize, ds.RasterXSize  # NY, NX


def write_2d_float32(file_handle, value, ny, nx):
    """NY*NX uniform float32 array, appended to the open binary file."""
    arr = np.full((ny, nx), value, dtype=np.float32)
    arr.tofile(file_handle)


# ==============================================================
# 1. Grid states (replaces MakeModelStateBin)
# ==============================================================
def generate_grid_states():
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    ny, nx = get_grid_dimensions()

    date_str = MODEL_START.replace('/', '.').replace('-', '.') + ".00.00"
    print(f"\n[step] grid states for {date_str} (NY={ny}, NX={nx})")

    # --- Interception state (.bin required for DHSVM ReadModelState) ---
    int_file = STATE_DIR / f"Interception.State.{date_str}.bin"
    with open(int_file, "wb") as f:
        for val in RAIN_INTERCEPTION: write_2d_float32(f, val, ny, nx)
        for _ in range(N_VEG_LAYERS): write_2d_float32(f, SNOW_INTERCEPTION_TOP, ny, nx)
        write_2d_float32(f, TEMP_INT_STORAGE, ny, nx)
    print(f"  -> {int_file.name}")

    # --- Snow state (.bin) ---
    snow_file = STATE_DIR / f"Snow.State.{date_str}.bin"
    with open(snow_file, "wb") as f:
        write_2d_float32(f, SNOW_MASK, ny, nx)
        write_2d_float32(f, DAYS_SINCE_LAST_SNOW, ny, nx)
        write_2d_float32(f, SNOW_WATER_EQUIV, ny, nx)
        write_2d_float32(f, LW_BOTTOM, ny, nx)
        write_2d_float32(f, T_BOTTOM, ny, nx)
        write_2d_float32(f, LW_TOP, ny, nx)
        write_2d_float32(f, T_TOP, ny, nx)
        write_2d_float32(f, COLD_CONTENT, ny, nx)
    print(f"  -> {snow_file.name}")

    # --- Soil state (.bin) ---
    soil_file = STATE_DIR / f"Soil.State.{date_str}.bin"
    with open(soil_file, "wb") as f:
        for val in SOIL_MOISTURE: write_2d_float32(f, val, ny, nx)
        write_2d_float32(f, SOIL_TEMP_SURF, ny, nx)
        for val in SOIL_TEMP_LAYERS: write_2d_float32(f, val, ny, nx)
        write_2d_float32(f, GROUND_HEAT, ny, nx)
        write_2d_float32(f, INITIAL_RUNOFF, ny, nx)
    print(f"  -> {soil_file.name}")


# ==============================================================
# 2. Channel states (replaces MakeChannelState.sh)
# ==============================================================
def generate_channel_state():
    date_str = MODEL_START.replace('/', '.').replace('-', '.') + ".00.00"

    # Channel state is a TEXT file and does NOT get .bin appended.
    out_file = STATE_DIR / f"Channel.State.{date_str}"
    print(f"[step] channel state -> {out_file.name}")

    class_widths = {}
    with open(CLASS_FILE, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 2:
                class_widths[int(parts[0])] = float(parts[1])

    with open(NETWORK_FILE, 'r') as f_in, open(out_file, 'w') as f_out:
        for line in f_in:
            if not line.strip() or line.startswith('#'): continue
            parts = line.split()
            seg_id   = int(parts[0])
            length_m = float(parts[3])
            class_id = int(parts[4])

            width_m = class_widths.get(class_id, 1.0)
            vol_m3  = width_m * length_m * INITIAL_DEPTH
            f_out.write(f"{seg_id} {vol_m3:.6f}\n")


# ==============================================================
# Entry
# ==============================================================
def run():
    print("\n=======================================================")
    print("  DHSVM STATE INITIALIZER (standalone)")
    print("=======================================================")
    generate_grid_states()
    generate_channel_state()
    print("=======================================================\n")


if __name__ == "__main__":
    run()
