# -*- coding: utf-8 -*-
# =====================================================================
# generate_dhsvm_states.py (Pure Python Unified Version - Final)
# PURPOSE:
#   100% Python replacement for the legacy MakeModelStateBin.c and 
#   MakeChannelState.sh. Generates physically consistent DHSVM initial 
#   states directly as flat binary matrices and routed text files.
#   
#   *FIX*: Explicitly appends '.bin' to grid state files to satisfy 
#   DHSVM's internal ReadModelState hardcoded path concatenation.
#
# AUTHOR: Yiyun Song
# DATE: 2026-04-06
# =====================================================================

import os
import numpy as np
from pathlib import Path
from osgeo import gdal

# ==============================================================
# Configuration & Paths
# ==============================================================
SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR.parent

DIR_INTERMED = WS / "Intermediate_GIS"
DIR_STREAMS  = WS / "DHSVM_input_streams"
STATE_DIR    = WS / "modelstate"

DEM_TIF = DIR_INTERMED / "elev_clipped.tif"

# Time & Channel Init
MODEL_START   = "01/01/2017-00"     # Format: MM/DD/YYYY-HH
INITIAL_DEPTH = 0.05                # Adjusted down to 5cm to prevent initial shocks

# Uniform State Configuration 
# Note: Lowered SOIL_MOISTURE slightly to reduce the "spin-up" overestimation
N_VEG_LAYERS          = 2
RAIN_INTERCEPTION     = [0.0, 0.0]
SNOW_INTERCEPTION_TOP = 0.0
TEMP_INT_STORAGE      = 0.0

SNOW_MASK             = 0.0         # C-engine casts all snow vars to float32
DAYS_SINCE_LAST_SNOW  = 0.0
SNOW_WATER_EQUIV      = 0.0
LW_BOTTOM             = 0.0
T_BOTTOM              = 0.0
LW_TOP                = 0.0
T_TOP                 = 0.0
COLD_CONTENT          = 0.0

N_SOIL_LAYERS         = 3
SOIL_MOISTURE         = [0.25, 0.25, 0.25, 0.25] # Tuned down from 0.30
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
    """Generates a NY*NX uniform float32 array and appends it to the binary file."""
    arr = np.full((ny, nx), value, dtype=np.float32)
    arr.tofile(file_handle)

# ==============================================================
# 1. Grid States (Replaces MakeModelStateBin)
# ==============================================================
def generate_grid_states():
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    ny, nx = get_grid_dimensions()
    
    date_str = MODEL_START.replace('/', '.').replace('-', '.') + ".00.00"
    print(f"\n[step] Generating binary grid states for {date_str} (NY={ny}, NX={nx})")
    
    # --- Interception State (Requires .bin for DHSVM ReadModelState) ---
    int_file = STATE_DIR / f"Interception.State.{date_str}.bin"
    with open(int_file, "wb") as f:
        for val in RAIN_INTERCEPTION: write_2d_float32(f, val, ny, nx)
        for _ in range(N_VEG_LAYERS): write_2d_float32(f, SNOW_INTERCEPTION_TOP, ny, nx)
        write_2d_float32(f, TEMP_INT_STORAGE, ny, nx)
    print(f"  -> {int_file.name} (pure Python binary)")

    # --- Snow State (Requires .bin) ---
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
    print(f"  -> {snow_file.name} (pure Python binary)")

    # --- Soil State (Requires .bin) ---
    soil_file = STATE_DIR / f"Soil.State.{date_str}.bin"
    with open(soil_file, "wb") as f:
        for val in SOIL_MOISTURE: write_2d_float32(f, val, ny, nx)
        write_2d_float32(f, SOIL_TEMP_SURF, ny, nx)
        for val in SOIL_TEMP_LAYERS: write_2d_float32(f, val, ny, nx)
        write_2d_float32(f, GROUND_HEAT, ny, nx)
        write_2d_float32(f, INITIAL_RUNOFF, ny, nx)
    print(f"  -> {soil_file.name} (pure Python binary)")

# ==============================================================
# 2. Channel States (Replaces MakeChannelState.sh)
# ==============================================================
def generate_channel_state():
    class_file   = DIR_STREAMS / "stream.class.dat"
    network_file = DIR_STREAMS / "stream.network.dat"
    date_str     = MODEL_START.replace('/', '.').replace('-', '.') + ".00.00"
    
    # Note: Channel State is a TEXT file and does NOT need .bin appended!
    out_file     = STATE_DIR / f"Channel.State.{date_str}"
    
    print(f"[step] Generating physics-based channel states -> {out_file.name}")
    
    class_widths = {}
    with open(class_file, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 2:
                class_widths[int(parts[0])] = float(parts[1])

    with open(network_file, 'r') as f_in, open(out_file, 'w') as f_out:
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
# Trigger (Optimized for QGIS & importlib Pipeline)
# ==============================================================
print("\n=======================================================")
print("  DHSVM PURE PYTHON STATE INITIALIZER")
print("=======================================================")
generate_grid_states()
generate_channel_state()
print("=======================================================\n")