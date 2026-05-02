# -*- coding: utf-8 -*-
# =====================================================================
# dem_to_dhsvm_bins.py
#
# PURPOSE:
#   Directly generates DHSVM binary input grids (DEM, Mask, Soil, Veg)
#   from the foundational clipped DEM. This script completely bypasses 
#   all legacy ESRI ASCII exports and C-engine (myconvert) steps, 
#   significantly accelerating and securing the preprocessing pipeline.
#
# EXPERIMENTAL DESIGN (Soil Depth Baselines):
#   Generates a suite of uniform soil depth baselines (e.g., 2m, 2.5m) 
#   for model sensitivity analysis. These utilize the exact same spatial 
#   mask to guarantee strict C-engine pixel alignment, serving as control 
#   groups against the dynamically routed topographic soil depth.
#
# PHYSICAL ASSUMPTION:
#   Follows the "One Soil / One Veg" uniformity simplification, which 
#   is optimal and standard for small experimental watersheds.
#
# AUTHOR:
#   Yiyun Song
#
# DATE:
#   2026-04-06
# =====================================================================

import os
import numpy as np
from osgeo import gdal
from pathlib import Path

# ==============================================================
# Configuration & Spatial Paths
# ==============================================================
# Automatically resolve the working directory
SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR.parent

# Input DEM (The foundational spatial framework)
DEM_IN  = WS / "Reprojected_DEM" / "elev_clipped.tif"

# Output Directory (Where all binaries will be stored for DHSVM)
OUT_DIR = WS / "DHSVM_input_binaries"

# Standard NoData value required by the DHSVM C-engine
DHSVM_NODATA = -9999.0

# --------------------------------------------------------------
# EXPERIMENTAL DESIGN: Uniform Soil Depth Baselines (in meters)
# --------------------------------------------------------------
# These uniform depths will be generated alongside the base maps
# to serve as baselines for subsurface stormflow sensitivity tests.
UNIFORM_SOIL_DEPTHS = [2.0, 2.5, 3.0, 3.5, 4.0]

# ==============================================================
# Core Generation Function
# ==============================================================
def generate_basemaps():
    print("\n--- Initializing DHSVM Base Map Generation Pipeline ---")
    
    if not DEM_IN.exists():
        print(f"[ERROR] Foundational DEM not found at: {DEM_IN}")
        print("        Please ensure upstream clipping tools completed successfully.")
        return
        
    os.makedirs(OUT_DIR, exist_ok=True)
    
    # 1. Read the foundational DEM using GDAL
    print(f"[step] Reading master topography: {DEM_IN.name}")
    ds = gdal.Open(str(DEM_IN))
    band = ds.GetRasterBand(1)
    dem_arr = band.ReadAsArray().astype(np.float32)
    dem_nd = band.GetNoDataValue()
    
    # 2. Define the Absolute Master Mask
    # Valid pixels are those that are strictly numeric and not NoData
    # This single mask enforces perfect pixel alignment across all layers.
    valid_mask = (dem_arr != dem_nd) & (~np.isnan(dem_arr))
    
    # 3. Memory Allocation for Standard DHSVM Arrays
    print("[step] Structuring physical memory arrays (Float32 and Int8)...")
    
    # DEM Array: Float32, strict NoData masking
    dem_out = np.full(dem_arr.shape, DHSVM_NODATA, dtype=np.float32)
    dem_out[valid_mask] = dem_arr[valid_mask]
    
    # MASK Array: Int8 (1 inside basin, 0 outside)
    mask_out = np.zeros(dem_arr.shape, dtype=np.int8)
    mask_out[valid_mask] = 1
    
    # SOIL Array: Int8 (Uniform type 1 inside basin, 0 outside)
    soil_out = np.zeros(dem_arr.shape, dtype=np.int8)
    soil_out[valid_mask] = 1
    
    # VEG Array: Int8 (Uniform type 1 inside basin, 0 outside)
    veg_out = np.zeros(dem_arr.shape, dtype=np.int8)
    veg_out[valid_mask] = 1
    
    # 4. Direct Flat Binary (.bin) Serialization (Base Maps)
    print(f"[step] Serializing standard base maps to: {OUT_DIR.name}/")
    
    dem_out.tofile(OUT_DIR / "dem.bin")
    mask_out.tofile(OUT_DIR / "mask.bin")
    soil_out.tofile(OUT_DIR / "soil.bin")
    veg_out.tofile(OUT_DIR / "veg.bin")
    
    print("  -> Generated: dem.bin  (Float32)")
    print("  -> Generated: mask.bin (Int8)")
    print("  -> Generated: soil.bin (Int8)")
    print("  -> Generated: veg.bin  (Int8)")

    # 5. Experimental Suites: Uniform Soil Depths
    print(f"\n[step] Generating uniform soil depth baselines...")
    for depth in UNIFORM_SOIL_DEPTHS:
        depth_out = np.full(dem_arr.shape, DHSVM_NODATA, dtype=np.float32)
        depth_out[valid_mask] = depth
        
        filename = f"soildepth_uniform_{depth:.1f}m.bin"
        depth_out.tofile(OUT_DIR / filename)
        print(f"  -> Generated: {filename} (Float32)")
    
    # Safely release GDAL file handlers from memory
    ds = None
    print("\n[SUCCESS] DHSVM spatial base maps and depth baselines successfully generated!")

# ==============================================================
# Execution Trigger
# ==============================================================
# Note: When executed via `exec()` in QGIS Python Console, the standard 
# `__name__ == '__main__'` check evaluates to False. Calling the function 
# directly ensures execution across all environments and master pipelines.
generate_basemaps()