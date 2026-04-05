# -*- coding: utf-8 -*-
# =====================================================================
# dem_to_dhsvm_bins.py
#
# PURPOSE:
#   Directly generates DHSVM binary input files (DEM, Mask, Soil, Veg)
#   from the foundational clipped DEM. This script bypasses all legacy 
#   intermediate ESRI ASCII exports and C-engine (myconvert) steps, 
#   significantly optimizing the preprocessing pipeline.
#
#   Physical Assumption: Follows the "One Soil / One Veg" uniformity 
#   simplification appropriate for small experimental watersheds.
#
# AUTHOR:
#   Yiyun Song
#
# DATE:
#   2026-04-04
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
    valid_mask = (dem_arr != dem_nd) & (~np.isnan(dem_arr))
    
    # 3. Memory Allocation for DHSVM Arrays
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
    
    # 4. Direct Flat Binary (.bin) Serialization
    print(f"[step] Serializing flat binary data to: {OUT_DIR.name}/")
    
    dem_bin_path  = OUT_DIR / "dem.bin"
    mask_bin_path = OUT_DIR / "mask.bin"
    soil_bin_path = OUT_DIR / "soil.bin"
    veg_bin_path  = OUT_DIR / "veg.bin"
    
    dem_out.tofile(dem_bin_path)
    mask_out.tofile(mask_bin_path)
    soil_out.tofile(soil_bin_path)
    veg_out.tofile(veg_bin_path)
    
    print(f"  -> Generated: {dem_bin_path.name}  (Float32)")
    print(f"  -> Generated: {mask_bin_path.name} (Int8)")
    print(f"  -> Generated: {soil_bin_path.name} (Int8)")
    print(f"  -> Generated: {veg_bin_path.name}  (Int8)")
    
    # Safely release GDAL file handlers from memory
    ds = None
    print("\n[SUCCESS] DHSVM spatial base maps successfully generated and aligned!")

# ==============================================================
# Execution Trigger
# ==============================================================
# Note: When executed via `exec()` in QGIS Python Console, the standard 
# `__name__ == '__main__'` check evaluates to False. Calling the function 
# directly ensures execution across all environments.
generate_basemaps()