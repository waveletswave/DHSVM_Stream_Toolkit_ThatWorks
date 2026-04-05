# -*- coding: utf-8 -*-
# =====================================================================
# soildepthscript.py (NumPy / DHSVM Binary Version - Boundary Aligned)
# PURPOSE:
#   Computes soil depth based on slope, source area (flow acc), and elevation.
#   Uses the exact PNNL weighting formula. Forces the DEM mask to prevent 
#   boundary shrinkage from upstream 3x3 raster operations, ensuring 
#   perfect pixel alignment for the DHSVM C-engine.
#
# AUTHOR: Yiyun Song
# DATE: 2026-04-01
# =====================================================================

import os
import numpy as np
from osgeo import gdal
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR.parent

def p(*xs) -> str:
    return str(WS.joinpath(*xs))

# =====================================================================
# PNNL Soil Depth Parameters (Southern Appalachians Calibration)
# =====================================================================

# Hydrologically active depth calibration based on Variable Source Area theory.
MIN_DEPTH = 1.6   # Must be > 1.5m (total root zone depth in .dhs)
MAX_DEPTH = 5.0   # Active subsurface stormflow zone

# # Hydrologically active depth calibration based on Variable Source Area (VSA) theory 
# # and Coweeta Hydrologic Laboratory field observations.

# # Minimum depth constrained by the basin-wide average soil thickness (~3.0m), 
# # allowing for localized geomorphic thinning on steep, erosional ridges (Swank & Crossley, 1988).
# MIN_DEPTH = 2.0   

# # Maximum active depth corresponds to the average depth of unweathered bedrock 
# # (saprolite thickness) observed in the Coweeta basin valleys (Swank & Douglass, 1975).
# MAX_DEPTH = 6.0

WT_SLOPE  = 0.7
WT_SOURCE = 0.0
WT_ELEV   = 0.3

MAX_SLOPE  = 30.0     
MAX_SOURCE = 100000.0 
MAX_ELEV   = 1500.0   

POW_SLOPE  = 0.25
POW_SOURCE = 1.0
POW_ELEV   = 0.75

DHSVM_NODATA = -9999.0

# =====================================================================
# Helper & Main Functions
# =====================================================================
def read_raster(filepath):
    ds = gdal.Open(filepath)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float32)
    return arr, band.GetNoDataValue(), ds.GetGeoTransform(), ds.GetProjection(), ds.RasterXSize, ds.RasterYSize

def generate_soildepth(elev_path, slope_path, flowacc_path, out_bin, out_tif):
    print("  [soildepth] Reading spatial inputs into memory...")
    elev_arr, elev_nd, gt, proj, cols, rows = read_raster(elev_path)
    slope_arr, slope_nd, _, _, _, _ = read_raster(slope_path)
    fac_arr, fac_nd, _, _, _, _ = read_raster(flowacc_path)

    # 1. Absolute Master Mask: Derived purely from the DEM
    valid_mask = (elev_arr != elev_nd) & (~np.isnan(elev_arr))

    # 2. Gap-filling: Force missing boundary pixels in slope/fac to 0.0 to match DEM mask
    print("  [soildepth] Aligning boundaries and sanitizing edge pixels...")
    elev_clean = np.where(valid_mask & (elev_arr >= 0), elev_arr, 0.0)
    slope_clean = np.where(valid_mask & (slope_arr != slope_nd) & (~np.isnan(slope_arr)) & (slope_arr >= 0), slope_arr, 0.0)
    fac_clean = np.where(valid_mask & (fac_arr != fac_nd) & (~np.isnan(fac_arr)), np.abs(fac_arr), 0.0)

    # 3. Limit processing bounds
    elev_lim = np.clip(elev_clean, 0.0, MAX_ELEV)
    slope_lim = np.clip(slope_clean, 0.0, MAX_SLOPE)
    fac_lim = np.clip(fac_clean, 0.0, MAX_SOURCE)

    # 4. Physical Calculation
    print("  [soildepth] Computing weighted topograpic functions...")
    term_slope  = WT_SLOPE * (1.0 - (slope_lim / MAX_SLOPE) ** POW_SLOPE)
    term_source = WT_SOURCE * ((fac_lim / MAX_SOURCE) ** POW_SOURCE)
    term_elev   = WT_ELEV * (1.0 - (elev_lim / MAX_ELEV) ** POW_ELEV)
    
    depth_calc = MIN_DEPTH + (MAX_DEPTH - MIN_DEPTH) * (term_slope + term_source + term_elev)
    
    # 5. Output array construction
    final_depth = np.full(elev_arr.shape, DHSVM_NODATA, dtype=np.float32)
    final_depth[valid_mask] = np.clip(depth_calc[valid_mask], MIN_DEPTH, MAX_DEPTH)

    print(f"  [soildepth] Writing flat binary -> {os.path.basename(out_bin)}")
    final_depth.tofile(out_bin)
    
    print(f"  [soildepth] Writing GeoTIFF -> {os.path.basename(out_tif)}")
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(out_tif, cols, rows, 1, gdal.GDT_Float32)
    out_ds.SetGeoTransform(gt)
    out_ds.SetProjection(proj)
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(final_depth)
    out_band.SetNoDataValue(DHSVM_NODATA)
    out_band.FlushCache()
    out_ds = None  # Safe memory release
    print("  [soildepth] Complete.")

# =====================================================================
# Execution Block
# =====================================================================
if __name__ != 'soildepthscript':
    print("\n--- Initializing Dynamic Soil Depth Module ---")
    dem_candidates = [p("Reprojected_DEM", "elev_clipped.tif"), p("Reprojected_DEM", "Cropped_DEM.tif"), p("dem.tif")]
    slope_candidates = [p("stream_slope.tif")]
    fac_candidates = [p("flow_acc.tif")]

    dem_in = next((c for c in dem_candidates if os.path.exists(c)), None)
    slope_in = next((c for c in slope_candidates if os.path.exists(c)), None)
    fac_in = next((c for c in fac_candidates if os.path.exists(c)), None)

    if dem_in and slope_in and fac_in:
        generate_soildepth(dem_in, slope_in, fac_in, p("soildepth.bin"), p("soildepth.tif"))
    else:
        print("[error] Missing input rasters. Ensure upstream tools completed successfully.")