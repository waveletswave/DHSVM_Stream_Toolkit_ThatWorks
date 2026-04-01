# -*- coding: utf-8 -*-
# =====================================================================
# soildepthscript.py (NumPy / DHSVM Binary Version - Auto Discovery)
# PURPOSE:
#   Computes soil depth based on slope, source area (flow acc), and elevation.
#   Uses the exact PNNL weighting formula.
#   Outputs a pure flat binary (.bin) file required by the DHSVM C-engine,
#   plus a (.tif) file for visual inspection in QGIS.
#
# RUN INSIDE: QGIS Python Console
#
# INPUT (Auto-discovered):
#   Workspace (WS)    = parent folder of this script
#   DEM Raster        = WS/Reprojected_DEM/elev_clipped.tif
#   Slope Raster      = WS/stream_slope.tif (from createstreamnetwork)
#   Flow Acc Raster   = WS/flow_acc.tif (from createstreamnetwork)
#
# OUTPUT:
#   soildepth.bin (Strict 32-bit float flat binary for DHSVM input)
#   soildepth.tif (For QGIS visualization)
#
# AUTHOR: Yiyun Song
# DATE: 2026-04-01
# =====================================================================

import os
import numpy as np
from osgeo import gdal
from pathlib import Path

# =====================================================================
#  Portable path configuration
# =====================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR.parent

def p(*xs) -> str:
    return str(WS.joinpath(*xs))

# =====================================================================
#  PNNL Soil Depth Parameters
# =====================================================================
MIN_DEPTH = 0.5
MAX_DEPTH = 3.0

WT_SLOPE  = 0.7
WT_SOURCE = 0.0
WT_ELEV   = 0.3

MAX_SLOPE  = 30.0     # Degrees
MAX_SOURCE = 100000.0 # Flow Accumulation limit
MAX_ELEV   = 1500.0   # Meters

POW_SLOPE  = 0.25
POW_SOURCE = 1.0
POW_ELEV   = 0.75

DHSVM_NODATA = -9999.0

# =====================================================================
#  Helper Functions
# =====================================================================
def read_raster(filepath):
    ds = gdal.Open(filepath)
    if not ds:
        raise FileNotFoundError(f"Could not open {filepath}")
    
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float32)
    nodata = band.GetNoDataValue()
    
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    
    return arr, nodata, gt, proj, cols, rows

# =====================================================================
#  Main Calculation Logic
# =====================================================================
def generate_soildepth(elev_path, slope_path, flowacc_path, out_bin, out_tif):
    print("[step] Reading input rasters into memory...")
    try:
        elev_arr, elev_nd, gt, proj, cols, rows = read_raster(elev_path)
        slope_arr, slope_nd, _, _, _, _ = read_raster(slope_path)
        fac_arr, fac_nd, _, _, _, _ = read_raster(flowacc_path)
    except Exception as e:
        print(f"[error] {e}")
        return

    print(f"[info] Parameters -> Min: {MIN_DEPTH}m, Max: {MAX_DEPTH}m")
    
    # Create valid data mask based on DEM NoData
    valid_mask = (elev_arr != elev_nd) & (~np.isnan(elev_arr))

    # Clean arrays and apply NoData thresholds for calculation
    elev_clean = np.where(valid_mask, elev_arr, 0.0)
    slope_clean = np.where(valid_mask, slope_arr, 0.0)
    fac_clean = np.where(valid_mask, np.abs(fac_arr), 0.0) # Ensure FlowAcc is positive

    print("[step] Applying max limits to parameters...")
    # Apply threshold limits A*(A<=max) + max*(A>max)
    elev_lim = np.clip(elev_clean, None, MAX_ELEV)
    slope_lim = np.clip(slope_clean, None, MAX_SLOPE)
    fac_lim = np.clip(fac_clean, None, MAX_SOURCE)

    print("[step] Computing PNNL weighted soil depth formula...")
    term_slope  = WT_SLOPE * (1.0 - (slope_lim / MAX_SLOPE) ** POW_SLOPE)
    term_source = WT_SOURCE * ((fac_lim / MAX_SOURCE) ** POW_SOURCE)
    term_elev   = WT_ELEV * (1.0 - (elev_lim / MAX_ELEV) ** POW_ELEV)
    
    depth_calc = MIN_DEPTH + (MAX_DEPTH - MIN_DEPTH) * (term_slope + term_source + term_elev)
    
    # Initialize final array with DHSVM NoData value
    final_depth = np.full(elev_arr.shape, DHSVM_NODATA, dtype=np.float32)
    
    # Inject calculated depths only into valid mask areas
    final_depth[valid_mask] = depth_calc[valid_mask]
    
    # Strict fallback clip for safety
    final_depth[valid_mask] = np.clip(final_depth[valid_mask], MIN_DEPTH, MAX_DEPTH)

    # EXPORT 1: DHSVM Flat Binary (.bin)
    print("[step] Writing Flat Binary (.bin) for DHSVM...")
    final_depth.tofile(out_bin)
    print(f"[ok] Binary file saved: {out_bin}")

    # EXPORT 2: GeoTIFF (.tif) for QGIS Visualization
    print("[step] Writing GeoTIFF (.tif) for visualization...")
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(out_tif, cols, rows, 1, gdal.GDT_Float32)
    out_ds.SetGeoTransform(gt)
    out_ds.SetProjection(proj)
    
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(final_depth)
    out_band.SetNoDataValue(DHSVM_NODATA)
    out_band.FlushCache()
    out_ds = None
    print(f"[ok] TIFF file saved: {out_tif}")

# =====================================================================
#  Execution Block (No __main__ lock for QGIS compatibility)
# =====================================================================
print("--- Starting Soil Depth Generation ---")

# 1. Auto-discover inputs generated by createstreamnetwork.py
dem_candidates = [p("Reprojected_DEM", "elev_clipped.tif"), p("Reprojected_DEM", "Cropped_DEM.tif"), p("dem.tif")]
slope_candidates = [p("stream_slope.tif")]
fac_candidates = [p("flow_acc.tif")]

dem_in = next((c for c in dem_candidates if os.path.exists(c)), None)
slope_in = next((c for c in slope_candidates if os.path.exists(c)), None)
fac_in = next((c for c in fac_candidates if os.path.exists(c)), None)

bin_out  = p("soildepth.bin")
tif_out  = p("soildepth.tif")

if dem_in and slope_in and fac_in:
    print(f"[ok] Found DEM: {dem_in}")
    print(f"[ok] Found Slope: {slope_in}")
    print(f"[ok] Found Flow Acc: {fac_in}")
    generate_soildepth(dem_in, slope_in, fac_in, bin_out, tif_out)
else:
    print("[error] Missing input rasters. Ensure createstreamnetwork.py has run successfully.")