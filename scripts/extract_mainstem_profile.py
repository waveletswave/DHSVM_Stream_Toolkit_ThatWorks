# -*- coding: utf-8 -*-
# =====================================================================
# extract_mainstem_profile.py
# PURPOSE:
#   Extract the main stem from the stream network, plot its 
#   longitudinal profile, detect knickpoints, and export data.
#
# RUN INSIDE: QGIS Python Console
#
# INPUT (auto-discovered):
#   Workspace (WS)    = parent folder of this script
#   DEM (projected)   = WS/Reprojected_DEM/* or WS/*
#   Stream Vector     = WS/streamfile.shp or WS/stream_order.shp
#
# OUTPUT:
#   stream_profile.png (saved to Workspace)
#   stream_profile_data.csv (saved to Workspace)
#
# AUTHOR: Y. Song
# =====================================================================

import os
import math
import csv
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

from qgis.core import (
    QgsVectorLayer, QgsRasterLayer, QgsPointXY, 
    QgsWkbTypes, QgsSpatialIndex, QgsRaster
)

# =====================================================================
#  Portable path configuration
# =====================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR.parent

def p(*xs) -> str:
    return str(WS.joinpath(*xs))

# =====================================================================
#  Geometry & Raster Helpers
# =====================================================================
def _line_coords(geom):
    """Extract coordinate sequence from geometry."""
    if geom.isMultipart():
        parts = geom.asMultiPolyline()
        if not parts: return []
        return max(parts, key=lambda L: len(L))
    else:
        return geom.asPolyline()

def _sample_raster_val(rl, x, y):
    """Sample band 1 value at a specific point with robust fallbacks."""
    provider = rl.dataProvider()
    pt = QgsPointXY(x, y)
    
    # Method 1: Standard sample
    val, ok = provider.sample(pt, 1)
    if ok and val is not None:
        try:
            fv = float(val)
            if math.isfinite(fv) and fv > -9000: 
                return fv
        except Exception:
            pass
            
    # Method 2: Fallback to identify tool
    ident = provider.identify(pt, QgsRaster.IdentifyFormatValue)
    if ident.isValid():
        res = ident.results()
        if 1 in res and res[1] is not None:
            try:
                fv = float(res[1])
                if math.isfinite(fv) and fv > -9000:
                    return fv
            except Exception:
                pass
                
    return None

# =====================================================================
#  Main Profiling Logic
# =====================================================================
def plot_stream_profile(stream_vec_path, dem_path, output_png, output_csv):
    vl = QgsVectorLayer(str(stream_vec_path), "streams", "ogr")
    rl = QgsRasterLayer(str(dem_path), "dem")
    
    if not vl.isValid() or not rl.isValid():
        print("[error] Failed to load vector or raster layers.")
        return

    features = list(vl.getFeatures())
    if not features:
        print("[warn] Stream vector is empty.")
        return

    print("[step] Building spatial index and topology...")
    sidx = QgsSpatialIndex()
    fid_to_feat = {}
    up_endpoints = {}
    down_endpoints = {}
    
    for ft in features:
        geom = ft.geometry()
        if not geom or geom.isEmpty():
            continue
            
        line = _line_coords(geom)
        if len(line) < 2:
            continue
            
        sidx.addFeature(ft)
        fid_to_feat[ft.id()] = ft
        up_endpoints[ft.id()] = line[0]
        down_endpoints[ft.id()] = line[-1]

    print("[step] Identifying main stem via maximum contributing area...")
    outlet_feat = max(features, key=lambda f: float(f.attribute("meanmsq")) if f.attribute("meanmsq") else 0.0)
    
    main_stem_fids = [outlet_feat.id()]
    current_fid = outlet_feat.id()
    
    px = abs(rl.rasterUnitsPerPixelX())
    py = abs(rl.rasterUnitsPerPixelY())
    search_tol = math.hypot(px, py) * 1.5

    while True:
        cur_start = up_endpoints[current_fid]
        candidate_ids = sidx.nearestNeighbor(cur_start, 10)
        upstream_cands = []
        
        for cid in candidate_ids:
            if cid == current_fid or cid in main_stem_fids:
                continue
            cand_down = down_endpoints.get(cid)
            if cand_down:
                dist = math.hypot(cand_down.x() - cur_start.x(), cand_down.y() - cur_start.y())
                if dist <= search_tol:
                    upstream_cands.append(cid)
        
        if not upstream_cands:
            break
            
        best_cand = max(upstream_cands, key=lambda cid: float(fid_to_feat[cid].attribute("meanmsq") or 0.0))
        main_stem_fids.append(best_cand)
        current_fid = best_cand

    main_stem_fids.reverse()
    
    print("[step] Extracting elevations along the main stem...")
    distances = []
    elevations = []
    cumulative_dist = 0.0
    missing_points = 0
    
    for fid in main_stem_fids:
        geom = fid_to_feat[fid].geometry()
        line = _line_coords(geom)
            
        for i, pt in enumerate(line):
            if i > 0:
                prev_pt = line[i-1]
                step = math.hypot(pt.x() - prev_pt.x(), pt.y() - prev_pt.y())
                cumulative_dist += step
                
            z = _sample_raster_val(rl, pt.x(), pt.y())
            
            if z is not None:
                distances.append(cumulative_dist)
                elevations.append(z)
            else:
                missing_points += 1

    if missing_points > 0:
        print(f"[warn] {missing_points} points were skipped because they fell on NoData or outside the DEM.")

    if len(elevations) < 5:
        print("[error] Not enough valid elevations read to compute profile.")
        return

    # =================================================================
    # Knickpoint Detection Math
    # =================================================================
    print("[step] Smoothing data and detecting knickpoints...")
    dist_arr = np.array(distances)
    elev_arr = np.array(elevations)
    
    # 1. Apply a moving average to elevation to filter raster noise
    window = 5 
    pad_width = window // 2
    elev_padded = np.pad(elev_arr, (pad_width, pad_width), mode='edge')
    elev_smooth = np.convolve(elev_padded, np.ones(window)/window, mode='valid')
    
    # 2. Calculate local slope (-dz/dx)
    dx = np.gradient(dist_arr)
    dz = np.gradient(elev_smooth)
    dx[dx == 0] = 1e-6 
    slope = -dz / dx 
    
    # 3. Identify statistically significant peaks in slope
    mean_slope = np.mean(slope)
    std_slope = np.std(slope)
    # Define threshold as 1.5 standard deviations above the mean slope
    threshold = mean_slope + 1.5 * std_slope 
    
    knickpoints_idx = []
    is_knickpoint = [False] * len(dist_arr)
    
    for i in range(1, len(slope) - 1):
        if slope[i] > threshold and slope[i] > slope[i-1] and slope[i] > slope[i+1]:
            knickpoints_idx.append(i)
            is_knickpoint[i] = True
            
    print(f"[info] Detected {len(knickpoints_idx)} potential knickpoints.")

    # =================================================================
    # Export to CSV
    # =================================================================
    print("[step] Exporting profile data to CSV...")
    with open(output_csv, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Distance_m', 'Elevation_m', 'Smoothed_Slope', 'Is_Knickpoint'])
        for d, e, s, k in zip(dist_arr, elev_arr, slope, is_knickpoint):
            writer.writerow([f"{d:.3f}", f"{e:.3f}", f"{s:.5f}", k])
    print(f"[ok] Data successfully saved to: {output_csv}")
    
    # =================================================================
    # Plotting
    # =================================================================
    print("[step] Generating plot...")
    plt.figure(figsize=(10, 5))
    
    # 1. Plot a faint base line to show the overall stream path
    plt.plot(dist_arr, elev_arr, color='gray', linewidth=1, alpha=0.5, zorder=1, label='Main Stem Path')
    
    # 2. Plot scatter points with a color gradient based on slope (coolwarm: blue=flat, red=steep)
    sc = plt.scatter(dist_arr, elev_arr, c=slope, cmap='coolwarm', s=25, zorder=3, alpha=0.9)
    
    # 3. Add a colorbar to explain the slope mapping
    cbar = plt.colorbar(sc)
    cbar.set_label('Local Slope (-dz/dx)', fontsize=11)
    
    # 4. Mark knickpoints (using large black empty circles to contrast with red steep points)
    if knickpoints_idx:
        plt.scatter(dist_arr[knickpoints_idx], elev_arr[knickpoints_idx], 
                    facecolors='none', edgecolors='black', s=80, linewidth=1.5, 
                    zorder=5, label='Knickpoints')
    
    plt.title('Main Stem Longitudinal Profile', fontsize=14, fontweight='bold')
    plt.xlabel('Distance from Headwater (m)', fontsize=12)
    plt.ylabel('Elevation (m)', fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.7)
    
    # Fill the area below the profile for visual depth
    bottom_limit = min(elev_arr) - 10
    plt.fill_between(dist_arr, elev_arr, bottom_limit, color='#e0e0e0', alpha=0.3, zorder=0)
    plt.ylim(bottom=bottom_limit)
    
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(str(output_png), dpi=300)
    plt.show()
    print(f"[ok] Profile plot successfully saved to: {output_png}")

# =====================================================================
#  Execution Configuration
# =====================================================================
stream_candidates = [
    p("streamfile.shp"),
    p("stream_order.shp")
]

stream_shp = None
for cand in stream_candidates:
    if os.path.exists(cand):
        temp_layer = QgsVectorLayer(cand, "test", "ogr")
        if temp_layer.isValid() and temp_layer.featureCount() > 0:
            stream_shp = cand
            break

dem_candidates = [
    p("Reprojected_DEM", "elev_clipped.tif"),
    p("Reprojected_DEM", "Cropped_DEM.tif"),
    p("Reprojected_DEM", "dem.tif"),
    p("dem.tif"),
    p("elev.tif")
]
dem_tif = next((c for c in dem_candidates if os.path.exists(c)), None)

output_png = p("stream_profile.png")
output_csv = p("stream_profile_data.csv")

if stream_shp and dem_tif:
    print(f"[info] Stream Vector: {stream_shp}")
    print(f"[info] DEM Raster: {dem_tif}")
    plot_stream_profile(stream_shp, dem_tif, output_png, output_csv)
else:
    print("[error] Could not locate a valid stream shapefile (with features) or DEM.")