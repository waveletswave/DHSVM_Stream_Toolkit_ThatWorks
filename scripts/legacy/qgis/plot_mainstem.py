# -*- coding: utf-8 -*-
# =====================================================================
# plot_mainstem.py
# PURPOSE:
#   Extract the main stem using DEM elevation to enforce physical flow 
#   direction. Calculates Upstream Accumulated Length as a proxy for 
#   Drainage Area. Plots Distance Upstream from Outlet and ensures 
#   clean, publication-ready figures without legend overlap.
#
# RUN INSIDE: QGIS Python Console
#
# INPUT:
#   Workspace (WS)    = parent folder of this script
#   DEM (projected)   = WS/Reprojected_DEM/* or WS/*
#   Stream Vector     = WS/streamfile.shp or WS/stream_order.shp
#
# OUTPUT:
#   stream_profile.png
#   stream_profile_data.csv
#
# AUTHOR: Y. Song
# DATE: 2026-02-27
# =====================================================================

import os
import math
import csv
import numpy as np
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt

from qgis.core import (
    QgsVectorLayer, QgsRasterLayer, QgsPointXY, 
    QgsSpatialIndex, QgsRaster
)

# =====================================================================
#  Workspace Configuration
# =====================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR.parent

def p(*xs) -> str:
    return str(WS.joinpath(*xs))

# =====================================================================
#  Geometry & Raster Helper Functions
# =====================================================================
def _line_coords(geom):
    """Extracts the coordinate sequence from a QGIS geometry object."""
    if geom.isMultipart():
        parts = geom.asMultiPolyline()
        if not parts: return []
        return max(parts, key=lambda L: len(L))
    else:
        return geom.asPolyline()

def _sample_raster_val(rl, x, y):
    """Robustly samples a raster value at a specific point."""
    provider = rl.dataProvider()
    pt = QgsPointXY(x, y)
    
    val, ok = provider.sample(pt, 1)
    if ok and val is not None:
        try:
            fv = float(val)
            if math.isfinite(fv) and fv > -9000: 
                return fv
        except Exception:
            pass
            
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
#  Main Analytical Workflow
# =====================================================================
def plot_stream_profile(stream_vec_path, dem_path, output_png, output_csv):
    """Executes the main stem extraction and plotting logic."""
    vl = QgsVectorLayer(str(stream_vec_path), "streams", "ogr")
    rl = QgsRasterLayer(str(dem_path), "dem")
    
    if not vl.isValid() or not rl.isValid():
        print("[error] Failed to load layers.")
        return

    features = list(vl.getFeatures())
    if not features:
        print("[warn] Stream vector is empty.")
        return

    print("[step] Extracting coordinates and establishing elevation routing...")
    sidx = QgsSpatialIndex()
    G = nx.DiGraph()
    endpoints_data = {}
    
    for ft in features:
        geom = ft.geometry()
        if not geom or geom.isEmpty(): continue
        line = _line_coords(geom)
        if len(line) < 2: continue
            
        pt_start = line[0]
        pt_end = line[-1]
        
        z_start = _sample_raster_val(rl, pt_start.x(), pt_start.y())
        z_end = _sample_raster_val(rl, pt_end.x(), pt_end.y())
        
        # Enforce physical routing: Water flows from higher Z to lower Z
        if z_start is not None and z_end is not None and z_start < z_end:
            up_pt = pt_end
            down_pt = pt_start
            line = list(reversed(line))
        else:
            up_pt = pt_start
            down_pt = pt_end
            
        seg_length = geom.length()
        G.add_node(ft.id(), length=seg_length, line=line)
        endpoints_data[ft.id()] = {'up': up_pt, 'down': down_pt}
        sidx.addFeature(ft)

    px = abs(rl.rasterUnitsPerPixelX())
    py = abs(rl.rasterUnitsPerPixelY())
    search_tol = math.hypot(px, py) * 1.5

    print("[step] Building network topology...")
    for fidA, dataA in endpoints_data.items():
        down_pt = dataA['down']
        candidate_ids = sidx.nearestNeighbor(down_pt, 10)
        
        for cid in candidate_ids:
            if cid == fidA: continue
            up_pt = endpoints_data[cid]['up']
            
            dist = math.hypot(down_pt.x() - up_pt.x(), down_pt.y() - up_pt.y())
            if dist <= search_tol:
                G.add_edge(fidA, cid)

    print("[step] Calculating Upstream Accumulated Length...")
    for n in G.nodes():
        ancestors = nx.ancestors(G, n)
        acc_len = sum(G.nodes[anc]['length'] for anc in ancestors) + G.nodes[n]['length']
        G.nodes[n]['acc_length'] = acc_len

    outlets = [n for n, d in G.out_degree() if d == 0]
    if not outlets:
        print("[error] No valid outlet found in the network topology.")
        return

    basin_outlet = max(outlets, key=lambda n: G.nodes[n]['acc_length'])
    headwaters = [n for n, d in G.in_degree() if d == 0]
    
    main_stem_fids = []
    max_dist = -1.0

    print("[step] Tracing the true physical main stem...")
    for hw in headwaters:
        if nx.has_path(G, hw, basin_outlet):
            path = nx.shortest_path(G, source=hw, target=basin_outlet)
            path_dist = sum(G.nodes[n]['length'] for n in path)
            
            if path_dist > max_dist:
                max_dist = path_dist
                main_stem_fids = path

    print(f"[info] Main stem identified: {len(main_stem_fids)} segments.")
    
    ordered_points_data = []
    for fid in main_stem_fids:
        line = G.nodes[fid]['line']
        acc_len = G.nodes[fid]['acc_length']
        
        for pt in line:
            if not ordered_points_data:
                ordered_points_data.append((pt, acc_len))
            else:
                last_pt = ordered_points_data[-1][0]
                dist_prev = math.hypot(pt.x() - last_pt.x(), pt.y() - last_pt.y())
                if dist_prev > 1e-3: 
                    ordered_points_data.append((pt, acc_len))

    print("[step] Sampling elevations...")
    downstream_distances = []
    elevations = []
    proxy_areas = []
    cumulative_dist = 0.0
    
    for i, (pt, acc_len) in enumerate(ordered_points_data):
        if i > 0:
            prev_pt = ordered_points_data[i-1][0]
            cumulative_dist += math.hypot(pt.x() - prev_pt.x(), pt.y() - prev_pt.y())
            
        z = _sample_raster_val(rl, pt.x(), pt.y())
        if z is not None:
            downstream_distances.append(cumulative_dist)
            elevations.append(z)
            proxy_areas.append(acc_len)

    if len(elevations) < 5:
        print("[error] Not enough valid elevations read.")
        return

    total_length = downstream_distances[-1]
    upstream_distances = [total_length - d for d in downstream_distances]

    print("[step] Smoothing data and analyzing topographic signals...")
    dist_arr = np.array(upstream_distances)
    elev_arr = np.array(elevations)
    area_arr = np.array(proxy_areas)
    
    window = 5 
    pad_width = window // 2
    elev_padded = np.pad(elev_arr, (pad_width, pad_width), mode='edge')
    elev_smooth = np.convolve(elev_padded, np.ones(window)/window, mode='valid')
    
    dx = np.gradient(dist_arr)
    dz = np.gradient(elev_smooth)
    dx[dx == 0] = 1e-6 
    slope = np.abs(dz / dx)
    
    mean_slope = np.mean(slope)
    std_slope = np.std(slope)
    threshold = mean_slope + 1.5 * std_slope 
    
    knickpoints_idx = []
    is_knickpoint = [False] * len(dist_arr)
    
    for i in range(1, len(slope) - 1):
        if slope[i] > threshold and slope[i] > slope[i-1] and slope[i] > slope[i+1]:
            knickpoints_idx.append(i)
            is_knickpoint[i] = True
            
    print(f"[info] Detected {len(knickpoints_idx)} robust knickpoints.")
    
    # =================================================================
    # Exporting & Visualization
    # =================================================================
    print("[step] Exporting and Plotting...")
    with open(output_csv, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Distance_Upstream_m', 'Elevation_m', 'Upstream_Length_m', 'Smoothed_Slope', 'Is_Knickpoint'])
        for d, e, a, s, k in zip(dist_arr, elev_arr, area_arr, slope, is_knickpoint):
            writer.writerow([f"{d:.3f}", f"{e:.3f}", f"{a:.1f}", f"{s:.5f}", k])
            
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    ax1.plot(dist_arr, elev_arr, color='gray', linewidth=1, alpha=0.5, zorder=1, label='Main Stem Path')
    sc = ax1.scatter(dist_arr, elev_arr, c=slope, cmap='coolwarm', s=25, zorder=3, alpha=0.9)
    
    if knickpoints_idx:
        ax1.scatter(dist_arr[knickpoints_idx], elev_arr[knickpoints_idx], 
                    facecolors='none', edgecolors='black', s=80, linewidth=1.5, 
                    zorder=5, label='Knickpoints')
    
    ax1.set_xlabel('Distance Upstream from Outlet (m)', fontsize=12)
    ax1.set_ylabel('Elevation (m)', fontsize=12, color='black')
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.grid(True, linestyle=':', alpha=0.7)
    
    bottom_limit = min(elev_arr) - 10
    ax1.fill_between(dist_arr, elev_arr, bottom_limit, color='#e0e0e0', alpha=0.3, zorder=0)
    ax1.set_ylim(bottom=bottom_limit)
    
    ax2 = ax1.twinx()
    ax2.plot(dist_arr, area_arr, color='forestgreen', linestyle=':', linewidth=2, alpha=0.8, label='Accumulated Upstream Length')
    ax2.set_ylabel('Upstream Accumulated Length (m)', fontsize=12, color='forestgreen')
    ax2.tick_params(axis='y', labelcolor='forestgreen')
    
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper center')
    
    cbar = fig.colorbar(sc, ax=ax2, pad=0.15)
    cbar.set_label('Local Slope (|dz/dx|)', fontsize=11)
    
    plt.title('Longitudinal Profile (Upstream from Outlet)', fontsize=14, fontweight='bold')
    
    fig.tight_layout()
    plt.savefig(str(output_png), dpi=300)
    plt.show()
    print(f"[ok] Outputs saved successfully.")

# =====================================================================
#  Execution Block
# =====================================================================
stream_candidates = [p("streamfile.shp"), p("stream_order.shp")]
stream_shp = next((c for c in stream_candidates if os.path.exists(c)), None)

dem_candidates = [
    p("Reprojected_DEM", "elev_clipped.tif"),
    p("Reprojected_DEM", "Cropped_DEM.tif"),
    p("Reprojected_DEM", "dem.tif"),
    p("dem.tif"), p("elev.tif")
]
dem_tif = next((c for c in dem_candidates if os.path.exists(c)), None)

output_png = p("stream_profile.png")
output_csv = p("stream_profile_data.csv")

if stream_shp and dem_tif:
    plot_stream_profile(stream_shp, dem_tif, output_png, output_csv)
else:
    print("[error] Could not locate a valid stream shapefile or DEM.")