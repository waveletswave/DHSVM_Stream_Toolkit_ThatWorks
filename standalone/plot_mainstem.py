# ⚠️ ============================================================
# ⚠️  STATUS: DRAFT / WORK-IN-PROGRESS — NOT YET VERIFIED
# ⚠️ ============================================================
# ⚠️  This standalone port has been written but has NOT been
# ⚠️  verified against qgis/plot_mainstem.py output. Do not use
# ⚠️  for production analyses until the byte-level CSV diff has
# ⚠️  been completed (see commit history for the verification
# ⚠️  protocol used for channelclass and rowcolmap).
# ⚠️
# ⚠️  TODO before declaring stable:
# ⚠️    1. Run qgis/plot_mainstem.py via QGIS Console -> /tmp/qgis_mainstem/
# ⚠️    2. Run this file via terminal -> /tmp/sa_mainstem/
# ⚠️    3. Diff stream_profile_data.csv (rows + max abs diff per column)
# ⚠️    4. Visual diff stream_profile.png
# ⚠️ ============================================================

# -*- coding: utf-8 -*-
# =====================================================================
# plot_mainstem.py (standalone / DCC-ready)
# =====================================================================
# Standalone (non-QGIS) sibling of qgis/plot_mainstem.py.
#
# Identifies the watershed outlet, traces the physically routed main stem
# upstream using DEM elevation, computes accumulated upstream length as a
# proxy for drainage area, detects knickpoints, and produces:
#   - stream_profile.png       (longitudinal profile figure)
#   - stream_profile_data.csv  (per-vertex distance, elevation, slope)
#
# Replaces QGIS dependencies with rasterio + geopandas + shapely:
#   QgsRasterLayer / provider.sample()  ->  rasterio.open() + src.sample()
#   QgsVectorLayer.getFeatures()        ->  geopandas.read_file()
#   QgsSpatialIndex                      ->  shapely.strtree.STRtree
#   QgsPointXY / QgsGeometry             ->  shapely Point / LineString
#
# Headless safety:
#   matplotlib uses the 'Agg' backend, which renders to file without
#   needing an X server. Required for HPC compute nodes (e.g. Duke DCC).
#
# CLI usage:
#   python plot_mainstem.py <stream.shp> <dem.tif> [--output-dir DIR]
#                                                  [--smoothing-window N]
#                                                  [--knickpoint-sigma F]
#
# Library usage:
#   from plot_mainstem import plot_stream_profile
#
# AUTHOR: Y. Song
# DATE:   2026-05-02
# =====================================================================

import os
import sys
import math
import csv
import argparse
from pathlib import Path

import numpy as np
import networkx as nx
import geopandas as gpd
import rasterio

# IMPORTANT: 'Agg' backend MUST be set before pyplot is imported.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from shapely.geometry import Point  # noqa: E402
from shapely.strtree import STRtree  # noqa: E402


# =====================================================================
# Geometry & raster helpers
# =====================================================================
def _line_coords(geom):
    """
    Extract a list of (x, y) coordinates from a (Multi)LineString.
    For MultiLineString, returns the longest part (matches QGIS version's
    `max(parts, key=lambda L: len(L))` behaviour).
    """
    if geom is None or geom.is_empty:
        return []
    if geom.geom_type == 'MultiLineString':
        parts = list(geom.geoms)
        if not parts:
            return []
        longest = max(parts, key=lambda p: len(p.coords))
        return [(c[0], c[1]) for c in longest.coords]
    if geom.geom_type == 'LineString':
        return [(c[0], c[1]) for c in geom.coords]
    return []


def _sample_raster_val(src, x, y, nodata):
    """
    Sample a single-band raster at world coords (x, y).
    Returns float or None for invalid / out-of-bounds / nodata.

    Mirrors the defensive logic of the QGIS version:
      - reject NaN
      - reject explicit nodata
      - reject sentinel values <= -9000 (catches unset nodata in some DEMs)
    """
    try:
        result = next(src.sample([(x, y)]))
        val = float(result[0])
    except (StopIteration, IndexError, ValueError, TypeError):
        return None

    if not math.isfinite(val):
        return None
    if nodata is not None and val == nodata:
        return None
    if val <= -9000:
        return None
    return val


# =====================================================================
# Main analytical workflow
# =====================================================================
def plot_stream_profile(stream_vec_path, dem_path,
                         output_png, output_csv,
                         smoothing_window=5,
                         knickpoint_sigma=1.5):
    """
    Identify the main stem, sample elevations along it, detect knickpoints,
    and produce a longitudinal profile PNG + CSV of per-vertex measurements.

    Parameters
    ----------
    stream_vec_path : str | Path
        Stream line layer (shapefile / GeoPackage).
    dem_path : str | Path
        DEM raster aligned with the stream layer.
    output_png : str | Path
        Output figure path.
    output_csv : str | Path
        Output CSV path.
    smoothing_window : int, default 5
        Moving-average window size used for slope/knickpoint detection.
    knickpoint_sigma : float, default 1.5
        Knickpoint threshold = mean(slope) + sigma * std(slope).
    """
    stream_vec_path = str(stream_vec_path)
    dem_path = str(dem_path)
    output_png = str(output_png)
    output_csv = str(output_csv)

    if not os.path.exists(stream_vec_path):
        print(f"[error] Cannot find stream layer: {stream_vec_path}")
        return
    if not os.path.exists(dem_path):
        print(f"[error] Cannot find DEM: {dem_path}")
        return

    # --- Load vector ---------------------------------------------------
    gdf = gpd.read_file(stream_vec_path)
    if len(gdf) == 0:
        print("[warn] Stream vector is empty.")
        return

    # All raster work happens within the rasterio context manager so file
    # handles close cleanly even if the function returns mid-way.
    with rasterio.open(dem_path) as src:
        nodata = src.nodata
        transform = src.transform

        # --- CRS handling: actively reproject (not just warn) ----------
        # `plot_mainstem.py` performs raster sampling, so a CRS mismatch
        # would silently produce garbage elevations. Reproject before
        # any sampling occurs.
        if gdf.crs is not None and src.crs is not None:
            try:
                same = (gdf.crs.to_wkt() == src.crs.to_wkt())
            except Exception:
                same = (str(gdf.crs) == str(src.crs))
            if not same:
                print("[reproject] Stream CRS != DEM CRS; "
                      "reprojecting stream layer to DEM CRS.")
                gdf = gdf.to_crs(src.crs)

        # --- Build directed graph with elevation routing ---------------
        print("[step] Establishing elevation-routed network topology...")
        G = nx.DiGraph()
        endpoints_data = {}

        for fid, row in gdf.iterrows():
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue
            line_coords = _line_coords(geom)
            if len(line_coords) < 2:
                continue

            pt_start = line_coords[0]
            pt_end   = line_coords[-1]
            z_start = _sample_raster_val(src, pt_start[0], pt_start[1], nodata)
            z_end   = _sample_raster_val(src, pt_end[0],   pt_end[1],   nodata)

            # Enforce physical routing: water flows from higher z to lower z
            if (z_start is not None and z_end is not None
                    and z_start < z_end):
                up_pt   = pt_end
                down_pt = pt_start
                line_coords = list(reversed(line_coords))
            else:
                up_pt   = pt_start
                down_pt = pt_end

            seg_length = geom.length
            G.add_node(fid, length=seg_length, line=line_coords)
            endpoints_data[fid] = {'up': up_pt, 'down': down_pt}

        if not endpoints_data:
            print("[error] No usable line segments found.")
            return

        # --- Build edges via STRtree on upstream endpoints -------------
        print("[step] Building network topology (spatial index on endpoints)...")
        px = abs(transform.a)
        py = abs(transform.e)
        search_tol = math.hypot(px, py) * 1.5

        fid_list  = list(endpoints_data.keys())
        up_points = [Point(endpoints_data[fid]['up']) for fid in fid_list]
        tree = STRtree(up_points)

        for fidA in fid_list:
            down_pt  = Point(endpoints_data[fidA]['down'])
            buffered = down_pt.buffer(search_tol)
            # shapely 2.x: STRtree.query() returns numpy array of indices
            cand_idx = tree.query(buffered)
            for idx in np.asarray(cand_idx).ravel():
                cid = fid_list[int(idx)]
                if cid == fidA:
                    continue
                if down_pt.distance(up_points[int(idx)]) <= search_tol:
                    G.add_edge(fidA, cid)

        # --- Compute upstream accumulated length per node --------------
        print("[step] Calculating upstream accumulated length...")
        for n in G.nodes():
            ancestors = nx.ancestors(G, n)
            acc_len = sum(G.nodes[anc]['length'] for anc in ancestors) \
                      + G.nodes[n]['length']
            G.nodes[n]['acc_length'] = acc_len

        # --- Identify outlet & main stem -------------------------------
        outlets = [n for n, d in G.out_degree() if d == 0]
        if not outlets:
            print("[error] No valid outlet found in the network topology.")
            return

        basin_outlet = max(outlets, key=lambda n: G.nodes[n]['acc_length'])
        headwaters   = [n for n, d in G.in_degree() if d == 0]

        main_stem_fids = []
        max_dist = -1.0
        print("[step] Tracing the physical main stem...")
        for hw in headwaters:
            if nx.has_path(G, hw, basin_outlet):
                path = nx.shortest_path(G, source=hw, target=basin_outlet)
                path_dist = sum(G.nodes[n]['length'] for n in path)
                if path_dist > max_dist:
                    max_dist = path_dist
                    main_stem_fids = path

        if not main_stem_fids:
            print("[error] Could not trace a path to the outlet.")
            return

        print(f"[info] Main stem identified: {len(main_stem_fids)} segments.")

        # --- Order vertices along the main stem ------------------------
        ordered_points_data = []
        for fid in main_stem_fids:
            line    = G.nodes[fid]['line']
            acc_len = G.nodes[fid]['acc_length']
            for pt in line:
                if not ordered_points_data:
                    ordered_points_data.append((pt, acc_len))
                else:
                    last_pt = ordered_points_data[-1][0]
                    dist_prev = math.hypot(pt[0] - last_pt[0],
                                           pt[1] - last_pt[1])
                    if dist_prev > 1e-3:
                        ordered_points_data.append((pt, acc_len))

        # --- Sample elevations along the main stem ---------------------
        print("[step] Sampling elevations along the main stem...")
        downstream_distances = []
        elevations           = []
        proxy_areas          = []
        cumulative_dist      = 0.0

        for i, (pt, acc_len) in enumerate(ordered_points_data):
            if i > 0:
                prev_pt = ordered_points_data[i - 1][0]
                cumulative_dist += math.hypot(pt[0] - prev_pt[0],
                                              pt[1] - prev_pt[1])
            z = _sample_raster_val(src, pt[0], pt[1], nodata)
            if z is not None:
                downstream_distances.append(cumulative_dist)
                elevations.append(z)
                proxy_areas.append(acc_len)

    # --- (rasterio src is closed past this point) ---------------------

    if len(elevations) < 5:
        print("[error] Not enough valid elevations were read.")
        return

    # --- Convert to upstream distance (matches QGIS version) -----------
    total_length = downstream_distances[-1]
    upstream_distances = [total_length - d for d in downstream_distances]

    # --- Smoothing & slope calculation ---------------------------------
    print("[step] Smoothing & analyzing topographic signals...")
    dist_arr = np.array(upstream_distances)
    elev_arr = np.array(elevations)
    area_arr = np.array(proxy_areas)

    window = int(smoothing_window)
    pad = window // 2
    elev_padded = np.pad(elev_arr, (pad, pad), mode='edge')
    elev_smooth = np.convolve(elev_padded, np.ones(window) / window,
                              mode='valid')

    dx = np.gradient(dist_arr)
    dz = np.gradient(elev_smooth)
    dx = np.where(dx == 0, 1e-6, dx)
    slope = np.abs(dz / dx)

    # --- Knickpoint detection ------------------------------------------
    mean_slope = float(np.mean(slope))
    std_slope  = float(np.std(slope))
    threshold  = mean_slope + knickpoint_sigma * std_slope

    knickpoints_idx = []
    is_knickpoint = [False] * len(dist_arr)
    for i in range(1, len(slope) - 1):
        if (slope[i] > threshold and slope[i] > slope[i - 1]
                and slope[i] > slope[i + 1]):
            knickpoints_idx.append(i)
            is_knickpoint[i] = True

    print(f"[info] Detected {len(knickpoints_idx)} robust knickpoints.")

    # --- Export CSV ----------------------------------------------------
    print(f"[step] Writing CSV  -> {output_csv}")
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)) or ".",
                exist_ok=True)
    with open(output_csv, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Distance_Upstream_m', 'Elevation_m',
                         'Upstream_Length_m', 'Smoothed_Slope',
                         'Is_Knickpoint'])
        for d, e, a, s, k in zip(dist_arr, elev_arr, area_arr,
                                 slope, is_knickpoint):
            writer.writerow([f"{d:.3f}", f"{e:.3f}",
                             f"{a:.1f}", f"{s:.5f}", k])

    # --- Plot ----------------------------------------------------------
    print(f"[step] Writing PNG  -> {output_png}")
    os.makedirs(os.path.dirname(os.path.abspath(output_png)) or ".",
                exist_ok=True)

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(dist_arr, elev_arr, color='gray', linewidth=1, alpha=0.5,
             zorder=1, label='Main Stem Path')
    sc = ax1.scatter(dist_arr, elev_arr, c=slope, cmap='coolwarm',
                     s=25, zorder=3, alpha=0.9)

    if knickpoints_idx:
        ax1.scatter(dist_arr[knickpoints_idx], elev_arr[knickpoints_idx],
                    facecolors='none', edgecolors='black',
                    s=80, linewidth=1.5, zorder=5, label='Knickpoints')

    ax1.set_xlabel('Distance Upstream from Outlet (m)', fontsize=12)
    ax1.set_ylabel('Elevation (m)', fontsize=12, color='black')
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.grid(True, linestyle=':', alpha=0.7)

    bottom_limit = float(np.min(elev_arr) - 10)
    ax1.fill_between(dist_arr, elev_arr, bottom_limit,
                     color='#e0e0e0', alpha=0.3, zorder=0)
    ax1.set_ylim(bottom=bottom_limit)

    ax2 = ax1.twinx()
    ax2.plot(dist_arr, area_arr, color='forestgreen', linestyle=':',
             linewidth=2, alpha=0.8,
             label='Accumulated Upstream Length')
    ax2.set_ylabel('Upstream Accumulated Length (m)', fontsize=12,
                   color='forestgreen')
    ax2.tick_params(axis='y', labelcolor='forestgreen')

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper center')

    cbar = fig.colorbar(sc, ax=ax2, pad=0.15)
    cbar.set_label('Local Slope (|dz/dx|)', fontsize=11)

    plt.title('Longitudinal Profile (Upstream from Outlet)',
              fontsize=14, fontweight='bold')

    fig.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close(fig)

    print(f"[ok] Outputs saved successfully.")


# =====================================================================
# CLI
# =====================================================================
def _build_parser():
    p = argparse.ArgumentParser(
        description=("Trace the physical main stem from a stream network "
                     "and DEM, detect knickpoints, and produce a "
                     "longitudinal profile (PNG) + per-vertex CSV.")
    )
    p.add_argument("stream",
                   help="Stream line layer (.shp / .gpkg / .geojson)")
    p.add_argument("dem",
                   help="DEM raster (.tif), should cover the stream extent")
    p.add_argument("--output-dir", default=None,
                   help="Directory for outputs (default: current directory)")
    p.add_argument("--png", default=None,
                   help="Output PNG path (overrides --output-dir for PNG)")
    p.add_argument("--csv", default=None,
                   help="Output CSV path (overrides --output-dir for CSV)")
    p.add_argument("--smoothing-window", type=int, default=5,
                   help="Moving-average window size (default: 5)")
    p.add_argument("--knickpoint-sigma", type=float, default=1.5,
                   help="Knickpoint sigma threshold (default: 1.5)")
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()

    out_dir = args.output_dir or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)

    out_png = args.png or os.path.join(out_dir, "stream_profile.png")
    out_csv = args.csv or os.path.join(out_dir, "stream_profile_data.csv")

    plot_stream_profile(
        stream_vec_path=args.stream,
        dem_path=args.dem,
        output_png=out_png,
        output_csv=out_csv,
        smoothing_window=args.smoothing_window,
        knickpoint_sigma=args.knickpoint_sigma,
    )