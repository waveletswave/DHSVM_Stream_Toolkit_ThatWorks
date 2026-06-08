# SUMMARY:      channelclass.py (standalone / DCC-ready, geopandas-based)
# USAGE:        Classify stream channels & write unique stream.class.dat
# ORG:          Pacific Northwest National Laboratory & Duke University
# AUTHOR:       Adapted by Y. Song
# LAST UPDATE:  2026-05-02
#
# Standalone (non-QGIS) sibling of qgis/channelclass.py.
# Uses geopandas + shapely instead of QgsVectorLayer / QVariant.
# Produces stream.class.dat that is BYTE-FOR-BYTE identical to the QGIS
# version when given the same input shapefile.
#
# Run from terminal:
#     python channelclass.py path/to/streamnet.shp path/to/output_dir
#
# Or import as a library:
#     from channelclass import channelclassfun
#     channelclassfun(streamnet_path, output_dir, mannings_n=0.045)
#
# Dependencies: geopandas, shapely, pandas, numpy  (all conda-forge)

import os
import math
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd


def channelclassfun(streamnet_path,
                    output_dir,
                    pixel_area_m2=None,            # kept for compatibility when acc_units != 'm2'
                    acc_units="m2",
                    mannings_n=0.045,
                    write_header=True,
                    # Candidate field names (case-insensitive); same defaults as QGIS version
                    slope_field_candidates=("slope", "slope_deg", "slp_mean", "slp", "slope_mean"),
                    area_field_candidates=("meanmsq", "mean_m2", "acc_mean"),
                    # Optional: export a debug CSV (segment -> assigned class)
                    debug_csv=False):
    """
    Classify by contributing area (m²) and slope (auto-detected -> tan(theta)).

    Reads a stream-network vector file (shapefile / GeoPackage / GeoJSON),
    classifies every segment into one of 18 (slope_band x area_bin) classes,
    writes new attribute fields back to the file, and emits stream.class.dat.

    Output behaviour (chanclass / Class / hydwidth / hyddepth fields, the
    table of class IDs, the stream.class.dat byte content) is preserved
    relative to qgis/channelclass.py for the same input.
    """

    # ---------------------------
    # 0) Load layer & locate fields
    # ---------------------------
    streamnet_path = str(streamnet_path)
    output_dir = str(output_dir)

    if not os.path.exists(streamnet_path):
        raise FileNotFoundError(f"Cannot find stream network: {streamnet_path}")

    gdf = gpd.read_file(streamnet_path)
    if len(gdf) == 0:
        raise RuntimeError(f"Stream network is empty: {streamnet_path}")

    def _field_name(columns, candidates):
        names = {c.lower(): c for c in columns}
        for cand in candidates:
            if cand.lower() in names:
                return names[cand.lower()]
        return None

    slope_name = _field_name(gdf.columns, slope_field_candidates)
    area_name  = _field_name(gdf.columns, area_field_candidates)
    if not slope_name:
        raise RuntimeError(
            f"Missing slope field; tried {slope_field_candidates}. "
            f"Available columns: {list(gdf.columns)}"
        )
    if not area_name:
        raise RuntimeError(
            f"Missing area field; tried {area_field_candidates}. "
            f"Available columns: {list(gdf.columns)}"
        )

    # ---------------------------
    # 1) Helpers: unit normalisation
    # ---------------------------
    def _to_float(x):
        """Coerce to float, returning None for NaN / non-numeric."""
        if pd.isna(x):
            return None
        try:
            return float(x)
        except (TypeError, ValueError):
            return None

    # Sample first 200 valid slope values to infer units (matches QGIS heuristic)
    sample = []
    for v in gdf[slope_name].head(200):
        fv = _to_float(v)
        if fv is not None:
            sample.append(fv)
    if not sample:
        raise RuntimeError(f"Slope field '{slope_name}' has no numeric values.")

    smax = max(sample)

    # Heuristic classification of slope-unit type (identical to QGIS version):
    #  - max > 6.0   -> degrees
    #  - max > 1.5   -> percent (0-100)
    #  - otherwise   -> tan(theta)
    if smax > 6.0:
        slope_mode = "deg"
    elif smax > 1.5:
        slope_mode = "pct"
    else:
        slope_mode = "tan"

    def slope_to_tan(val):
        v = _to_float(val)
        if v is None:
            return None
        if slope_mode == "deg":
            return math.tan(math.radians(v))
        elif slope_mode == "pct":
            # 10 percent -> 0.10 ~= tan(theta) (small-angle approximation,
            # sufficient for most hydrologic slopes)
            return v / 100.0
        return v  # already tan

    def area_to_m2(val):
        a = _to_float(val)
        if a is None:
            return None
        if acc_units.lower() == "m2":
            return a
        if pixel_area_m2:
            return a * float(pixel_area_m2)
        return a  # fallback: keep raw value

    # ---------------------------
    # 2) Classification table (identical to QGIS version)
    #
    #    tan(slope) bands:
    #      gentle:   tan <= 0.002
    #      moderate: 0.002 < tan <= 0.1
    #      steep:    tan > 0.1
    #
    #    Within each band, contributing area (m²) is split into six classes:
    #      <=1e6, <=1e7, <=2e7, <=3e7, <=4e7, >4e7
    #
    #    gentle band   -> class IDs  1..6
    #    moderate band -> class IDs  7..12
    #    steep band    -> class IDs 13..18
    # ---------------------------
    area_edges = [1e6, 1e7, 2e7, 3e7, 4e7]
    CLASS_TABLE = {
        "gentle": [
            (1, 0.5, 0.03), (2, 1.0, 0.03), (3, 2.0, 0.03),
            (4, 3.0, 0.03), (5, 4.0, 0.03), (6, 4.5, 0.03),
        ],
        "moderate": [
            (7,  0.5, 0.05), (8,  1.0, 0.05), (9,  2.0, 0.05),
            (10, 3.0, 0.05), (11, 4.0, 0.05), (12, 4.5, 0.05),
        ],
        "steep": [
            (13, 0.5, 0.10), (14, 1.0, 0.10), (15, 2.0, 0.10),
            (16, 3.0, 0.10), (17, 4.0, 0.10), (18, 4.5, 0.10),
        ],
    }

    def slope_band(tan_s):
        if tan_s is None:
            # Conservative default: treat missing slope as steep
            return "steep"
        if tan_s <= 0.002:
            return "gentle"
        if tan_s <= 0.1:
            return "moderate"
        return "steep"

    def area_bin(a_m2):
        if a_m2 is None:
            return 5
        for i, edge in enumerate(area_edges):
            if a_m2 <= edge:
                return i
        return 5

    # ---------------------------
    # 3) Per-segment classification
    #     (loop kept row-by-row for 1:1 parity with QGIS version; the data
    #     is small enough that vectorising buys nothing meaningful)
    # ---------------------------
    n = len(gdf)
    chanclass = np.zeros(n, dtype=np.int32)
    hydwidth  = np.zeros(n, dtype=np.float64)
    hyddepth  = np.zeros(n, dtype=np.float64)

    class_defs   = {}
    class_counts = {}
    assignments  = []  # for optional debug CSV

    slope_vals = gdf[slope_name].tolist()
    area_vals  = gdf[area_name].tolist()

    for i in range(n):
        tan_s = slope_to_tan(slope_vals[i])
        a_m2  = area_to_m2(area_vals[i])

        # PATCH B: safety guards (identical to QGIS version)
        if (tan_s is None) or (not math.isfinite(tan_s)) or (tan_s < 0):
            tan_s = 0.01
        if (a_m2 is None) or (not math.isfinite(a_m2)) or (a_m2 <= 0):
            a_m2 = max(1.0, float(pixel_area_m2) if pixel_area_m2 else 1.0)

        band = slope_band(tan_s)
        abin = area_bin(a_m2)
        cid, w, d = CLASS_TABLE[band][abin]

        chanclass[i] = int(cid)
        hydwidth[i]  = float(w)
        hyddepth[i]  = float(d)

        class_defs.setdefault(cid, (w, d))
        class_counts[cid] = class_counts.get(cid, 0) + 1

        if debug_csv:
            assignments.append((i, float(tan_s), float(a_m2), int(cid)))

    # Write fields back to the GeoDataFrame
    gdf["chanclass"] = chanclass
    gdf["Class"]     = chanclass        # mirror, matches QGIS field convention
    gdf["hydwidth"]  = hydwidth
    gdf["hyddepth"]  = hyddepth

    # ---------------------------
    # 4) Save back to the same vector file (preserve original driver)
    # ---------------------------
    ext = os.path.splitext(streamnet_path)[1].lower()
    driver_map = {
        ".shp":     "ESRI Shapefile",
        ".gpkg":    "GPKG",
        ".geojson": "GeoJSON",
        ".json":    "GeoJSON",
    }
    driver = driver_map.get(ext)
    if driver:
        gdf.to_file(streamnet_path, driver=driver)
    else:
        gdf.to_file(streamnet_path)  # let geopandas guess

    # ---------------------------
    # 5) Write stream.class.dat (byte-identical format to QGIS version)
    # ---------------------------
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "stream.class.dat")
    with open(out_path, "w") as f:
        if write_header:
            f.write("#ID W  D   n    inf\n")
        for cid in sorted(class_defs.keys()):
            w, d = class_defs[cid]
            f.write(f"{cid:<3d} {w:4.1f} {d:5.3f} {mannings_n:0.4f} 0.0\n")

    # Optional debug CSV
    if debug_csv:
        dbg = os.path.join(output_dir, "channelclass_assignments.csv")
        with open(dbg, "w") as f:
            f.write("fid,tan_slope,area_m2,class\n")
            for row in assignments:
                f.write("{},{:.6f},{:.3f},{}\n".format(*row))
        print(f"[debug] Wrote {dbg}")

    # ---------------------------
    # 6) Console summary (mirrors QGIS version)
    # ---------------------------
    print(f"stream.class.dat written: {out_path}")
    print("Class counts:", {k: class_counts[k] for k in sorted(class_counts)})
    if acc_units.lower() != "m2":
        print("[WARN] acc_units != 'm2'. For area-units pipeline, prefer acc_units='m2'.")
    print(f"[info] Slope mode auto-detected as: {slope_mode} (converted to tan for binning)")


# ==============================================================
# CLI entry point
# ==============================================================
def _build_parser():
    p = argparse.ArgumentParser(
        description=(
            "Classify DHSVM stream segments by slope band x contributing-area bin. "
            "Writes chanclass / Class / hydwidth / hyddepth back to the input vector "
            "and emits stream.class.dat in the output directory."
        )
    )
    p.add_argument("streamnet",
                   help="Path to stream network (.shp / .gpkg / .geojson)")
    p.add_argument("output_dir",
                   help="Directory to write stream.class.dat (created if missing)")
    p.add_argument("--pixel-area-m2", type=float, default=None,
                   help="Pixel area in m² (only needed when --acc-units != m2)")
    p.add_argument("--acc-units", default="m2", choices=["m2", "cells"],
                   help="Units of the contributing-area field (default: m2)")
    p.add_argument("--mannings-n", type=float, default=0.045,
                   help="Manning's n written to stream.class.dat (default: 0.045)")
    p.add_argument("--no-header", action="store_true",
                   help="Suppress the leading '#ID W D n inf' header line")
    p.add_argument("--debug-csv", action="store_true",
                   help="Also emit channelclass_assignments.csv for QC")
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()
    channelclassfun(
        streamnet_path=args.streamnet,
        output_dir=args.output_dir,
        pixel_area_m2=args.pixel_area_m2,
        acc_units=args.acc_units,
        mannings_n=args.mannings_n,
        write_header=(not args.no_header),
        debug_csv=args.debug_csv,
    )