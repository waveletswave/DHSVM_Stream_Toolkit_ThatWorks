# -*- coding: utf-8 -*-
# =====================================================================
# channelclass_standalone.py  -  geopandas rewrite of channelclass.py
#
# Sub-step B (part 2) of the #6 vector-IO rebuild. Same classification logic
# as the QGIS channelclass.py (slope->tan, area bins, CLASS_TABLE, PATCH B
# guards, slope-mode auto-detect) copied verbatim; only the QgsVectorLayer
# field IO is swapped for geopandas. Writes stream.class.dat.
#
# prep calls channelclassfun(vl_lines.source(), DIR_STREAMS) with all defaults:
#   acc_units="m2", mannings_n=0.045, pixel_area_m2=None, write_header=True.
#
# Run:  python3 channelclass_standalone.py
# =====================================================================

import os
import math
import statistics
import geopandas as gpd
from pathlib import Path

OUT_DIR    = Path("/work/ys451/dhsvm_ca/standalone_dev/outputs")
STREAM_IN  = OUT_DIR / "streamfile_attr.shp"          # from vector_attrs.py
STREAMS_OUT_DIR = OUT_DIR / "DHSVM_input_streams"     # where stream.class.dat lands


def channelclassfun(streamnet_path,
                    output_dir,
                    pixel_area_m2=None,
                    acc_units="m2",
                    mannings_n=0.045,
                    write_header=True,
                    slope_field_candidates=("slope", "slope_deg", "slp_mean", "slp", "slope_mean"),
                    area_field_candidates=("meanmsq", "mean_m2", "acc_mean"),
                    debug_csv=False,
                    write_back=True):
    """Classify by contributing area (m2) and slope (converted to tan). Verbatim
    logic from the QGIS version; geopandas IO.

    write_back=True also writes the per-segment chanclass/Class/hydwidth/hyddepth
    columns back onto the attribute shapefile (mirrors prep's in-place edit), so
    the downstream stream.network.dat step can read per-segment class. This is the
    option-(ii) equivalent of prep editing vl_lines in place.
    """

    gdf = gpd.read_file(str(streamnet_path))

    def _field_name(cols, candidates):
        names = {c.lower(): c for c in cols}
        for c in candidates:
            if c.lower() in names:
                return names[c.lower()]
        return None

    slope_name = _field_name(gdf.columns, slope_field_candidates)
    area_name = _field_name(gdf.columns, area_field_candidates)
    if not slope_name:
        raise RuntimeError(f"Missing slope field; tried {slope_field_candidates}")
    if not area_name:
        raise RuntimeError(f"Missing area field; tried {area_field_candidates}")

    def _to_float(x, default=None):
        if x is None:
            return default
        try:
            return float(x)
        except Exception:
            return default

    # --- infer slope units from first <=200 records ---
    sample = []
    for i, v in enumerate(gdf[slope_name].tolist()):
        if i > 199:
            break
        fv = _to_float(v)
        if fv is not None:
            sample.append(fv)
    if not sample:
        raise RuntimeError("Slope field has no numeric values.")

    smax = max(sample)
    smed = statistics.median(sample)

    if smax > 6.0:
        slope_mode = "deg"
    elif smax > 1.5:
        slope_mode = "pct"
    else:
        slope_mode = "tan"

    def slope_to_tan(val):
        v = _to_float(val, default=None)
        if v is None:
            return None
        if slope_mode == "deg":
            return math.tan(math.radians(v))
        elif slope_mode == "pct":
            return v / 100.0
        else:
            return v

    def area_to_m2(val):
        a = _to_float(val, default=None)
        if a is None:
            return None
        if acc_units.lower() == "m2":
            return a
        if pixel_area_m2:
            return a * float(pixel_area_m2)
        return a

    # --- bin edges + class table (verbatim) ---
    area_edges = [1e6, 1e7, 2e7, 3e7, 4e7]
    CLASS_TABLE = {
        "gentle": [
            (1, 0.5, 0.03), (2, 1.0, 0.03), (3, 2.0, 0.03),
            (4, 3.0, 0.03), (5, 4.0, 0.03), (6, 4.5, 0.03)
        ],
        "moderate": [
            (7, 0.5, 0.05), (8, 1.0, 0.05), (9, 2.0, 0.05),
            (10, 3.0, 0.05), (11, 4.0, 0.05), (12, 4.5, 0.05)
        ],
        "steep": [
            (13, 0.5, 0.10), (14, 1.0, 0.10), (15, 2.0, 0.10),
            (16, 3.0, 0.10), (17, 4.0, 0.10), (18, 4.5, 0.10)
        ],
    }

    def slope_band(tan_s):
        if tan_s is None:
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

    class_defs, class_counts = {}, {}
    assignments = []
    seg_chanclass = []   # per-segment, for write-back (mirrors prep ft["chanclass"])
    seg_width = []
    seg_depth = []

    # --- per-segment classification (verbatim guards) ---
    for idx in range(len(gdf)):
        tan_s = slope_to_tan(gdf.iloc[idx][slope_name])
        a_m2 = area_to_m2(gdf.iloc[idx][area_name])

        if (tan_s is None) or (not math.isfinite(tan_s)) or (tan_s < 0):
            tan_s = 0.01
        if (a_m2 is None) or (not math.isfinite(a_m2)) or (a_m2 <= 0):
            a_m2 = max(1.0, float(pixel_area_m2) if pixel_area_m2 else 1.0)

        band = slope_band(tan_s)
        abin = area_bin(a_m2)
        cid, w, d = CLASS_TABLE[band][abin]

        seg_chanclass.append(int(cid))
        seg_width.append(float(w))
        seg_depth.append(float(d))

        class_defs.setdefault(cid, (w, d))
        class_counts[cid] = class_counts.get(cid, 0) + 1

        if debug_csv:
            assignments.append((idx, float(tan_s), float(a_m2), int(cid)))

    # --- write stream.class.dat (unique class rows, verbatim format) ---
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(str(output_dir), "stream.class.dat")
    with open(out_path, "w") as f:
        if write_header:
            f.write("#ID W  D   n    inf\n")
        for cid in sorted(class_defs.keys()):
            w, d = class_defs[cid]
            f.write(f"{cid:<3d} {w:4.1f} {d:5.3f} {mannings_n:0.4f} 0.0\n")

    print(f"stream.class.dat written: {out_path}")
    print("Class counts:", {k: class_counts[k] for k in sorted(class_counts)})
    print(f"[info] Slope mode auto-detected as: {slope_mode} (converted to tan for binning)")

    # --- write-back per-segment class columns onto the attr shapefile ---
    if write_back:
        gdf["chanclass"] = seg_chanclass
        gdf["Class"] = seg_chanclass          # mirror, as in the QGIS version
        gdf["hydwidth"] = seg_width
        gdf["hyddepth"] = seg_depth
        gdf.to_file(str(streamnet_path))
        print(f"[info] wrote per-segment chanclass/Class/hydwidth/hyddepth back to "
              f"{os.path.basename(str(streamnet_path))}")

    return out_path


if __name__ == "__main__":
    channelclassfun(STREAM_IN, STREAMS_OUT_DIR)
