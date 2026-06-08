# SUMMARY:      channelclass.py (QGIS, PNNL-compatible, area-units, robust)
# USAGE:        Classify stream channels & write unique stream.class.dat
# ORG:          Pacific Northwest National Laboratory & Duke University
# AUTHOR:       Adapted by Y. Song
# LAST UPDATE:  2025-10-30

import os, math, statistics
from qgis.core import QgsVectorLayer, QgsField
from qgis.PyQt.QtCore import QVariant


def channelclassfun(streamnet_path,
                    output_dir,
                    pixel_area_m2=None,     # kept for compatibility when acc_units != 'm2'
                    acc_units="m2",
                    mannings_n=0.045,
                    write_header=True,
                    # Candidate field names (case-insensitive)
                    slope_field_candidates=("slope", "slope_deg", "slp_mean", "slp", "slope_mean"),
                    area_field_candidates=("meanmsq", "mean_m2", "acc_mean"),
                    # Optional: export a debug CSV (segment → assigned class)
                    debug_csv=False):
    """
    Classify by contributing area (m²) & slope (internally converted to tan(slope)).

    Requires segment-level attributes already populated by upstream steps:
      - slope raster zonal mean in degrees / percent / tan (auto-detected)
      - mean contributing area in m²
        (or in "number of cells" if acc_units != 'm2' and pixel_area_m2 is set)

    Writes:
      - Fields on stream layer: chanclass (int), Class (int), hydwidth (m), hyddepth (m)
      - File: stream.class.dat with unique classes sorted by ID.
    """

    # ---------------------------
    # 0) Load layer & locate fields
    # ---------------------------
    vl = QgsVectorLayer(streamnet_path, "stream_network", "ogr")
    if not vl.isValid():
        raise RuntimeError(f"Cannot load stream network: {streamnet_path}")

    def _field_name(layer, candidates):
        names = {f.name().lower(): f.name() for f in layer.fields()}
        for c in candidates:
            if c.lower() in names:
                return names[c.lower()]
        return None

    slope_name = _field_name(vl, slope_field_candidates)
    area_name  = _field_name(vl, area_field_candidates)
    if not slope_name:
        raise RuntimeError(f"Missing slope field; tried {slope_field_candidates}")
    if not area_name:
        raise RuntimeError(f"Missing area field; tried {area_field_candidates}")

    # ---------------------------
    # 1) Helpers: unit normalization
    # ---------------------------
    def _to_float(x, default=None):
        if x is None:
            return default
        if isinstance(x, QVariant) and x.isNull():
            return default
        try:
            return float(x)
        except Exception:
            return default

    # Sample a subset of records to infer slope units
    sample = []
    for i, ft in enumerate(vl.getFeatures()):
        if i > 199:
            break
        v = _to_float(ft[slope_name])
        if v is not None:
            sample.append(v)
    if not sample:
        raise RuntimeError("Slope field has no numeric values.")

    smax = max(sample)
    smed = statistics.median(sample)

    # Heuristic classification of slope-unit type:
    #  - If max > ~6  → treat as degrees (values > 343° are unrealistic, but we use
    #    a generous threshold to avoid edge cases).
    #  - If 1.5 < max ≤ 100 → commonly percent (0–100).
    #  - Otherwise treat as tan(θ).
    if smax > 6.0:
        slope_mode = "deg"
    elif smax > 1.5:
        slope_mode = "pct"
    else:
        slope_mode = "tan"

    def slope_to_tan(val):
        """Convert slope to tan(slope) according to the inferred unit."""
        v = _to_float(val, default=None)
        if v is None:
            return None
        if slope_mode == "deg":
            return math.tan(math.radians(v))
        elif slope_mode == "pct":
            # 10 percent → 0.10 ≈ tan(θ) (small-angle approximation,
            # sufficient for most hydrologic slopes)
            return v / 100.0
        else:
            return v  # already tan

    # Convert contributing area to m²
    def area_to_m2(val):
        a = _to_float(val, default=None)
        if a is None:
            return None
        if acc_units.lower() == "m2":
            return a
        # Also support "number of cells" as the unit when pixel_area_m2 is given
        if pixel_area_m2:
            return a * float(pixel_area_m2)
        return a  # If no conversion is possible, keep the raw value

    # ---------------------------
    # 2) Configurable bin edges (default matches your existing logic)
    #
    #    tan(slope) is split into three bands:
    #       gentle:   tan ≤ 0.002
    #       moderate: 0.002 < tan ≤ 0.1
    #       steep:    tan > 0.1
    #
    #    Within each band, contributing area (m²) is split into six classes:
    #       ≤1e6, ≤1e7, ≤2e7, ≤3e7, ≤4e7, >4e7
    #
    #    gentle band → class IDs  1.. 6
    #    moderate    → class IDs  7..12
    #    steep       → class IDs 13..18
    #
    #    Each class is associated with (width, depth) as defined below.
    # ---------------------------
    area_edges = [1e6, 1e7, 2e7, 3e7, 4e7]
    CLASS_TABLE = {
        "gentle": [
            (1, 0.5, 0.03), (2, 1.0, 0.03), (3, 2.0, 0.03),
            (4, 3.0, 0.03), (5, 4.0, 0.03), (6, 4.5, 0.03)
        ],
        "moderate": [
            (7, 0.5, 0.05), (8, 1.0, 0.05), (9, 2.0, 0.05),
            (10,3.0, 0.05), (11,4.0, 0.05), (12,4.5, 0.05)
        ],
        "steep": [
            (13,0.5, 0.10), (14,1.0, 0.10), (15,2.0, 0.10),
            (16,3.0, 0.10), (17,4.0, 0.10), (18,4.5, 0.10)
        ],
    }

    def slope_band(tan_s):
        if tan_s is None:
            # Conservative choice: treat as "steep" if slope is missing
            return "steep"
        if tan_s <= 0.002:
            return "gentle"
        if tan_s <= 0.1:
            return "moderate"
        return "steep"

    def area_bin(a_m2):
        """Return a 0..5 index for the six area classes."""
        if a_m2 is None:
            return 5
        for i, edge in enumerate(area_edges):
            if a_m2 <= edge:
                return i
        return 5

    # ---------------------------
    # 3) Prepare output fields
    # ---------------------------
    vl.startEditing()
    for name, qtype in (("chanclass", QVariant.Int),
                        ("Class",     QVariant.Int),
                        ("hydwidth",  QVariant.Double),
                        ("hyddepth",  QVariant.Double)):
        if vl.fields().indexOf(name) < 0:
            vl.addAttribute(QgsField(name, qtype))
    vl.updateFields()

    class_defs, class_counts = {}, {}
    assignments = []  # Optional: for debug CSV

    # ---------------------------
    # 4) Per-segment classification and attribute update
    # ---------------------------
    for ft in vl.getFeatures():
        tan_s = slope_to_tan(ft[slope_name])
        a_m2  = area_to_m2(ft[area_name])

        # --- PATCH B: safety guards for classification inputs ---
        if (tan_s is None) or (not math.isfinite(tan_s)) or (tan_s < 0):
            # Use a modest default slope if value is missing or invalid
            tan_s = 0.01

        if (a_m2 is None) or (not math.isfinite(a_m2)) or (a_m2 <= 0):
            # Use pixel area as minimum plausible contributing area if available;
            # otherwise fall back to 1.0 m² to avoid zero/negative values.
            a_m2 = max(1.0, float(pixel_area_m2) if pixel_area_m2 else 1.0)

        band = slope_band(tan_s)
        abin = area_bin(a_m2)

        cid, w, d = CLASS_TABLE[band][abin]

        ft["chanclass"] = int(cid)
        ft["Class"]     = int(cid)   # Mirror into a commonly used field name
        ft["hydwidth"]  = float(w)
        ft["hyddepth"]  = float(d)
        vl.updateFeature(ft)

        class_defs.setdefault(cid, (w, d))
        class_counts[cid] = class_counts.get(cid, 0) + 1

        if debug_csv:
            assignments.append((
                int(ft.id()),
                float(tan_s if tan_s is not None else -1),
                float(a_m2 if a_m2 is not None else -1),
                int(cid)
            ))

    vl.commitChanges()

    # ---------------------------
    # 5) Write stream.class.dat (unique class rows)
    # ---------------------------
    out_path = os.path.join(output_dir, "stream.class.dat")
    with open(out_path, "w") as f:
        if write_header:
            f.write("#ID W  D   n    inf\n")
        for cid in sorted(class_defs.keys()):
            w, d = class_defs[cid]
            f.write(f"{cid:<3d} {w:4.1f} {d:5.3f} {mannings_n:0.4f} 0.0\n")

    if debug_csv:
        dbg = os.path.join(output_dir, "channelclass_assignments.csv")
        with open(dbg, "w") as f:
            f.write("fid,tan_slope,area_m2,class\n")
            for row in assignments:
                f.write("{},{:.6f},{:.3f},{}\n".format(*row))
        print(f"[debug] Wrote {dbg}")

    print(f"stream.class.dat written: {out_path}")
    print("Class counts:", {k: class_counts[k] for k in sorted(class_counts)})
    if acc_units.lower() != "m2":
        print("[WARN] acc_units != 'm2'. For area-units pipeline, prefer acc_units='m2'.")
    print(f"[info] Slope mode auto-detected as: {slope_mode} (converted to tan for binning)")