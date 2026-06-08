# -*- coding: utf-8 -*-
# =====================================================================
# make_4class_veg_map.py
#
# EXPERIMENT R (RtltS4h) - Step 1: Generate the 4-class vegetation map
#
# PURPOSE:
#   Build a 4-class burn-severity vegetation map by classifying the
#   RdNBR raster directly with the standard Caldwell 2020 Table 2
#   thresholds, instead of starting from the existing 3-class map.
#
#   This produces a finer spatial framework than veg_bs.bin (3-class)
#   or veg_bs_2c.bin (2-class), and aligns one-to-one with the burn
#   classes used in the Sentinel-2 RS feasibility test (so Cohen's d
#   diagnostics in the RS work directly correspond to DHSVM classes).
#
# CLASS DEFINITION (RdNBR thresholds, Caldwell et al. 2020 Table 2):
#     Class 1: Unburned-Moderate    RdNBR <  62
#     Class 2: Moderate             62  <= RdNBR < 181
#     Class 3: Moderate-High        181 <= RdNBR < 541
#     Class 4: High                 RdNBR >= 541
#
# RATIONALE:
#   - Unlike make_2class_veg_map.py, this script does NOT collapse
#     existing classes; it re-classifies the RdNBR raster from scratch
#     onto the DHSVM grid. This guarantees every class comes from the
#     same underlying data product.
#   - The RS feasibility test (Sentinel-2 NBR Cohen's d analysis)
#     showed strong, monotonic separation between these four classes
#     in CA_TO (d = 6.59 in 2017 NBR for High vs Unburned-Moderate),
#     supporting their use as DHSVM input classes.
#
# INPUTS:
#   1. The DHSVM grid template (use veg_bs.bin's GeoTIFF for georeference)
#        veg_bs.tif  (74 x 82, EPSG:??, DHSVM grid)
#   2. The RdNBR raster (from Pete's GIS data)
#        RdNBR_20160609_20170726.tif  (10 m, EPSG:32617, much finer)
#
# OUTPUTS:
#   veg_bs_4c.bin  (int8, 74 x 82, row-major; for DHSVM)
#   veg_bs_4c.tif  (GeoTIFF, for visual verification in QGIS)
#   make_4class_veg_map_summary.txt
#
# VERIFICATION:
#   - Total active cell count must equal veg_bs.bin's active cell count
#     (4334 cells; matches CA watershed)
#   - Cells outside watershed stay = 0
#   - Class fractions should approximately match Caldwell 2020 Table 1
#     for CA basin: High ~ 21%, Moderate-High ~ 17%, Moderate ~ 34%,
#     Unburned-Moderate ~ 28%
#   - Class 4 cluster should be in the NW headwater corner (CA-TO area)
#
# AUTHOR: Yiyun Song
# DATE:   2026-05-04
# =====================================================================

import sys
from pathlib import Path
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling


# ==============================================================
# Configuration
# ==============================================================

# DHSVM input directory (where veg_bs.bin lives)
DHSVM_INPUT_DIR = Path(
    "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CA_0406/DHSVM_input_binaries"
)

# Existing veg_bs files (used as the DHSVM grid template + active mask)
TEMPLATE_BIN = DHSVM_INPUT_DIR / "veg_bs.bin"
TEMPLATE_TIF = DHSVM_INPUT_DIR / "veg_bs.tif"

# Source RdNBR raster
RDNBR_TIF = Path(
    "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/"
    "GIS/burn_severity/RdNBR_20160609_20170726.tif"
)

# Outputs
OUTPUT_BIN = DHSVM_INPUT_DIR / "veg_bs_4c.bin"
OUTPUT_TIF = DHSVM_INPUT_DIR / "veg_bs_4c.tif"

SUMMARY_TXT = (
    Path(__file__).resolve().parent / "make_4class_veg_map_summary.txt"
    if "__file__" in dir() else Path.cwd() / "make_4class_veg_map_summary.txt"
)

# Expected DHSVM grid dimensions (from veg_bs.bin metadata)
EXPECTED_NROWS = 74
EXPECTED_NCOLS = 82

# Caldwell 2020 Table 2 RdNBR thresholds for 4-class classification
# Class 1: Unburned-Moderate  RdNBR < 62
# Class 2: Moderate           62 <= RdNBR < 181
# Class 3: Moderate-High      181 <= RdNBR < 541
# Class 4: High               RdNBR >= 541
CALDWELL_THRESHOLDS = [62.0, 181.0, 541.0]
CLASS_NAMES = {
    1: "Unburned-Moderate",
    2: "Moderate",
    3: "Moderate-High",
    4: "High",
}

# Cross-check: published CA watershed fractions from Caldwell 2020 Table 1
# (used only for sanity, not for hard validation)
EXPECTED_CA_FRACTIONS = {
    1: 0.28,   # Unburned-Moderate ~ 28%
    2: 0.34,   # Moderate          ~ 34%
    3: 0.17,   # Moderate-High     ~ 17%
    4: 0.21,   # High              ~ 21%
}


# ==============================================================
# Core routines
# ==============================================================

def reproject_rdnbr_to_template(rdnbr_tif, template_tif):
    """Reproject the RdNBR raster onto the DHSVM grid (template_tif)."""
    with rasterio.open(template_tif) as tpl:
        dst_meta = tpl.meta.copy()
        dst_shape = (tpl.height, tpl.width)
        dst_transform = tpl.transform
        dst_crs = tpl.crs

    with rasterio.open(rdnbr_tif) as src:
        src_data = src.read(1).astype(np.float32)
        src_transform = src.transform
        src_crs = src.crs
        src_nodata = src.nodata

    # Use 'average' resampling: each DHSVM coarse cell averages many
    # 10 m RdNBR pixels. This is the right choice for a continuous
    # severity index that we will threshold afterwards.
    dst = np.full(dst_shape, np.nan, dtype=np.float32)
    reproject(
        source=src_data,
        destination=dst,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        resampling=Resampling.average,
        src_nodata=src_nodata,
        dst_nodata=np.nan,
    )
    # Clean: extreme negatives are NoData markers, not real values
    dst[dst < -1e30] = np.nan
    return dst, dst_meta


def classify_rdnbr_to_4class(rdnbr_arr, active_mask):
    """Apply Caldwell thresholds, restricted to active (in-watershed) cells."""
    out = np.zeros_like(rdnbr_arr, dtype=np.int8)  # 0 = NoData / outside

    # Only classify within active mask AND where RdNBR is finite
    valid = active_mask & ~np.isnan(rdnbr_arr)

    out[valid & (rdnbr_arr <  CALDWELL_THRESHOLDS[0])] = 1
    out[valid & (rdnbr_arr >= CALDWELL_THRESHOLDS[0])
              & (rdnbr_arr <  CALDWELL_THRESHOLDS[1])] = 2
    out[valid & (rdnbr_arr >= CALDWELL_THRESHOLDS[1])
              & (rdnbr_arr <  CALDWELL_THRESHOLDS[2])] = 3
    out[valid & (rdnbr_arr >= CALDWELL_THRESHOLDS[2])] = 4

    # Edge case: active cells where RdNBR is NaN after reprojection.
    # Fall back to the existing veg_bs class for that cell? Or set to
    # Class 1 (no detected burn)? Choose the latter as a conservative
    # default; report how many such cells exist.
    fallback_cells = active_mask & np.isnan(rdnbr_arr)
    n_fallback = int(fallback_cells.sum())
    out[fallback_cells] = 1  # default to Unburned-Moderate

    return out, n_fallback


# ==============================================================
# Main
# ==============================================================

def main():
    print("=" * 72)
    print("EXPERIMENT R - Step 1: Generate 4-class vegetation map")
    print("=" * 72)
    print()

    # Pre-flight checks
    if not TEMPLATE_BIN.exists():
        sys.exit(f"[ERROR] Template binary not found: {TEMPLATE_BIN}")
    if not TEMPLATE_TIF.exists():
        sys.exit(f"[ERROR] Template GeoTIFF not found: {TEMPLATE_TIF}")
    if not RDNBR_TIF.exists():
        sys.exit(f"[ERROR] RdNBR raster not found: {RDNBR_TIF}")

    print(f"  DHSVM template (.bin): {TEMPLATE_BIN}")
    print(f"  DHSVM template (.tif): {TEMPLATE_TIF}")
    print(f"  RdNBR source:          {RDNBR_TIF}")
    print(f"  Output:                {OUTPUT_BIN}")
    print()

    # ---- STEP 1: read template, derive active mask ----
    print("[STEP 1] Reading 3-class template to derive active mask...")
    arr_3c = np.fromfile(TEMPLATE_BIN, dtype=np.int8)
    expected_n = EXPECTED_NROWS * EXPECTED_NCOLS
    if arr_3c.size != expected_n:
        sys.exit(f"[ERROR] Expected {expected_n} cells, got {arr_3c.size}")
    arr_3c = arr_3c.reshape(EXPECTED_NROWS, EXPECTED_NCOLS)
    active_mask = (arr_3c != 0)
    n_active = int(active_mask.sum())
    print(f"  Grid: {EXPECTED_NROWS} x {EXPECTED_NCOLS}")
    print(f"  Active (in-watershed) cells: {n_active}")
    print()

    # ---- STEP 2: reproject RdNBR onto DHSVM grid ----
    print("[STEP 2] Reprojecting RdNBR onto DHSVM grid...")
    rdnbr_dhsvm, dhsvm_meta = reproject_rdnbr_to_template(RDNBR_TIF, TEMPLATE_TIF)
    n_finite = int(np.isfinite(rdnbr_dhsvm).sum())
    print(f"  Resampled shape:    {rdnbr_dhsvm.shape}")
    print(f"  Finite RdNBR cells: {n_finite}")
    if n_finite > 0:
        print(f"  RdNBR range:        "
              f"{np.nanmin(rdnbr_dhsvm):.0f} .. {np.nanmax(rdnbr_dhsvm):.0f}")
    print()

    # ---- STEP 3: classify into 4 classes ----
    print("[STEP 3] Classifying with Caldwell 2020 Table 2 thresholds...")
    print(f"  Class 1 (Unburned-Moderate): RdNBR <  {CALDWELL_THRESHOLDS[0]:.0f}")
    print(f"  Class 2 (Moderate):          {CALDWELL_THRESHOLDS[0]:.0f} "
          f"<= RdNBR < {CALDWELL_THRESHOLDS[1]:.0f}")
    print(f"  Class 3 (Moderate-High):     {CALDWELL_THRESHOLDS[1]:.0f} "
          f"<= RdNBR < {CALDWELL_THRESHOLDS[2]:.0f}")
    print(f"  Class 4 (High):              RdNBR >= {CALDWELL_THRESHOLDS[2]:.0f}")
    print()

    arr_4c, n_fallback = classify_rdnbr_to_4class(rdnbr_dhsvm, active_mask)
    if n_fallback > 0:
        print(f"  [INFO] {n_fallback} active cells had no RdNBR coverage; "
              f"defaulted to Class 1 (Unburned-Moderate).")

    # ---- STEP 4: verify ----
    print()
    print("[STEP 4] Verifying 4-class output...")

    # Class distribution
    class_counts = {}
    for c in [1, 2, 3, 4]:
        n = int((arr_4c == c).sum())
        class_counts[c] = n
    n_zero = int((arr_4c == 0).sum())
    n_active_out = sum(class_counts.values())

    print(f"  Class distribution (out of {n_active_out} active cells):")
    print(f"  {'Class':<6}{'Name':<22}{'Count':<8}{'Fraction':<10}"
          f"{'Caldwell expect.':<18}")
    for c in [1, 2, 3, 4]:
        n = class_counts[c]
        frac = n / n_active_out if n_active_out > 0 else 0
        expected = EXPECTED_CA_FRACTIONS[c]
        delta = frac - expected
        marker = "OK" if abs(delta) < 0.10 else "CHECK"
        print(f"  {c:<6}{CLASS_NAMES[c]:<22}{n:<8}{frac:<10.2%}"
              f"~{expected:.0%}  diff {delta:+.2%}  [{marker}]")

    # Conservation check
    if n_active_out != n_active:
        print(f"  [ERROR] Active cell count changed: "
              f"{n_active} -> {n_active_out}")
        sys.exit(1)
    if n_zero != (expected_n - n_active):
        print(f"  [ERROR] NoData count changed: "
              f"{expected_n - n_active} -> {n_zero}")
        sys.exit(1)
    print(f"  Conservation: total active cells unchanged ({n_active}). OK")
    print()

    # ---- STEP 5: write outputs ----
    print("[STEP 5] Writing 4-class outputs...")

    arr_4c.astype(np.int8).tofile(OUTPUT_BIN)
    print(f"  Wrote: {OUTPUT_BIN}  ({OUTPUT_BIN.stat().st_size:,} bytes)")

    profile = dhsvm_meta.copy()
    profile.update(dtype="int8", nodata=0, count=1)
    with rasterio.open(OUTPUT_TIF, "w", **profile) as dst:
        dst.write(arr_4c.astype(np.int8), 1)
    print(f"  Wrote: {OUTPUT_TIF}")
    print()

    # ---- STEP 6: summary ----
    print("[STEP 6] Writing summary...")
    lines = []
    lines.append("=" * 72)
    lines.append("Experiment R (RtltS4h) - 4-class veg map summary")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Source data:")
    lines.append(f"  DHSVM template (.bin): {TEMPLATE_BIN}")
    lines.append(f"  DHSVM template (.tif): {TEMPLATE_TIF}")
    lines.append(f"  RdNBR raster:          {RDNBR_TIF}")
    lines.append("")
    lines.append("Output:")
    lines.append(f"  veg_bs_4c.bin: {OUTPUT_BIN}")
    lines.append(f"  veg_bs_4c.tif: {OUTPUT_TIF}")
    lines.append("")
    lines.append("Grid:")
    lines.append(f"  Dimensions: {EXPECTED_NROWS} x {EXPECTED_NCOLS} "
                 f"= {expected_n} cells")
    lines.append(f"  Active cells (in CA watershed): {n_active}")
    if n_fallback > 0:
        lines.append(f"  Cells without RdNBR coverage (fallback to Class 1): "
                     f"{n_fallback}")
    lines.append("")
    lines.append("Classification (Caldwell 2020 Table 2 RdNBR thresholds):")
    lines.append(f"  Class 1 Unburned-Moderate: RdNBR <  62")
    lines.append(f"  Class 2 Moderate:           62 <= RdNBR < 181")
    lines.append(f"  Class 3 Moderate-High:     181 <= RdNBR < 541")
    lines.append(f"  Class 4 High:              RdNBR >= 541")
    lines.append("")
    lines.append("Class distribution:")
    for c in [1, 2, 3, 4]:
        n = class_counts[c]
        frac = n / n_active_out if n_active_out > 0 else 0
        expected = EXPECTED_CA_FRACTIONS[c]
        lines.append(
            f"  Class {c} {CLASS_NAMES[c]:<22}: {n:5d} cells "
            f"({frac:.2%}, Caldwell ~{expected:.0%})"
        )
    lines.append("")
    lines.append(f"Conservation: total active cells unchanged ({n_active}).")
    lines.append("")
    lines.append("Notes for downstream use in DHSVM:")
    lines.append("  - This map mirrors the burn-severity classes used in the")
    lines.append("    Sentinel-2 RS feasibility test (CA_TO 2017 NBR Cohen's d")
    lines.append("    = 6.59 between Class 4 High and Class 1 Unburned-Moderate),")
    lines.append("    so the spatial framework in DHSVM corresponds directly to")
    lines.append("    a measured RS gradient.")
    lines.append("  - Per-class LAI multipliers in RtltS4h will be inverted by")
    lines.append("    DHSVM (grid search or optimization on PBIAS / NSE), not")
    lines.append("    derived from RS via empirical NBR-to-LAI relationships.")
    lines.append("    The RS data's role is to define WHICH cells belong to")
    lines.append("    which class, not to assign LAI values.")
    SUMMARY_TXT.write_text("\n".join(lines))
    print(f"  Saved: {SUMMARY_TXT}")
    print()
    print("=" * 72)
    print("DONE. veg_bs_4c.bin is ready for use in RtltS4h .dhs templates.")
    print("=" * 72)


if __name__ == "__main__":
    main()
else:
    main()