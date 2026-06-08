# -*- coding: utf-8 -*-
# =====================================================================
# make_2class_veg_map.py
#
# EXPERIMENT D - Step 1: Generate the 2-class vegetation map
#
# PURPOSE:
#   Collapse the 3-class burn-severity veg map (veg_bs.bin from Step 1
#   of Experiment B) into a 2-class map for Experiment D:
#
#     veg_bs.bin          ->   veg_bs_2c.bin
#     Class 1 (unburned)       Class 1 (unburned/moderate merged)
#     Class 2 (moderate)       Class 1
#     Class 3 (severe)         Class 2 (severe burn isolated)
#
#   Rationale (from Brad's discussion):
#     - 'Severe burn' (~21% of basin, NW headwater cluster) is the
#       physical analog of CA-TO sub-basin (~100% severe burn).
#     - Merging Class 1 + 2 reduces a degree of freedom and lets us
#       directly compare 1-class (uniform) -> 2-class (Exp D) ->
#       3-class (Exp B) frameworks.
#
# INPUTS:
#   /Users/.../DEM_CA_0406/DHSVM_input_binaries/veg_bs.bin  (int8, 74x82)
#
# OUTPUTS:
#   /Users/.../DEM_CA_0406/DHSVM_input_binaries/veg_bs_2c.bin  (int8, 74x82)
#   /Users/.../DEM_CA_0406/DHSVM_input_binaries/veg_bs_2c.tif  (for visual check)
#   make_2class_veg_map_summary.txt
#
# VERIFICATION:
#   - The new Class 2 fraction must match the old Class 3 fraction (~21%)
#   - The new Class 1 fraction must match the old (Class 1 + Class 2)
#     fraction (~79%)
#   - Total active cell count must be unchanged (4334 cells)
#
# AUTHOR: Yiyun Song
# DATE:   2026-04-30
# =====================================================================

import sys
from pathlib import Path
import numpy as np
import rasterio


# ==============================================================
# Configuration
# ==============================================================

INPUT_BIN  = Path("/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CA_0406/DHSVM_input_binaries/veg_bs.bin")
INPUT_TIF  = Path("/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CA_0406/DHSVM_input_binaries/veg_bs.tif")

OUTPUT_BIN = Path("/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CA_0406/DHSVM_input_binaries/veg_bs_2c.bin")
OUTPUT_TIF = Path("/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CA_0406/DHSVM_input_binaries/veg_bs_2c.tif")

# Where to put the summary
SUMMARY_TXT = Path(__file__).resolve().parent / "make_2class_veg_map_summary.txt" \
              if "__file__" in dir() else Path.cwd() / "make_2class_veg_map_summary.txt"

# Expected dimensions (from veg_bs.bin metadata)
EXPECTED_NROWS = 74
EXPECTED_NCOLS = 82

# Class remapping: old class -> new class
CLASS_REMAP = {
    1: 1,   # unburned-low      -> unburned/moderate merged (new class 1)
    2: 1,   # moderate          -> unburned/moderate merged (new class 1)
    3: 2,   # severe burn       -> severe burn isolated     (new class 2)
}

# Cross-check expectations from Step 1 results
EXPECTED_OLD = {1: 1115, 2: 2303, 3: 916}   # cell counts
EXPECTED_NEW = {1: 3418, 2: 916}            # cell counts (1115+2303, 916)


# ==============================================================
# Main
# ==============================================================

def main():
    print("=" * 72)
    print("EXPERIMENT D - Step 1: Generate 2-class vegetation map")
    print("=" * 72)
    print()

    # Pre-flight
    if not INPUT_BIN.exists():
        sys.exit(f"[ERROR] Input not found: {INPUT_BIN}")
    if not INPUT_TIF.exists():
        print(f"[WARN] Input TIFF not found at {INPUT_TIF}")
        print(f"       Will skip TIFF output but still produce .bin")
    print(f"  Input  .bin: {INPUT_BIN}")
    print(f"  Output .bin: {OUTPUT_BIN}")
    print()

    # Read the 3-class binary (int8, row-major)
    print("[STEP 1] Reading 3-class binary...")
    arr_3c = np.fromfile(INPUT_BIN, dtype=np.int8)
    n_cells = arr_3c.size
    expected = EXPECTED_NROWS * EXPECTED_NCOLS
    if n_cells != expected:
        sys.exit(f"[ERROR] Expected {expected} cells "
                 f"({EXPECTED_NROWS}x{EXPECTED_NCOLS}), got {n_cells}")
    arr_3c = arr_3c.reshape(EXPECTED_NROWS, EXPECTED_NCOLS)

    # Class distribution in the input
    print()
    print("  3-class input distribution:")
    total_active = 0
    for c in [1, 2, 3]:
        n = int((arr_3c == c).sum())
        total_active += n
        marker = ""
        if c in EXPECTED_OLD:
            diff = n - EXPECTED_OLD[c]
            marker = f"  (Step 1 reported {EXPECTED_OLD[c]}, diff={diff:+d})"
        print(f"    Class {c}: {n:5d} cells{marker}")
    n_zero = int((arr_3c == 0).sum())
    print(f"    Class 0/nodata: {n_zero} cells")
    print(f"    Total active: {total_active}")
    print()

    # Apply remap
    print("[STEP 2] Applying class remap...")
    print("  Old class 1 (unburned-low)  -> new class 1")
    print("  Old class 2 (moderate)      -> new class 1")
    print("  Old class 3 (severe)        -> new class 2")
    print("  (Outside-watershed 0 stays 0)")
    print()

    arr_2c = np.zeros_like(arr_3c)
    for old_c, new_c in CLASS_REMAP.items():
        mask = (arr_3c == old_c)
        arr_2c[mask] = new_c

    # Verify: 0 stays 0
    assert ((arr_3c == 0) == (arr_2c == 0)).all(), "nodata mask changed!"

    # Distribution check
    print("[STEP 3] Verifying 2-class output...")
    n_new_1 = int((arr_2c == 1).sum())
    n_new_2 = int((arr_2c == 2).sum())

    print(f"  2-class output distribution:")
    print(f"    Class 1 (unburned/moderate merged): {n_new_1:5d} cells "
          f"({100*n_new_1/total_active:5.2f}%)")
    print(f"    Class 2 (severe burn isolated):     {n_new_2:5d} cells "
          f"({100*n_new_2/total_active:5.2f}%)")
    print()

    # Cross-check
    cross_check_pass = True
    if n_new_1 != EXPECTED_NEW[1]:
        print(f"  [WARN] Class 1 count: expected {EXPECTED_NEW[1]}, got {n_new_1}")
        cross_check_pass = False
    if n_new_2 != EXPECTED_NEW[2]:
        print(f"  [WARN] Class 2 count: expected {EXPECTED_NEW[2]}, got {n_new_2}")
        cross_check_pass = False
    # Conservation: total active must be unchanged
    if (n_new_1 + n_new_2) != total_active:
        print(f"  [ERROR] Active cell count changed: "
              f"{total_active} -> {n_new_1 + n_new_2}")
        sys.exit(1)
    if cross_check_pass:
        print("  ✓ All cross-checks passed.")
    print()

    # Write 2-class binary
    print("[STEP 4] Writing 2-class binary...")
    arr_2c.astype(np.int8).tofile(OUTPUT_BIN)
    print(f"  Wrote: {OUTPUT_BIN}  ({OUTPUT_BIN.stat().st_size:,} bytes)")
    print()

    # Write 2-class GeoTIFF (for visual verification, copy georef from input)
    if INPUT_TIF.exists():
        print("[STEP 5] Writing 2-class GeoTIFF (for visual verification)...")
        with rasterio.open(INPUT_TIF) as src:
            profile = src.profile.copy()
            profile.update(dtype='int8', nodata=0, count=1)
            with rasterio.open(OUTPUT_TIF, 'w', **profile) as dst:
                dst.write(arr_2c.astype(np.int8), 1)
        print(f"  Wrote: {OUTPUT_TIF}")
        print()

    # Write summary
    print("[STEP 6] Writing summary...")
    lines = []
    lines.append("=" * 72)
    lines.append("Experiment D - 2-class veg map summary")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Source: 3-class veg_bs.bin (Step 1 of Experiment B)")
    lines.append(f"  Input file:  {INPUT_BIN}")
    lines.append(f"  Output file: {OUTPUT_BIN}")
    lines.append(f"  Dimensions:  {EXPECTED_NROWS} x {EXPECTED_NCOLS} = {expected} cells")
    lines.append(f"  Active cells (in watershed): {total_active}")
    lines.append("")
    lines.append("Class remap:")
    lines.append("  3-class -> 2-class")
    lines.append("    Class 1 (unburned-low)  -> Class 1 (unburned/moderate merged)")
    lines.append("    Class 2 (moderate)      -> Class 1 (unburned/moderate merged)")
    lines.append("    Class 3 (severe burn)   -> Class 2 (severe burn isolated)")
    lines.append("")
    lines.append("3-class input distribution:")
    for c in [1, 2, 3]:
        n = int((arr_3c == c).sum())
        lines.append(f"  Class {c}: {n:5d} cells "
                     f"({100*n/total_active:5.2f}%)")
    lines.append("")
    lines.append("2-class output distribution:")
    lines.append(f"  Class 1 (unburned/moderate merged): {n_new_1:5d} cells "
                 f"({100*n_new_1/total_active:5.2f}%)")
    lines.append(f"  Class 2 (severe burn isolated):     {n_new_2:5d} cells "
                 f"({100*n_new_2/total_active:5.2f}%)")
    lines.append("")
    lines.append("Conservation: total active cells unchanged "
                 f"({total_active}).")
    lines.append("")
    lines.append("Spatial structure (from inspection of veg_bs.tif):")
    lines.append("  - Class 2 forms a contiguous cluster in the NW headwater corner")
    lines.append("    of the basin, geographically corresponding to the CA-TO sub-basin.")
    lines.append("  - This makes Exp D's burned class a natural physical analog to")
    lines.append("    the CA-TO watershed and supports later mechanistic comparison.")
    SUMMARY_TXT.write_text("\n".join(lines))
    print(f"  Saved: {SUMMARY_TXT}")
    print()
    print("=" * 72)
    print("DONE. veg_bs_2c.bin is ready for use in Exp D .dhs templates.")
    print("=" * 72)


if __name__ == "__main__":
    main()
else:
    # This branch runs when the script is exec()'d in a QGIS Python console
    # (or any other environment where __name__ != "__main__"). We still want
    # the script to do its job in that case.
    main()