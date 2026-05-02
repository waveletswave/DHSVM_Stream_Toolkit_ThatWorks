# -*- coding: utf-8 -*-
# =====================================================================
# rdnbr_to_veg_bs.py
#
# PURPOSE:
#   Generates a 3-class burn-severity vegetation map (`veg_bs.bin`) for
#   DHSVM Experiment B from Sentinel-2 RdNBR data. This replaces the
#   uniform `veg.bin` (all type 1) with a spatially heterogeneous map
#   that assigns each grid cell to one of three burn-severity classes
#   based on the Camp Branch Fire (Caldwell et al. 2020) thresholds.
#
# CLASSIFICATION (RdNBR thresholds from Caldwell et al. 2020 Table 2,
# 4-class scheme, with their lower two classes already merged into
# "Unburned-Moderate" because RdNBR cannot reliably discern them
# one year post-fire):
#
#   Class 1 (Unburned-Low):  RdNBR < 62
#       Corresponds to Caldwell's "Unburned-Moderate". Renamed here
#       to avoid confusion with our Class 2.
#   Class 2 (Moderate):      62 <= RdNBR <= 541
#       Merges Caldwell's "Moderate" (62-181) and "Moderate-High"
#       (182-541). Both have likely intact understory.
#   Class 3 (High):          RdNBR > 541
#       Corresponds to Caldwell's "High". Stand-replacement burn,
#       understory destroyed.
#
# DHSVM IMPLICATIONS:
#   The .dhs vegetation block must define these three classes. The key
#   per-class flag is `Understory Present`:
#       Class 1 -> TRUE  (intact forest, understory present)
#       Class 2 -> TRUE  (overstory damaged but understory survives)
#       Class 3 -> FALSE (stand-replacement, bare-soil evap re-enters)
#   This addresses the SoilEvap=0 mechanism documented in
#   DHSVM_SoilEvap_Trace.md (MassEnergyBalance.c line 391).
#
# RESAMPLING:
#   RdNBR is at 10m (Sentinel-2 native), DEM is at 30m. Aggregation is
#   done via gdal.Warp with resampleAlg='average', which correctly
#   handles arbitrary spatial alignment between the two grids. (The
#   RdNBR origin is offset from a natural 30m grid by 20m, so naive
#   numpy 3x3 block aggregation would produce incorrect results.)
#
# REFERENCE:
#   Caldwell, P.V., et al. 2020. Watershed-scale vegetation, water
#   quantity, and water quality responses to wildfire in the southern
#   Appalachian mountain region. Hydrological Processes 34: 5188-5209.
#
# AUTHOR: Yiyun Song
# DATE:   2026-04-29
# =====================================================================

import os
import sys
import numpy as np
from osgeo import gdal, osr
from pathlib import Path

# ==============================================================
# Configuration: input paths and classification thresholds
# ==============================================================

# Auto-resolve workspace following the existing toolkit convention
SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR.parent

# Input DEM (provides target grid: extent, transform, projection, mask)
DEM_IN  = WS / "Reprojected_DEM" / "elev_clipped.tif"

# Input RdNBR (Sentinel-2 derived, fixed location outside WS)
RDNBR_IN = Path("/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/burn_severity/RdNBR_20160609_20170726.tif")

# Output directory (matches existing toolkit convention)
OUT_DIR = WS / "DHSVM_input_binaries"

# DHSVM standard NoData value
DHSVM_NODATA = -9999.0

# Classification thresholds (from Caldwell 2020, see header)
THR_LOW_MOD  = 62.0    # Class 1 / Class 2 boundary
THR_MOD_HIGH = 541.0   # Class 2 / Class 3 boundary

# Class names for reporting
CLASS_NAMES = {
    1: "Unburned-Low",
    2: "Moderate",
    3: "High",
}

# Caldwell 2020 Table 1 reference values for CA watershed (for cross-check)
CA_REFERENCE = {
    "watershed_mean_RdNBR": 292.0,
    "class_pct_expected": {
        1: 27.0,   # Caldwell's "Unburned-Moderate"
        2: 52.0,   # Caldwell's "Moderate" (35%) + "Moderate-High" (17%)
        3: 21.0,   # Caldwell's "High"
    },
}


# ==============================================================
# Helper: read a raster's CRS as a normalized WKT string
# ==============================================================
def get_crs_info(raster_path):
    """
    Returns a dict with CRS information for a raster:
      - 'wkt':         the projection WKT string
      - 'epsg':        the EPSG code if identifiable, else None
      - 'auth_name':   typically 'EPSG'
      - 'projcs_name': human-readable projected CS name
    """
    ds = gdal.Open(str(raster_path))
    if ds is None:
        raise RuntimeError(f"Cannot open raster: {raster_path}")
    wkt = ds.GetProjection()
    if not wkt:
        raise RuntimeError(f"Raster has no CRS defined: {raster_path}")

    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    srs.AutoIdentifyEPSG()  # try to attach EPSG code
    epsg = srs.GetAuthorityCode(None)
    auth = srs.GetAuthorityName(None)
    projcs = srs.GetAttrValue("PROJCS") or srs.GetAttrValue("GEOGCS") or "unknown"

    ds = None
    return {
        "wkt": wkt,
        "epsg": epsg,
        "auth_name": auth,
        "projcs_name": projcs,
    }


def verify_crs_match(crs_a, crs_b, name_a="DEM", name_b="RdNBR"):
    """
    Raises if CRSs don't match. Uses EPSG when available, falls back to
    SpatialReference.IsSame() comparison for robustness.
    """
    print(f"\n[CRS check]")
    print(f"  {name_a:6s} CRS: {crs_a['projcs_name']}  (EPSG:{crs_a['epsg']})")
    print(f"  {name_b:6s} CRS: {crs_b['projcs_name']}  (EPSG:{crs_b['epsg']})")

    # Best case: both have identifiable EPSG and they match
    if crs_a["epsg"] and crs_b["epsg"] and crs_a["epsg"] == crs_b["epsg"]:
        print(f"  -> Match (EPSG:{crs_a['epsg']})")
        return

    # Fallback: use OGR SpatialReference equality
    srs_a = osr.SpatialReference(); srs_a.ImportFromWkt(crs_a["wkt"])
    srs_b = osr.SpatialReference(); srs_b.ImportFromWkt(crs_b["wkt"])
    if srs_a.IsSame(srs_b):
        print(f"  -> Match (WKT-equivalent, EPSG codes differ or absent)")
        return

    # Mismatch - raise with diagnostic info
    raise RuntimeError(
        f"CRS MISMATCH between {name_a} and {name_b}.\n"
        f"  {name_a}: {crs_a['projcs_name']} (EPSG:{crs_a['epsg']})\n"
        f"  {name_b}: {crs_b['projcs_name']} (EPSG:{crs_b['epsg']})\n"
        f"  Both should be UTM Zone 17N WGS 84 (EPSG:32617) for CA.\n"
        f"  Reproject one of them with gdalwarp before running this script."
    )


# ==============================================================
# Helper: read DEM and derive target grid + master mask
# ==============================================================
def read_dem_grid(dem_path):
    """
    Reads the DEM and returns the canonical "target grid" parameters
    plus the master valid mask. All downstream rasters must align to
    this grid.
    """
    print(f"\n[step] Reading DEM: {dem_path.name}")
    ds = gdal.Open(str(dem_path))
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float32)
    nd = band.GetNoDataValue()
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    nx = ds.RasterXSize
    ny = ds.RasterYSize

    # Master mask: same definition as dem_to_dhsvm_bins.py
    valid_mask = (arr != nd) & ~np.isnan(arr)
    n_valid = int(valid_mask.sum())

    # Compute extent for diagnostic and Warp output bounds
    px = gt[1]
    py = abs(gt[5])
    x_min = gt[0]
    y_max = gt[3]
    x_max = x_min + nx * px
    y_min = y_max - ny * py

    print(f"  Grid:    {ny} rows x {nx} cols, pixel = {px}m x {py}m")
    print(f"  Origin:  ({x_min:.1f}, {y_max:.1f})  [UL]")
    print(f"  Extent:  x = [{x_min:.1f}, {x_max:.1f}], "
          f"y = [{y_min:.1f}, {y_max:.1f}]")
    print(f"  Valid:   {n_valid} cells inside watershed (out of {ny*nx})")

    ds = None
    return {
        "shape": (ny, nx),
        "geotransform": gt,
        "projection": proj,
        "valid_mask": valid_mask,
        "extent": (x_min, y_min, x_max, y_max),
        "pixel_size": (px, py),
    }


# ==============================================================
# Core: aggregate RdNBR onto the DEM target grid via GDAL Warp
# ==============================================================
def aggregate_rdnbr_to_dem_grid(rdnbr_path, dem_grid, work_dir):
    """
    Resamples the 10m RdNBR raster onto the DEM 30m grid using GDAL's
    'average' algorithm. Output is exactly aligned to the DEM grid.

    Returns the aggregated RdNBR as a numpy array (NaN where source
    coverage is incomplete), same shape as the DEM.
    """
    print(f"\n[step] Aggregating RdNBR (10m) -> DEM grid (30m) via gdal.Warp")

    x_min, y_min, x_max, y_max = dem_grid["extent"]
    px, py = dem_grid["pixel_size"]

    # Warp to a temp GeoTIFF first (more transparent than warping to
    # memory; the temp file can be inspected for debugging)
    tmp_warped = work_dir / "_tmp_rdnbr_warped_to_DEM_grid.tif"

    warp_options = gdal.WarpOptions(
        format="GTiff",
        outputBounds=(x_min, y_min, x_max, y_max),  # match DEM extent
        xRes=px, yRes=py,                           # match DEM pixel size
        targetAlignedPixels=False,                  # bounds already aligned
        resampleAlg="average",                      # nodata-aware mean
        srcNodata=-3.4e38,                          # RdNBR's NoData
        dstNodata=np.nan,                           # represent gaps as NaN
        dstSRS=dem_grid["projection"],              # force CRS to DEM's
        outputType=gdal.GDT_Float32,
        multithread=True,
    )

    result = gdal.Warp(str(tmp_warped), str(rdnbr_path), options=warp_options)
    if result is None:
        raise RuntimeError("gdal.Warp failed")
    result = None  # close

    # Read back the warped result
    ds = gdal.Open(str(tmp_warped))
    arr = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
    ds = None

    # Sanity: shape must match DEM grid
    expected_shape = dem_grid["shape"]
    if arr.shape != expected_shape:
        raise RuntimeError(
            f"Warped raster shape {arr.shape} != DEM shape {expected_shape}. "
            f"This should not happen; check Warp parameters."
        )

    # Clean up temp file (keep commented out if you want to inspect it)
    try:
        tmp_warped.unlink()
    except OSError:
        pass

    n_valid = int(np.sum(~np.isnan(arr)))
    print(f"  Aggregated raster: {arr.shape}, {n_valid} cells with valid RdNBR")
    return arr


# ==============================================================
# Core: classify aggregated RdNBR into 3 burn-severity classes
# ==============================================================
def classify(rdnbr_agg, valid_mask):
    """
    Returns:
      - veg_bs: int8 array, same shape, with 0 outside basin and 1/2/3 inside
      - n_gap:  count of cells inside basin where RdNBR is NaN (gap cells,
                defaulted to class 1)
    """
    print(f"\n[step] Classifying into 3 burn-severity classes")
    print(f"  Class 1 (Unburned-Low): RdNBR < {THR_LOW_MOD}")
    print(f"  Class 2 (Moderate):     {THR_LOW_MOD} <= RdNBR <= {THR_MOD_HIGH}")
    print(f"  Class 3 (High):         RdNBR > {THR_MOD_HIGH}")

    veg_bs = np.zeros(rdnbr_agg.shape, dtype=np.int8)

    classifiable = valid_mask & ~np.isnan(rdnbr_agg)
    veg_bs[classifiable & (rdnbr_agg <  THR_LOW_MOD)]                                   = 1
    veg_bs[classifiable & (rdnbr_agg >= THR_LOW_MOD) & (rdnbr_agg <= THR_MOD_HIGH)]    = 2
    veg_bs[classifiable & (rdnbr_agg >  THR_MOD_HIGH)]                                  = 3

    # Edge case: cells inside watershed but no RdNBR coverage.
    # Default to class 1 (assume no fire effect) and report.
    gap_cells = valid_mask & np.isnan(rdnbr_agg)
    n_gap = int(gap_cells.sum())
    if n_gap > 0:
        print(f"  WARNING: {n_gap} watershed cells have no RdNBR coverage")
        print(f"           Defaulting these to class 1 (Unburned-Low).")
        print(f"           If n_gap is large, inspect _tmp_rdnbr_warped*.tif")
        veg_bs[gap_cells] = 1

    return veg_bs, n_gap


# ==============================================================
# Output: write .bin, .tif, and .txt
# ==============================================================
def write_outputs(veg_bs, rdnbr_agg, dem_grid, n_gap, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. veg_bs.bin (DHSVM-ready flat binary, int8) ---
    bin_path = out_dir / "veg_bs.bin"
    veg_bs.tofile(bin_path)
    print(f"\n[step] Wrote {bin_path.name}  ({veg_bs.nbytes} bytes, int8)")

    # --- 2. veg_bs.tif (diagnostic GeoTIFF for QGIS visualization) ---
    tif_path = out_dir / "veg_bs.tif"
    ny, nx = dem_grid["shape"]
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(str(tif_path), nx, ny, 1, gdal.GDT_Byte,
                       options=["COMPRESS=LZW"])
    ds.SetGeoTransform(dem_grid["geotransform"])
    ds.SetProjection(dem_grid["projection"])
    band = ds.GetRasterBand(1)
    band.WriteArray(veg_bs.astype(np.uint8))
    band.SetNoDataValue(0)
    band.FlushCache()
    ds = None
    print(f"[step] Wrote {tif_path.name}  (GeoTIFF for QGIS inspection)")

    # --- 3. veg_bs_summary.txt (statistics + cross-check) ---
    txt_path = out_dir / "veg_bs_summary.txt"
    valid_mask = dem_grid["valid_mask"]
    n_inside = int(valid_mask.sum())

    # Basin-mean aggregated RdNBR (for cross-check vs Caldwell 2020 = 292)
    basin_rdnbr = rdnbr_agg[valid_mask]
    basin_rdnbr_finite = basin_rdnbr[~np.isnan(basin_rdnbr)]
    if len(basin_rdnbr_finite) > 0:
        basin_mean = float(np.mean(basin_rdnbr_finite))
        basin_min  = float(np.min(basin_rdnbr_finite))
        basin_max  = float(np.max(basin_rdnbr_finite))
    else:
        basin_mean = basin_min = basin_max = float("nan")

    with open(txt_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("VEG_BS.BIN GENERATION SUMMARY\n")
        f.write("DHSVM Experiment B - 3-class burn severity vegetation map\n")
        f.write("=" * 70 + "\n\n")

        f.write("Source: Sentinel-2 RdNBR (Miller & Thode 2007 formulation)\n")
        f.write(f"        Pre-fire image:  June 9, 2016\n")
        f.write(f"        Post-fire image: July 26, 2017\n\n")

        f.write("Classification thresholds (Caldwell et al. 2020, Table 2):\n")
        f.write(f"  Class 1 (Unburned-Low):  RdNBR < {THR_LOW_MOD}\n")
        f.write(f"  Class 2 (Moderate):      {THR_LOW_MOD} <= RdNBR <= {THR_MOD_HIGH}\n")
        f.write(f"  Class 3 (High):          RdNBR > {THR_MOD_HIGH}\n\n")
        f.write("Note: Class 1 corresponds to Caldwell's 'Unburned-Moderate' but is\n")
        f.write("renamed here to avoid confusion with our Class 2. Class 2 merges\n")
        f.write("Caldwell's Moderate and Moderate-High classes.\n\n")

        f.write("Grid:\n")
        f.write(f"  Shape:           {dem_grid['shape'][0]} rows x {dem_grid['shape'][1]} cols\n")
        f.write(f"  Pixel size:      {dem_grid['pixel_size'][0]}m\n")
        f.write(f"  Cells inside basin: {n_inside}\n")
        f.write(f"  Watershed area:  {n_inside * dem_grid['pixel_size'][0]**2 / 1e6:.3f} km^2\n\n")

        f.write("Class distribution (inside watershed):\n")
        for cls in (1, 2, 3):
            cnt = int((veg_bs == cls).sum())
            pct = 100 * cnt / n_inside if n_inside > 0 else 0
            expect = CA_REFERENCE["class_pct_expected"][cls]
            diff = pct - expect
            f.write(f"  Class {cls} ({CLASS_NAMES[cls]:14s}): "
                    f"{cnt:>6} cells ({pct:5.1f}%)   "
                    f"Caldwell 2020 CA: ~{expect:.0f}%   diff: {diff:+5.1f}%\n")

        if n_gap > 0:
            f.write(f"\nGap cells (basin cells with no RdNBR coverage): {n_gap}\n")
            f.write(f"  These were defaulted to class 1.\n")

        f.write("\nCross-check: basin-mean aggregated RdNBR\n")
        f.write(f"  This run:        {basin_mean:.1f}\n")
        f.write(f"  Caldwell 2020:   292 (CA watershed)\n")
        f.write(f"  Difference:      {basin_mean - 292.0:+.1f}\n")
        f.write(f"  RdNBR range in basin: [{basin_min:.1f}, {basin_max:.1f}]\n\n")

        f.write("Reference:\n")
        f.write("  Caldwell, P.V., et al. 2020. Watershed-scale vegetation, water\n")
        f.write("  quantity, and water quality responses to wildfire in the southern\n")
        f.write("  Appalachian mountain region, United States. Hydrological Processes\n")
        f.write("  34: 5188-5209. https://doi.org/10.1002/hyp.13922\n")

    print(f"[step] Wrote {txt_path.name}")

    # Echo cross-check to console for immediate feedback
    print(f"\n" + "=" * 70)
    print("CROSS-CHECK against Caldwell 2020 Table 1 (CA watershed)")
    print("=" * 70)
    print(f"  Basin-mean RdNBR: {basin_mean:.1f}   (Caldwell: 292,  diff: {basin_mean-292:+.1f})")
    print(f"  Class distribution:")
    for cls in (1, 2, 3):
        cnt = int((veg_bs == cls).sum())
        pct = 100 * cnt / n_inside if n_inside > 0 else 0
        expect = CA_REFERENCE["class_pct_expected"][cls]
        marker = "  OK" if abs(pct - expect) < 10 else "  CHECK"
        print(f"    Class {cls} ({CLASS_NAMES[cls]:14s}): "
              f"{pct:5.1f}%   (Caldwell: ~{expect:.0f}%)  {marker}")


# ==============================================================
# Main pipeline
# ==============================================================
def main():
    print("\n" + "=" * 70)
    print("  STEP 1: GENERATE 3-CLASS BURN-SEVERITY VEGETATION MAP")
    print("=" * 70)

    # --- Pre-flight: input file existence ---
    if not DEM_IN.exists():
        sys.exit(f"[ERROR] DEM not found: {DEM_IN}")
    if not RDNBR_IN.exists():
        sys.exit(f"[ERROR] RdNBR not found: {RDNBR_IN}")

    print(f"\nInputs:")
    print(f"  DEM:    {DEM_IN}")
    print(f"  RdNBR:  {RDNBR_IN}")
    print(f"Output:")
    print(f"  Dir:    {OUT_DIR}")

    # --- CRS check ---
    crs_dem   = get_crs_info(DEM_IN)
    crs_rdnbr = get_crs_info(RDNBR_IN)
    verify_crs_match(crs_dem, crs_rdnbr, name_a="DEM", name_b="RdNBR")

    # --- Read DEM grid ---
    dem_grid = read_dem_grid(DEM_IN)

    # --- Aggregate RdNBR onto DEM grid ---
    rdnbr_agg = aggregate_rdnbr_to_dem_grid(RDNBR_IN, dem_grid, OUT_DIR)

    # --- Classify ---
    veg_bs, n_gap = classify(rdnbr_agg, dem_grid["valid_mask"])

    # --- Write outputs ---
    write_outputs(veg_bs, rdnbr_agg, dem_grid, n_gap, OUT_DIR)

    print("\n" + "=" * 70)
    print("  STEP 1 COMPLETE")
    print("=" * 70)
    print(f"\nNext steps:")
    print(f"  1. Open veg_bs.tif in QGIS, overlay on DEM, sanity-check pattern")
    print(f"     against Caldwell 2020 Figure 1b (CA watershed burn severity)")
    print(f"  2. Compare class percentages above against Caldwell Table 1")
    print(f"  3. If both pass, proceed to Step 2 (per-class LAI from MODIS)")
    print()


# ==============================================================
# Trigger (works in both QGIS Python Console and standalone)
# ==============================================================
main()