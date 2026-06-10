# -*- coding: utf-8 -*-
# =====================================================================
# soildepth.py  -  standalone (no-QGIS, DCC) port of soildepthscript.py
#
# Stage: soil depth (PNNL weighting formula). Pure NumPy + GDAL.
# Computes depth from elevation, slope, and flow accumulation, writes the
# DHSVM flat binary (soildepth.bin) plus a diagnostic GeoTIFF (soildepth.tif).
#
# Byte-parity port. The PNNL parameter block and the generate_soildepth math
# are copied verbatim from soildepthscript.py. Only changed: path resolution
# (now via paths.py) and the run-on-import guard -> run().
#
# IMPORTANT (slope input): prep_dhsvm_inputs.py passes `stream_slope` to
# generate_soildepth, but slope_fill ran earlier and overwrote stream_slope.tif
# in place, so the file soildepth actually reads is the FILLED slope. The
# standalone feeds SLOPE_FILLED (the slope stage output), or the 318 boundary
# cells slope_fill repaired would diverge -- the exact slope=0 -> max-depth ring
# artifact Tier B fixed.
#
# Inputs are all validated stage outputs:
#   elev  -> clip stage (ELEV_CLIPPED)
#   slope -> slope stage (SLOPE_FILLED)
#   facc  -> hydrology stage (FLOW_ACC)
# soildepth.bin goes in BIN_DIR (shared with the base-map binaries, decision A).
# =====================================================================

import os
import numpy as np
from osgeo import gdal

from paths import ELEV_CLIPPED, SLOPE_FILLED, FLOW_ACC, BIN_DIR, SOILDEPTH_TIF

ELEV_TIF  = ELEV_CLIPPED
SLOPE_TIF = SLOPE_FILLED          # FILLED -- see header
FACC_TIF  = FLOW_ACC
OUT_BIN   = BIN_DIR / "soildepth.bin"
OUT_TIF   = SOILDEPTH_TIF

# =====================================================================
# PNNL Soil Depth Parameters  (copied verbatim from soildepthscript.py)
# Southern Appalachians calibration (VSA theory + Coweeta field obs).
# =====================================================================
MIN_DEPTH = 2.0
MAX_DEPTH = 6.0

WT_SLOPE  = 0.7
WT_SOURCE = 0.0
WT_ELEV   = 0.3

MAX_SLOPE  = 30.0
MAX_SOURCE = 100000.0
MAX_ELEV   = 1500.0

POW_SLOPE  = 0.25
POW_SOURCE = 1.0
POW_ELEV   = 0.75

DHSVM_NODATA = -9999.0

# =====================================================================
# Helper & main  (verbatim math; only paths/IO wiring differs)
# =====================================================================
def read_raster(filepath):
    ds = gdal.Open(str(filepath))
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float32)
    return arr, band.GetNoDataValue(), ds.GetGeoTransform(), ds.GetProjection(), ds.RasterXSize, ds.RasterYSize


def generate_soildepth(elev_path, slope_path, flowacc_path, out_bin, out_tif):
    print("  [soildepth] Reading spatial inputs into memory...")
    elev_arr, elev_nd, gt, proj, cols, rows = read_raster(elev_path)
    slope_arr, slope_nd, _, _, _, _ = read_raster(slope_path)
    fac_arr, fac_nd, _, _, _, _ = read_raster(flowacc_path)

    # 1. Absolute master mask: derived purely from the DEM
    valid_mask = (elev_arr != elev_nd) & (~np.isnan(elev_arr))

    # 2. Gap-fill: force missing boundary pixels in slope/fac to 0.0 to match DEM mask
    print("  [soildepth] Aligning boundaries and sanitizing edge pixels...")
    elev_clean = np.where(valid_mask & (elev_arr >= 0), elev_arr, 0.0)
    slope_clean = np.where(valid_mask & (slope_arr != slope_nd) & (~np.isnan(slope_arr)) & (slope_arr >= 0), slope_arr, 0.0)
    fac_clean = np.where(valid_mask & (fac_arr != fac_nd) & (~np.isnan(fac_arr)), np.abs(fac_arr), 0.0)

    # 3. Limit processing bounds
    elev_lim = np.clip(elev_clean, 0.0, MAX_ELEV)
    slope_lim = np.clip(slope_clean, 0.0, MAX_SLOPE)
    fac_lim = np.clip(fac_clean, 0.0, MAX_SOURCE)

    # 4. Physical calculation
    print("  [soildepth] Computing weighted topographic functions...")
    term_slope  = WT_SLOPE * (1.0 - (slope_lim / MAX_SLOPE) ** POW_SLOPE)
    term_source = WT_SOURCE * ((fac_lim / MAX_SOURCE) ** POW_SOURCE)
    term_elev   = WT_ELEV * (1.0 - (elev_lim / MAX_ELEV) ** POW_ELEV)

    depth_calc = MIN_DEPTH + (MAX_DEPTH - MIN_DEPTH) * (term_slope + term_source + term_elev)

    # 5. Output array construction
    final_depth = np.full(elev_arr.shape, DHSVM_NODATA, dtype=np.float32)
    final_depth[valid_mask] = np.clip(depth_calc[valid_mask], MIN_DEPTH, MAX_DEPTH)

    print(f"  [soildepth] Writing flat binary -> {os.path.basename(str(out_bin))}")
    final_depth.tofile(str(out_bin))

    print(f"  [soildepth] Writing GeoTIFF -> {os.path.basename(str(out_tif))}")
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(str(out_tif), cols, rows, 1, gdal.GDT_Float32)
    out_ds.SetGeoTransform(gt)
    out_ds.SetProjection(proj)
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(final_depth)
    out_band.SetNoDataValue(DHSVM_NODATA)
    out_band.FlushCache()
    out_ds = None
    print("  [soildepth] Complete.")


def run():
    print("\n--- Dynamic Soil Depth Module (standalone) ---")
    for tag, pth in [("elev", ELEV_TIF), ("slope(filled)", SLOPE_TIF), ("flow_acc", FACC_TIF)]:
        if not pth.exists():
            raise FileNotFoundError(f"[error] missing {tag} input: {pth}")
    BIN_DIR.mkdir(parents=True, exist_ok=True)
    generate_soildepth(ELEV_TIF, SLOPE_TIF, FACC_TIF, OUT_BIN, OUT_TIF)


if __name__ == "__main__":
    run()
