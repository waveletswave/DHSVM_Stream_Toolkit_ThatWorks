# -*- coding: utf-8 -*-
# =====================================================================
# slope_fill.py
# PURPOSE:
#   Condition the slope raster before any consumer reads it. r.slope.aspect
#   leaves NoData on the outermost ring of the clipped DEM (no full 3x3
#   neighbourhood), and on cells ringing former interior DEM voids. Those
#   cells are valid in the DEM. soildepthscript previously filled them with
#   slope=0, and slope=0 receives the MAXIMUM slope bonus (term_slope=WT_SLOPE),
#   producing a spurious deep ring. This module fills the missing-slope cells
#   with the mean of their valid 8-neighbours instead, recovering a physically
#   reasonable slope and leaving the depth formula untouched.
#
#   Portable: pure GDAL + NumPy, no QGIS. Callable from the QGIS orchestrator
#   (prep_dhsvm_inputs.py) and from the standalone/cluster pipeline.
#
# AUTHOR: Y. Song (slope-conditioning step)
# =====================================================================

import numpy as np
from osgeo import gdal


def _neighbour_mean_fill(arr, fillable, max_iters=8):
    """Fill `fillable` cells with the mean of valid 8-neighbours, iterating
    inward. `arr` uses NaN for invalid cells. Returns (filled_array,
    filled_mask, still_unfilled_mask)."""
    work = arr.copy()
    remaining = fillable.copy()
    filled_any = np.zeros_like(fillable)
    shifts = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    for _ in range(max_iters):
        if not remaining.any():
            break
        nb_sum = np.zeros_like(work)
        nb_cnt = np.zeros(work.shape, dtype=float)
        for dy, dx in shifts:
            sh = np.full_like(work, np.nan)
            ys = slice(max(0, dy), work.shape[0] + min(0, dy))
            xs = slice(max(0, dx), work.shape[1] + min(0, dx))
            yd = slice(max(0, -dy), work.shape[0] + min(0, -dy))
            xd = slice(max(0, -dx), work.shape[1] + min(0, -dx))
            sh[yd, xd] = work[ys, xs]
            ok = ~np.isnan(sh)
            nb_sum[ok] += sh[ok]
            nb_cnt[ok] += 1.0
        can = remaining & (nb_cnt > 0)
        work[can] = nb_sum[can] / nb_cnt[can]
        filled_any |= can
        remaining &= ~can
    return work, filled_any, remaining


def fill_slope_nodata(slope_path, dem_path, out_path=None, max_iters=8, verbose=True):
    """Read slope + DEM, fill (DEM-valid AND slope-NoData) cells by neighbour
    interpolation, write the result. If out_path is None, overwrites slope_path.
    Returns a dict of counts. Pure GDAL + NumPy.

    The DEM passed here must be the same grid the slope was derived from
    (in prep_dhsvm_inputs.py that is `elev`)."""
    sds = gdal.Open(slope_path)
    if sds is None:
        raise RuntimeError(f"cannot open slope raster: {slope_path}")
    sb = sds.GetRasterBand(1)
    slope = sb.ReadAsArray().astype(np.float64)
    s_nd = sb.GetNoDataValue()
    gt, proj = sds.GetGeoTransform(), sds.GetProjection()
    cols, rows = sds.RasterXSize, sds.RasterYSize

    dds = gdal.Open(dem_path)
    if dds is None:
        raise RuntimeError(f"cannot open DEM: {dem_path}")
    db = dds.GetRasterBand(1)
    elev = db.ReadAsArray().astype(np.float64)
    e_nd = db.GetNoDataValue()
    if elev.shape != slope.shape:
        raise RuntimeError(f"grid mismatch: slope {slope.shape} vs DEM {elev.shape}")

    valid_dem = np.isfinite(elev) & ((elev != e_nd) if e_nd is not None else True)
    slope_valid = np.isfinite(slope) & ((slope != s_nd) if s_nd is not None else True) & (slope >= 0)
    missing = valid_dem & ~slope_valid
    n_missing = int(missing.sum())

    nd_out = s_nd if s_nd is not None else -9999.0
    if n_missing == 0:
        if verbose:
            print("  [slope_fill] no missing-slope cells; raster already complete.")
        # still (re)write if a different out_path was requested
        if out_path and out_path != slope_path:
            _write(out_path, np.where(valid_dem, slope, nd_out), gt, proj, cols, rows, nd_out)
        sds = dds = None
        return {"missing": 0, "filled": 0, "unfilled": 0}

    work = np.where(slope_valid, slope, np.nan)
    filled, filled_mask, still = _neighbour_mean_fill(work, missing, max_iters)
    n_filled = int(filled_mask.sum())
    n_left = int((still & missing).sum())

    out_arr = np.where(valid_dem, filled, nd_out)
    out_arr = np.where(np.isnan(out_arr), nd_out, out_arr).astype(np.float32)

    target = out_path or slope_path
    _write(target, out_arr, gt, proj, cols, rows, nd_out)
    sds = dds = None

    if verbose:
        fv = filled[filled_mask]
        print(f"  [slope_fill] filled {n_filled} missing-slope cells "
              f"(of {n_missing}); {n_left} unfilled. "
              f"filled slope mean {np.nanmean(fv):.2f} deg -> {target}")
    return {"missing": n_missing, "filled": n_filled, "unfilled": n_left}


def _write(path, arr, gt, proj, cols, rows, nodata):
    drv = gdal.GetDriverByName('GTiff')
    ds = drv.Create(path, cols, rows, 1, gdal.GDT_Float32)
    ds.SetGeoTransform(gt)
    ds.SetProjection(proj)
    bnd = ds.GetRasterBand(1)
    bnd.WriteArray(arr.astype(np.float32))
    bnd.SetNoDataValue(float(nodata))
    bnd.FlushCache()
    ds = None


# Optional standalone use (mirrors fix_slope_nodata.py): writes _filled by default
if __name__ == "__main__":
    import argparse
    from pathlib import Path
    ap = argparse.ArgumentParser()
    ap.add_argument('--slope', required=True)
    ap.add_argument('--dem', required=True)
    ap.add_argument('--out', default=None,
                    help="default: <slope>_filled.tif (does NOT overwrite)")
    a = ap.parse_args()
    out = a.out or str(Path(a.slope).with_name(Path(a.slope).stem + "_filled.tif"))
    fill_slope_nodata(a.slope, a.dem, out_path=out)