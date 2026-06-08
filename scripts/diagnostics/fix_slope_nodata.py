#!/usr/bin/env python3
"""
fix_slope_nodata.py  --  repair the boundary-ring NoData in stream_slope.tif

Problem
-------
r.slope.aspect cannot compute slope on the outermost ring of the clipped DEM
(no full 3x3 neighbourhood), so those cells are NoData in stream_slope.tif
while the DEM is valid. soildepthscript.py fills that NoData with 0, and
slope=0 gets the MAXIMUM slope bonus (term_slope = WT_SLOPE), producing a
spurious ~5 m depth ring. This is the 3.3 artifact.

Fix (neighbour interpolation)
-----------------------------
Treat the ring as real terrain whose slope was simply not computable at the
edge: fill each missing-slope cell with the mean of its valid 8-neighbours
(iterated inward for any multi-cell gaps). This is deterministic, transparent,
and leaves the depth formula and its parameters untouched, so it is orthogonal
to the uniform-vs-variable decision. Writes a NEW raster (stream_slope_filled.tif);
the original is not modified.

Usage
-----
    python fix_slope_nodata.py \
        --slope ../Intermediate_GIS/stream_slope.tif \
        --dem   ../Intermediate_GIS/elev_clipped.tif

  (zero-arg also works: auto-resolves the same paths soildepthscript.py uses,
   WS = this script's parent-of-parent dir.)

Options
    --out <path>        output raster (default: <slope dir>/stream_slope_filled.tif)
    --max-iters N       inward fill passes (default 6; a 1-cell ring needs 1)
    --no-depth-check    skip the soildepth before/after report

Read/write uses rasterio (present in the dhsvm_rs env).
"""

import argparse
import sys
from pathlib import Path
import numpy as np

# ----- soildepth formula constants (match qgis_CA/soildepthscript.py) --------
MIN_DEPTH, MAX_DEPTH = 2.0, 6.0
WT_SLOPE, WT_ELEV    = 0.7, 0.3          # WT_SOURCE = 0 at current params
MAX_SLOPE, POW_SLOPE = 30.0, 0.25
MAX_ELEV,  POW_ELEV  = 1500.0, 0.75


def autoresolve():
    WS = Path(__file__).resolve().parent.parent
    dem_c   = [WS / "Intermediate_GIS" / "elev_clipped.tif",
               WS / "Reprojected_DEM" / "elev_clipped.tif",
               WS / "dem.tif"]
    slope_c = [WS / "Intermediate_GIS" / "stream_slope.tif",
               WS / "stream_slope.tif"]
    dem   = next((str(c) for c in dem_c   if c.exists()), None)
    slope = next((str(c) for c in slope_c if c.exists()), None)
    return dem, slope


def neighbour_mean_fill(arr, fillable, max_iters):
    """Fill `fillable` cells with the mean of their valid 8-neighbours,
    iterating inward. `arr` uses NaN for all invalid cells. Returns the
    filled array and the boolean mask of cells that were actually filled."""
    work = arr.copy()
    remaining = fillable.copy()
    filled_any = np.zeros_like(fillable)
    shifts = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
    for _ in range(max_iters):
        if not remaining.any():
            break
        nb_sum = np.zeros_like(work)
        nb_cnt = np.zeros(work.shape, dtype=float)
        for dy, dx in shifts:
            sh = np.full_like(work, np.nan)
            ys = slice(max(0,dy), work.shape[0] + min(0,dy))
            xs = slice(max(0,dx), work.shape[1] + min(0,dx))
            yd = slice(max(0,-dy), work.shape[0] + min(0,-dy))
            xd = slice(max(0,-dx), work.shape[1] + min(0,-dx))
            sh[yd, xd] = work[ys, xs]
            valid = ~np.isnan(sh)
            nb_sum[valid] += sh[valid]
            nb_cnt[valid] += 1.0
        can = remaining & (nb_cnt > 0)
        work[can] = nb_sum[can] / nb_cnt[can]
        filled_any |= can
        remaining &= ~can
    return work, filled_any, remaining


def soildepth(slope_arr, elev_arr, valid):
    """Faithful current-parameter depth (WT_SOURCE=0). NaN outside `valid`."""
    el = np.clip(elev_arr, 0.0, MAX_ELEV)
    sl = np.clip(slope_arr, 0.0, MAX_SLOPE)
    ts = WT_SLOPE * (1.0 - (sl / MAX_SLOPE) ** POW_SLOPE)
    te = WT_ELEV  * (1.0 - (el / MAX_ELEV) ** POW_ELEV)
    d  = MIN_DEPTH + (MAX_DEPTH - MIN_DEPTH) * (ts + te)
    out = np.full(slope_arr.shape, np.nan)
    out[valid] = np.clip(d[valid], MIN_DEPTH, MAX_DEPTH)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slope', default=None)
    ap.add_argument('--dem',   default=None)
    ap.add_argument('--out',   default=None)
    ap.add_argument('--max-iters', type=int, default=6)
    ap.add_argument('--no-depth-check', action='store_true')
    a = ap.parse_args()

    try:
        import rasterio
    except ModuleNotFoundError:
        sys.exit("[error] rasterio not importable. Run in the env that has it:\n"
                 "    conda activate dhsvm_rs")

    slope_fp, dem_fp = a.slope, a.dem
    if not (slope_fp and dem_fp):
        dem_r, slope_r = autoresolve()
        slope_fp = slope_fp or slope_r
        dem_fp   = dem_fp or dem_r
        if not (slope_fp and dem_fp):
            sys.exit("[error] could not resolve --slope/--dem; pass them explicitly.")
        print(f">>> auto-resolved\n    slope: {slope_fp}\n    dem  : {dem_fp}\n")

    out_fp = a.out or str(Path(slope_fp).with_name("stream_slope_filled.tif"))

    with rasterio.open(slope_fp) as src:
        profile = src.profile
        slope = src.read(1).astype(np.float64)
        s_nd = src.nodata
    with rasterio.open(dem_fp) as src:
        elev = src.read(1).astype(np.float64)
        e_nd = src.nodata
    if slope.shape != elev.shape:
        sys.exit(f"[error] grid mismatch: slope {slope.shape} vs dem {elev.shape}. "
                 f"Use the grid-matched Intermediate_GIS rasters.")

    # masks
    valid_dem = np.isfinite(elev)
    if e_nd is not None: valid_dem &= (elev != e_nd)
    slope_valid = np.isfinite(slope)
    if s_nd is not None: slope_valid &= (slope != s_nd)
    missing = valid_dem & ~slope_valid          # the 318 cells

    # --- diagnose ------------------------------------------------------------
    n_missing = int(missing.sum())
    print(f"DEM-valid cells           : {int(valid_dem.sum())}")
    print(f"slope NoData & DEM-valid  : {n_missing}  "
          f"({100*n_missing/max(1,valid_dem.sum()):.2f}% of DEM-valid)")
    if n_missing == 0:
        print("Nothing to fill. stream_slope has no boundary-ring NoData.")
        return

    # boundary vs interior: a missing cell touching an outside-basin neighbour
    outside = ~valid_dem
    touch_out = np.zeros_like(missing)
    for dy in (-1,0,1):
        for dx in (-1,0,1):
            if dy==0 and dx==0: continue
            sh = np.zeros_like(outside)
            ys = slice(max(0,dy), outside.shape[0]+min(0,dy))
            xs = slice(max(0,dx), outside.shape[1]+min(0,dx))
            yd = slice(max(0,-dy), outside.shape[0]+min(0,-dy))
            xd = slice(max(0,-dx), outside.shape[1]+min(0,-dx))
            sh[yd, xd] = outside[ys, xs]
            touch_out |= sh
    on_boundary = int((missing & touch_out).sum())
    interior    = int((missing & ~touch_out).sum())
    print(f"  on basin boundary       : {on_boundary}  (expected: ~all of them)")
    print(f"  interior holes          : {interior}  (expected: 0)\n")

    # --- fill ----------------------------------------------------------------
    work = np.where(slope_valid, slope, np.nan)
    filled, filled_mask, still = neighbour_mean_fill(work, missing, a.max_iters)
    n_filled = int(filled_mask.sum()); n_left = int((still & missing).sum())
    print(f"filled {n_filled} cells by neighbour mean "
          f"({a.max_iters} passes); {n_left} still unfilled")
    if n_left:
        print("  (unfilled cells had no valid neighbour within range; "
              "raise --max-iters or check the slope footprint)")
    fv = filled[filled_mask]
    iv = slope[slope_valid]
    print(f"  filled-cell slope (deg) : mean {fv.mean():.2f}  min {fv.min():.2f}  "
          f"max {fv.max():.2f}   | interior slope mean {iv.mean():.2f}\n")

    # --- write corrected raster ---------------------------------------------
    out_arr = np.where(valid_dem, filled, (s_nd if s_nd is not None else -9999.0))
    out_arr = np.where(np.isnan(out_arr),
                       (s_nd if s_nd is not None else -9999.0), out_arr)
    profile.update(dtype='float32', nodata=(s_nd if s_nd is not None else -9999.0),
                   count=1)
    with rasterio.open(out_fp, 'w', **profile) as dst:
        dst.write(out_arr.astype('float32'), 1)
    print(f"[written] {out_fp}\n")

    # --- soildepth before/after ---------------------------------------------
    if a.no_depth_check:
        return
    slope_before = np.where(slope_valid, slope, 0.0)   # pipeline fills NoData->0
    slope_after  = filled
    d_before = soildepth(slope_before, elev, valid_dem)
    d_after  = soildepth(slope_after,  elev, valid_dem)

    others = valid_dem & ~missing
    max_other = float(np.nanmax(np.abs(d_before[others] - d_after[others])))
    db_m, da_m = d_before[missing], d_after[missing]
    print("soildepth before/after (current params MS30 POW0.25 WS0):")
    print(f"  the {n_missing} ring cells : depth {db_m.mean():.2f} m (max {db_m.max():.2f}) "
          f"-> {da_m.mean():.2f} m (max {da_m.max():.2f})")
    print(f"  other {int(others.sum())} cells : max|change| {max_other:.3g}  "
          f"(expected 0 -> byte-identical)")
    bv, av = d_before[valid_dem], d_after[valid_dem]
    print(f"  basin mean   {bv.mean():.2f} -> {av.mean():.2f} m")
    print(f"  basin median {np.median(bv):.2f} -> {np.median(av):.2f} m")
    print(f"  basin p95    {np.percentile(bv,95):.2f} -> {np.percentile(av,95):.2f} m")
    print(f"  basin max    {bv.max():.2f} -> {av.max():.2f} m")
    # 3.3 artifact count after: deepest cells should no longer be ring cells
    topN = 20
    flat = np.where(valid_dem, d_after, -1).ravel()
    idx = np.argsort(flat)[-topN:]
    n_ring_top = int(missing.ravel()[idx].sum())
    print(f"  of the {topN} deepest cells after fill, {n_ring_top} are ring cells "
          f"(was 20)\n")
    print("Next: if the before/after looks right, regenerate soildepth from the "
          "filled slope (or fold the fill into createstreamnetwork right after "
          "r.slope.aspect). The regression check is that only these "
          f"{n_missing} cells change.")


if __name__ == "__main__":
    main()
