#!/usr/bin/env python3
"""
compare_soildepth.py -- verify the slope-conditioning regression.

Compares OLD vs NEW soildepth.tif and confirms the changes are confined to
the boundary/interior cells the conditioning filled. Cross-checks the changed
cells against the OLD slope's NoData mask: this proves the soildepth changes
come from the conditioning step, not from any non-determinism elsewhere in
the rerun. Read-only.

    python compare_soildepth.py \
        --old       ../Intermediate_GIS/soildepth_OLD.tif \
        --new       ../Intermediate_GIS/soildepth.tif \
        --old-slope ../Intermediate_GIS/stream_slope_OLD.tif \
        --dem       ../Intermediate_GIS/elev_clipped.tif
"""
import argparse, sys
import numpy as np


def read(fp):
    import rasterio
    with rasterio.open(fp) as s:
        return s.read(1).astype(np.float64), s.nodata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--old', required=True)
    ap.add_argument('--new', required=True)
    ap.add_argument('--old-slope', default=None)
    ap.add_argument('--dem', default=None)
    ap.add_argument('--tol', type=float, default=1e-3)
    a = ap.parse_args()
    try:
        import rasterio  # noqa
    except ModuleNotFoundError:
        sys.exit("[error] rasterio not importable. conda activate dhsvm_rs")

    old, old_nd = read(a.old)
    new, new_nd = read(a.new)
    if old.shape != new.shape:
        sys.exit(f"[error] shape mismatch: old {old.shape} vs new {new.shape}")

    valid = np.isfinite(old) & np.isfinite(new)
    if old_nd is not None: valid &= (old != old_nd)
    if new_nd is not None: valid &= (new != new_nd)
    if a.dem:
        dem, dem_nd = read(a.dem)
        vd = np.isfinite(dem) & ((dem != dem_nd) if dem_nd is not None else True)
        valid &= vd

    diff = np.abs(old - new)
    changed = valid & (diff > a.tol)
    n_changed = int(changed.sum())
    n_valid = int(valid.sum())
    print(f"valid cells compared : {n_valid}")
    print(f"cells changed (>{a.tol:g}) : {n_changed}")
    print(f"cells identical      : {n_valid - n_changed}")
    print(f"max |diff| over valid: {float(diff[valid].max()):.6g}")
    if n_changed:
        print(f"  at changed cells: old depth mean {old[changed].mean():.2f} "
              f"(max {old[changed].max():.2f}) -> new mean {new[changed].mean():.2f} "
              f"(max {new[changed].max():.2f})")

    # cross-check: changed cells should equal the OLD slope's NoData-in-basin set
    if a.old_slope:
        slope_old, s_nd = read(a.old_slope)
        s_missing = valid & (~np.isfinite(slope_old) | ((slope_old == s_nd) if s_nd is not None else False))
        n_missing = int(s_missing.sum())
        extra   = int((changed & ~s_missing).sum())   # changed but slope was NOT missing -> bad
        missed  = int((s_missing & ~changed).sum())    # slope missing but depth didn't change
        print(f"\ncross-check vs OLD slope NoData:")
        print(f"  old-slope NoData cells (in basin): {n_missing}")
        print(f"  changed AND slope-was-missing    : {int((changed & s_missing).sum())}")
        print(f"  changed but slope was NOT missing: {extra}  (must be 0)")
        print(f"  slope-missing but depth unchanged: {missed}")
        if extra == 0 and n_changed == n_missing:
            print("  -> CLEAN: every changed cell is a conditioned cell, and only "
                  "those changed. Regression confirmed.")
        elif extra == 0:
            print("  -> changed cells are all conditioned cells (subset); counts "
                  "differ slightly (a conditioned cell may map to the same depth).")
        else:
            print("  -> WARNING: some changed cells were NOT slope-NoData. That "
                  "points to non-determinism upstream (DEM clip / r.watershed / "
                  "r.slope.aspect), not the conditioning. Investigate those.")


if __name__ == "__main__":
    main()
