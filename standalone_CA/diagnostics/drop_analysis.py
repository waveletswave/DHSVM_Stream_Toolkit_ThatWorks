# -*- coding: utf-8 -*-
# =====================================================================
# drop_analysis.py  -  objective stream threshold by constant stream drop
#
# Determines the stream initiation threshold on a principled basis instead of
# by visual inspection, then expresses it as a physical support area so it
# transfers across resolutions and watersheds.
#
# Method (Broscoe 1959; Tarboton, Bras, and Rodriguez-Iturbe 1991):
#   The constant stream drop law holds that the mean elevation drop along
#   streams of each Strahler order is statistically the same. We sweep a range
#   of support-area thresholds. For each threshold we extract the stream
#   network, assign Strahler order, and measure the drop of every stream. We
#   test whether the mean drop of first-order streams differs from the higher
#   orders by a t-test. The smallest threshold for which the difference is not
#   significant (absolute t below 2) is the finest network consistent with the
#   law, and is taken as the objective threshold.
#
# Flow routing here is pyflwdir's D8, built from the DEM. This is independent
# of the GRASS MFD routing in the pipeline by design: the analysis yields a
# physical support area A_c, an algorithm-independent geomorphic quantity. The
# pipeline consumes that area and computes its own cell threshold with its own
# routing and cell size. A_c is the bridge between the two.
#
# Output:
#   - a table to stdout: threshold, cells, area, drainage density, |t|, pass
#   - the objective threshold and its A_c, against the current visual value
#   - a CSV of the full sweep for plotting (quick-look of |t| vs threshold)
#
# Input is one DEM (the clipped DEM from the clip stage). Nothing else; this
# does not touch the pipeline outputs.
#
# Run:  python3 drop_analysis.py [elev_clipped.tif] [--csv out.csv]
# =====================================================================

import sys
import argparse
import numpy as np
import rasterio
from scipy import stats
import pyflwdir

# Resolve the default DEM and the current visual threshold from paths.py when
# importable, so this script tracks the pipeline's settings. Fall back to the
# CA defaults when run standalone outside the package.
try:
    from paths import ELEV_CLIPPED, STREAM_SOURCE_AREA_M2
    DEFAULT_DEM = str(ELEV_CLIPPED)
    CURRENT_AREA_M2 = float(STREAM_SOURCE_AREA_M2)
except Exception:
    DEFAULT_DEM = "/work/ys451/dhsvm_ca/standalone_dev/outputs/elev_clipped.tif"
    CURRENT_AREA_M2 = 47571.5


def load_dem(path):
    with rasterio.open(path) as src:
        elev = src.read(1).astype("float64")
        nodata = src.nodata
        transform = src.transform
    cell_area = abs(transform.a * transform.e)
    return elev, nodata, transform, cell_area


def build_flow(elev, nodata, transform):
    """D8 flow direction from the DEM. pyflwdir conditions internally."""
    flw = pyflwdir.from_dem(
        data=elev, nodata=nodata, transform=transform, latlon=False
    )
    return flw


def stream_drops(flw, elev_flat, stream_mask, strord):
    """Drop of every stream segment: elev at head minus elev at outlet.

    streams() returns one feature per stream segment with properties idx (the
    head cell, linear index) and idx_ds (the segment outlet's downstream cell).
    The segment runs from idx down to the cell whose downstream pointer leaves
    the segment. We take head and tail elevation from the segment's own cells.
    Returns a list of (strahler_order, drop).
    """
    feats = flw.streams(mask=stream_mask, strord=strord)
    out = []
    for f in feats:
        props = f["properties"]
        order = int(props["strord"])
        coords_idx = f["geometry"]  # LineString in map coords; use cell idxs
        # The segment's cells in order: head is props['idx']; the polyline
        # vertices correspond to the segment path. Head and tail elevation come
        # from the first and last vertex cell. We recover cell row/col from the
        # vertex coordinates via the inverse transform.
        line = coords_idx["coordinates"]
        # head = first vertex, tail = last vertex
        head_xy = line[0]
        tail_xy = line[-1]
        out.append((order, head_xy, tail_xy))
    return out


def drops_from_coords(segments, elev, transform):
    """Convert (order, head_xy, tail_xy) to (order, drop) using DEM at the
    head and tail cells. drop = elev[head] - elev[tail]."""
    inv = ~transform
    rows = []
    for order, head_xy, tail_xy in segments:
        hc, hr = inv * (head_xy[0], head_xy[1])
        tc, tr = inv * (tail_xy[0], tail_xy[1])
        hr, hc = int(hr), int(hc)
        tr, tc = int(tr), int(tc)
        if not (0 <= hr < elev.shape[0] and 0 <= hc < elev.shape[1]):
            continue
        if not (0 <= tr < elev.shape[0] and 0 <= tc < elev.shape[1]):
            continue
        drop = elev[hr, hc] - elev[tr, tc]
        rows.append((order, drop))
    return rows


def evaluate_threshold(flw, elev, transform, uparea_cells, cell_area,
                       thresh_cells, basin_cells):
    """One threshold: extract streams, Strahler order, drops, t-test.

    Returns dict with cells, area, n_streams, n_order1, n_higher, drainage
    density, t statistic, pass flag, and negative-drop count.
    """
    stream_mask = uparea_cells > thresh_cells
    n_stream_cells = int(stream_mask.sum())
    if n_stream_cells < 2:
        return None
    strord = flw.stream_order(type="strahler", mask=stream_mask)
    segs = stream_drops(flw, elev, stream_mask, strord)
    rows = drops_from_coords(segs, elev, transform)
    if not rows:
        return None
    arr = np.array(rows, dtype="float64")
    orders = arr[:, 0].astype(int)
    drops = arr[:, 1]

    o1 = drops[orders == 1]
    hi = drops[orders >= 2]
    n_neg = int((drops <= 0).sum())

    # drainage density: total stream length over basin area. Approximate stream
    # length as stream-cell count times cell size (a standard cell proxy).
    cell_len = np.sqrt(cell_area)
    total_len_m = n_stream_cells * cell_len
    basin_area_m2 = basin_cells * cell_area
    dd = (total_len_m / basin_area_m2) * 1000.0  # km / km2 = 1/km

    # Welch t-test, first order vs all higher orders.
    if o1.size >= 2 and hi.size >= 2:
        t, p = stats.ttest_ind(o1, hi, equal_var=False)
        t_abs = abs(t)
    else:
        t_abs = np.nan
        p = np.nan

    return dict(
        cells=thresh_cells,
        area_m2=thresh_cells * cell_area,
        area_km2=thresh_cells * cell_area / 1e6,
        n_streams=len(rows),
        n_order1=int(o1.size),
        n_higher=int(hi.size),
        max_order=int(orders.max()),
        dd=dd,
        t_abs=t_abs,
        p=p,
        n_neg=n_neg,
        passes=(not np.isnan(t_abs)) and (t_abs < 2.0),
    )


def main():
    ap = argparse.ArgumentParser(description="Constant stream drop analysis")
    ap.add_argument("dem", nargs="?", default=DEFAULT_DEM)
    ap.add_argument("--tmin", type=int, default=10, help="min threshold (cells)")
    ap.add_argument("--tmax", type=int, default=300, help="max threshold (cells)")
    ap.add_argument("--step", type=int, default=10, help="threshold step (cells)")
    ap.add_argument("--csv", default=None, help="write the sweep to this CSV")
    args = ap.parse_args()

    print("=======================================================")
    print("  CONSTANT STREAM DROP ANALYSIS")
    print("=======================================================")
    print(f"  DEM   : {args.dem}")

    elev, nodata, transform, cell_area = load_dem(args.dem)
    valid = np.isfinite(elev) & (elev != nodata)
    basin_cells = int(valid.sum())
    print(f"  cell area {cell_area:.2f} m2   valid cells {basin_cells}")
    print(f"  current visual threshold: {CURRENT_AREA_M2:.1f} m2 "
          f"= {CURRENT_AREA_M2/cell_area:.1f} cells = {CURRENT_AREA_M2/1e6:.5f} km2")

    flw = build_flow(elev, nodata, transform)
    uparea_cells = np.asarray(flw.upstream_area(unit="cell"))

    print(f"\n  sweep {args.tmin}..{args.tmax} step {args.step} cells")
    print(f"  {'cells':>6} {'km2':>9} {'n_str':>6} {'o1':>4} {'hi':>4} "
          f"{'maxO':>4} {'Dd':>6} {'|t|':>7} {'neg':>4} {'pass':>5}")

    results = []
    for thr in range(args.tmin, args.tmax + 1, args.step):
        r = evaluate_threshold(flw, elev, transform, uparea_cells,
                               cell_area, thr, basin_cells)
        if r is None:
            continue
        results.append(r)
        t_str = f"{r['t_abs']:.2f}" if not np.isnan(r['t_abs']) else "  nan"
        print(f"  {r['cells']:>6} {r['area_km2']:>9.5f} {r['n_streams']:>6} "
              f"{r['n_order1']:>4} {r['n_higher']:>4} {r['max_order']:>4} "
              f"{r['dd']:>6.2f} {t_str:>7} {r['n_neg']:>4} "
              f"{'yes' if r['passes'] else 'no':>5}")

    # Objective threshold: smallest threshold that passes (|t| < 2).
    passing = [r for r in results if r["passes"]]
    print("\n=======================================================")
    if passing:
        obj = min(passing, key=lambda r: r["cells"])
        print(f"  objective threshold (smallest |t|<2):")
        print(f"    {obj['cells']} cells = {obj['area_m2']:.1f} m2 "
              f"= {obj['area_km2']:.5f} km2   (|t|={obj['t_abs']:.2f})")
        print(f"  current visual: {CURRENT_AREA_M2/cell_area:.0f} cells "
              f"= {CURRENT_AREA_M2/1e6:.5f} km2")
        ratio = obj["area_m2"] / CURRENT_AREA_M2
        print(f"  objective / visual = {ratio:.2f}x")
        if ratio > 1.3:
            print(f"  -> visual threshold is too dense; objective wants a larger area")
        elif ratio < 0.77:
            print(f"  -> visual threshold is too sparse; objective wants a smaller area")
        else:
            print(f"  -> visual threshold is within ~30% of the objective value")
    else:
        print("  no threshold in the swept range satisfies |t|<2.")
        print("  widen the range (--tmin/--tmax) or inspect the drop table above.")
    print("=======================================================")

    if args.csv:
        import csv
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["cells", "area_m2", "area_km2", "n_streams",
                        "n_order1", "n_higher", "max_order", "drainage_density",
                        "t_abs", "p", "n_neg", "passes"])
            for r in results:
                w.writerow([r["cells"], f"{r['area_m2']:.1f}",
                            f"{r['area_km2']:.6f}", r["n_streams"],
                            r["n_order1"], r["n_higher"], r["max_order"],
                            f"{r['dd']:.4f}", f"{r['t_abs']:.4f}",
                            f"{r['p']:.4f}", r["n_neg"], int(r["passes"])])
        print(f"  sweep written to {args.csv}")


if __name__ == "__main__":
    main()
