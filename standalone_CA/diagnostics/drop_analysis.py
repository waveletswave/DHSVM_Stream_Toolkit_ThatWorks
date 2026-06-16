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
#   orders by a t-test. The objective threshold is the start of the first
#   sustained band of consecutive thresholds for which the difference is not
#   significant (absolute t below 2). Requiring a sustained band, rather than
#   the single smallest passing threshold, rejects isolated noise dips where
#   |t| falls below 2 at one threshold by chance. These are common on small
#   basins, where the t-test samples are small and |t| is jumpy.
#
# Flow routing here is pyflwdir's D8, built from the DEM. This is independent
# of the GRASS MFD routing in the pipeline by design: the analysis yields a
# physical support area A_c, an algorithm-independent geomorphic quantity. The
# pipeline consumes that area and computes its own cell threshold with its own
# routing and cell size. A_c is the bridge between the two.
#
# Output:
#   - a table to stdout: threshold, cells, area, drainage density, |t|, pass
#   - the objective threshold and its A_c. If a reference area is given (via
#     --compare-area, or DHSVM_STREAM_SOURCE_AREA_M2 set in the env), the
#     objective is compared against it; otherwise only the objective is shown,
#     with a hint to set DHSVM_STREAM_SOURCE_AREA_M2 from it
#   - a CSV of the full sweep for plotting (quick-look of |t| vs threshold)
#   - with --emit-area: only the objective A_c in m2 to stdout, nothing else,
#     so a script can set DHSVM_STREAM_SOURCE_AREA_M2 from it. Exits nonzero if
#     no threshold in the swept range passes.
#
# Input is one DEM (the clipped DEM from the clip stage). Nothing else; this
# does not touch the pipeline outputs.
#
# Run:  python3 drop_analysis.py [elev_clipped.tif] [--csv out.csv]
#       python3 drop_analysis.py elev_clipped.tif --emit-area --tmin 50 --tmax 800
# =====================================================================

import os
import sys
import argparse
import numpy as np
import rasterio
from scipy import stats
import pyflwdir

# paths.py is the pipeline's single source of truth. It lives in pipeline/, one
# level up from diagnostics/, so add that directory to the path explicitly. The
# import then resolves from any working directory, and a missing paths.py fails
# loudly instead of silently falling back to a stale CA constant.
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
from paths import ELEV_CLIPPED, STREAM_SOURCE_AREA_M2

DEFAULT_DEM = str(ELEV_CLIPPED)
# Whether DHSVM_STREAM_SOURCE_AREA_M2 was set by the user (a basin-specific
# pinned value) versus left at the paths.py default, which is CA's value. The
# default is meaningful only for CA, so for any other basin we do not compare
# the objective against it. That stale-CA comparison confused new users.
ENV_AREA_SET = "DHSVM_STREAM_SOURCE_AREA_M2" in os.environ


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


def first_sustained_band(results, min_band):
    """Start of the first sustained passing band in the sweep.

    A sustained band is a run of at least min_band consecutive passing
    thresholds (consecutive entries in the swept sequence, so it tracks the
    --step spacing). Returns (start_result, band_results); the objective is the
    band's smallest threshold. Returns (None, []) if no run reaches min_band.

    This rejects isolated single-threshold passes, and any run shorter than
    min_band, which are noise on small basins where the t-test is jumpy. The
    plain smallest-passing rule would take those noise dips. With min_band=1
    this reduces to the smallest passing threshold (the old behavior).
    """
    n = len(results)
    i = 0
    while i < n:
        if not results[i]["passes"]:
            i += 1
            continue
        j = i
        while j < n and results[j]["passes"]:
            j += 1
        run = results[i:j]                       # maximal consecutive pass run
        if len(run) >= min_band:
            return run[0], run
        i = j
    return None, []


def main():
    ap = argparse.ArgumentParser(description="Constant stream drop analysis")
    ap.add_argument("dem", nargs="?", default=DEFAULT_DEM)
    ap.add_argument("--tmin", type=int, default=10, help="min threshold (cells)")
    ap.add_argument("--tmax", type=int, default=300, help="max threshold (cells)")
    ap.add_argument("--step", type=int, default=10, help="threshold step (cells)")
    ap.add_argument("--min-band", type=int, default=3,
                    help="min consecutive passing thresholds to count as a "
                         "sustained band; the objective is the band's smallest "
                         "threshold. Rejects isolated single-threshold passes "
                         "(noise on small basins). Use 2 to reject only single "
                         "points, or 1 for the old smallest-passing behavior")
    ap.add_argument("--csv", default=None, help="write the sweep to this CSV")
    ap.add_argument("--compare-area", type=float, default=None,
                    help="reference support area in m2 to compare the objective "
                         "against, e.g. a previously chosen threshold. Default: "
                         "compare only if DHSVM_STREAM_SOURCE_AREA_M2 is set in "
                         "the environment, otherwise just report the objective")
    ap.add_argument("--emit-area", action="store_true",
                    help="print only the objective support area in m2 to stdout "
                         "(machine readable, for DHSVM_STREAM_SOURCE_AREA_M2); "
                         "suppresses the table and exits nonzero if none pass")
    args = ap.parse_args()
    emit = args.emit_area

    # Optional reference area for the "objective vs reference" comparison:
    # --compare-area wins; else the configured value only if the user set it in
    # the env (basin-specific); else None, so a new basin gets no comparison
    # against CA's default value.
    if args.compare_area is not None:
        compare = args.compare_area
    elif ENV_AREA_SET:
        compare = float(STREAM_SOURCE_AREA_M2)
    else:
        compare = None

    if not emit:
        print("=======================================================")
        print("  CONSTANT STREAM DROP ANALYSIS")
        print("=======================================================")
        print(f"  DEM   : {args.dem}")

    elev, nodata, transform, cell_area = load_dem(args.dem)
    valid = np.isfinite(elev) & (elev != nodata)
    basin_cells = int(valid.sum())
    if not emit:
        print(f"  cell area {cell_area:.2f} m2   valid cells {basin_cells}")
        if compare is not None:
            print(f"  reference threshold: {compare:.1f} m2 "
                  f"= {compare/cell_area:.1f} cells = {compare/1e6:.5f} km2")

    flw = build_flow(elev, nodata, transform)
    uparea_cells = np.asarray(flw.upstream_area(unit="cell"))

    if not emit:
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
        if not emit:
            t_str = f"{r['t_abs']:.2f}" if not np.isnan(r['t_abs']) else "  nan"
            print(f"  {r['cells']:>6} {r['area_km2']:>9.5f} {r['n_streams']:>6} "
                  f"{r['n_order1']:>4} {r['n_higher']:>4} {r['max_order']:>4} "
                  f"{r['dd']:>6.2f} {t_str:>7} {r['n_neg']:>4} "
                  f"{'yes' if r['passes'] else 'no':>5}")

    # Objective threshold: start of the first sustained passing band, the first
    # run of at least --min-band consecutive thresholds with |t|<2. This rejects
    # isolated single-threshold passes, noise on small basins where the t-test
    # is jumpy. Computed once, used by both --emit-area and the human summary.
    passing = [r for r in results if r["passes"]]
    obj, band = first_sustained_band(results, args.min_band)
    # passing thresholds below the chosen band (runs shorter than --min-band),
    # reported for transparency so the user sees what was rejected and why.
    below = [r for r in passing
             if obj is not None and r["cells"] < obj["cells"]]

    if emit:
        if obj is None:
            print(f"drop_analysis: no sustained band (>= {args.min_band} "
                  f"consecutive |t|<2) in {args.tmin}..{args.tmax} cells",
                  file=sys.stderr)
            sys.exit(1)
        print(f"{obj['area_m2']:.1f}")
        return

    print("\n=======================================================")
    if obj is not None:
        band_lo, band_hi = band[0]["cells"], band[-1]["cells"]
        print(f"  objective threshold (start of first sustained band, "
              f">= {args.min_band} consecutive |t|<2):")
        print(f"    {obj['cells']} cells = {obj['area_m2']:.1f} m2 "
              f"= {obj['area_km2']:.5f} km2   (|t|={obj['t_abs']:.2f})")
        print(f"    band: {band_lo} to {band_hi} cells "
              f"({len(band)} consecutive thresholds)")
        if below:
            sk = ", ".join(str(r["cells"]) for r in below)
            print(f"    passes below the band (runs < {args.min_band}, "
                  f"rejected as noise): {sk} cells")
        if compare is not None:
            print(f"  reference: {compare/cell_area:.0f} cells "
                  f"= {compare/1e6:.5f} km2")
            ratio = obj["area_m2"] / compare
            print(f"  objective / reference = {ratio:.2f}x")
            if ratio > 1.3:
                print(f"  -> reference picks a denser network than the objective")
            elif ratio < 0.77:
                print(f"  -> reference picks a sparser network than the objective")
            else:
                print(f"  -> objective is within ~30% of the reference value")
        else:
            print(f"  set DHSVM_STREAM_SOURCE_AREA_M2={obj['area_m2']:.1f} "
                  f"to use this threshold for the pipeline.")
    else:
        print(f"  no sustained band of >= {args.min_band} consecutive passes "
              f"(|t|<2) in the swept range.")
        if passing:
            iso = ", ".join(str(r["cells"]) for r in passing)
            print(f"  only short/isolated passes at: {iso} cells "
                  f"(rejected; lower --min-band to accept them).")
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
