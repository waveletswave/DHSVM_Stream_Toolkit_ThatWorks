#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
soildepth_sensitivity.py  --  B-Q2 desk sensitivity (READ-ONLY)

Reimplements the EXACT soil-depth formula from qgis_CA/soildepthscript.py
and sweeps the B-Q2 candidate parameters (MAX_SLOPE, POW_SLOPE, WT_SOURCE,
MAX_SOURCE) over the real CA rasters.

It only READS rasters and PRINTS / writes a CSV of diagnostics. It does NOT
write any soildepth.bin / .tif and does NOT run DHSVM. This is the
static-analysis step; no pipeline state changes until you decide which
values to commit.

Real data:
    python soildepth_sensitivity.py \
        --dem   <WS>/Reprojected_DEM/elev_clipped.tif \
        --slope <WS>/stream_slope.tif \
        --fac   <WS>/flow_acc.tif \
        --csv   bq2_sensitivity_CA.csv

  (use whatever the soildepthscript.py path candidates resolve to on your tree;
   <WS> is the parent dir of qgis_CA/, i.e. where stream_slope.tif lives.)

Preview the output format without real data (synthetic CA-like rasters):
    python soildepth_sensitivity.py --demo

Options:
    (no args)                          auto-resolve rasters the same way
                                       qgis_CA/soildepthscript.py does (WS =
                                       script's parent-of-parent dir); explicit
                                       --dem/--slope/--fac override it
    --verify <soildepth.tif>           compare the baseline reconstruction
                                       against an existing soildepth.tif
                                       cell-by-cell (regression / faithfulness
                                       check); MATCH => sweep rows trustworthy
    --fixgap {none|zero|exclude}       gap-fill artifact handling for the sweep.
                                       none=faithful (current pipeline); zero=
                                       slope->0 cells get no slope bonus; exclude
                                       =drop them. Use to see the de-artifacted
                                       distribution (verify auto-skips here).
    --maxsource {max|p99|p95|<float>}  fac normalization used for WT_SOURCE>0
                                       rows (default: max = raw |max| of fac)
    --renorm                           when WT_SOURCE>0, rescale slope+elev so
                                       all three weights sum to 1 (default OFF
                                       = faithful additive, matches current code)
"""

import argparse
import sys
from pathlib import Path
import numpy as np

# ----- fixed parameters (match qgis_CA/soildepthscript.py; NOT swept) --------
MIN_DEPTH  = 2.0
MAX_DEPTH  = 6.0
WT_SLOPE_0 = 0.7        # baseline slope weight
WT_ELEV_0  = 0.3        # baseline elev weight
MAX_ELEV   = 1500.0
POW_ELEV   = 0.75
POW_SOURCE = 1.0
DHSVM_NODATA = -9999.0

# ----- candidate sweep values (from v2 handoff B-Q2 table) -------------------
SWEEP_MAX_SLOPE = [30.0, 45.0, 60.0]      # 30  = current
SWEEP_POW_SLOPE = [0.25, 0.5, 1.0]        # 0.25 = current
SWEEP_WT_SOURCE = [0.0, 0.1, 0.2, 0.3]    # 0.0  = current


# =============================================================================
# Core formula  (faithful to qgis_CA/soildepthscript.py, step 3-4)
# =============================================================================
def gapfill(elev, slope, fac, elev_nd, slope_nd, fac_nd):
    """Replicate the exact masking/gap-fill of soildepthscript.py."""
    valid = (elev != elev_nd) & (~np.isnan(elev))
    elev_clean  = np.where(valid & (elev >= 0), elev, 0.0)
    slope_clean = np.where(valid & (slope != slope_nd) & (~np.isnan(slope)) & (slope >= 0), slope, 0.0)
    fac_clean   = np.where(valid & (fac != fac_nd) & (~np.isnan(fac)), np.abs(fac), 0.0)
    # cells that are valid in the DEM but had NO valid slope -> gap-filled to 0
    # (these become spuriously DEEP because term_slope peaks at slope=0; cf. 3.3)
    slope_gapfilled = valid & ~((slope != slope_nd) & (~np.isnan(slope)) & (slope >= 0))
    return valid, elev_clean, slope_clean, fac_clean, slope_gapfilled


def compute_depth(elev_clean, slope_clean, fac_clean, valid,
                  MAX_SLOPE, POW_SLOPE, WT_SOURCE, MAX_SOURCE, renorm=False,
                  noslope=None):
    if renorm and WT_SOURCE > 0:
        scale = (1.0 - WT_SOURCE) / (WT_SLOPE_0 + WT_ELEV_0)
        WT_SLOPE, WT_ELEV = WT_SLOPE_0 * scale, WT_ELEV_0 * scale
    else:
        WT_SLOPE, WT_ELEV = WT_SLOPE_0, WT_ELEV_0

    elev_lim  = np.clip(elev_clean,  0.0, MAX_ELEV)
    slope_lim = np.clip(slope_clean, 0.0, MAX_SLOPE)
    fac_lim   = np.clip(fac_clean,   0.0, MAX_SOURCE)

    term_slope = WT_SLOPE * (1.0 - (slope_lim / MAX_SLOPE) ** POW_SLOPE)
    if noslope is not None:
        # gap-fill fix (zero mode): cells with no real slope get no slope bonus,
        # instead of the spurious MAX bonus they get when slope is filled to 0.
        term_slope = np.where(noslope, 0.0, term_slope)
    if WT_SOURCE > 0:
        term_source = WT_SOURCE * ((fac_lim / MAX_SOURCE) ** POW_SOURCE)
    else:
        term_source = np.zeros_like(fac_lim)
    term_elev = WT_ELEV * (1.0 - (elev_lim / MAX_ELEV) ** POW_ELEV)

    depth_calc = MIN_DEPTH + (MAX_DEPTH - MIN_DEPTH) * (term_slope + term_source + term_elev)
    depth = np.full(elev_clean.shape, np.nan, dtype=np.float64)
    depth[valid] = np.clip(depth_calc[valid], MIN_DEPTH, MAX_DEPTH)
    return depth


# =============================================================================
# Diagnostics
# =============================================================================
def diag(depth, slope_clean, fac_clean, valid, label):
    d = depth[valid]; s = slope_clean[valid]; f = fac_clean[valid]
    q1, q3 = np.percentile(s, 25), np.percentile(s, 75)
    f90 = np.percentile(f, 90)
    return {
        'label': label,
        'mean': d.mean(), 'std': d.std(),
        'min': d.min(), 'max': d.max(),
        'p5':  np.percentile(d, 5),  'p25': np.percentile(d, 25),
        'p50': np.percentile(d, 50), 'p75': np.percentile(d, 75),
        'p95': np.percentile(d, 95),
        'pct_at_min': 100.0 * np.mean(d <= MIN_DEPTH + 1e-6),
        'pct_at_max': 100.0 * np.mean(d >= MAX_DEPTH - 1e-6),
        'd_gentle': d[s <= q1].mean(),    # bottom slope quartile (valley-like)
        'd_steep':  d[s >= q3].mean(),    # top slope quartile (ridge-like)
        'contrast': d[s <= q1].mean() - d[s >= q3].mean(),
        'd_topfac': d[f >= f90].mean(),   # channel-proximal (top fac decile)
    }


def print_table(rows):
    hdr = (f"{'config':<34}{'mean':>6}{'std':>6}{'p5':>6}{'p50':>6}{'p95':>6}"
           f"{'%@2m':>7}{'%@6m':>7}{'gentle':>8}{'steep':>7}{'contr':>7}{'topfac':>8}")
    print(hdr); print('-' * len(hdr))
    for r in rows:
        print(f"{r['label']:<34}{r['mean']:>6.2f}{r['std']:>6.2f}{r['p5']:>6.2f}"
              f"{r['p50']:>6.2f}{r['p95']:>6.2f}{r['pct_at_min']:>7.1f}{r['pct_at_max']:>7.1f}"
              f"{r['d_gentle']:>8.2f}{r['d_steep']:>7.2f}{r['contrast']:>7.2f}{r['d_topfac']:>8.2f}")


def write_csv(rows, path):
    cols = ['label','mean','std','min','max','p5','p25','p50','p75','p95',
            'pct_at_min','pct_at_max','d_gentle','d_steep','contrast','d_topfac']
    with open(path, 'w') as fh:
        fh.write(','.join(cols) + '\n')
        for r in rows:
            fh.write(','.join(str(r[c]) if c == 'label' else f"{r[c]:.4f}" for c in cols) + '\n')
    print(f"\n[csv] wrote {len(rows)} rows -> {path}")


# =============================================================================
# Input loading
# =============================================================================
def _read_band(fp):
    """Read band-1 array + nodata. Use osgeo.gdal if importable, else rasterio.
    The numbers are identical either way; only the binding differs."""
    try:
        from osgeo import gdal
    except ModuleNotFoundError:
        gdal = None
    if gdal is not None:
        ds = gdal.Open(fp)
        if ds is None:
            sys.exit(f"[error] could not open raster: {fp}")
        b = ds.GetRasterBand(1)
        return b.ReadAsArray().astype(np.float64), b.GetNoDataValue()
    try:
        import rasterio
    except ModuleNotFoundError:
        sys.exit("[error] neither osgeo (gdal) nor rasterio is importable here. "
                 "Activate the env that has them and rerun, e.g.:\n"
                 "    conda activate dhsvm_rs")
    with rasterio.open(fp) as ds:
        arr = ds.read(1).astype(np.float64)
        nd = ds.nodata
    return arr, nd


def load_real(dem, slope, fac):
    e, e_nd = _read_band(dem); s, s_nd = _read_band(slope); f, f_nd = _read_band(fac)
    if e_nd is None: e_nd = DHSVM_NODATA
    if s_nd is None: s_nd = DHSVM_NODATA
    if f_nd is None: f_nd = DHSVM_NODATA
    return e, s, f, e_nd, s_nd, f_nd


def make_demo(n=140, seed=1):
    """Synthetic CA-like rasters for format preview + self-validation.
    Slope mean ~26 deg / max ~45; heavy-tailed fac with a central channel;
    elev 650-1150 m. A 2-px ring is valid-in-DEM but NoData-in-slope so the
    boundary gap-fill artifact (3.3) is exercised."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n]
    elev = np.clip(700 + 400 * (yy / n) + 30 * rng.standard_normal((n, n)), 650, 1150)
    slope = 45.0 * rng.beta(5.0, 3.0, size=(n, n))           # mean ~28, max ~45
    fac = rng.gamma(0.6, 8.0, size=(n, n))
    col = np.clip((n // 2 + 8 * np.sin(yy[:, 0] / 8)).astype(int), 1, n - 2)
    ramp = np.linspace(50, 3388, n)
    for r in range(n):
        c = col[r]
        fac[r, c-1:c+2] += ramp[r] / 3.0
        slope[r, c-1:c+2] *= 0.3
        elev[r, c-1:c+2] -= 40
    fac = np.clip(fac, 0, 3388)
    nd = DHSVM_NODATA
    # DEM NoData border = 3 px; slope NoData border = 5 px -> 2-px gap-fill ring
    elev_nd = elev.copy()
    eb = np.zeros((n, n), bool); eb[:3] = eb[-3:] = eb[:, :3] = eb[:, -3:] = True
    elev_nd[eb] = nd
    slope_nd_arr = slope.copy()
    sb = np.zeros((n, n), bool); sb[:5] = sb[-5:] = sb[:, :5] = sb[:, -5:] = True
    slope_nd_arr[sb] = nd
    return elev_nd, slope_nd_arr, fac, nd, nd, nd


# =============================================================================
# Auto-resolve raster paths the SAME way qgis_CA/soildepthscript.py does:
#   WS = SCRIPT_DIR.parent  (so a script in <WS>/scripts/ resolves WS=<WS>),
#   same candidate filenames. Lets a zero-arg run read exactly the rasters
#   soildepthscript.py reads. Explicit --dem/--slope/--fac still override.
# =============================================================================
def autoresolve():
    WS = Path(__file__).resolve().parent.parent
    dem_c   = [WS / "Reprojected_DEM" / "elev_clipped.tif",
               WS / "Reprojected_DEM" / "Cropped_DEM.tif",
               WS / "dem.tif"]
    slope_c = [WS / "stream_slope.tif"]
    fac_c   = [WS / "flow_acc.tif"]
    dem   = next((str(c) for c in dem_c   if c.exists()), None)
    slope = next((str(c) for c in slope_c if c.exists()), None)
    fac   = next((str(c) for c in fac_c   if c.exists()), None)
    return WS, dem, slope, fac, dem_c, slope_c, fac_c


# =============================================================================
# Main
# =============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dem'); ap.add_argument('--slope'); ap.add_argument('--fac')
    ap.add_argument('--demo', action='store_true')
    ap.add_argument('--csv', default=None)
    ap.add_argument('--maxsource', default='max',
                    help="fac normalization for WT_SOURCE>0 rows: max|p99|p95|<float>")
    ap.add_argument('--renorm', action='store_true')
    ap.add_argument('--verify', default=None,
                    help="path to existing soildepth.tif; compares baseline "
                         "reconstruction against it cell-by-cell (regression check)")
    ap.add_argument('--fixgap', choices=['none', 'zero', 'exclude'], default='none',
                    help="gap-fill artifact handling for the sweep. none=faithful "
                         "(current pipeline, slope->0 cells get spurious MAX depth); "
                         "zero=those cells get no slope bonus; exclude=drop them from "
                         "the domain. Use to see the de-artifacted distribution.")
    a = ap.parse_args()

    if a.demo:
        print(">>> DEMO MODE: synthetic CA-like rasters (NOT your real data)\n")
        e, s, f, e_nd, s_nd, f_nd = make_demo()
    elif a.dem and a.slope and a.fac:
        e, s, f, e_nd, s_nd, f_nd = load_real(a.dem, a.slope, a.fac)
    else:
        WS, dem, slope, fac, dem_c, slope_c, fac_c = autoresolve()
        if dem and slope and fac:
            print(f">>> AUTO-RESOLVED rasters (WS = {WS})")
            print(f"    dem  : {dem}")
            print(f"    slope: {slope}")
            print(f"    fac  : {fac}\n")
            e, s, f, e_nd, s_nd, f_nd = load_real(dem, slope, fac)
        else:
            msg = ["[usage] no --dem/--slope/--fac given and auto-resolve failed.",
                   f"  WS (script's parent-of-parent dir) = {WS}",
                   "  looked for:"]
            for c in dem_c + slope_c + fac_c:
                msg.append(f"    {'FOUND  ' if c.exists() else 'missing'}  {c}")
            msg.append("  fix: pass the paths explicitly with --dem/--slope/--fac,")
            msg.append("       or run --demo for a synthetic preview.")
            sys.exit("\n".join(msg))

    valid, e_c, s_c, f_c, s_gap = gapfill(e, s, f, e_nd, s_nd, f_nd)

    # --- gap-fill artifact handling (affects the sweep / stats / 3.3) --------
    if a.fixgap == 'exclude':
        valid_eff = valid & ~s_gap
        noslope = None
        print(f">>> GAP-FILL FIX = exclude: dropping {int(s_gap.sum())} "
              f"slope-gapfilled cells from the domain.")
        if a.verify:
            print("    (--verify skipped: the existing soildepth.tif still "
                  "contains the artifact, so it won't match the fixed field.)")
            a.verify = None
        print()
    elif a.fixgap == 'zero':
        valid_eff = valid
        noslope = s_gap
        print(f">>> GAP-FILL FIX = zero: {int(s_gap.sum())} slope-gapfilled "
              f"cells get no slope bonus (term_slope=0).")
        if a.verify:
            print("    (--verify skipped: the existing soildepth.tif still "
                  "contains the artifact, so it won't match the fixed field.)")
            a.verify = None
        print()
    else:
        valid_eff = valid
        noslope = None

    nvalid = int(valid_eff.sum())

    # --- input context -------------------------------------------------------
    sv = s_c[valid_eff]; fv = f_c[valid_eff]
    print(f"Valid cells: {nvalid}")
    print(f"slope (deg) : mean {sv.mean():.2f}  max {sv.max():.2f}  "
          f">30deg {100*np.mean(sv>30):.1f}%  >45deg {100*np.mean(sv>45):.1f}%")
    print(f"flow_acc    : max {fv.max():.0f}  p99 {np.percentile(fv,99):.0f}  "
          f"p95 {np.percentile(fv,95):.0f}  p90 {np.percentile(fv,90):.0f}  median {np.median(fv):.0f}")
    print(f"slope-gapfilled valid cells (3.3 artifact source): "
          f"{int((s_gap).sum())} ({100*np.mean(s_gap[valid]):.2f}% of valid)\n")

    # resolve MAX_SOURCE for WT_SOURCE>0 rows
    if a.maxsource == 'max':   max_src = float(fv.max())
    elif a.maxsource == 'p99': max_src = float(np.percentile(fv, 99))
    elif a.maxsource == 'p95': max_src = float(np.percentile(fv, 95))
    else:                      max_src = float(a.maxsource)
    print(f"MAX_SOURCE for source-term rows = {max_src:.0f}  "
          f"(renorm weights = {a.renorm})\n")

    # --- baseline (current pipeline values) ----------------------------------
    base = compute_depth(e_c, s_c, f_c, valid_eff, 30.0, 0.25, 0.0, 100000.0,
                         noslope=noslope)
    blabel = "BASELINE MS30 POW0.25 WS0" + ("" if a.fixgap == 'none'
                                            else f" [fixgap={a.fixgap}]")
    rows = [diag(base, s_c, f_c, valid_eff, blabel)]

    # --- optional faithfulness check vs existing soildepth.tif ---------------
    if a.verify:
        ref, ref_nd = _read_band(a.verify)
        if ref.shape != base.shape:
            print(f"[verify] SHAPE MISMATCH: soildepth.tif {ref.shape} vs "
                  f"computed {base.shape}. Likely a different DEM/grid than the "
                  f"slope/fac you passed. Check that all three inputs are the "
                  f"grid-matched Intermediate_GIS rasters.\n")
        else:
            if ref_nd is None: ref_nd = DHSVM_NODATA
            both = valid_eff & (ref != ref_nd) & (~np.isnan(ref))
            d = np.abs(base[both] - ref[both])
            r = np.corrcoef(base[both], ref[both])[0, 1] if both.sum() > 1 else float('nan')
            n_off = int(np.sum(d > 1e-3))
            print(f"[verify] vs {a.verify}")
            print(f"         compared {int(both.sum())} cells | max|diff| {d.max():.6g} | "
                  f"mean|diff| {d.mean():.6g} | corr {r:.6f} | cells>1e-3: {n_off}")
            if d.max() < 1e-3:
                print("         -> MATCH: baseline reconstruction reproduces "
                      "soildepth.tif. Sweep rows are trustworthy.\n")
            else:
                print("         -> DIVERGENCE: reconstruction does not match. "
                      "Either wrong rasters, or the formula/params on your tree "
                      "differ from qgis_CA/soildepthscript.py. Resolve before "
                      "trusting the sweep.\n")

    # --- full sweep ----------------------------------------------------------
    for MS in SWEEP_MAX_SLOPE:
        for PW in SWEEP_POW_SLOPE:
            for WS in SWEEP_WT_SOURCE:
                ms_src = max_src if WS > 0 else 100000.0
                d = compute_depth(e_c, s_c, f_c, valid_eff, MS, PW, WS, ms_src,
                                  renorm=a.renorm, noslope=noslope)
                lbl = f"MS{MS:.0f} POW{PW:g} WS{WS:g}"
                rows.append(diag(d, s_c, f_c, valid_eff, lbl))

    print("SWEEP (depth in m). gentle/steep = bottom/top slope quartile; "
          "contr = gentle-steep; topfac = top fac decile\n")
    print_table(rows)

    # --- MAX_SOURCE leverage at a fixed slope config -------------------------
    print("\nMAX_SOURCE leverage (at MS45 POW1.0 WS0.2), shows how fac "
          "normalization changes valley deepening:")
    msrc_rows = []
    for tag, ms in [('100000(current)', 100000.0), ('max', float(fv.max())),
                    ('p99', float(np.percentile(fv,99))), ('p95', float(np.percentile(fv,95)))]:
        d = compute_depth(e_c, s_c, f_c, valid_eff, 45.0, 1.0, 0.2, ms,
                          renorm=a.renorm, noslope=noslope)
        msrc_rows.append(diag(d, s_c, f_c, valid_eff, f"MAXSRC={tag}"))
    print_table(msrc_rows)

    # --- boundary / gap-fill artifact (3.3) ----------------------------------
    if a.fixgap == 'none' and s_gap.any():
        dg = base[s_gap]
        print(f"\n[3.3 artifact] {s_gap.sum()} valid cells had slope gap-filled -> "
              f"depth mean {dg.mean():.2f} m (vs basin mean {base[valid_eff].mean():.2f} m), "
              f"max {dg.max():.2f} m. These are mask-edge cells, not real coves.")
        topN = 20
        flat = np.where(valid_eff, base, -1).ravel()
        idx = np.argsort(flat)[-topN:]
        n_gap_in_top = int(s_gap.ravel()[idx].sum())
        print(f"[3.3 artifact] of the {topN} DEEPEST valid cells, "
              f"{n_gap_in_top} are slope-gapfilled edge cells.")
    elif a.fixgap != 'none':
        print(f"\n[3.3] gap-fill fix '{a.fixgap}' applied to {int(s_gap.sum())} "
              f"cells; the artifact tail above is removed from the rows shown.")

    if a.csv:
        write_csv(rows, a.csv)


if __name__ == '__main__':
    main()