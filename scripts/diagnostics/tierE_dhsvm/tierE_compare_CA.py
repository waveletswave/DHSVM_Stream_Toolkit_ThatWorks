# -*- coding: utf-8 -*-
# =====================================================================
# tierE_compare_CA.py  -  manuscript runs, control reruns and Tier E
#                         reruns of the CA configurations, compared
#
# For each run (CA_S4h, CA_LAI70, CA_LAI20) three output sets are read
# from TestCase/CA/output, each a DHSVM prefix:
#   manuscript  CA_S4h, CA_LAI70, CA_LAI20 (April 2026; for CA_LAI20 the
#               May 13 set 512_LAI20 is checked as well)
#   control     old_<tag>    old network, run today (tag S4h, LAI70,
#               LAI20; the prefixes are short because this DHSVM build
#               accepts an Output Directory of at most 78 characters)
#   Tier E      new_<tag>    new network, run today
#   --base apr  oA_<tag> and nA_<tag>: the same pair on the April
#               inputs (DEM_CA_apr); oA_* should be identical to the
#               manuscript outputs, nA_* is the strict R11 rerun
#
# Reported per run:
#   1. the Mass.Final.Balance totals of every set side by side
#   2. reproducibility: manuscript versus control. Identical output at
#      the printed precision means today's inputs reproduce the April
#      runs; then manuscript versus Tier E is a clean network comparison
#   3. network effect: control versus Tier E, column by column for
#      Aggregated.Values and Mass.Balance (sums, max single-step
#      difference, that difference relative to the column's own peak),
#      and the routed basin outflow from the Stream.Flow "Totals" rows
#      (total_outflow over every segment without a downstream segment,
#      which exists for the old runs too, which had no SAVE rows)
#   4. R12: outlet_1 + outlet_2 in Streamflow.Only (Tier E) must equal
#      the Stream.Flow total outflow step by step
#   5. cumulative ChannelInt versus cumulative routed outflow
#
# Usage:  python3 tierE_compare_CA.py [--runs CA_S4h CA_LAI70 CA_LAI20]
#         python3 tierE_compare_CA.py --base apr
#         python3 tierE_compare_CA.py --smoke     (the nine-day sets
#             oldsm_S4h and newsm_S4h; the manuscript CA_S4h is
#             compared over the same nine days)
#         python3 tierE_compare_CA.py --adhoc CA_LAI20 512_LAI20
#             (any two prefixes)
# Writes output/tierE_compare_<run>_<pair>.csv (Aggregated.Values).
# =====================================================================

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

CASE = Path(os.environ.get(
    "TIERE_CASE",
    "/Users/benthosyy/Desktop/CodeBits/DHSVM-PNNL-2025/TestCase/CA"))
RUNS = ["CA_S4h", "CA_LAI70", "CA_LAI20"]
MANUSCRIPT = {"CA_S4h": ["CA_S4h"],
              "CA_LAI70": ["CA_LAI70"],
              "CA_LAI20": ["CA_LAI20", "512_LAI20"]}
TAGS = {"CA_S4h": "S4h", "CA_LAI70": "LAI70", "CA_LAI20": "LAI20"}
WORDS = {"jun": {"ctrl": "old", "tierE": "new", "ww": "ww"},
         "apr": {"ctrl": "oA", "tierE": "nA", "ww": "wA"}}


def prefix_for(kind, run, smoke, base="jun"):
    """Same rule as tierE_rerun_CA.py: old_S4h, new_LAI70, oldsm_S4h,
    oA_S4h, nA_LAI20"""
    return f"{WORDS[base][kind]}{'sm' if smoke else ''}_{TAGS[run]}"


CELL_M = 28.1577401
N_CELLS = 4334
BASIN_AREA_M2 = N_CELLS * CELL_M * CELL_M       # 3.436e6 m2, as in the log
MM_PER_M3 = 1000.0 / BASIN_AREA_M2
FINAL_KEYS = ["Initial Storage", "Precip/Inflow", "ET", "ChannelInt",
              "Final Storage", "Mass Error"]
TOL = 1e-9          # relative to the column peak: "identical"


def out_path(prefix, name):
    return CASE / "output" / (prefix + name)


def parse_dates(s):
    s = s.astype(str).str.strip().str.replace(".", "/", regex=False)
    return pd.to_datetime(s, format="%m/%d/%Y-%H:%M:%S")


def read_table(path):
    """Aggregated.Values / Mass.Balance: header line, whitespace fields."""
    with open(path) as f:
        header = f.readline()
    cols = [c for c in re.split(r"\s+", header.strip()) if c]
    df = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None)
    if df.shape[1] > len(cols):
        df = df.iloc[:, :len(cols)]
    elif df.shape[1] < len(cols):
        cols = cols[:df.shape[1]]
    df.columns = cols
    df.index = parse_dates(df[cols[0]])
    df = df.drop(columns=[cols[0]]).apply(pd.to_numeric, errors="coerce")
    return df[~df.index.duplicated(keep="first")]


def read_stream_totals(path):
    """Stream.Flow rows ending in "Totals" (channel.c: date, 0,
    total_lateral_inflow, total_outflow, total_storage,
    total_storage_change, total_error; m3 per timestep)."""
    rows = []
    with open(path) as f:
        for line in f:
            if not line.rstrip().endswith('"Totals"'):
                continue
            p = line.split()
            rows.append((p[0], float(p[2]), float(p[3]), float(p[4]),
                         float(p[5]), float(p[6])))
    df = pd.DataFrame(rows, columns=["date", "lateral_in", "outflow",
                                     "storage", "dstorage", "error"])
    df.index = parse_dates(df["date"])
    return df.drop(columns=["date"])


def read_streamflow_only(path):
    with open(path) as f:
        header = f.readline()
    cols = [c for c in re.split(r"\s+", header.strip()) if c]
    if len(cols) < 2:
        return None
    df = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None, names=cols)
    df.index = parse_dates(df[cols[0]])
    df = df.drop(columns=[cols[0]]).apply(pd.to_numeric, errors="coerce")
    return df.dropna(how="all")


def read_final(prefix):
    p = out_path(prefix, "Mass.Final.Balance")
    vals = {}
    if not p.exists():
        return vals
    for ln in p.read_text().splitlines():
        for key in FINAL_KEYS:
            if re.search(rf"\b{re.escape(key)} ", ln):
                vals[key] = float(ln.split()[-1])
    return vals


def compare_tables(old, new, label, csv=None):
    common = [c for c in old.columns if c in new.columns]
    idx = old.index.intersection(new.index)
    o, n = old.loc[idx, common], new.loc[idx, common]
    rows = []
    for c in common:
        d = (n[c] - o[c]).abs()
        peak = float(o[c].abs().max())
        rows.append(dict(column=c, sum_old=float(o[c].sum()),
                         sum_new=float(n[c].sum()),
                         max_step_diff=float(d.max()),
                         rel_to_peak=(float(d.max()) / peak if peak > 0
                                      else 0.0),
                         when=(d.idxmax() if d.max() > 0 else "")))
    tab = pd.DataFrame(rows).sort_values("rel_to_peak", ascending=False)
    changed = tab[tab["rel_to_peak"] > TOL]
    print(f"  {label}: {len(idx)} common timesteps (old {len(old)}, new "
          f"{len(new)} rows), {len(common)} columns; "
          f"{len(changed)} differ (rel > {TOL:g})"
          + ("" if len(changed) else ": IDENTICAL at printed precision"))
    if len(changed):
        with pd.option_context("display.width", 140, "display.max_rows", 60,
                               "display.float_format", "{:.6g}".format):
            print(changed.head(15).to_string(index=False))
    if csv is not None:
        tab.to_csv(csv, index=False)
        print(f"  table written: {csv}")
    return tab


def compare_pair(a, b, label, csv_stem=None):
    """Full comparison of two prefixes; returns True if everything read
    is identical at printed precision."""
    print("-" * 72)
    print(f"{label}: {a}* versus {b}*")
    same = True
    for name in ("Aggregated.Values", "Mass.Balance"):
        pa, pb = out_path(a, name), out_path(b, name)
        if not (pa.exists() and pb.exists()):
            print(f"  {name}: missing ({pa.exists()}, {pb.exists()})")
            same = False
            continue
        csv = (CASE / "output" / f"{csv_stem}.csv"
               if csv_stem and name == "Aggregated.Values" else None)
        tab = compare_tables(read_table(pa), read_table(pb), name, csv)
        same = same and not (tab["rel_to_peak"] > TOL).any()
    pa, pb = out_path(a, "Stream.Flow"), out_path(b, "Stream.Flow")
    if pa.exists() and pb.exists():
        ta, tb = read_stream_totals(pa), read_stream_totals(pb)
        idx = ta.index.intersection(tb.index)
        ta, tb = ta.loc[idx], tb.loc[idx]
        qa, qb = ta["outflow"], tb["outflow"]
        d = qb - qa
        print(f"  Stream.Flow totals: {len(idx)} common timesteps")
        print(f"    routed outflow {qa.sum() * MM_PER_M3:.3f} mm versus "
              f"{qb.sum() * MM_PER_M3:.3f} mm, difference "
              f"{d.sum() * MM_PER_M3:+.3f} mm "
              f"({100 * d.sum() / qa.sum():+.4f} %)"
              if qa.sum() > 0 else "    routed outflow: zero in the first")
        print(f"    max |step difference| {d.abs().max():.4g} m3/step "
              f"({100 * d.abs().max() / qa.max():.3f} % of peak "
              f"{qa.max():.4g}), RMSE {np.sqrt((d ** 2).mean()):.4g} "
              f"m3/step, r {qa.corr(qb):.6f}"
              if qa.max() > 0 else "")
        print(f"    channel storage at the end {ta['storage'].iloc[-1]:.4g}"
              f" versus {tb['storage'].iloc[-1]:.4g} m3; max |error| "
              f"{ta['error'].abs().max():.3g} versus "
              f"{tb['error'].abs().max():.3g} m3/step")
        # the manuscript's simulated Q: total_lateral_inflow (column 3 of
        # the Totals row, DHSVM_CA_Stitch_Optimal_LAI_V2.py parts[2]),
        # summed to daily mm
        la = ta["lateral_in"].resample("1D").sum() * MM_PER_M3
        lb = tb["lateral_in"].resample("1D").sum() * MM_PER_M3
        dl = lb - la
        print(f"    channel lateral inflow (manuscript Q): "
              f"{la.sum():.3f} versus {lb.sum():.3f} mm, difference "
              f"{dl.sum():+.3f} mm ({100 * dl.sum() / la.sum():+.4f} %); "
              f"daily max |diff| {dl.abs().max():.4g} mm/d, RMSE "
              f"{np.sqrt((dl ** 2).mean()):.4g} mm/d, r {la.corr(lb):.6f}"
              if la.sum() > 0 else
              "    channel lateral inflow: zero in the first")
        same = same and d.abs().max() == 0 and dl.abs().max() == 0
    else:
        print(f"  Stream.Flow: missing ({pa.exists()}, {pb.exists()})")
        same = False
    return same


def check_r12(prefix):
    so = out_path(prefix, "Streamflow.Only")
    sf = out_path(prefix, "Stream.Flow")
    if not (so.exists() and sf.exists()):
        print(f"  R12 {prefix}: files missing")
        return
    only = read_streamflow_only(so)
    if only is None:
        print(f"  R12 {prefix}: Streamflow.Only has no SAVE columns")
        return
    tot = read_stream_totals(sf)
    idx = only.index.intersection(tot.index)
    dd = (only.loc[idx].sum(axis=1) - tot.loc[idx, "outflow"]).abs()
    peak = tot.loc[idx, "outflow"].max()
    print(f"  R12 {prefix}: columns {list(only.columns)}; "
          f"|sum of SAVE outflows - Stream.Flow total outflow| max "
          f"{dd.max():.4g} m3/step over {len(idx)} steps, "
          f"{dd.max() / peak:.2e} of the peak total {peak:.4g} "
          f"(both files print 5 significant digits, so up to about "
          f"1e-4 of the peak is rounding)")
    for c in only.columns:
        print(f"    {c}: {only[c].sum() * MM_PER_M3:.3f} mm, peak "
              f"{only[c].max():.4g} m3/step, share "
              f"{100 * only[c].sum() / only.sum().sum():.2f} %")


def channelint_vs_routed(prefix):
    mb, sf = out_path(prefix, "Mass.Balance"), out_path(prefix, "Stream.Flow")
    if not (mb.exists() and sf.exists()):
        return
    t = read_table(mb)
    col = [c for c in t.columns if c.startswith("ChannelInt")]
    if not col:
        return
    ci = t[col[0]].sum() * 1000.0
    q = read_stream_totals(sf)["outflow"].sum() * MM_PER_M3
    print(f"  {prefix}: cumulative ChannelInt {ci:.3f} mm, routed outflow "
          f"{q:.3f} mm, ratio {q / ci:.4f}" if ci > 0 else
          f"  {prefix}: cumulative ChannelInt is zero")


def final_table(prefixes):
    rows = []
    for p in prefixes:
        v = read_final(p)
        if v:
            mt = out_path(p, "Mass.Final.Balance").stat().st_mtime
            v = dict(prefix=p, written=pd.Timestamp(mt, unit="s")
                     .strftime("%Y-%m-%d %H:%M"), **v)
            rows.append(v)
    if rows:
        with pd.option_context("display.width", 160,
                               "display.float_format", "{:.3f}".format):
            print(pd.DataFrame(rows).to_string(index=False))


def compare_run(run, smoke, base="jun"):
    ctrl = prefix_for("ctrl", run, smoke, base)
    tier = prefix_for("tierE", run, smoke, base)
    manu = MANUSCRIPT[run]
    print("=" * 72)
    print(f"{run}: manuscript {manu}, control {ctrl}, Tier E {tier}")
    print("=" * 72)
    print("1. Mass.Final.Balance (mm)" + (
        "; the smoke sets cover nine days, so only Initial Storage is "
        "comparable with the manuscript" if smoke else ""))
    final_table(manu + [ctrl, tier])
    print("2. reproducibility: manuscript versus control")
    for m in manu:
        if out_path(m, "Aggregated.Values").exists() and \
                out_path(ctrl, "Aggregated.Values").exists():
            same = compare_pair(m, ctrl, f"manuscript {m} vs control",
                                f"tierE_compare_{run}_{m}_vs_{ctrl}")
            print(f"  => {m} {'IS' if same else 'is NOT'} reproduced by the "
                  f"control run")
        else:
            print(f"  {m} or {ctrl}: Aggregated.Values missing")
    print("3. network effect: control versus Tier E")
    if out_path(ctrl, "Aggregated.Values").exists() and \
            out_path(tier, "Aggregated.Values").exists():
        compare_pair(ctrl, tier, "control vs Tier E",
                     f"tierE_compare_{run}_{ctrl}_vs_{tier}")
    else:
        print(f"  {ctrl} or {tier}: Aggregated.Values missing")
    print("4. R12, SAVE outlets versus the routed total")
    check_r12(tier)
    print("5. cumulative ChannelInt versus routed outflow")
    for p in manu[:1] + [ctrl, tier]:
        channelint_vs_routed(p)
    # the WW-DHSVM network on the same inputs (cross-engine comparison,
    # 2026-09-23), when its outputs exist
    ww = prefix_for("ww", run, smoke, base)
    if out_path(ww, "Aggregated.Values").exists():
        print("6. cross-engine: control versus WW-DHSVM network, and Tier E "
              "versus WW-DHSVM")
        final_table([ctrl, tier, ww])
        compare_pair(ctrl, ww, "control vs WW-DHSVM",
                     f"tierE_compare_{run}_{ctrl}_vs_{ww}")
        if out_path(tier, "Aggregated.Values").exists():
            compare_pair(tier, ww, "Tier E vs WW-DHSVM",
                         f"tierE_compare_{run}_{tier}_vs_{ww}")
        check_r12(ww)
        channelint_vs_routed(ww)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", default=RUNS, choices=RUNS)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--base", default="jun", choices=list(WORDS),
                    help="jun: old_/new_ prefixes; apr: oA_/nA_")
    ap.add_argument("--adhoc", nargs=2, metavar=("PREFIX_A", "PREFIX_B"),
                    help="compare any two output prefixes and stop")
    args = ap.parse_args()
    if args.adhoc:
        a, b = args.adhoc
        final_table([a, b])
        compare_pair(a, b, "adhoc", f"tierE_compare_adhoc_{a}_vs_{b}")
        check_r12(b)
        for p in (a, b):
            channelint_vs_routed(p)
        return
    for run in (["CA_S4h"] if args.smoke else args.runs):
        compare_run(run, args.smoke, args.base)


if __name__ == "__main__":
    main()
