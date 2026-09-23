# -*- coding: utf-8 -*-
# =====================================================================
# tierE_eval_CA.py  -  the manuscript's streamflow metrics, recomputed
#                      with the manuscript's own code, for the reruns
#
# Imports DHSVM_CA_Stitch_Optimal_LAI_V2.py (SpongeBurn/10_CA_Calib) by
# path and uses its parse_streamflow, load_observed, calc_metrics,
# get_file_paths and WATERSHED_AREA_M2 unchanged. For each base (a
# triple of output prefixes for S4h, LAI70 and LAI20) it prints the
# per-year metrics of each run (calendar years, as in
# find_best_run_per_year) and the stitched series (2017 from the LAI20
# run, 2018 from the LAI70 run; windows 2017-02-01 to 2018-12-31, as in
# evaluate_stitched), then the differences between bases.
#
# Bases:
#   manuscript  CA_S4h   CA_LAI70   CA_LAI20    April 2026 outputs
#   old         old_S4h  old_LAI70  old_LAI20   today's inputs, old network
#   new         new_S4h  new_LAI70  new_LAI20   today's inputs, Tier E
#   oA / nA     oA_*     nA_*                   April inputs, old / Tier E
# Missing bases are skipped. Writes output/tierE_eval_CA.csv.
#
# Run with the Python environment used for the manuscript scripts
# (pandas, numpy, scipy, scikit-learn):
#   python3 tierE_eval_CA.py
#   python3 tierE_eval_CA.py --bases manuscript old new
# =====================================================================

import argparse
import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd

CASE = Path(os.environ.get(
    "TIERE_CASE",
    "/Users/benthosyy/Desktop/CodeBits/DHSVM-PNNL-2025/TestCase/CA"))
EVAL_SCRIPT = Path(os.environ.get(
    "TIERE_EVAL_SCRIPT",
    "/Users/benthosyy/Desktop/CodeBits/DHSVM-PNNL-2025/TestCase/LAI/"
    "SpongeBurn/10_CA_Calib/DHSVM_CA_Stitch_Optimal_LAI_V2.py"))

BASES = {
    "manuscript": ("CA_S4h", "CA_LAI70", "CA_LAI20"),
    "old": ("old_S4h", "old_LAI70", "old_LAI20"),
    "new": ("new_S4h", "new_LAI70", "new_LAI20"),
    "oA": ("oA_S4h", "oA_LAI70", "oA_LAI20"),
    "nA": ("nA_S4h", "nA_LAI70", "nA_LAI20"),
}
YEARS = {"2017": ("2017-01-01", "2017-12-31"),
         "2018": ("2018-01-01", "2018-12-31")}
STITCH = {"all": ("2017-02-01", "2018-12-31"),
          "2017": ("2017-02-01", "2017-12-31"),
          "2018": ("2018-01-01", "2018-12-31")}
METRICS = ("NSE", "r", "RMSE", "PBIAS")


def load_manuscript_module():
    assert EVAL_SCRIPT.exists(), f"evaluation script missing: {EVAL_SCRIPT}"
    spec = importlib.util.spec_from_file_location("stitch_v2", EVAL_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def flow_path(mod, tag):
    p = mod.get_file_paths(tag)
    if not p.exists():
        p = CASE / "output" / f"{tag}Stream.Flow"
    return p


def run_metrics(mod, tag, df_obs):
    p = flow_path(mod, tag)
    if not p.exists():
        return None, None
    df = mod.parse_streamflow(p, mod.WATERSHED_AREA_M2)
    df = df.merge(df_obs, left_index=True, right_index=True, how="left")
    rows = {}
    for year, (a, b) in YEARS.items():
        rows[year] = mod.calc_metrics(df.loc[a:b, "Q_mm_sim"],
                                      df.loc[a:b, "Q_mm_obs"])
    return df, rows


def stitched_metrics(mod, s20, s70):
    st = pd.concat([s20.loc["2017-01-01":"2017-12-31"],
                    s70.loc["2018-01-01":"2018-12-31"]])
    return {k: mod.calc_metrics(st.loc[a:b, "Q_mm_sim"],
                                st.loc[a:b, "Q_mm_obs"])
            for k, (a, b) in STITCH.items()}


def fmt(m):
    return (f"NSE {m['NSE']:6.3f}  r {m['r']:.3f}  RMSE {m['RMSE']:.3f}  "
            f"PBIAS {m['PBIAS']:+5.1f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bases", nargs="+", default=list(BASES),
                    choices=list(BASES))
    args = ap.parse_args()
    mod = load_manuscript_module()
    print(f"manuscript code: {EVAL_SCRIPT}")
    print(f"observed: {mod.OBS_CSV}; basin area {mod.WATERSHED_AREA_M2:g} m2")
    df_obs = mod.load_observed()
    records = []
    stitched = {}
    for base in args.bases:
        tags = BASES[base]
        series = {}
        print("=" * 72)
        print(f"base {base}: {tags}")
        for tag in tags:
            df, rows = run_metrics(mod, tag, df_obs)
            if df is None:
                print(f"  {tag}: Stream.Flow missing, base skipped")
                break
            series[tag] = df
            for year, m in rows.items():
                print(f"  {tag:10s} {year}: {fmt(m)}")
                records.append(dict(base=base, run=tag, window=year, **m))
        if len(series) < 3:
            continue
        st = stitched_metrics(mod, series[tags[2]], series[tags[1]])
        stitched[base] = st
        print(f"  stitched (2017 {tags[2]}, 2018 {tags[1]}):")
        for k, m in st.items():
            print(f"    {k:5s} {fmt(m)}")
            records.append(dict(base=base, run="stitched", window=k, **m))
    if "manuscript" in stitched:
        print("=" * 72)
        print("stitched overall, difference from the manuscript "
              "(NSE, r, RMSE, PBIAS as printed by calc_metrics)")
        ref = stitched["manuscript"]["all"]
        for base, st in stitched.items():
            if base == "manuscript":
                continue
            d = {k: st["all"][k] - ref[k] for k in METRICS}
            print(f"  {base:10s} " + "  ".join(
                f"{k} {d[k]:+.3f}" for k in METRICS))
    out = CASE / "output" / "tierE_eval_CA.csv"
    pd.DataFrame(records).to_csv(out, index=False)
    print(f"table written: {out}")


if __name__ == "__main__":
    main()
