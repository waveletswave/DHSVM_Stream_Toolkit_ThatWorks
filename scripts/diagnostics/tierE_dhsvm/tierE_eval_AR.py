# -*- coding: utf-8 -*-
# =====================================================================
# tierE_eval_AR.py  -  the manuscript's AR streamflow metrics (Fig 2a),
#                      recomputed with the manuscript's own code for the
#                      reruns
#
# Imports SpongeBurn/24_FIG2_DHSVM/10_fig5_CA_AR_unburned.py by path and
# uses its parse_streamflow, load_ar_observed, window_metrics,
# AR_AREA_M2, AR_OBS_CSV, WIN_2017, WIN_2018, EVAL_START and EVAL_END
# unchanged. For each prefix (S4h = the manuscript run, oA_S4h = the
# control rerun, nA_S4h = April inputs with the Tier E network, new_S4h
# = the current pipeline outputs) it prints the 2017, 2018 and overall
# metrics as the figure script computes them (2017 window 2017-02-01 to
# 2017-12-31, overall 2017-02-01 to 2018-12-31) and the differences from
# the manuscript. Missing prefixes are skipped. Writes
# output/tierE_eval_AR.csv.
#
# Run with the Python environment used for the manuscript scripts:
#   python3 tierE_eval_AR.py
#   python3 tierE_eval_AR.py --prefixes S4h oA_S4h nA_S4h new_S4h
# =====================================================================

import argparse
import importlib.util
import os
from pathlib import Path

import pandas as pd

CASE = Path(os.environ.get(
    "TIERE_CASE",
    "/Users/benthosyy/Desktop/CodeBits/DHSVM-PNNL-2025/TestCase/AR"))
FIG_SCRIPT = Path(os.environ.get(
    "TIERE_FIG_SCRIPT",
    "/Users/benthosyy/Desktop/CodeBits/DHSVM-PNNL-2025/TestCase/LAI/"
    "SpongeBurn/24_FIG2_DHSVM/10_fig5_CA_AR_unburned.py"))
PREFIXES = ["S4h", "oA_S4h", "nA_S4h", "new_S4h", "wA_S4h"]
METRICS = ("NSE", "r", "RMSE", "PBIAS")


def load_fig_module():
    assert FIG_SCRIPT.exists(), f"figure script missing: {FIG_SCRIPT}"
    spec = importlib.util.spec_from_file_location("fig2_ar", FIG_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fmt(m):
    return (f"NSE {m['NSE']:6.3f}  r {m['r']:.3f}  RMSE {m['RMSE']:.3f}  "
            f"PBIAS {m['PBIAS']:+5.1f}  n {m['n']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefixes", nargs="+", default=PREFIXES)
    args = ap.parse_args()
    mod = load_fig_module()
    windows = {"2017": mod.WIN_2017, "2018": mod.WIN_2018,
               "all": (mod.EVAL_START, mod.EVAL_END)}
    print(f"manuscript code: {FIG_SCRIPT}")
    print(f"observed: {mod.AR_OBS_CSV}; basin area {mod.AR_AREA_M2:g} m2")
    df_obs = mod.load_ar_observed(mod.AR_OBS_CSV)
    records = []
    results = {}
    for prefix in args.prefixes:
        flow = CASE / "output" / f"{prefix}Stream.Flow"
        if not flow.exists():
            print(f"  {prefix}: Stream.Flow missing, skipped")
            continue
        df = mod.parse_streamflow(flow, mod.AR_AREA_M2)
        df = df.merge(df_obs, left_index=True, right_index=True, how="left")
        df = df.loc["2017-01-01":"2018-12-31"]
        results[prefix] = {}
        for name, (a, b) in windows.items():
            m = mod.window_metrics(df, a, b)
            results[prefix][name] = m
            print(f"  {prefix:8s} {name:5s} {fmt(m)}")
            records.append(dict(prefix=prefix, window=name, **m))
    if "S4h" in results:
        print("=" * 72)
        print("overall (2017-02 to 2018-12), difference from the manuscript "
              "run S4h, as printed by calc_metrics")
        ref = results["S4h"]["all"]
        for prefix, r in results.items():
            if prefix == "S4h":
                continue
            d = {k: r["all"][k] - ref[k] for k in METRICS}
            print(f"  {prefix:8s} " + "  ".join(
                f"{k} {d[k]:+.3f}" for k in METRICS))
    out = CASE / "output" / "tierE_eval_AR.csv"
    pd.DataFrame(records).to_csv(out, index=False)
    print(f"table written: {out}")


if __name__ == "__main__":
    main()
