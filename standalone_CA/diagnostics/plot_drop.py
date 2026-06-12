# -*- coding: utf-8 -*-
# =====================================================================
# plot_drop.py  -  plot the constant stream drop sweep
#
# Reads the CSV written by drop_analysis.py and plots the absolute t-statistic
# against support area. The constant-drop criterion is the line |t| = 2; the
# threshold is the smallest support area whose |t| drops below it. The passing
# band (every threshold with |t| < 2) is shaded. The objective threshold and the
# current visual threshold are marked, so the figure shows at a glance that the
# visual choice sits inside the band the geomorphic law allows.
#
# Headless: uses the Agg backend, writes a PNG, no display needed.
#
# Run:  python3 plot_drop.py [drop_sweep.csv] [--out drop_sweep.png]
# =====================================================================

import sys
import csv
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Current visual threshold, for the marker. Resolve from paths.py when the
# script runs inside the package; fall back to the CA value otherwise.
try:
    from paths import STREAM_SOURCE_AREA_M2
    VISUAL_AREA_KM2 = float(STREAM_SOURCE_AREA_M2) / 1e6
except Exception:
    VISUAL_AREA_KM2 = 47571.5 / 1e6


def load_sweep(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            def f(key):
                v = r[key]
                try:
                    return float(v)
                except (ValueError, TypeError):
                    return float("nan")
            rows.append(dict(
                cells=int(r["cells"]),
                area_km2=f("area_km2"),
                t_abs=f("t_abs"),
                passes=int(r["passes"]),
            ))
    return rows


def main():
    ap = argparse.ArgumentParser(description="Plot the stream drop sweep")
    ap.add_argument("csv", nargs="?", default="drop_sweep.csv")
    ap.add_argument("--out", default="drop_sweep.png")
    args = ap.parse_args()

    rows = load_sweep(args.csv)

    # Series for the curve. Keep only finite |t| so the line breaks where the
    # test becomes undefined (no first-order streams left).
    finite = [r for r in rows if r["t_abs"] == r["t_abs"]]  # drop NaN
    x = [r["area_km2"] for r in finite]
    y = [r["t_abs"] for r in finite]

    passing = [r for r in rows if r["passes"] == 1]
    obj_area = min((r["area_km2"] for r in passing), default=None)
    band_lo = min((r["area_km2"] for r in passing), default=None)
    band_hi = max((r["area_km2"] for r in passing), default=None)

    fig, ax = plt.subplots(figsize=(7.0, 4.4))

    # passing band
    if band_lo is not None:
        ax.axvspan(band_lo, band_hi, color="#cfe8d4", alpha=0.6, zorder=0,
                   label="passing band (|t| < 2)")

    # the curve
    ax.plot(x, y, "-o", color="#1f4e79", markersize=4, linewidth=1.4,
            zorder=3, label="|t| first vs higher order")

    # criterion line
    ax.axhline(2.0, color="#c0392b", linestyle="--", linewidth=1.2, zorder=2,
               label="constant-drop criterion (|t| = 2)")

    # objective threshold (annotation placed high, left of its line)
    if obj_area is not None:
        ax.axvline(obj_area, color="#27772f", linewidth=1.4, zorder=2)
        ax.annotate(f"objective\n{obj_area:.4f} km$^2$",
                    xy=(obj_area, ax.get_ylim()[1] * 0.88),
                    xytext=(-8, 0), textcoords="offset points",
                    color="#27772f", fontsize=9, va="top", ha="right")

    # visual threshold (annotation placed lower, right of its line) so the two
    # labels do not collide when the thresholds are close together
    ax.axvline(VISUAL_AREA_KM2, color="#d68910", linewidth=1.4, zorder=2)
    ax.annotate(f"visual\n{VISUAL_AREA_KM2:.4f} km$^2$",
                xy=(VISUAL_AREA_KM2, ax.get_ylim()[1] * 0.62),
                xytext=(8, 0), textcoords="offset points",
                color="#d68910", fontsize=9, va="top", ha="left")

    ax.set_xlabel("support area  (km$^2$)")
    ax.set_ylabel("|t|  (first order vs higher order stream drop)")
    ax.set_title("Constant stream drop analysis, CA")
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
