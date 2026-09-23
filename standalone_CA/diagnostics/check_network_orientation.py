# -*- coding: utf-8 -*-
# =====================================================================
# check_network_orientation.py  -  Tier E read-only diagnostic
#
# Question: is the stream network written by stream_network.py oriented
# with the flow?  r.watershed (no -a) writes negative accumulation for
# cells that may receive flow from outside the region; on a DEM clipped
# to the basin the whole network is negative.  stream_network.py compares
# raw values, so the end with the SMALLER magnitude is taken as
# downstream.  This script measures the consequences on the real
# pipeline outputs.  It writes nothing into the pipeline output
# directory; a text report goes to stdout and a per-segment CSV to
# --report-dir.
#
# Part A reads only the written files (rasters + stream.*.dat).
# Part B re-runs the unmodified build_directed_by_FA on
# streamfile_attr.shp and compares the orientation it chooses with the
# orientation implied by |accumulation| and by elevation, per segment.
# Part C walks the D8 direction raster (flow_dir.tif, GRASS encoding
# 1..8 counter-clockwise from NE, negative = drains out of the region)
# from both ends of every segment: the end from which the walk reaches
# the other end is upstream.  It also lists the stream-raster outlets
# (stream cells whose D8 successor is not a stream cell).
#
# Invariants (PASS/FAIL; exit code 1 only with --strict):
#   I1  the outlet segment (down = 0) holds the max-|acc| channel cell
#   I2  the outlet segment's downstream endpoint is the lowest channel
#       cell
#   I3  every segment: |acc| at the downstream end >= at the upstream end
#   I4  every segment: z at the upstream end >= at the downstream end
#       (informational: r.watershed does not fill sinks)
#   I5  every segment: in-degree <= 3
#   I6  exactly one sink before the sink merge
#
# Run with the geo env active, from anywhere:
#   python3 check_network_orientation.py --out-dir <pipeline OUT> \
#       --pipeline-dir <repo>/standalone_CA/pipeline \
#       --label CA_28m --report-dir <where the CSV goes>
# --pipeline-dir defaults to ../pipeline relative to this file, so the
# script also runs unchanged from standalone_CA/diagnostics/.
# Overrides: --network, --map, --class-file, --streamfile-attr, --dem,
#   --acc, --stream-raster, --flow-dir, --no-recompute (skip Part B),
#   --strict.
# =====================================================================

import argparse
import csv
import hashlib
import math
import os
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


# ----------------------------- helpers --------------------------------
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_raster(path):
    import rasterio
    with rasterio.open(path) as ds:
        arr = ds.read(1).astype(np.float64)
        return arr, ds.nodata, ds.transform, ds.crs, ds.shape


def valid_mask(arr, nodata):
    m = np.isfinite(arr)
    if nodata is not None and np.isfinite(nodata):
        m &= arr != nodata
    return m


def boundary_mask(valid):
    """Valid cells with an invalid or off-grid 8-neighbour."""
    padded = np.pad(valid, 1, constant_values=False)
    nb_ok = np.ones_like(valid, dtype=bool)
    nr, nc = valid.shape
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            nb_ok &= padded[1 + dr:1 + dr + nr, 1 + dc:1 + dc + nc]
    return valid & ~nb_ok


def parse_network(path):
    """stream.network.dat: segid nin slope length class down [SAVE]."""
    rows = {}
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            save = None
            if len(p) >= 7 and p[6].upper().startswith("SAVE"):
                save = " ".join(p[6:])
            rows[int(p[0])] = dict(nin=int(p[1]), slope=float(p[2]),
                                   length=float(p[3]), cls=int(p[4]),
                                   down=int(p[5]), save=save)
    return rows


def parse_map(path):
    """stream.map.dat: col row id length height width aspect."""
    recs = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            recs.append(dict(col=int(p[0]), row=int(p[1]), id=int(p[2]),
                             length=float(p[3]), height=float(p[4]),
                             width=float(p[5]), aspect=float(p[6])))
    return recs


def parse_classes(path):
    classes = {}
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            classes[int(p[0])] = dict(width=float(p[1]), depth=float(p[2]),
                                      n=float(p[3]), inf=float(p[4]))
    return classes


def fmt(v, nd=3):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def git_head(path):
    try:
        r = subprocess.run(["git", "-C", str(path), "rev-parse",
                            "--short", "HEAD"],
                           capture_output=True, text=True)
        return r.stdout.strip() or "NA"
    except Exception:
        return "NA"


# ----------------------------- arguments ------------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="Tier E stream-network orientation diagnostic")
    ap.add_argument("--out-dir", required=True,
                    help="pipeline output directory (paths.OUT)")
    ap.add_argument("--label", default="run")
    ap.add_argument("--report-dir", default=".")
    ap.add_argument("--pipeline-dir", default=str(HERE.parent / "pipeline"),
                    help="standalone_CA/pipeline of the Toolkit clone")
    ap.add_argument("--dem", default=None)
    ap.add_argument("--acc", default=None)
    ap.add_argument("--stream-raster", default=None)
    ap.add_argument("--flow-dir", default=None)
    ap.add_argument("--streamfile-attr", default=None)
    ap.add_argument("--network", default=None)
    ap.add_argument("--map", default=None)
    ap.add_argument("--class-file", default=None)
    ap.add_argument("--no-recompute", action="store_true",
                    help="skip Part B")
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 if a hard invariant fails")
    return ap.parse_args()


def resolve_inputs(args):
    out = Path(args.out_dir).resolve()
    streams = out / "DHSVM_input_streams"

    def pick(v, default):
        return Path(v) if v else default
    return dict(
        out=out,
        dem=pick(args.dem, out / "elev_clipped.tif"),
        acc=pick(args.acc, out / "flow_acc.tif"),
        srast=pick(args.stream_raster, out / "stream_raster.tif"),
        fdir=pick(args.flow_dir, out / "flow_dir.tif"),
        attr=pick(args.streamfile_attr, out / "streamfile_attr.shp"),
        net=pick(args.network, streams / "stream.network.dat"),
        map=pick(args.map, streams / "stream.map.dat"),
        cls=pick(args.class_file, streams / "stream.class.dat"),
    )


# ----------------------------- part A ---------------------------------
def part_a(inp, args):
    """Written files only.  Returns the state Part B needs."""
    dem, dem_nd, T, crs, shape = read_raster(inp["dem"])
    acc, acc_nd, Ta, _, shape_a = read_raster(inp["acc"])
    srast, sr_nd, _, _, shape_s = read_raster(inp["srast"])
    if shape != shape_a or shape != shape_s:
        raise SystemExit(f"[error] raster shapes differ: dem {shape} "
                         f"acc {shape_a} stream {shape_s}")
    if [round(v, 6) for v in T][:6] != [round(v, 6) for v in Ta][:6]:
        print("[warn] dem and acc transforms differ")
    nrows, ncols = shape
    px, py = abs(T.a), abs(T.e)
    cell_area = px * py

    vdem = valid_mask(dem, dem_nd)
    vacc = valid_mask(acc, acc_nd) & vdem
    stream = valid_mask(srast, sr_nd) & (srast != 0) & vdem
    bnd = boundary_mask(vdem)
    absacc = np.abs(acc)

    print()
    print(f"grid: {nrows} rows x {ncols} cols, cell {px:.3f} x {py:.3f} m "
          f"(area {cell_area:.2f} m2), CRS {crs}")
    print(f"valid DEM cells {int(vdem.sum())}, valid acc cells "
          f"{int(vacc.sum())}, stream-raster cells {int(stream.sum())}, "
          f"boundary cells {int(bnd.sum())}")

    # D1
    a = acc[vacc]
    s = acc[stream & vacc]
    print()
    print("D1  accumulation sign")
    print(f"    all valid cells: min {a.min():.1f}  max {a.max():.1f}  "
          f"negative {int((a < 0).sum())}  zero {int((a == 0).sum())}  "
          f"positive {int((a > 0).sum())}")
    print(f"    stream-raster cells: {int((s < 0).sum())} of {s.size} "
          f"negative ({100.0 * (s < 0).sum() / max(s.size, 1):.1f} %)")

    # D2
    am = np.where(vacc, absacc, -np.inf)
    r_max, c_max = np.unravel_index(int(np.argmax(am)), shape)
    acc_max = float(absacc[r_max, c_max])
    z_max = float(dem[r_max, c_max])
    ds = np.where(stream, dem, np.inf)
    r_low, c_low = np.unravel_index(int(np.argmin(ds)), shape)
    z_low = float(dem[r_low, c_low])
    dv = np.where(vdem, dem, np.inf)
    r_dm, c_dm = np.unravel_index(int(np.argmin(dv)), shape)
    print()
    print("D2  pour point candidates (0-based row, col; row 0 = north, "
          "as in stream.map.dat)")
    print(f"    max |acc| cell: (row {r_max}, col {c_max})  "
          f"|acc| {acc_max:.1f} cells  z {z_max:.2f} m  "
          f"on mask boundary: {bool(bnd[r_max, c_max])}  "
          f"in stream raster: {bool(stream[r_max, c_max])}")
    print(f"    lowest stream-raster cell: (row {r_low}, col {c_low})  "
          f"z {z_low:.2f} m  |acc| {absacc[r_low, c_low]:.1f}")
    print(f"    lowest valid DEM cell: (row {r_dm}, col {c_dm})  "
          f"z {float(dem[r_dm, c_dm]):.2f} m  "
          f"in stream raster: {bool(stream[r_dm, c_dm])}")
    print(f"    distance max|acc| cell -> lowest stream cell: "
          f"{math.hypot(r_max - r_low, c_max - c_low):.1f} cells")

    # D3 written network
    net = parse_network(inp["net"])
    recs = parse_map(inp["map"])
    classes = parse_classes(inp["cls"]) if inp["cls"].exists() else {}
    cells_by_id = defaultdict(list)
    owner = defaultdict(list)
    for r in recs:
        cells_by_id[r["id"]].append((r["row"], r["col"]))
        owner[(r["row"], r["col"])].append(r["id"])
    indeg = Counter()
    for sid, row in net.items():
        if row["down"] != 0:
            indeg[row["down"]] += 1
    outlets = [sid for sid, row in net.items() if row["down"] == 0]
    saves = [sid for sid, row in net.items() if row["save"]]
    missing_map = [sid for sid in net if sid not in cells_by_id]
    missing_net = sorted(set(cells_by_id) - set(net))
    multi = sum(1 for v in owner.values() if len(v) > 1)
    off = sum(1 for (rr, cc) in owner
              if not (0 <= rr < nrows and 0 <= cc < ncols
                      and stream[rr, cc]))
    print()
    print("D3  written network (stream.network.dat / stream.map.dat)")
    print(f"    segments {len(net)}, map records {len(recs)}, distinct "
          f"map cells {len(owner)}, cells under 2+ segments {multi}, "
          f"map cells not in stream raster {off}")
    if missing_map or missing_net:
        print(f"    [warn] ids in network not in map: {missing_map}; "
              f"in map not in network: {missing_net}")
    print(f"    outlet segments (down = 0): {outlets}   "
          f"SAVE flags: {saves if saves else 'none'}")
    nin_hist = dict(sorted(Counter(r["nin"] for r in net.values()).items()))
    ind_hist = dict(sorted(Counter(indeg[s] for s in net).items()))
    print(f"    order (nin) histogram: {nin_hist}")
    print(f"    in-degree histogram: {ind_hist}")
    if indeg:
        print(f"    highest in-degree segments: {indeg.most_common(3)}")
    if classes:
        used = dict(sorted(Counter(r["cls"] for r in net.values()).items()))
        table = ", ".join(f"{k}: W {v['width']} D {v['depth']} n {v['n']}"
                          for k, v in sorted(classes.items()))
        print(f"    classes used: {used}   class table: {table}")

    print()
    print("D3  outlet segment(s) as written")
    for sid in outlets:
        row = net[sid]
        cells = [(rr, cc) for (rr, cc) in cells_by_id.get(sid, [])
                 if 0 <= rr < nrows and 0 <= cc < ncols]
        print(f"    segment {sid}: nin {row['nin']}  "
              f"slope {row['slope']:.5f}  length {row['length']:.3f}  "
              f"class {row['cls']}  in-degree {indeg[sid]}")
        if cells:
            av = [absacc[rr, cc] for rr, cc in cells if vacc[rr, cc]]
            zv = [dem[rr, cc] for rr, cc in cells if vdem[rr, cc]]
            dist = min(math.hypot(rr - r_max, cc - c_max)
                       for rr, cc in cells)
            print(f"      map cells {len(cells)}  |acc| on cells "
                  f"{fmt(min(av) if av else None, 1)}.."
                  f"{fmt(max(av) if av else None, 1)} "
                  f"(basin max {acc_max:.1f})  z "
                  f"{fmt(min(zv) if zv else None, 2)}.."
                  f"{fmt(max(zv) if zv else None, 2)} m "
                  f"(lowest stream cell {z_low:.2f})")
            print(f"      mean (row {np.mean([c[0] for c in cells]):.1f}, "
                  f"col {np.mean([c[1] for c in cells]):.1f})  "
                  f"nearest cell to max|acc| cell: {dist:.1f} cells")

    print()
    print("D4  segment(s) holding the max-|acc| cell and the lowest "
          "stream cell")
    owners_max = owner.get((r_max, c_max), [])
    owners_low = owner.get((r_low, c_low), [])
    for tag, owners in (("max |acc| cell", owners_max),
                        ("lowest stream cell", owners_low)):
        if not owners:
            print(f"    {tag}: not in stream.map.dat")
        for sid in owners:
            row = net.get(sid)
            if row is None:
                print(f"    {tag}: segment {sid} not in network")
                continue
            verdict = ("IS the written outlet" if row["down"] == 0
                       else "NOT the written outlet")
            print(f"    {tag}: segment {sid}  down {row['down']}  "
                  f"nin {row['nin']}  in-degree {indeg[sid]}  "
                  f"<-- {verdict}")

    max_indeg = max(indeg.values()) if indeg else 0
    print()
    print(f"D5  max in-degree {max_indeg}  "
          "(dendritic expectation: 2, at most 3)")

    results = {}
    results["I1 outlet segment holds the max-|acc| cell (informational: "
            "MFD accumulation is not monotone)"] = (
        bool(outlets) and any(o in owners_max for o in outlets))
    results["I2 lowest stream cell lies in an outlet segment"] = (
        bool(outlets) and any(o in owners_low for o in outlets))
    results["I5 max in-degree <= 3"] = max_indeg <= 3
    st_outlets = outlets

    return dict(dem=dem, acc=acc, absacc=absacc, vacc=vacc, vdem=vdem,
                stream=stream, T=T, px=px, py=py, nrows=nrows,
                ncols=ncols, cell_area=cell_area, net=net,
                acc_max=acc_max, z_low=z_low, results=results,
                outlets=st_outlets)


# ----------------------------- part B ---------------------------------
def part_b(inp, args, st, gdf):
    try:
        import stream_network as sn
        sn.build_directed_by_FA
    except (ImportError, AttributeError):
        print()
        print("B   skipped: stream_network.build_directed_by_FA is not in "
              "the pipeline directory (retired in Tier E)")
        return None, None
    print()
    print("B   recompute with the unmodified build_directed_by_FA")
    facc_arr, facc_inv = sn._open_raster(str(inp["acc"]))
    elev_arr, elev_inv = sn._open_raster(str(inp["dem"]))
    px, py = st["px"], st["py"]
    diag = math.hypot(px, py)
    (geoms, down_map, indeg_b, shreve, prop, az, topo, new_id,
     up_end, down_end) = sn.build_directed_by_FA(gdf, facc_arr, facc_inv,
                                                 px, py, diag)
    fids = list(range(len(gdf)))
    tol = sn.FA_TOL_MULT * diag
    net = st["net"]
    results = st["results"]

    # B2 written vs recomputed
    mism_down = mism_nin = 0
    for fid in fids:
        w = net.get(new_id[fid])
        if w is None:
            mism_down += 1
            continue
        d = down_map[fid]
        if w["down"] != (0 if d == -1 else new_id[d]):
            mism_down += 1
        if w["nin"] != max(1, prop[fid]):
            mism_nin += 1
    print(f"    B2 written vs recomputed: down mismatches {mism_down}, "
          f"nin mismatches {mism_nin} (0 and 0: the written files come "
          "from this code on these inputs)")

    # B3 per-segment orientation
    rows = []
    ends = {}
    n_ca = n_cz = n_az = n_at = n_zt = n_i3 = n_i4 = 0
    for fid in fids:
        line = sn._line_coords(geoms[fid])
        if len(line) < 2:
            continue
        a, b = line[0], line[-1]
        ends[fid] = (a, b, len(line))
        acc_a = sn._sample(facc_arr, facc_inv, a[0], a[1], 0.0)
        acc_b = sn._sample(facc_arr, facc_inv, b[0], b[1], 0.0)
        z_a = sn._sample(elev_arr, elev_inv, a[0], a[1], float("nan"))
        z_b = sn._sample(elev_arr, elev_inv, b[0], b[1], float("nan"))
        code = "b" if down_end[fid] == (b[0], b[1]) else "a"
        if abs(acc_a) == abs(acc_b):
            by_acc = "tie"
            n_at += 1
        else:
            by_acc = "b" if abs(acc_b) > abs(acc_a) else "a"
        if not (np.isfinite(z_a) and np.isfinite(z_b)) or z_a == z_b:
            by_z = "tie"
            n_zt += 1
        else:
            by_z = "b" if z_b < z_a else "a"
        n_ca += int(by_acc != "tie" and code == by_acc)
        n_cz += int(by_z != "tie" and code == by_z)
        n_az += int(by_acc != "tie" and by_z != "tie" and by_acc == by_z)
        acc_dn = abs(sn._sample(facc_arr, facc_inv, *down_end[fid], 0.0))
        acc_up = abs(sn._sample(facc_arr, facc_inv, *up_end[fid], 0.0))
        z_dn = sn._sample(elev_arr, elev_inv, *down_end[fid], float("nan"))
        z_up = sn._sample(elev_arr, elev_inv, *up_end[fid], float("nan"))
        i3 = acc_dn >= acc_up
        i4 = (not (np.isfinite(z_dn) and np.isfinite(z_up))) or z_up >= z_dn
        n_i3 += int(i3)
        n_i4 += int(i4)
        d = down_map[fid]
        rows.append(dict(
            fid=fid, segid=new_id[fid],
            down_segid=(0 if d == -1 else new_id[d]),
            nin=max(1, prop[fid]), indeg=indeg_b[fid],
            a_x=round(a[0], 3), a_y=round(a[1], 3),
            b_x=round(b[0], 3), b_y=round(b[1], 3),
            acc_a_raw=acc_a, acc_b_raw=acc_b, z_a=z_a, z_b=z_b,
            code_down_end=code, down_by_abs_acc=by_acc,
            down_by_elevation=by_z,
            code_eq_acc=(code == by_acc), code_eq_z=(code == by_z),
            abs_acc_down_end=acc_dn, abs_acc_up_end=acc_up,
            z_down_end=z_dn, z_up_end=z_up,
            I3_acc_monotone=i3, I4_z_monotone=i4,
            azimuth_deg=round(az[fid], 1),
            length_m=round(float(geoms[fid].length), 3)))
    n = len(rows)
    print(f"    B3 segments {n}: code agrees with |acc| orientation on "
          f"{n_ca}, with elevation orientation on {n_cz}; |acc| and "
          f"elevation agree on {n_az}; |acc| ties {n_at}, z ties {n_zt}")
    print(f"       I3 (|acc| down >= up) holds on {n_i3}/{n};  "
          f"I4 (z up >= down) holds on {n_i4}/{n}")

    # B4 sinks before the merge
    outlet_fids = [fid for fid in fids if down_map[fid] == -1]
    merged = []
    if len(outlet_fids) == 1:
        o = outlet_fids[0]
        ux, uy = up_end[o]
        for fid in fids:
            if fid == o or down_map[fid] != o:
                continue
            x, y = down_end[fid]
            probes = [(x, y)]
            for s in sn.AHEAD_STEPS:
                probes.append(sn._step_along_az(x, y, az[fid],
                                                s * max(px, py)))
            if min(math.hypot(qx - ux, qy - uy) for qx, qy in probes) > tol:
                merged.append(fid)
    sinks_before = len(outlet_fids) + len(merged)
    print(f"    B4 outlet fid {outlet_fids} (segid "
          f"{[new_id[f] for f in outlet_fids]}), links to it classified "
          f"as merged sinks: {len(merged)} -> sinks before merge = "
          f"{sinks_before}")
    if merged:
        print(f"       merged (fid, segid): "
              f"{[(f, new_id[f]) for f in merged]}")
    if len(outlet_fids) == 1:
        o = outlet_fids[0]
        z_o = sn._sample(elev_arr, elev_inv, *down_end[o], float("nan"))
        acc_o = abs(sn._sample(facc_arr, facc_inv, *down_end[o], 0.0))
        T = st["T"]
        col_o = int(math.floor((down_end[o][0] - T.c) / px))
        row_o = int(math.floor((T.f - down_end[o][1]) / py))
        print(f"       outlet downstream endpoint: (row {row_o}, "
              f"col {col_o})  z {z_o:.2f} m  |acc| {acc_o:.1f}  "
              f"(lowest stream cell z {st['z_low']:.2f}, basin max "
              f"|acc| {st['acc_max']:.1f})")
    results["I3 |acc| monotone on every segment (informational)"] = (
        n_i3 == n)
    results["I4 elevation monotone on every segment (informational)"] = (
        n_i4 == n)
    results["I6 exactly one sink before merge (old code path)"] = (
        sinks_before == 1)

    # B5 coincident endpoints at downstream ends
    hist = Counter()
    for fid in fids:
        x, y = down_end[fid]
        k = sum(1 for g in fids if g != fid
                and math.hypot(up_end[g][0] - x, up_end[g][1] - y) < 1e-6)
        hist[k] += 1
    print(f"    B5 other segments' UP endpoints coincident with each "
          f"segment's DOWN endpoint: {dict(sorted(hist.items()))}  "
          "(2+ = the d=0 tie; expected none with correct orientation)")

    return rows, ends


# ----------------------------- part C ---------------------------------
D8 = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1),
      5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}


def walk_d8(fdir, valid, start, target, max_steps):
    """Follow the D8 raster from start; True if target is reached."""
    nr, nc = fdir.shape
    r, c = start
    for step in range(max_steps):
        if (r, c) == target:
            return True, step, "reached"
        d = fdir[r, c]
        if not np.isfinite(d) or int(abs(d)) not in D8:
            return False, step, "pit"
        dr, dc = D8[int(abs(d))]
        r2, c2 = r + dr, c + dc
        if not (0 <= r2 < nr and 0 <= c2 < nc) or not valid[r2, c2]:
            return False, step, "exit"
        r, c = r2, c2
    return False, max_steps, "maxsteps"


def part_c(inp, st, rows, ends):
    print()
    print("C   D8 direction raster (flow_dir.tif) walked from both ends")
    fdir, fd_nd, Tf, _, shape_f = read_raster(inp["fdir"])
    if shape_f != st["dem"].shape:
        print(f"    [warn] flow_dir shape {shape_f} differs; Part C skipped")
        return
    vdem, stream, T = st["vdem"], st["stream"], st["T"]
    px, py = st["px"], st["py"]
    vals = np.abs(fdir[vdem & np.isfinite(fdir)])
    codes = dict(sorted(Counter(vals.astype(int).tolist()).items()))
    n_neg = int((fdir[vdem & np.isfinite(fdir)] < 0).sum())
    print(f"    codes on valid cells: {codes}  negative (drain out) {n_neg}")

    def rc(x, y):
        return (int(math.floor((T.f - y) / py)),
                int(math.floor((x - T.c) / px)))

    # stream-raster outlets: stream cells whose successor is not a stream
    nr, nc = fdir.shape
    outlets = []
    pits = 0
    for r, c in zip(*np.nonzero(stream)):
        d = fdir[r, c]
        if not np.isfinite(d) or int(abs(d)) not in D8:
            pits += 1
            outlets.append((r, c, "pit"))
            continue
        dr, dc = D8[int(abs(d))]
        r2, c2 = r + dr, c + dc
        inside = 0 <= r2 < nr and 0 <= c2 < nc and vdem[r2, c2]
        if not inside:
            outlets.append((r, c, "exit"))
        elif not stream[r2, c2]:
            outlets.append((r, c, "leaves stream raster"))
    outlets.sort(key=lambda t: -st["absacc"][t[0], t[1]])
    print(f"    stream-raster outlet cells by direction: {len(outlets)} "
          f"(pits {pits})")
    for r, c, why in outlets[:12]:
        print(f"      (row {r}, col {c})  |acc| {st['absacc'][r, c]:.1f}  "
              f"z {st['dem'][r, c]:.2f}  {why}")
    st["results"]["I7 stream raster has exactly one outlet cell "
                  "(informational)"] = (len(outlets) == 1)
    st["results"]["I9 outlet rows in the network == outlet cells in "
                  "the raster"] = (len(st["outlets"]) == len(outlets))

    n_det = n_da = n_dz = n_dc = n_amb = 0
    by_fid = {row["fid"]: row for row in rows}
    for fid, (a, b, nvert) in ends.items():
        row = by_fid.get(fid)
        if row is None:
            continue
        ra, rb = rc(*a), rc(*b)
        steps = 4 * nvert + 10
        ok_ab, s_ab, why_ab = walk_d8(fdir, vdem, ra, rb, steps)
        ok_ba, s_ba, why_ba = walk_d8(fdir, vdem, rb, ra, steps)
        if ok_ab and not ok_ba:
            by_dir = "b"
        elif ok_ba and not ok_ab:
            by_dir = "a"
        else:
            by_dir = "ambiguous"
        row["walk_a_to_b"] = f"{ok_ab}:{why_ab}:{s_ab}"
        row["walk_b_to_a"] = f"{ok_ba}:{why_ba}:{s_ba}"
        row["down_by_direction"] = by_dir
        if by_dir == "ambiguous":
            n_amb += 1
            continue
        n_det += 1
        n_da += int(row["down_by_abs_acc"] == by_dir)
        n_dz += int(row["down_by_elevation"] == by_dir)
        n_dc += int(row["code_down_end"] == by_dir)
    print(f"    segments oriented by the walk: {n_det}, ambiguous {n_amb}")
    print(f"    direction agrees with |acc| on {n_da}, with elevation on "
          f"{n_dz}, with the written orientation on {n_dc} "
          "(the written one is the digitised direction when Part B "
          "was skipped)")
    st["results"]["I8 every segment oriented by the D8 walk"] = (
        n_amb == 0 and n_det == len(rows))


def rows_from_geometry(gdf, st):
    """Per-segment rows when Part B cannot run: endpoints, |acc| and
    elevation at both ends, and the orientation each implies."""
    from shapely.geometry import MultiLineString
    T, dem, absacc = st["T"], st["dem"], st["absacc"]
    px, py = st["px"], st["py"]
    nr, nc = dem.shape

    def sample(arr, x, y):
        r = int(math.floor((T.f - y) / py))
        c = int(math.floor((x - T.c) / px))
        if 0 <= r < nr and 0 <= c < nc and np.isfinite(arr[r, c]):
            return float(arr[r, c])
        return float("nan")
    rows, ends = [], {}
    for fid, geom in enumerate(gdf.geometry):
        if geom is None or geom.is_empty:
            continue
        if isinstance(geom, MultiLineString):
            geom = max(geom.geoms, key=lambda g: g.length)
        line = list(geom.coords)
        a, b = line[0], line[-1]
        ends[fid] = (a, b, len(line))
        acc_a, acc_b = sample(absacc, *a), sample(absacc, *b)
        z_a, z_b = sample(dem, *a), sample(dem, *b)
        by_acc = ("tie" if acc_a == acc_b else
                  ("b" if acc_b > acc_a else "a"))
        by_z = ("tie" if z_a == z_b else ("b" if z_b < z_a else "a"))
        seg = int(gdf.iloc[fid]["segid"]) if "segid" in gdf.columns else -1
        rows.append(dict(fid=fid, segid=seg,
                         a_x=round(a[0], 3), a_y=round(a[1], 3),
                         b_x=round(b[0], 3), b_y=round(b[1], 3),
                         abs_acc_a=acc_a, abs_acc_b=acc_b, z_a=z_a, z_b=z_b,
                         code_down_end="b",
                         down_by_abs_acc=by_acc, down_by_elevation=by_z,
                         length_m=round(float(geom.length), 3)))
    return rows, ends


def write_csv(args, rows):
    report_dir = Path(args.report_dir).resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    csv_path = report_dir / f"tierE_segments_{args.label}.csv"
    keys = []
    for row in rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"    per-segment table: {csv_path}")


# ----------------------------- main -----------------------------------
def main():
    args = parse_args()
    inp = resolve_inputs(args)
    pipeline = Path(args.pipeline_dir).resolve()
    if not (pipeline / "paths.py").exists():
        raise SystemExit(f"[error] paths.py not found in {pipeline}; "
                         "pass --pipeline-dir")
    # keep paths.py from creating a stray output directory on import
    os.environ.setdefault("DHSVM_OUT", str(inp["out"]))
    sys.path.insert(0, str(pipeline))

    import rasterio
    import geopandas as gpd
    print("=" * 70)
    print(f"Tier E network orientation diagnostic   label={args.label}")
    print("=" * 70)
    print(f"pipeline dir {pipeline}")
    print(f"repo HEAD {git_head(pipeline)}   rasterio {rasterio.__version__}"
          f"   geopandas {gpd.__version__}   numpy {np.__version__}")
    print("inputs:")
    for tag in ("dem", "acc", "srast", "fdir", "attr", "net", "map", "cls"):
        p = inp[tag]
        state = f"sha256={sha256(p)[:16]}" if p.exists() else "MISSING"
        print(f"  {tag:6s} {p}  {state}")
    for tag in ("dem", "acc", "srast", "net", "map"):
        if not inp[tag].exists():
            raise SystemExit(f"[error] required input missing: "
                             f"{inp[tag]}")

    st = part_a(inp, args)

    gdf = None
    if inp["attr"].exists():
        gdf = gpd.read_file(inp["attr"])
        print()
        print("D7  streamfile_attr.shp")
        print(f"    features {len(gdf)}  columns {list(gdf.columns)}")
        if "meanmsq" in gdf.columns:
            mm = gdf["meanmsq"].to_numpy(dtype=float)
            uniq = np.unique(np.round(mm, 3))
            print(f"    meanmsq: {len(uniq)} unique value(s) "
                  f"{uniq[:10]}  min {mm.min():.2f}  max {mm.max():.2f}"
                  f"  cell_area {st['cell_area']:.2f}")
        for col in ("chanclass", "hydwidth", "hyddepth"):
            if col in gdf.columns:
                print(f"    {col}: {dict(Counter(gdf[col].tolist()))}")
    else:
        print()
        print(f"D7  streamfile_attr.shp missing at {inp['attr']}; "
              "Part B skipped")

    if gdf is not None and not args.no_recompute:
        rows, ends = part_b(inp, args, st, gdf)
        if rows is None:
            rows, ends = rows_from_geometry(gdf, st)
        if inp["fdir"].exists():
            part_c(inp, st, rows, ends)
        else:
            print()
            print(f"C   flow_dir.tif missing at {inp['fdir']}; skipped")
        write_csv(args, rows)

    print()
    print("invariants")
    failed = 0
    for k in sorted(st["results"]):
        ok = st["results"][k]
        print(f"    {'PASS' if ok else 'FAIL'}  {k}")
        failed += int((not ok) and ("informational" not in k))
    print()
    print(f"{args.label}: {failed} hard invariant(s) failed")
    if args.strict and failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
